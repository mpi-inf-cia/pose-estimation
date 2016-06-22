#pragma once
#ifndef POSE_SOLVER_CALLBACK_HXX
#define POSE_SOLVER_CALLBACK_HXX

#include <cstddef>
#include <limits>
#include <iostream>
#include <vector>
#include <string>

#include <andres/graph/complete-graph.hxx>
#include <andres/graph/components.hxx>
#include <andres/timer.hxx>

#include "problem.hxx"
#include "solution.hxx"

namespace pose {

template<class S = std::size_t>
struct FeasibleSolutionCallbackEmpty {
    typedef S size_type;
    typedef Solution<size_type> SolutionType;

    void operator()(const SolutionType& solution) const
        {}
};

template<
    class ILP_SOLVER,
    class FEASBILE_SOLUTION_CALLBACK = FeasibleSolutionCallbackEmpty<>
>
class PoseEstimator {
public:
    typedef ILP_SOLVER IlpSolver;
    typedef FEASBILE_SOLUTION_CALLBACK FeasibleSolutionCallback;
    typedef std::size_t size_type;
    typedef Problem<double, size_type> ProblemType;
    typedef Solution<size_type> SolutionType;
    typedef andres::Timer<double> TimerType;

    PoseEstimator(ProblemType&, SolutionType&, const SolutionType &, FeasibleSolutionCallback = FeasibleSolutionCallback(), const bool = false, const int = 36000, const double = 1.0 / 255.0);

private:
    typedef typename ProblemType::JoinIndexType JoinIndexType;
    typedef typename ProblemType::JoinDataType JoinDataType;
    typedef typename ProblemType::JoinMapType JoinMapType;
    typedef typename andres::graph::CompleteGraph<> CompleteGraphType;
    typedef andres::graph::ComponentsBySearch<CompleteGraphType> ComponentsType;

    struct SubgraphMask {
        typedef PoseEstimator<IlpSolver, FeasibleSolutionCallback> PoseEstimatorType;

        SubgraphMask(const PoseEstimatorType& poseEstimator)
            : poseEstimator_(poseEstimator) {}
        bool vertex(const std::size_t v) const
            { return true; }
        bool edge(const std::size_t e) const
            {
                const size_type v0 = poseEstimator_.detectionGraph_.vertexOfEdge(e, 0);
                const size_type v1 = poseEstimator_.detectionGraph_.vertexOfEdge(e, 1);
                return poseEstimator_.ilpSolver_.label(poseEstimator_.y(v0, v1)) == 1;
            }
        const PoseEstimatorType& poseEstimator_;
    };

    class Callback
    :   public IlpSolver::Callback
    {
    public:        
        typedef typename IlpSolver::Callback IlpSolverCallback;
        typedef PoseEstimator<IlpSolver, FeasibleSolutionCallback> PoseEstimatorType;
        typedef typename PoseEstimatorType::size_type size_type;
        typedef typename PoseEstimatorType::ComponentsType ComponentsType;

        struct SubgraphMask {
            typedef Callback CallbackType;

            SubgraphMask(CallbackType& callback)
                : callback_(callback) {}
            bool vertex(const std::size_t v) const
                { return true; }
            bool edge(const std::size_t e) const
                {
                    const size_type v0 = callback_.poseEstimator_.detectionGraph_.vertexOfEdge(e, 0);
                    const size_type v1 = callback_.poseEstimator_.detectionGraph_.vertexOfEdge(e, 1);
                    return callback_.label(callback_.poseEstimator_.y(v0, v1)) == 1;
                }
            Callback& callback_;
        };

        Callback(PoseEstimatorType&);
        void separateAndAddLazyConstraints();

    private:
        size_type separateAndAddViolatedUniquenessConstraints();
        size_type separateAndAddViolatedSingleClusterConstraints();
        size_type separateAndAddViolatedCouplingConstraints();
        size_type separateAndAddViolated3CycleConstraints();
        size_type separateAndAddViolatedLinearizationConstraints();
        size_type separateAndAddViolatedMustSelectClassConstraints();

        PoseEstimatorType& poseEstimator_;
    };

    // ilp
    void setObjectiveFunction();
    size_type addAllImpossiblePartClassConstraints();
    size_type addAllLinearizationConstraints();
    size_type addAllCouplingConstraints();
    size_type addAllUniquenessConstraints();
    size_type addPartialSolutionConstraints(const SolutionType &partial);


    // variable indexing
    size_type x(const size_type, const size_type) const;
    size_type y(const size_type, const size_type) const;

    size_type numberOfVariablesX_;
    size_type numberOfVariablesY_;
    size_type numberOfVariablesZ_;
    ProblemType& problem_;
    double epsilonProbability_;
    CompleteGraphType detectionGraph_;
    TimerType timer_;
    FeasibleSolutionCallback feasibleSolutionCallback_;
    IlpSolver ilpSolver_;
    bool withSingleClusterConstraints_;
};

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::PoseEstimator(
    ProblemType& problem,
    SolutionType& solution,
    const SolutionType& partialSolution,
    FeasibleSolutionCallback feasibleSolutionCallback,
    const bool withSingleClusterConstraints,
    const int timeLimit,
    const double epsilonProbability
)
:   numberOfVariablesX_(),
    numberOfVariablesY_(),
    numberOfVariablesZ_(),
    problem_(problem),
    epsilonProbability_(epsilonProbability),
    detectionGraph_(problem.numberOfDetections()),
    timer_(),
    feasibleSolutionCallback_(feasibleSolutionCallback),
    ilpSolver_(),
    withSingleClusterConstraints_(withSingleClusterConstraints)
{
    ilpSolver_.setVerbosity(false);
    ilpSolver_.setRelativeGap(0.01);
    ilpSolver_.setNumberOfThreads(1);
    ilpSolver_.setTimeLimit(timeLimit);

    timer_.start();

    setObjectiveFunction();
    if(partialSolution.size() != 0)
        addPartialSolutionConstraints(partialSolution);
    addAllImpossiblePartClassConstraints();
    addAllLinearizationConstraints();
    addAllCouplingConstraints();
    addAllUniquenessConstraints();

    Callback callback(*this);
    ilpSolver_.setCallback(callback);

    ilpSolver_.optimize();

    timer_.stop();

    // save solution
    ComponentsType components;
    components.build(detectionGraph_, SubgraphMask(*this));
    solution.resize(problem_.numberOfDetections());
    for(size_type d = 0; d < problem.numberOfDetections(); ++d) {
       solution[d].clusterIndex_ = components.labels_[d];
       solution[d].partClass_ = std::numeric_limits<size_type>::max(); // suppressed
       for(size_type c = 0; c < problem.numberOfPartClasses(); ++c) {
            if(ilpSolver_.label(x(d, c)) == 1) {
                solution[d].partClass_ = c;
                break;
            }
       }
    }
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
void
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::setObjectiveFunction() {
    std::cout << "setting up objective function:" << std::endl;

    numberOfVariablesX_ = problem_.numberOfDetections() * problem_.numberOfPartClasses();
    numberOfVariablesY_ = detectionGraph_.numberOfEdges();

    std::vector<double> coefficients(numberOfVariablesX_ + numberOfVariablesY_);

    // set coefficients of variables x
    for(size_type d = 0; d < problem_.numberOfDetections(); ++d)
    for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
        const double p = problem_.getPartClassProbability(d, c);
        if(p > epsilonProbability_) {
            const size_type vi = x(d, c);
            coefficients[vi] = std::log( (1.0 - p) / p );
        }
    }

    // introduce variables z and set their coefficients:
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
        for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
            const JoinIndexType joinIndex(d0, d1, c0, c1);
            auto it = problem_.joinMap().find(joinIndex);
            if(it != problem_.joinMap().end()) {
                const double p = it->second.getProbability();
                if(p > epsilonProbability_) {
                    const double c = std::log( (1.0 - p) / p );
                    it->second.setVariableIndex(coefficients.size());
                    coefficients.push_back(c);
                    ++numberOfVariablesZ_;
                }
            }
        }
    }

    ilpSolver_.initModel(coefficients.size(), coefficients.data());

    std::cout << "   " << numberOfVariablesX_ << " variables x" << std::endl
        << "   " << numberOfVariablesY_ << " variables y" << std::endl
        << "   " << numberOfVariablesZ_ << " variables z not known to be zero" << std::endl;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAllImpossiblePartClassConstraints() {
    std::cout << "adding all impossible part class constraints: " << std::flush;

    size_type n = 0;
    for(size_type d = 0; d < problem_.numberOfDetections(); ++d)
    for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
        const double p = problem_.getPartClassProbability(d, c);
        if(p <= epsilonProbability_) {
            const size_type vi[] = {x(d, c)};
            const double c[] = {1.0};
            const double lowerBound = 0.0;
            const double upperBound = 0.0;
            ilpSolver_.addConstraint(vi, vi + 1, c, lowerBound, upperBound); // x_{d, c} = 0
	    n++;
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAllLinearizationConstraints() {
    std::cout << "adding all linearization constraints: " << std::flush;

    std::size_t n = 0;
    size_type vi[] = {0, 0, 0, 0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 2.0;
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        vi[0] = y(d0, d1);
        { // if y_{d0, d1} = 1
            const double c[] = {1.0, 1.0, 1.0, -1.0};
            for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0) {
                vi[1] = x(d0, c0);
                for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
                    vi[2] = x(d1, c1);
                    const JoinIndexType joinIndex(d0, d1, c0, c1);
                    auto it = problem_.joinMap().find(joinIndex);
                    if(it == problem_.joinMap().end()
                    || it->second.getProbability() <= epsilonProbability_) { // if z_{d0, d1, c0, c1} is not explicitly in the ILP
                        // z_{d0, d1, c0, c1} = 0. Thus:
                        ilpSolver_.addConstraint(vi, vi + 3, c, lowerBound, upperBound);
                        ++n;
                    }
                    else { // if z_{d0, d1, c0, c1} is explicitly in the ILP
                        vi[3] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}
                        ilpSolver_.addConstraint(vi, vi + 4, c, lowerBound, upperBound);
                        ++n;
                    }
                }
            }
        }
        { // if y_{d0, d1} = 0
            const double c[] = {-1.0, 1.0};
            for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0)
            for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
                const JoinIndexType joinIndex(d0, d1, c0, c1);
                auto it = problem_.joinMap().find(joinIndex);
                if (it != problem_.joinMap().end()) {// if z_{d0, d1, c0, c1} is explicitly in the ILP
                    vi[1] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}
                    ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, 0.0);
                    ++n;
                }
             }
         }
    }

    const double c[] = {-1.0, 1.0};
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0) {
        vi[0] = x(d0, c0);
        for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1)
        for(size_type c1 = 0; c1 < problem_.numberOfPartClasses(); ++c1) {
            {
                const JoinIndexType joinIndex(d0, d1, c0, c1);
                auto it = problem_.joinMap().find(joinIndex);
                if (it != problem_.joinMap().end()) { // if z_{d0, d1, c0, c1} is explicitly in the ILP
                    vi[1] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}
                    ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, 0.0);
                    ++n;
                }
            }{
                const JoinIndexType joinIndex(d1, d0, c1, c0); // note: different order
                auto it = problem_.joinMap().find(joinIndex);
                if (it != problem_.joinMap().end()) { // if z_{d1, d0, c1, c0} is explicitly in the ILP
                    vi[1] = it->second.getVariableIndex(); // z_{d1, d0, c1, c0}
                    ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, 0.0);
                    ++n;
                }
            }
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAllCouplingConstraints() {
    std::cout << "adding all coupling constraints: " << std::flush;

    size_type n = 0;
    std::vector<size_type> vi(1 + problem_.numberOfPartClasses());
    std::vector<double> c(1 + problem_.numberOfPartClasses(), 1.0);
    c.back() = -1.0;
    const double lowerBound = 0.0;
    const double upperBound = std::numeric_limits<double>::infinity();
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < problem_.numberOfDetections(); ++d1) {
        vi.back() = y(d0, d1);
        // d0
        {
            for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
                vi[c] = x(d0, c); // d0
            }
            ilpSolver_.addConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 0 <= (sum_{c \in C} x_{d0, c}) - y_{d0, d1}
            ++n;
        }
        // d1
        {
            for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
                vi[c] = x(d1, c); // d1
            }
            ilpSolver_.addConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 0 <= (sum_{c \in C} x_{d1, c}) - y_{d0, d1}
            ++n;
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addAllUniquenessConstraints() {
    std::cout << "adding all uniqueness constraints: " << std::flush;

    size_type n = 0;
    size_type vi[] = {0, 0};
    const double c[] = {1.0, 1.0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 1.0;
    for(size_type d = 0; d < problem_.numberOfDetections(); ++d) {
        for(size_type c0 = 0; c0 < problem_.numberOfPartClasses(); ++c0) {
            vi[0] = x(d, c0);
            for(size_type c1 = c0 + 1; c1 < problem_.numberOfPartClasses(); ++c1) {
                vi[1] = x(d, c1);
                ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, upperBound); // x_{d, c0} + x_{d, c1} <= 1
                ++n;
            }
        }
    }

    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::addPartialSolutionConstraints(const SolutionType &partial) {

    std::cout << "solution and problem size: " << partial.size() << " " << problem_.numberOfDetections() << std::endl;

    const size_type max_val = std::numeric_limits<size_type>::max();

    size_type n = 0;

    for(size_type d = 0; d < problem_.numberOfDetections(); ++d) { //problem_.numberOfDetections()
        auto det = partial[d];
        auto part_class_label = det.partClass_;
        auto cluster = det.clusterIndex_;
        bool suppressed = (part_class_label == max_val && cluster != max_val);
        //std::cout << "detection " << d << " suppressed " << suppressed << std::endl;
        const double bound = suppressed ? 0.0 : 1.0;
        for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
            if((part_class_label == c) || suppressed)
            {
                const size_type vi[] = {x(d, c)};
                const double coef[] = {1.0};
                ilpSolver_.addConstraint(vi, vi + 1, coef, bound, bound); // x_{d, c} = 1
                //std::cout << "setting constraint " << d << " c " << c << " " << bound << std::endl;
                n++;
            }
        }
    }

    std::cout << "adding partial solution constraints: " << std::flush;

    std::cout << n << " ";
    n = 0;

    for(size_type d0 = 0; d0 < problem_.numberOfDetections()-1; ++d0) {
        auto cluster0 = partial[d0].clusterIndex_;
        auto part_class_label0 = partial[d0].partClass_;
        if (cluster0 == max_val || part_class_label0 == max_val) {
            //std::cout << "detection " << d0 << " is unlabeled" << std::endl;
            continue;
        }
        for(size_type d1 = d0+1; d1 < problem_.numberOfDetections(); ++d1) {
            auto cluster1 = partial[d1].clusterIndex_;
            auto part_class_label1 = partial[d1].partClass_;
            if (cluster1 == max_val || part_class_label1 == max_val) {
                continue;
            }

            const size_type vi[] = {y(d0, d1)};
            const double c[] = {1.0};
            const double bound = (cluster0 == cluster1) ? 1.0 : 0.0;
            //std::cout << "yconstr " << d0 << " " << d1 << " " << bound << std::endl;
            ilpSolver_.addConstraint(vi, vi + 1, c, bound, bound); // y_{d0, d1} = 1 or 0
            n++;
        }
    }
/*
    std::cout << n << " ";
    n = 0;

    for(size_type d0 = 0; d0 < problem_.numberOfDetections()-1; ++d0) {
        auto cluster0 = partial[d0].clusterIndex_;
        auto c0 = partial[d0].partClass_;
        if (cluster0 == max_val || c0 == max_val) {
            //std::cout << "detection " << d0 << " is unlabeled" << std::endl;
            continue;
        }
        for(size_type d1 = d0+1; d1 < problem_.numberOfDetections(); ++d1) {
            auto cluster1 = partial[d1].clusterIndex_;
            auto c1 = partial[d1].partClass_;
            // only add this constraint if detection d1 is yet unlebeled
            if (cluster1 != max_val || c1 != max_val) {
                continue;
            }

            const size_type vi[] = {x(d1, c0), y(d0, d1)};
            const double c[] = {1.0, 1.0};
            const double lowerBound = -std::numeric_limits<double>::infinity();
            const double upperBound = 1.0;
            //std::cout << "yconstr " << d0 << " " << d1 << " " << bound << std::endl;
            ilpSolver_.addConstraint(vi, vi + 2, c, lowerBound, upperBound); // y_{d0, d1} = 1 or 0
            n++;
        }
    }
*/
    std::cout << n << std::endl;

    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::x(
    const size_type d,
    const size_type c
) const {
    assert(d < problem_.numberOfDetections());
    assert(c < problem_.numberOfPartClasses());
    return c + d * problem_.numberOfPartClasses();
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::y(
    const size_type d0,
    const size_type d1
) const {
    assert(d0 < problem_.numberOfDetections());
    assert(d1 < problem_.numberOfDetections());
    return numberOfVariablesX_ + detectionGraph_.findEdge(d0, d1).second;
}

// Callback

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
inline
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::Callback(
    PoseEstimatorType& poseEstimator
)
:   IlpSolverCallback(poseEstimator.ilpSolver_),
    poseEstimator_(poseEstimator)
{}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
inline void
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddLazyConstraints() {
    poseEstimator_.timer_.stop();
    const double elapsedSeconds = poseEstimator_.timer_.elapsedSeconds();
    poseEstimator_.timer_.start();

    std::cout
        << elapsedSeconds
        << '\t' << IlpSolverCallback::objectiveBound()
        << '\t' << IlpSolverCallback::objectiveBest()
        << std::flush;

    const size_type nLinearization = separateAndAddViolatedLinearizationConstraints();
    std::cout << '\t' << nLinearization << std::flush;

    const size_type nUniqueness = separateAndAddViolatedUniquenessConstraints();
    std::cout << '\t' << nUniqueness << std::flush;

    const size_type nCoupling = separateAndAddViolatedCouplingConstraints();
    std::cout << '\t' << nCoupling << std::flush;

    const size_type nCycle = separateAndAddViolated3CycleConstraints();
    std::cout << '\t' << nCycle << std::flush;

    //const size_type nMustSelectClass = separateAndAddViolatedMustSelectClassConstraints();
    //std::cout << '\t' << nMustSelectClass << std::flush;

    size_type nSingleCluster = 0;
    if(poseEstimator_.withSingleClusterConstraints_) {
        nSingleCluster = separateAndAddViolatedSingleClusterConstraints();
    }
    std::cout << '\t' << nSingleCluster << std::flush;

    std::cout << std::endl;

    if(nLinearization + nUniqueness + nCoupling + nCycle + nSingleCluster == 0) {
        // save intermediate feasible solution
        ComponentsType components;
        components.build(poseEstimator_.detectionGraph_, SubgraphMask(*this));
        typename PoseEstimatorType::SolutionType solution(poseEstimator_.problem_.numberOfDetections());
        for(size_type d = 0; d < poseEstimator_.problem_.numberOfDetections(); ++d) {
           solution[d].clusterIndex_ = components.labels_[d];
           solution[d].partClass_ = std::numeric_limits<size_type>::max(); // suppressed
           for(size_type c = 0; c < poseEstimator_.problem_.numberOfPartClasses(); ++c) {
                if(this->label(poseEstimator_.x(d, c)) == 1) {
                    solution[d].partClass_ = c;
                    break;
                }
           }
        }
        poseEstimator_.feasibleSolutionCallback_(solution);
    }
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolatedUniquenessConstraints() {
    size_type n = 0;
    size_type vi[] = {0, 0};
    const double c[] = {1.0, 1.0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 1.0;
    for(size_type d = 0; d < poseEstimator_.problem_.numberOfDetections(); ++d) {
        for(size_type c0 = 0; c0 < poseEstimator_.problem_.numberOfPartClasses(); ++c0) {
            vi[0] = poseEstimator_.x(d, c0);
            if(this->label(vi[0]) == 1) {
                for(size_type c1 = c0 + 1; c1 < poseEstimator_.problem_.numberOfPartClasses(); ++c1) {
                    vi[1] = poseEstimator_.x(d, c1);
                    if(this->label(vi[1]) == 1) {
                        this->addLazyConstraint(vi, vi + 2, c, lowerBound, upperBound); // x_{d, c0} + x_{d, c1} <= 1
                        ++n;
                    }
                }
            }
        }
    }
    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolatedMustSelectClassConstraints() {
    size_type n = 0;
    std::vector<size_type> vi(problem_.numberOfPartClasses());
    std::vector<double> c(problem_.numberOfPartClasses(), 1.0);
    const double lowerBound = 1;
    const double upperBound = std::numeric_limits<double>::infinity();
    for(size_type d0 = 0; d0 < problem_.numberOfDetections(); ++d0) {
        double sum = 0;
        for(size_type c = 0; c < problem_.numberOfPartClasses(); ++c) {
            vi[c] = x(d0, c); // d0
            sum += ilpSolver_.label(vi[c]);
        }
        if(sum == 0) {
            ilpSolver_.addConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 1 <= (sum_{c \in C} x_{d0, c})
            ++n;
        }
    }
    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolatedSingleClusterConstraints() {
    size_type n = 0;
    size_type vi[] = {0, 0, 0};
    const double c[] = {1.0, 1.0, -1.0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 1.0;
    for(size_type d0 = 0; d0 < poseEstimator_.problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < poseEstimator_.problem_.numberOfDetections(); ++d1) {
        vi[2] = poseEstimator_.y(d0, d1);
        if(this->label(vi[2]) == 0) { // y_{d0, d1} = 0
            for(size_type c0 = 0; c0 < poseEstimator_.problem_.numberOfPartClasses(); ++c0) {
                vi[0] = poseEstimator_.x(d0, c0);
                if(this->label(vi[0]) == 1) { // x_{d0, c0} = 1
                    for(size_type c1 = 0; c1 < poseEstimator_.problem_.numberOfPartClasses(); ++c1) {
                        vi[1] = poseEstimator_.x(d1, c1);
                        if(this->label(vi[1]) == 1) { // x_{d1, c0} = 1
                            this->addLazyConstraint(vi, vi + 3, c, lowerBound, upperBound); // x_{d0, c0} + x_{d1, c1} - y <= 1
                            ++n;
                        }
                    }
                }
            }
        }
    }
    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolatedCouplingConstraints() {
    size_type n = 0;
    std::vector<size_type> vi(1 + poseEstimator_.problem_.numberOfPartClasses());
    std::vector<double> c(1 + poseEstimator_.problem_.numberOfPartClasses(), 1.0);
    c.back() = -1.0;
    const double lowerBound = 0.0;
    const double upperBound = std::numeric_limits<double>::infinity();
    for(size_type d0 = 0; d0 < poseEstimator_.problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < poseEstimator_.problem_.numberOfDetections(); ++d1) {
        vi.back() = poseEstimator_.y(d0, d1);
        if(this->label(vi.back()) == 1) {
            // d0
            {
                double sum = 0;
                for(size_type c = 0; c < poseEstimator_.problem_.numberOfPartClasses(); ++c) {
                    vi[c] = poseEstimator_.x(d0, c); // d0
                    sum += this->label(vi[c]);
                }
                if(sum == 0) {
                    this->addLazyConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 0 <= (sum_{c \in C} x_{d0, c}) - y_{d0, d1}
                    ++n;
                }
            }
            // d1
            {
                double sum = 0;
                for(size_type c = 0; c < poseEstimator_.problem_.numberOfPartClasses(); ++c) {
                    vi[c] = poseEstimator_.x(d1, c); // d1
                    sum += this->label(vi[c]);
                }
                if(sum == 0) {
                    this->addLazyConstraint(vi.begin(), vi.end(), c.begin(), lowerBound, upperBound); // 0 <= (sum_{c \in C} x_{d1, c}) - y_{d0, d1}
                    ++n;
                }
            }
        }
    }

    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolated3CycleConstraints() {
    std::size_t n = 0;
    /* working Theta(n^3) implementation
    std::size_t vi[] = {0, 0, 0};
    for(std::size_t j = 0; j < poseEstimator_.problem_.numberOfDetections(); ++j) {
        for(std::size_t k = j + 1; k < poseEstimator_.problem_.numberOfDetections(); ++k) {
            vi[0] = poseEstimator_.y(j, k);
            for(std::size_t l = k + 1; l < poseEstimator_.problem_.numberOfDetections(); ++l) {
                vi[1] = poseEstimator_.y(k, l);
                vi[2] = poseEstimator_.y(j, l);
                const double lowerBound = -1.0;
                const double upperBound = 1.0;//std::numeric_limits<double>::infinity();
                if(this->label(vi[0]) == 1) {
                    if(this->label(vi[1]) == 1) {
                        if(this->label(vi[2]) == 0) {
                            const double coefficients[] = {1.0, 1.0, -1.0};
                            this->addLazyConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                            ++n;
                        }
                    }
                    else {
                        if(this->label(vi[2]) == 1) {
                            const double coefficients[] = {1.0, -1.0, 1.0};
                            this->addLazyConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                            ++n;
                        }
                    }
                }
                else {
                    if(this->label(vi[1]) == 1 && this->label(vi[2]) == 1) {
                        const double coefficients[] = {-1.0, 1.0, 1.0};
                        this->addLazyConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                        ++n;
                    }
                }
            }
        }
    }
    */

    // Theta(n^2) implementation
    ComponentsType components;
    components.build(poseEstimator_.detectionGraph_, SubgraphMask(*this));
    std::size_t vi[] = {0, 0, 0};
    for(std::size_t edge = 0; edge < poseEstimator_.detectionGraph_.numberOfEdges(); ++edge) {
        const std::size_t v0 = poseEstimator_.detectionGraph_.vertexOfEdge(edge, 0);
        const std::size_t v1 = poseEstimator_.detectionGraph_.vertexOfEdge(edge, 1);
        vi[0] = poseEstimator_.y(v0, v1);
        if(this->label(vi[0]) == 0) { // cut
            for(size_t v2 = 0; v2 < poseEstimator_.detectionGraph_.numberOfVertices(); ++v2) {
                vi[1] = poseEstimator_.y(v2, v0);
                vi[2] = poseEstimator_.y(v2, v1);
                if(this->label(vi[1]) == 1 && this->label(vi[2]) == 1) { // join, join
                    const double coefficients[] = {-1.0, 1.0, 1.0};
                    const double lowerBound = -1.0;
                    const double upperBound = 1.0;//std::numeric_limits<double>::infinity();
                    this->addLazyConstraint(vi, vi + 3, coefficients, lowerBound, upperBound);
                    ++n;
                }
            }
        }
    }
    return n;
}

template<class ILP_SOLVER, class FEASBILE_SOLUTION_CALLBACK>
typename PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::size_type
PoseEstimator<ILP_SOLVER, FEASBILE_SOLUTION_CALLBACK>::Callback::separateAndAddViolatedLinearizationConstraints() {
    std::size_t n = 0;
    size_type vi[] = {0, 0, 0, 0};
    const double lowerBound = -std::numeric_limits<double>::infinity();
    const double upperBound = 2.0;
    for(size_type d0 = 0; d0 < poseEstimator_.problem_.numberOfDetections(); ++d0)
    for(size_type d1 = d0 + 1; d1 < poseEstimator_.problem_.numberOfDetections(); ++d1) {
        vi[0] = poseEstimator_.y(d0, d1);
        if(this->label(vi[0]) == 1) { // if y_{d0, d1} = 1
            const double c[] = {1.0, 1.0, 1.0, -1.0};
            for(size_type c0 = 0; c0 < poseEstimator_.problem_.numberOfPartClasses(); ++c0) {
                vi[1] = poseEstimator_.x(d0, c0);
                if(this->label(vi[1]) == 1) { // if x_{d0, c0} = 1
                    for(size_type c1 = 0; c1 < poseEstimator_.problem_.numberOfPartClasses(); ++c1) {
                        vi[2] = poseEstimator_.x(d1, c1);
                        if(this->label(vi[2]) == 1) { // if x_{d1, c1} = 1
                            const JoinIndexType joinIndex(d0, d1, c0, c1);
                            auto it = poseEstimator_.problem_.joinMap().find(joinIndex);
                            if(it == poseEstimator_.problem_.joinMap().end()
                            || it->second.getProbability() <= poseEstimator_.epsilonProbability_) { // if z_{d0, d1, c0, c1} is not explicitly in the ILP
                                // z_{d0, d1, c0, c1} = 0. Thus:
                                this->addLazyConstraint(vi, vi + 3, c, lowerBound, upperBound);
                                ++n;
                            }
                            else { // if z_{d0, d1, c0, c1} is explicitly in the ILP
                                vi[3] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}
                                if(this->label(vi[3]) == 0) { // if z_{d0, d1, c0, c1} = 0
                                    this->addLazyConstraint(vi, vi + 4, c, lowerBound, upperBound);
                                    ++n;
                                }

                            }
                        }
                    }
                }
            }
        }
        else { // if y_{d0, d1} = 0
            const double c[] = {-1.0, 1.0};
            for(size_type c0 = 0; c0 < poseEstimator_.problem_.numberOfPartClasses(); ++c0)
            for(size_type c1 = 0; c1 < poseEstimator_.problem_.numberOfPartClasses(); ++c1) {
                const JoinIndexType joinIndex(d0, d1, c0, c1);
                auto it = poseEstimator_.problem_.joinMap().find(joinIndex);
                if (it != poseEstimator_.problem_.joinMap().end()) {// if z_{d0, d1, c0, c1} is explicitly in the ILP
                    vi[1] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}
                    if(this->label(vi[1]) == 1) { // if z_{d0, d1, c0, c1} = 1
                        this->addLazyConstraint(vi, vi + 2, c, lowerBound, 0.0);
                        ++n;
                    }
                }
             }
         }
    }

    const double c[] = {-1.0, 1.0};
    for(size_type d0 = 0; d0 < poseEstimator_.problem_.numberOfDetections(); ++d0)
    for(size_type c0 = 0; c0 < poseEstimator_.problem_.numberOfPartClasses(); ++c0) {
        vi[0] = poseEstimator_.x(d0, c0);
        if(this->label(vi[0]) == 0) { // if x_{d0, c0} = 0
            for(size_type d1 = d0 + 1; d1 < poseEstimator_.problem_.numberOfDetections(); ++d1)
            for(size_type c1 = 0; c1 < poseEstimator_.problem_.numberOfPartClasses(); ++c1) {
                {
                    const JoinIndexType joinIndex(d0, d1, c0, c1);
                    auto it = poseEstimator_.problem_.joinMap().find(joinIndex);
                    if (it != poseEstimator_.problem_.joinMap().end()) { // if z_{d0, d1, c0, c1} is explicitly in the ILP
                        vi[1] = it->second.getVariableIndex(); // z_{d0, d1, c0, c1}
                        if(this->label(vi[1]) == 1) { // if z_{d0, d1, c0, c1} = 1
                            this->addLazyConstraint(vi, vi + 2, c, lowerBound, 0.0);
                            ++n;
                        }
                    }
                }{
                    const JoinIndexType joinIndex(d1, d0, c1, c0); // note: different order
                    auto it = poseEstimator_.problem_.joinMap().find(joinIndex);
                    if (it != poseEstimator_.problem_.joinMap().end()) { // if z_{d1, d0, c1, c0} is explicitly in the ILP
                        vi[1] = it->second.getVariableIndex(); // z_{d1, d0, c1, c0}
                        if(this->label(vi[1]) == 1) { // if z_{d1, d0, c1, c0} = 1
                            this->addLazyConstraint(vi, vi + 2, c, lowerBound, 0.0);
                            ++n;
                        }
                    }
                }
            }
        }
    }

    return n;
}

} // namespace pose

#endif
