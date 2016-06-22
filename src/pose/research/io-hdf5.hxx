#include "andres/marray-hdf5.hxx"

#include "pose/problem.hxx"
#include "pose/solution.hxx"

typedef double value_type;
typedef std::size_t size_type;
typedef pose::Problem<value_type, size_type> Problem;
typedef pose::Solution<size_type> Solution;

void
loadProblem(
    const std::string problemFileName,
    Problem& problem
) {
    std::cout << "loading problem: " << std::flush;

    andres::Marray<value_type> partClassProbabilityTable;
    andres::Marray<value_type> joinProbabilityTable;
    {
        hid_t problemFile = andres::hdf5::openFile(problemFileName);
        andres::hdf5::load(problemFile, "part-class-probabilities", partClassProbabilityTable);
        andres::hdf5::load(problemFile, "join-probabilities", joinProbabilityTable);
        andres::hdf5::closeFile(problemFile);
    }

    if(partClassProbabilityTable.dimension() != 2) {
        throw std::runtime_error("partClassProbabilityTable.dimension() != 2");
    }

    const size_type numberOfDetections = partClassProbabilityTable.shape(0);
    const size_type numberOfPartClasses = partClassProbabilityTable.shape(1);

    if(joinProbabilityTable.dimension() != 2) {
        throw std::runtime_error("joinProbabilityTable.dimension() != 2");
    }

    if(joinProbabilityTable.shape(1) != 5) {
        throw std::runtime_error("joinProbabilityTable.shape(1) != 5\n \expecting a 2-dimensional matrix with 5 columns:\ndetection0, detection1, partClass0, partClass1, join probability");
    }

    problem.assign(numberOfPartClasses, numberOfDetections);

    for(size_type d0 = 0; d0 < numberOfDetections; ++d0)
    for(size_type c0 = 0; c0 < numberOfPartClasses; ++c0) {
        problem.setPartClassProbability(d0, c0, partClassProbabilityTable(d0, c0));
    }

    for(size_type j = 0; j < joinProbabilityTable.shape(0); ++j) {
        problem.setJoinProbability(
            joinProbabilityTable(j, 0),
            joinProbabilityTable(j, 1),
            joinProbabilityTable(j, 2),
            joinProbabilityTable(j, 3),
            joinProbabilityTable(j, 4)
        );
    }

    std::cout << numberOfDetections << " detections and "
        << numberOfPartClasses << " part classes" << std::endl;
}

void
saveSolution(
    const std::string solutionFileName,
    const Solution& solution
) {
    size_type shape[] = {solution.size(), 2};
    andres::Marray<size_type> m(shape, shape + 2);
    for(size_type d = 0; d < solution.size(); ++d) {
        m(d, 0) = solution[d].partClass_;
        m(d, 1) = solution[d].clusterIndex_;
    }

    /*
    std::cout << "solution:" << std::endl
        << std::cout << m.transposedView().asString() << std::endl;
    */

    hid_t file = andres::hdf5::createFile(solutionFileName);
    andres::hdf5::save(file, "detection-parts-and-clusters", m);
    andres::hdf5::closeFile(file);
}

void
loadSolution(
    const std::string solutionFileName,
    Solution& solution
) {
    typedef pose::Detection<size_type> Detection;
    andres::Marray<size_type> m;
    {
        hid_t file = andres::hdf5::openFile(solutionFileName);
        andres::hdf5::load(file, "detection-parts-and-clusters", m);
        andres::hdf5::closeFile(file);
    }

    const size_type numberOfDetections = m.shape(0);
    solution.resize(numberOfDetections);
    for(size_type d = 0; d < numberOfDetections; ++d) {
        solution[d].partClass_ = m(d, 0);
        solution[d].clusterIndex_ = m(d, 1);
        //std::cout << solution[d].partClass_ << " " << solution[d].clusterIndex_ << std::endl;
    }
}
