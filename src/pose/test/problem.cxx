#include <cstddef>
#include <stdexcept>

#include "pose/problem.hxx"

inline void test(const bool condition) {
    if(!condition) throw std::logic_error("test failed.");
}

int main() {
    typedef std::size_t size_type;
    typedef pose::Problem<> Problem;
    typedef Problem::JoinDataType JoinData;
    typedef Problem::JoinIndexType JoinIndex;

    const size_type numberOfPartClasses = 3;
    const size_type numberOfDetections = 5;

    Problem problem(numberOfPartClasses, numberOfDetections);

    // test initial state
    test(problem.numberOfPartClasses() == numberOfPartClasses);
    test(problem.numberOfDetections() == numberOfDetections);

    for(size_type d0 = 0; d0 < numberOfDetections; ++d0)
    for(size_type c0 = 0; c0 < numberOfPartClasses; ++c0) {
        test(problem.getPartClassProbability(d0, c0) == 0);
    }

    for(size_type d0 = 0; d0 < numberOfDetections; ++d0)
    for(size_type d1 = 0; d1 < numberOfDetections; ++d1)
    for(size_type c0 = 0; c0 < numberOfPartClasses; ++c0)
    for(size_type c1 = 0; c1 < numberOfPartClasses; ++c1) {
        const JoinIndex joinIndex(d0, d1, c0, c1);
        test(problem.joinMap().find(joinIndex) == problem.joinMap().end());
    }

    // test adding of probabilities
    problem.setPartClassProbability(2, 1, 0.5);
    test(problem.getPartClassProbability(2, 1) == 0.5);

    problem.setJoinProbability(2, 0, 2, 1, 0.5);
    problem.setJoinProbability(1, 0, 2, 1, 0.4);
    problem.setJoinProbability(2, 0, 1, 0, 0.3);
    test(problem.joinMap().find(JoinIndex(2, 0, 2, 1))->second.getProbability() == 0.5);
    test(problem.joinMap().find(JoinIndex(1, 0, 2, 1))->second.getProbability() == 0.4);
    test(problem.joinMap().find(JoinIndex(2, 0, 1, 0))->second.getProbability() == 0.3);
    return 0;
}
