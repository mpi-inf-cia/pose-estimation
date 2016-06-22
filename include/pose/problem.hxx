#pragma once
#ifndef POSE_PROBLEM_HXX
#define POSE_PROBLEM_HXX

#include <cassert>
#include <cstddef>

#include <andres/marray.hxx>

#include "join.hxx"

namespace pose {

template<class T = double, class S = std::size_t>
class Problem {
public:
    typedef T value_type;
    typedef S size_type;
    typedef JoinData<size_type> JoinDataType;
    typedef JoinIndex<size_type> JoinIndexType;
    typedef JoinMap<size_type> JoinMapType;

    Problem(const size_type = 0, const size_type = 0);
    void assign(const size_type, const size_type);
    void setPartClassProbability(const size_type, const size_type, const value_type);
    void setJoinProbability(const size_type, const size_type, const size_type, const size_type, const value_type);
    JoinMapType& joinMap();

    size_type numberOfPartClasses() const;
    size_type numberOfDetections() const;
    value_type getPartClassProbability(const size_type, const size_type) const;
    const JoinMapType& joinMap() const;

private:
    typedef andres::Marray<value_type> PartClassProbabilityMatrix;

    size_type numberOfPartClasses_;
    size_type numberOfDetections_;
    PartClassProbabilityMatrix partClassProbabilities_;
    JoinMapType joinMap_;
};

template<class T, class S>
Problem<T, S>::Problem(
    const size_type numberOfPartClasses,
    const size_type numberOfDetections
)
:   numberOfPartClasses_(numberOfPartClasses),
    numberOfDetections_(numberOfDetections),
    partClassProbabilities_(),
    joinMap_()
{    
    if(numberOfPartClasses_ != 0 && numberOfDetections_ != 0) {
        size_type shape[] = {numberOfDetections_, numberOfPartClasses_};
        partClassProbabilities_.resize(shape, shape + 2);
    }
}

template<class T, class S>
void
Problem<T, S>::assign(
    const size_type numberOfPartClasses,
    const size_type numberOfDetections
) {
    numberOfPartClasses_ = numberOfPartClasses;
    numberOfDetections_ = numberOfDetections;
    if(numberOfPartClasses_ != 0 && numberOfDetections_ != 0) {
        size_type shape[] = {numberOfDetections_, numberOfPartClasses_};
        partClassProbabilities_ = PartClassProbabilityMatrix(shape, shape + 2);
    }
    joinMap_.clear();
}

template<class T, class S>
inline void
Problem<T, S>::setPartClassProbability(
    const size_type detection,
    const size_type partClass,
    const value_type probability
) {
    assert(detection < numberOfDetections());
    assert(partClass < numberOfPartClasses());
    assert(probability >= 0 && probability <= 1);

    partClassProbabilities_(detection, partClass) = probability;
}

template<class T, class S>
inline void
Problem<T, S>::setJoinProbability(
    const size_type detection0,
    const size_type detection1,
    const size_type partClass0,
    const size_type partClass1,
    const value_type probability
) {
    JoinIndexType joinIndex(detection0, detection1, partClass0, partClass1);
    JoinDataType joinData(probability);
    joinMap_[joinIndex] = joinData;
}

template<class T, class S>
typename Problem<T, S>::JoinMapType&
Problem<T, S>::joinMap() {
    return joinMap_;
}

template<class T, class S>
inline typename Problem<T, S>::size_type
Problem<T, S>::numberOfPartClasses() const {
    return numberOfPartClasses_;
}

template<class T, class S>
inline typename Problem<T, S>::size_type
Problem<T, S>::numberOfDetections() const {
    return numberOfDetections_;
}

template<class T, class S>
inline typename Problem<T, S>::value_type
Problem<T, S>::getPartClassProbability(
    const size_type detection,
    const size_type partClass
) const {
    return partClassProbabilities_(detection, partClass);
}

template<class T, class S>
const typename Problem<T, S>::JoinMapType&
Problem<T, S>::joinMap() const {
    return joinMap_;
}

} // namespace pose

#endif
