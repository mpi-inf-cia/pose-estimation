#pragma once
#ifndef POSE_DETECTION_HXX
#define POSE_DETECTION_HXX

#include <cstddef>
namespace pose {

template<class S = std::size_t>
struct Detection {
    typedef S size_type;

    Detection()
        :   partClass_(),
            clusterIndex_()
        {}


    size_type partClass_;
    size_type clusterIndex_;
};


} // namespace pose

#endif
