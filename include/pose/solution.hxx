#pragma once
#ifndef POSE_SOLUTION_HXX
#define POSE_SOLUTION_HXX

#include <cstddef>
#include <vector>

#include "detection.hxx"

namespace pose {

template<class S = std::size_t>
using Solution = std::vector<Detection<S> >;

} // namespace pose

#endif
