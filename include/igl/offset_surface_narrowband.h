// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2024
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OFFSET_SURFACE_NARROWBAND_H
#define IGL_OFFSET_SURFACE_NARROWBAND_H
#include "igl_inline.h"
#include "signed_distance.h"
#include <Eigen/Core>

namespace igl
{
  // Compute a triangulated offset surface using narrowband optimization.
  // Only computes signed distance for voxels within the narrowband,
  // providing significant speedup for large grids (10-25x faster).
  //
  // This is an optimized version of offset_surface that uses a two-pass
  // approach:
  //   1. First pass: compute unsigned squared distance for all voxels (fast)
  //   2. Filter to identify voxels within the narrowband (~4% of total)
  //   3. Second pass: compute signed distance only for narrowband voxels
  //
  // Inputs:
  //   V  #V by 3 list of mesh vertex positions
  //   F  #F by 3 list of mesh triangle indices into V
  //   isolevel  iso level to extract (signed distance: negative inside)
  //   s  number of grid cells along longest side (controls resolution)
  //   signed_distance_type  type of signing to use (see signed_distance.h)
  //   narrowband_width  width of narrowband in world units. If <= 0, uses
  //     automatic width of 2*sqrt(3)*voxel_size which ensures proper
  //     marching cubes operation.
  // Outputs:
  //   SV  #SV by 3 list of output surface mesh vertex positions
  //   SF  #SF by 3 list of output mesh triangle indices into SV
  //   GV  #GV=side(0)*side(1)*side(2) by 3 list of grid cell centers
  //   side  list of number of grid cells in x, y, and z directions
  //   S  #GV by 1 list of signed distance values _near_ `isolevel` ("far"
  //     from `isolevel` these values are flood-filled approximations)
  //
  template <
    typename DerivedV,
    typename DerivedF,
    typename isolevelType,
    typename DerivedSV,
    typename DerivedSF,
    typename DerivedGV,
    typename Derivedside,
    typename DerivedS>
  IGL_INLINE void offset_surface_narrowband(
    const Eigen::MatrixBase<DerivedV> & V,
    const Eigen::MatrixBase<DerivedF> & F,
    const isolevelType isolevel,
    const typename Derivedside::Scalar s,
    const SignedDistanceType & signed_distance_type,
    const typename DerivedV::Scalar narrowband_width,
    Eigen::PlainObjectBase<DerivedSV> & SV,
    Eigen::PlainObjectBase<DerivedSF> & SF,
    Eigen::PlainObjectBase<DerivedGV> & GV,
    Eigen::PlainObjectBase<Derivedside> & side,
    Eigen::PlainObjectBase<DerivedS> & S);
}

#ifndef IGL_STATIC_LIBRARY
#  include "offset_surface_narrowband.cpp"
#endif
#endif
