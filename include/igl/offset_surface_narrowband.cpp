// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2024
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include "offset_surface_narrowband.h"
#include "voxel_grid.h"
#include "signed_distance.h"
#include "flood_fill.h"
#include "marching_cubes.h"
#include "parallel_for.h"
#include "AABB.h"
#include <vector>
#include <cmath>
#include <limits>

template <
  typename DerivedV,
  typename DerivedF,
  typename isolevelType,
  typename DerivedSV,
  typename DerivedSF,
  typename DerivedGV,
  typename Derivedside,
  typename DerivedS>
IGL_INLINE void igl::offset_surface_narrowband(
  const Eigen::MatrixBase<DerivedV> & V,
  const Eigen::MatrixBase<DerivedF> & F,
  const isolevelType isolevel,
  const typename Derivedside::Scalar s,
  const SignedDistanceType & signed_distance_type,
  const typename DerivedV::Scalar narrowband_width_in,
  Eigen::PlainObjectBase<DerivedSV> & SV,
  Eigen::PlainObjectBase<DerivedSF> & SF,
  Eigen::PlainObjectBase<DerivedGV> & GV,
  Eigen::PlainObjectBase<Derivedside> & side,
  Eigen::PlainObjectBase<DerivedS> & S)
{
  typedef typename DerivedV::Scalar Scalar;
  typedef typename DerivedF::Scalar Index;

  // Step 1: Create voxel grid (same as original offset_surface)
  igl::voxel_grid(V, isolevel, s, 1, GV, side);

  const int n_points = GV.rows();
  if (n_points == 0)
  {
    SV.resize(0, 3);
    SF.resize(0, 3);
    S.resize(0, 1);
    return;
  }

  // Compute voxel size (handle edge case of single voxel)
  const Scalar h = (side(0) > 1) ?
                   (GV.col(0).maxCoeff() - GV.col(0).minCoeff()) / ((Scalar)(side(0) - 1)) :
                   Scalar(1);

  // Marching cubes needs accurate values within sqrt(3)*h of isolevel
  // This is the voxel diagonal - the maximum distance a vertex can be from cell center
  const Scalar mc_margin = std::sqrt(3.0) * h;

  // Compute narrowband width (auto if <= 0)
  // Default: 2 * mc_margin ensures we capture all voxels that could affect marching cubes
  // plus some margin for the signed_distance bounds filtering
  const Scalar narrowband_width = (narrowband_width_in <= 0) ?
                                   (2.0 * mc_margin) :
                                   narrowband_width_in;

  // For filtering, we use UNSIGNED distance around |isolevel|
  // This is because sqrD gives unsigned distance, which is always >= 0
  const Scalar abs_isolevel = std::abs(static_cast<Scalar>(isolevel));
  const Scalar filter_lower = (abs_isolevel > narrowband_width) ?
                               (abs_isolevel - narrowband_width) : Scalar(0);
  const Scalar filter_upper = abs_isolevel + narrowband_width;

  // Bounds for signed_distance call (used for optimization/early termination)
  // These are signed distance bounds, matching the original offset_surface behavior
  const Scalar sd_lower_bound = isolevel - mc_margin;
  const Scalar sd_upper_bound = isolevel + mc_margin;

  // Step 2: Build AABB tree once for fast distance queries
  AABB<DerivedV, 3> tree;
  tree.init(V, F);

  // Step 3: FAST PASS - Compute unsigned squared distance for ALL points
  // This is much faster than signed_distance because it only uses the AABB tree
  // without any sign computation (winding number, pseudonormal, etc.)
  Eigen::Matrix<Scalar, Eigen::Dynamic, 1> sqrD(n_points);

  igl::parallel_for(n_points, [&](const int i) {
    Eigen::Matrix<Scalar, 1, 3> c;
    int idx;
    sqrD(i) = tree.squared_distance(V, F, GV.row(i), idx, c);
  }, 10000);

  // Step 4: Identify narrowband voxels
  // Filter based on UNSIGNED distance being near |isolevel|
  // This works correctly for both positive and negative isolevel values
  std::vector<int> band_indices;
  band_indices.reserve(n_points / 10);  // Estimate ~10% in band (usually ~4%)

  for (int i = 0; i < n_points; ++i)
  {
    const Scalar dist = std::sqrt(sqrD(i));
    // Check if unsigned distance is within the narrowband around |isolevel|
    if (dist >= filter_lower && dist <= filter_upper)
    {
      band_indices.push_back(i);
    }
  }

  const int n_band = static_cast<int>(band_indices.size());

  // Step 5: Initialize S with NaN (will be flood-filled later)
  S.resize(n_points, 1);
  S.setConstant(std::numeric_limits<Scalar>::quiet_NaN());

  // Step 6: Compute SIGNED distance only for narrowband voxels
  if (signed_distance_type == SIGNED_DISTANCE_TYPE_UNSIGNED)
  {
    // For unsigned type, just use sqrt of squared distance (already computed)
    igl::parallel_for(n_band, [&](const int j) {
      const int i = band_indices[j];
      S(i) = std::sqrt(sqrD(i));
    }, 10000);
  }
  else
  {
    // Extract narrowband points into a contiguous matrix
    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> GV_band(n_band, 3);
    for (int j = 0; j < n_band; ++j)
    {
      GV_band.row(j) = GV.row(band_indices[j]);
    }

    // Compute signed distance for narrowband only
    // This is the expensive operation that we're avoiding for ~96% of voxels
    Eigen::Matrix<Scalar, Eigen::Dynamic, 1> S_band;
    Eigen::Matrix<Index, Eigen::Dynamic, 1> I_band;
    Eigen::Matrix<Scalar, Eigen::Dynamic, 3> C_band, N_band;

    igl::signed_distance(
      GV_band, V, F,
      signed_distance_type,
      sd_lower_bound, sd_upper_bound,
      S_band, I_band, C_band, N_band);

    // Scatter results back to full grid
    for (int j = 0; j < n_band; ++j)
    {
      S(band_indices[j]) = S_band(j);
    }
  }

  // Step 7: Flood fill NaN regions
  // This propagates valid signed distance values to fill the NaN regions,
  // ensuring proper inside/outside determination for marching cubes
  igl::flood_fill(side, S);

  // Step 8: Extract isosurface via marching cubes
  DerivedS SS = S.array() - isolevel;
  igl::marching_cubes(SS, GV, side(0), side(1), side(2), 0, SV, SF);
}

#ifdef IGL_STATIC_LIBRARY
// Explicit template instantiation
template void igl::offset_surface_narrowband<
  Eigen::Matrix<double, -1, -1, 0, -1, -1>,
  Eigen::Matrix<int, -1, -1, 0, -1, -1>,
  double,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>,
  Eigen::Matrix<int, -1, -1, 0, -1, -1>,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>,
  Eigen::Matrix<int, 1, 3, 1, 1, 3>,
  Eigen::Matrix<double, -1, 1, 0, -1, 1>
>(
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&,
  double,
  Eigen::Matrix<int, 1, 3, 1, 1, 3>::Scalar,
  igl::SignedDistanceType const&,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<int, 1, 3, 1, 1, 3>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>&);

// Additional common instantiation
template void igl::offset_surface_narrowband<
  Eigen::Matrix<double, -1, -1, 0, -1, -1>,
  Eigen::Matrix<int, -1, -1, 0, -1, -1>,
  int,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>,
  Eigen::Matrix<int, -1, -1, 0, -1, -1>,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>,
  Eigen::Matrix<int, 1, 3, 1, 1, 3>,
  Eigen::Matrix<double, -1, 1, 0, -1, 1>
>(
  Eigen::MatrixBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>> const&,
  Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>> const&,
  int,
  Eigen::Matrix<int, 1, 3, 1, 1, 3>::Scalar,
  igl::SignedDistanceType const&,
  Eigen::Matrix<double, -1, -1, 0, -1, -1>::Scalar,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<int, -1, -1, 0, -1, -1>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, -1, 0, -1, -1>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<int, 1, 3, 1, 1, 3>>&,
  Eigen::PlainObjectBase<Eigen::Matrix<double, -1, 1, 0, -1, 1>>&);
#endif
