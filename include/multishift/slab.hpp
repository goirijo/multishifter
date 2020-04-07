#ifndef SLAB_HH
#define SLAB_HH

#include "casmutils/xtal/structure.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <multishift/definitions.hpp>

namespace casmutils
{
namespace xtal
{
/// Create a superlattice using the provided integer transformation matrix
xtal::Lattice make_superlattice(const xtal::Lattice& tiling_unit, const Eigen::Matrix3i col_transf_mat);

/// Given a lattice and a vector of integer Miller indices, return the smallest superlattice
/// that has the a and b vectors spanning the specified plane
xtal::Lattice make_sliced_lattice(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes);

/// Given a structure and a vector of integer Miller indices, return the smallest superstructure
/// that has the a and b vectors spanning the specified plane
xtal::Structure make_sliced_structure(const xtal::Structure& unit_structure, const Eigen::Vector3i& miller_indexes);
} // namespace xtal
} // namespace casmutils

namespace mush
{
} // namespace mush

#endif
