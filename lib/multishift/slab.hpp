#ifndef SLAB_HH
#define SLAB_HH

#include "./definitions.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>

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
namespace cu = casmutils;
/// Given a slab, where the surface plane has already been exposed to the a-b vectors, create a
/// superstructure by stacking units along the c direction
cu::xtal::Structure make_stacked_slab(const cu::xtal::Structure& slab_unit, int stacks);

} // namespace mush

#endif
