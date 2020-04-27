#ifndef SLAB_HH
#define SLAB_HH

#include "./definitions.hpp"
#include "casmutils/xtal/site.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>
#include <vector>

namespace casmutils
{
namespace xtal
{
        /// Invert the directions of the lattice vectors if necessary, such that
        /// the column vector matrix has a positive determinant
        Lattice make_right_handed(const Lattice& left_handed_lattice);

/// Create a superlattice using the provided integer transformation matrix
xtal::Lattice make_superlattice(const xtal::Lattice& tiling_unit, const Eigen::Matrix3i col_transf_mat);

/// Given a lattice and a vector of integer Miller indices, return the smallest superlattice
/// that has the a and b vectors spanning the specified plane
xtal::Lattice make_sliced_lattice(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes);

/// Given a structure and a vector of integer Miller indices, return the smallest superstructure
/// that has the a and b vectors spanning the specified plane
xtal::Structure make_sliced_structure(const xtal::Structure& unit_structure, const Eigen::Vector3i& miller_indexes);

/// Translate the given basis by the specified cartesian value
std::vector<xtal::Site> make_translated_basis(const std::vector<xtal::Site>& basis, const Eigen::Vector3d& shift);
} // namespace xtal
} // namespace casmutils

namespace mush
{
namespace cu = casmutils;
/// Given the primitive structure, create the smallest unit possible that exposes the surface
/// of the specified miller indices along the a-b vectors
constexpr auto make_slab_unit=cu::xtal::make_sliced_structure;

/// Given a slab, where the surface plane has already been exposed to the a-b vectors, create a
/// superstructure by stacking units along the c direction
cu::xtal::Structure make_stacked_slab(const cu::xtal::Structure& slab_unit, int stacks);

/// Tranlsate the basis in the given structure such that the basis at the specified
/// index ends up at the origin
cu::xtal::Structure make_floored_structure(const cu::xtal::Structure& shiftable_struc, int floor_atom_ix);

//TODO:
/// Given a stacked slab lattice, reduce angles as much as possible in the c direction but still
/// keep the a-b vectors invariant
cu::xtal::Lattice make_reduced_angles(const cu::xtal::Lattice& primitive_unit, const cu::xtal::Lattice& slab_unit_lattice);

} // namespace mush

#endif
