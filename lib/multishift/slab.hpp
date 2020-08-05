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
} // namespace xtal
} // namespace casmutils

namespace mush
{
namespace cu = casmutils;
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
