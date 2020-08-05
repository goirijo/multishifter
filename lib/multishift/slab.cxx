#include "./slab.hpp"
#include "casmutils/xtal/structure.hpp"
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/xtal/frankenstein.hpp>

// Extract to casm-utilities
#include <casm/crystallography/Lattice.hh>
#include <casm/crystallography/Superlattice.hh>
#include <vector>

namespace mush
{
cu::xtal::Structure make_stacked_slab(const cu::xtal::Structure& slab_unit, int stacks)
{
    // Always stack along c-direction
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, stacks;

    return cu::xtal::make_superstructure(slab_unit, stack_mat);
}

cu::xtal::Structure make_floored_structure(const cu::xtal::Structure& shiftable_struc, int floor_atom_ix)
{
    // Index 0 means no translation
    if (floor_atom_ix < 1)
    {
        return shiftable_struc;
    }

    const Eigen::Vector3d cart_shift = -shiftable_struc.basis_sites()[floor_atom_ix].cart();
    return cu::frankenstein::translate_basis(shiftable_struc,cart_shift);
}
} // namespace mush
