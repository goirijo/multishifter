#include "./slab.hpp"
#include "casmutils/xtal/structure.hpp"
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/site.hpp>
#include <casmutils/xtal/structure_tools.hpp>

// Extract to casm-utilities
#include <casm/crystallography/Lattice.hh>
#include <casm/crystallography/Superlattice.hh>
#include <vector>

namespace
{

template <typename T>
bool almost_equal(T lhs, T rhs, T tol = 1e-12)
{
    return std::abs(lhs - rhs) < tol;
}
namespace cu = casmutils;
/// Ensure that both lattices have parallel vectors axb
bool ab_plane_conserved(const cu::xtal::Lattice& lhs, const cu::xtal::Lattice& rhs)
{
    Eigen::Vector3d lhs_norm = lhs[0].cross(lhs[1]).normalized();
    Eigen::Vector3d rhs_norm = rhs[0].cross(rhs[1]).normalized();

    std::cout << "DEBUGGING: lhs_norm.transpose() is " << lhs_norm.transpose() << std::endl;
    std::cout << "DEBUGGING: rhs_norm.transpose() is " << rhs_norm.transpose() << std::endl;

    double norm_dot = lhs_norm.dot(rhs_norm);
    if (almost_equal(norm_dot, 1.0) || almost_equal(norm_dot, -1.0))
    {
        std::cout << "DEBUGGING: norm_dot is " << norm_dot << std::endl;

        return true;
    }

    return false;
}
} // namespace

namespace casmutils
{
namespace xtal
{
Lattice make_superlattice(const Lattice& tiling_unit, const Eigen::Matrix3i col_transf_mat)
{
    return Lattice(CASM::xtal::make_superlattice(tiling_unit.__get(), col_transf_mat));
}

Lattice make_sliced_lattice(const Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes)
{
    auto orthoscore = [=](const Lattice& l) { return std::abs(l.a().normalized().dot(l.b().normalized())); };

    // The 0 means "get the smallest cell possible", and I wish it was the default
    Lattice shift_lattice = unit_lattice.__get().lattice_in_plane(miller_indexes, 0);

    Lattice best_lattice=shift_lattice;
    double last_score = orthoscore(shift_lattice);
    for (const Eigen::Matrix3i& mat : CASM::unimodular_matrices())
    {
        // discard any transformation on the a or b vector that isn't a linear combination of a and b
        // discard any transformation that modifies c
        if (mat(0, 2) != 0 || mat(1, 2) != 0 || mat(2, 0) != 0 || mat(2, 1) != 0 || mat(2, 2) != 1)
        {
            continue;
        }

        /* auto candidate_lat_mat = shift_lattice.column_vector_matrix() * mat.cast<double>(); */
        Lattice new_lattice = make_superlattice(shift_lattice,mat);
        double new_score = orthoscore(new_lattice);
        if (new_score < last_score)
        {
            best_lattice=new_lattice;
            last_score = new_score;
        }
    }

    assert(::ab_plane_conserved(best_lattice, shift_lattice));
    return best_lattice;
}

Structure make_sliced_structure(const Structure& unit_structure, const Eigen::Vector3i& miller_indexes)
{
    const Lattice& unit_lattice = unit_structure.lattice();
    Lattice sliced_lattice = make_sliced_lattice(unit_lattice, miller_indexes);

    // TODO: Refactor out, so that it's available in casmutils
    CASM::xtal::Superlattice slice_superlattice(unit_lattice.__get(), sliced_lattice.__get());
    // TODO: Fix this annoying casting issue with Eigen
    Eigen::Matrix3i slice_transf_mat = slice_superlattice.transformation_matrix().cast<int>();

    Structure sliced_structure = make_superstructure(unit_structure, slice_transf_mat);
    sliced_structure.within();
    return sliced_structure;
}

std::vector<Site> make_translated_basis(const std::vector<Site>& basis, const Eigen::Vector3d& shift)
{
    std::vector<Site> translated_basis;
    for (const auto& site : basis)
    {
        translated_basis.emplace_back(Coordinate(site.cart() + shift), site.label());
    }

    return translated_basis;
}
} // namespace xtal
} // namespace casmutils

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
    auto translated_basis = cu::xtal::make_translated_basis(shiftable_struc.basis_sites(), cart_shift);

    cu::xtal::Structure shifted_struc(shiftable_struc.lattice(), translated_basis);
    shifted_struc.within();

    return shifted_struc;
}
} // namespace mush
