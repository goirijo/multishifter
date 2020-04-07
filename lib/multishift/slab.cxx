#include "casmutils/xtal/site.hpp"
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <multishift/slab.hpp>

// Extract to casm-utilities
#include <casm/crystallography/Lattice.hh>

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
    auto lhs_norm = lhs[0].cross(lhs[1]).normalized();
    auto rhs_norm = rhs[0].cross(rhs[1]).normalized();

    double norm_dot = lhs_norm.dot(rhs_norm);
    if (almost_equal(norm_dot, 1.0) || almost_equal(norm_dot, -1.0))
    {
        return true;
    }

    return false;
}
} // namespace

namespace casmutils
{
namespace xtal
{
xtal::Lattice make_superlattice(const xtal::Lattice& tiling_unit, const Eigen::Matrix3i col_transf_mat)
{
    return xtal::Lattice(CASM::xtal::make_superlattice(tiling_unit.__get(), col_transf_mat));
}

xtal::Lattice make_sliced_lattice(const xtal::Lattice& unit_lattice, const Eigen::Vector3i& miller_indexes) {
    // The 0 means "get the smallest cell possible", and I wish it was the default
    xtal::Lattice shift_lattice = unit_lattice.__get().lattice_in_plane(miller_indexes, 0);

    Eigen::Matrix3d best_lat_mat = shift_lattice.column_vector_matrix();
    double orthoscore=std::abs(best_lat_mat.col(0).normalized().dot(best_lat_mat.col(1).normalized()));
    auto uni_mats=CASM::unimodular_matrices();
    for(const auto& mat : uni_mats)
    {
        //discard any transformation on the a or b vector that isn't a linear combination of a and b
        //discard any transformation that modifies c
        if(mat(0,2)!=0 || mat(1,2)!=0 || mat(2,0)!=0 || mat(2,1)!=0 || mat(2,2)!=1)
        {
            continue;
        }

        auto candidate_lat_mat=shift_lattice.column_vector_matrix()*mat.cast<double>();
        double new_orthoscore=std::abs(candidate_lat_mat.col(0).normalized().dot(candidate_lat_mat.col(1).normalized()));

        if(new_orthoscore<orthoscore)
        {
            best_lat_mat=candidate_lat_mat;
            orthoscore=new_orthoscore;
        }
    }

    xtal::Lattice best_lattice(best_lat_mat);
    assert(::ab_plane_conserved(best_lattice,shift_lattice));

    return best_lattice;


}

xtal::Structure make_sliced_structure(const xtal::Structure& unit_structure, const Eigen::Vector3i& miller_indexes) {
    const Lattice& unit_lattice=unit_structure.lattice();
    xtal::Lattice sliced_lattice=make_sliced_lattice(unit_lattice,miller_indexes);

    //TODO: Refactor out, so that it's available in casmutils
    CASM::xtal::Superlattice slice_superlattice(unit_lattice.__get(),sliced_lattice.__get());
    //TODO: Fix this annoying casting issue with Eigen
    Eigen::Matrix3i slice_transf_mat=slice_superlattice.transformation_matrix().cast<int>();

    Structure sliced_structure=make_superstructure(unit_structure,slice_transf_mat);
    //TODO: The basis isn't within
    return sliced_structure;
}
} // namespace xtal
} // namespace casmutils

namespace mush
{
}
