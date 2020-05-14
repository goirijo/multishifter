#include "./slicer.hpp"
#include "casmutils/xtal/site.hpp"
#include "multishift/slab.hpp"
#include <algorithm>

namespace mush
{
Slicer::Slicer(const cu::xtal::Structure& prim, const Eigen::Vector3i& miller_indexes)
    : prim(prim), sliced_prim(cu::xtal::make_sliced_structure(prim, miller_indexes))
{
    this->generate_all_possible_floored_structures();
}

void Slicer::generate_all_possible_floored_structures()
{
    for (int i = 0; i < sliced_prim.basis_sites().size(); ++i)
    {
        cu::xtal::SiteEquals_f match_to_ith(sliced_prim.basis_sites()[i], 1e-8);
        const auto& prim_basis=prim.basis_sites();
        if (std::find_if(prim_basis.begin(), prim_basis.end(), match_to_ith) != prim_basis.end())
        {
            floored_sliced_prims.emplace_back(mush::make_floored_structure(sliced_prim,i));
        }
    }
    return;
}
} // namespace mush
