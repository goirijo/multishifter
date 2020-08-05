#include "./slicer.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "multishift/slab.hpp"
#include <algorithm>
#include <vector>

namespace mush
{
Slicer::Slicer(const cu::xtal::Structure& prim, const Eigen::Vector3i& miller_indexes)
    : prim(prim), sliced_prim(cu::xtal::slice_along_plane(prim, miller_indexes)), miller_indexes(miller_indexes)
{
    this->generate_all_possible_floored_structures();
}

void Slicer::generate_all_possible_floored_structures()
{
    for (int i = 0; i < prim.basis_sites().size(); ++i)
    {
        auto floored_prim=make_floored_structure(prim,i);
        floored_sliced_prims.emplace_back(slice_along_plane(floored_prim,miller_indexes));
    }
    return;
}
} // namespace mush
