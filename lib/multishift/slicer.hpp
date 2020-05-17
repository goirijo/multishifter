#ifndef SLICER_HH
#define SLICER_HH

#include "./slab.hpp"
#include "casmutils/xtal/structure.hpp"
#include <vector>

namespace mush
{
    /**
     * Given a primitive structure and miller indexes, generates
     * all structures relevant to the slicing step:
     * the sliced primitive, and a set of translationally
     * equivalent structures, with the basis translated
     * to expose every possible atom to the ab-plane
     */

    struct Slicer
    {
        using Structure=cu::xtal::Structure;

        Slicer(const Structure& prim, const Eigen::Vector3i& miller_indexes);
        Structure prim;
        Structure sliced_prim;
        std::vector<Structure> floored_sliced_prims;

        private:
        const Eigen::Vector3i miller_indexes;
        void generate_all_possible_floored_structures();
    };
}

#endif
