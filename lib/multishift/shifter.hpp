#ifndef SHIFTER_HH
#define SHIFTER_HH

#include "./slab.hpp"
#include "casmutils/xtal/structure.hpp"
#include "multishift/shift.hpp"
#include <bits/c++config.h>
#include <vector>

namespace mush
{
    /**
     * Given a slab and the density along the a and b
     * vectors, generate all shift structures, and a record
     * of which ones are symmetrically equivalent.
     */

    struct Shifter
    {
        using Structure=cu::xtal::Structure;

        Shifter(const Structure& slab, int a_max, int b_max);
        std::vector<Structure> shifted_structures;
        std::vector<Structure> wigner_seitz_shifted_structures;
        /// Minimal information to determine what shift has been applied to which structure 
        std::vector<ShiftRecord> shift_records; 
        /// For each index i, shows the indexes of the structures that are equivalent to i.
        std::vector<std::vector<std::size_t>> equivalence_map;
    };
}

#endif
