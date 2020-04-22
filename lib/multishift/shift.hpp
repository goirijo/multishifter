#ifndef SHIFT_HH
#define SHIFT_HH

#include<vector>

namespace casmutils
{
    namespace xtal
    {
        class Structure;
    }
}

namespace mush
{
    namespace cu=casmutils;
    /// Add empty space over the a-b plane of the given structure, keeping all the atoms in place
    /* std::vector<cu::xtal::Structure> make_cleaved_structure(const cu::xtal::Structure& slab, double cleavage_value); */

    ///For each cleavage value, create a new structure that has that amount of empty space over the a-b plane
    std::vector<cu::xtal::Structure> make_cleaved_structures(const cu::xtal::Structure& slab, const std::vector<double>& cleavage_values);
}

#endif
