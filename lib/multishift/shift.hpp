#ifndef SHIFT_HH
#define SHIFT_HH

#include "casmutils/xtal/lattice.hpp"
#include <bits/c++config.h>
#include <casmutils/xtal/structure.hpp>
#include <utility>
#include <vector>

namespace mush
{
/// Simple data structure to keep track of what index corresponds to a particular shift.
/// When generating a uniform grid of shifts, only a density along each a anb b direction
/// needs to be specified, and then all shifts are retuned in an unrolled vector.
/// This data container keeps track of the index along the a direction, b direction, and
/// position of the unrolled vector.
struct ShiftRecord
{
    ShiftRecord(int a, int b, int index) : a(a), b(b), index(index) {}

    int index;
    int a;
    int b;

    bool operator==(const ShiftRecord& lhs) const
    {
        return (this->a==lhs.a && this->b==lhs.b && this->index==lhs.index);
    }
};

inline std::ostream& operator<<(std::ostream& os, const ShiftRecord& sr)
{
    os<<sr.a<<", "<<sr.b<<": "<<sr.index;
    return os;
}

namespace cu = casmutils;
/// Add empty space over the a-b plane of the given structure, keeping all the atoms in place
/* std::vector<cu::xtal::Structure> make_cleaved_structure(const cu::xtal::Structure& slab, double cleavage_value); */

/// For each cleavage value, create a new structure that has that amount of empty space over the a-b plane
std::vector<cu::xtal::Structure> make_cleaved_structures(const cu::xtal::Structure& slab, const std::vector<double>& cleavage_values);

/// Given a density along each of the a and b directions, create an unrolled vector of shifts
/// along the ab-plane with a uniform distribution. The density specifies the number of unique
/// shifts, so periodic images are NOT counted in the density. This means that even numbers yield
/// shifts at midway points. For example, a 2x2 grid gives 4 points: origin, halfway along a,
/// halfway along b, and dead center of the ab-plane.
/// Two values are returned: a vector of the shifts themselves, and a vector of ShiftRecord, one
/// for each shift, which relate integer representations of the shifts to the correct index
/// of the vector of shifts. Following the previous example, the shifts and their corresponding
/// records might be:
///
/// origin : {0,0,0}
/// halfway along b: {0,1,1}
/// halfway along a: {1,0,2}
/// center of plane: {1,1,3}
std::pair<std::vector<Eigen::Vector3d>, std::vector<ShiftRecord>>
make_uniform_in_plane_shift_vectors(const cu::xtal::Lattice& slab_lattice, int a_max, int b_max);

/// Identical to make_uniform_in_plane_shift_vectors, except the sifts are kept within the
/// Wigner-Seitz cell.
std::pair<std::vector<Eigen::Vector3d>, std::vector<ShiftRecord>>
make_uniform_in_plane_wigner_seitz_shift_vectors(const cu::xtal::Lattice& slab_lattice, int a_max, int b_max);

/// Given a list of shift vectors parallel to the ab-plane of the structure, make a new structure
/// such that the periodic images along the vertical direction have been shifted by that amount.
/// The resulting structre will be a series of staggered slabs
std::vector<cu::xtal::Structure> make_shifted_structures(const cu::xtal::Structure& slab, std::vector<Eigen::Vector3d>& shifts);

/// Given a list of slab structures with different shifts applied, return a list of indexes
/// that describe which structures are equivalent to each other. The vector at index i contains
/// a the indexes of the structures that are equivalent to the structure at index i.
std::vector<std::vector<std::size_t>> categorize_equivalently_shifted_structures(const std::vector<cu::xtal::Structure>& shifted_structures);
} // namespace mush

#endif
