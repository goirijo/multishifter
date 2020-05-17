#include <cassert>
#include <unordered_map>
#include "./shifter.hpp"
#include "casmutils/xtal/lattice.hpp"

namespace mush
{
Shifter::Shifter(const Structure& slab, int a_max, int b_max)
{
    const cu::xtal::Lattice& slab_lat = slab.lattice();
    auto [shift_vectors, _shift_records] = make_uniform_in_plane_shift_vectors(slab_lat, a_max, b_max);
    std::swap(_shift_records,this->shift_records);  //Discard temporary variable _shift_records
    auto [wg_shift_vectors, wg_shift_records] = make_uniform_in_plane_wigner_seitz_shift_vectors(slab_lat, a_max, b_max);

    assert(shift_records==wg_shift_records);

    shifted_structures=make_shifted_structures(slab,shift_vectors);
    wigner_seitz_shifted_structures=make_shifted_structures(slab,wg_shift_vectors);

    equivalence_map=categorize_equivalently_shifted_structures(shifted_structures);
}
} // namespace mush
