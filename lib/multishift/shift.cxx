#include "./shift.hpp"
#include "casmutils/xtal/coordinate.hpp"
#include <casmutils/definitions.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <casmutils/xtal/lattice.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/mapping/structure_mapping.hpp>
#include <utility>
#include <cassert>

namespace mush
{
std::vector<cu::xtal::Structure> make_cleaved_structures(const cu::xtal::Structure& slab, const std::vector<double>& cleavage_values)
{
    // TODO: This could go astray if you get a left handed lattice
    std::vector<cu::xtal::Structure> cleaved_structures;

    Eigen::Vector3d unit_normal = slab.lattice().a().cross(slab.lattice().b()).normalized();
    for (double cleave : cleavage_values)
    {
        Eigen::Vector3d new_c_vector = slab.lattice().c() + cleave * unit_normal;
        cu::xtal::Lattice cleaved_lattice(slab.lattice().a(), slab.lattice().b(), new_c_vector);

        cleaved_structures.emplace_back(slab.set_lattice(cleaved_lattice, cu::xtal::CART));
    }
    return cleaved_structures;
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<ShiftRecord>>
make_uniform_in_plane_shift_vectors(const cu::xtal::Lattice& slab_lattice, int a_max, int b_max)
{
    std::vector<Eigen::Vector3d> shifts;
    std::vector<ShiftRecord> records;
    int ix = 0;
    for (int a = 0; a < a_max; ++a)
    {
        for (int b = 0; b < b_max; ++b)
        {
            Eigen::Vector3d frac_coord(static_cast<double>(a) / a_max, static_cast<double>(b) / b_max, 0);
            shifts.emplace_back(cu::xtal::make_fractional(frac_coord, slab_lattice));

            records.emplace_back(a, b, ix);
            ++ix;
        }
    }
    return std::make_pair(std::move(shifts), std::move(records));
}

std::pair<std::vector<Eigen::Vector3d>, std::vector<ShiftRecord>>
make_uniform_in_plane_wigner_seitz_shift_vectors(const cu::xtal::Lattice& slab_lattice, int a_max, int b_max)
{
    auto standard_shifts_and_records = make_uniform_in_plane_shift_vectors(slab_lattice, a_max, b_max);
    const auto& standard_shifts = standard_shifts_and_records.first;

    std::vector<Eigen::Vector3d> wigner_seitz_shifts;
    for (const auto& standard_shift : standard_shifts)
    {
        wigner_seitz_shifts.emplace_back(cu::xtal::bring_within_wigner_seitz(standard_shift,slab_lattice));
    }

    return std::make_pair(wigner_seitz_shifts,std::move(standard_shifts_and_records.second));
}

std::vector<cu::xtal::Structure> make_shifted_structures(const cu::xtal::Structure& slab, std::vector<Eigen::Vector3d>& shifts)
{
    //for debugging:
    Eigen::Vector3d plane_normal=slab.lattice().a().cross(slab.lattice().b());

    std::vector<cu::xtal::Structure> shifted_structures;
    for(const Eigen::Vector3d& shift : shifts)
    {
        assert(std::abs(shift.dot(plane_normal))<1e-10);
        
        cu::xtal::Lattice shifted_lat(slab.lattice().a(),slab.lattice().b(),slab.lattice().c()+shift);
        shifted_structures.emplace_back(slab.set_lattice(shifted_lat,cu::xtal::CART));
    }
    return shifted_structures;
}

cu::mapping::MappingInput make_shifted_structures_categorization_map_strategy()
{
    cu::mapping::MappingInput map_strategy;
    map_strategy.k_best_maps=0;
    map_strategy.min_cost=1e-8;
    map_strategy.use_crystal_symmetry=true;

    return map_strategy;
}

std::vector<std::vector<std::size_t>> categorize_equivalently_shifted_structures(const std::vector<cu::xtal::Structure>& shifted_structures)
{
    std::vector<std::vector<std::size_t>> index_map;
    cu::mapping::MappingInput map_strategy=make_shifted_structures_categorization_map_strategy();
    for(const auto shifted_structure : shifted_structures)
    {
        cu::mapping::StructureMapper_f map_to_shifted(shifted_structure,map_strategy);
        index_map.push_back({});
        for(int s_ix=0; s_ix<shifted_structures.size(); ++s_ix)
        {
            auto results=map_to_shifted(shifted_structures[s_ix]);
            if(results.size()>0)
            {
                index_map.back().emplace_back(s_ix);
            }
        }
    }
return index_map;
}
} // namespace mush
