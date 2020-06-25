#include "./slice.hpp"
#include "./misc.hpp"
#include "casmutils/xtal/structure.hpp"
#include "multishift/shift.hpp"
#include "multishift/twist.hpp"
#include <bits/c++config.h>
#include <casmutils/xtal/structure_tools.hpp>
#include <cassert>
#include <filesystem>
#include <string>
#include <unordered_map>
#include <unordered_set>

void write_slicer_structures(const mush::Slicer& slicer, const mush::fs::path& slices_path, std::ostream& log)
{
    log << "Back up prim used...\n";
    mush::cu::xtal::write_poscar(slicer.prim, slices_path / "prim.vasp");

    log << "Write sliced prim...\n";
    mush::cu::xtal::write_poscar(slicer.sliced_prim, slices_path / "sliced_prim.vasp");

    log << "Write aligned sliced prim...\n";
    auto aligned_sliced_lat=mush::make_aligned_lattice(slicer.sliced_prim.lattice());
    mush::cu::xtal::Structure aligned_sliced_prim=slicer.sliced_prim;
    //TODO: Is this really how it works?
    aligned_sliced_prim.set_lattice(aligned_sliced_lat,mush::cu::xtal::FRAC);
    mush::cu::xtal::write_poscar(aligned_sliced_prim, slices_path / "aligned_sliced_prim.vasp");

    log << "Write all possible basis translations...";
    for (int i = 0; i < slicer.floored_sliced_prims.size(); ++i)
    {
        log << "..";
        mush::fs::path target = slices_path / ("sliced_prim_floor." + std::to_string(i) + ".vasp");
        mush::cu::xtal::write_poscar(slicer.floored_sliced_prims[i], target);
    }
    log << std::endl;

    return;
}

std::string make_cleave_dirname(double cleavage)
{
    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << cleavage;
    return "cleave__" + cleavestream.str();
}

std::string make_shift_dirname(int a, int b)
{
    return "shift__" + std::to_string(a)+"."+std::to_string(b);
}

std::unordered_map<std::string, MultiRecord> write_cleaver_structures(const std::vector<Structure>& cleaved_structures,
                                                                      const std::vector<double>& cleavage_values,
                                                                      const mush::fs::path& cleaved_root,
                                                                      std::ostream& log,
                                                                      const MultiRecord& starting_state)
{
    std::unordered_map<std::string, MultiRecord> path_record;
    assert(cleaved_structures.size() == cleavage_values.size());
    log << "Write cleaved structures to..." << cleaved_root << std::endl;
    for (int i = 0; i < cleaved_structures.size(); ++i)
    {
        MultiRecord cleaved_state = starting_state;
        cleaved_state.cleavage = cleavage_values[i];

        //TODO: You have to update the cleavage values of the equivalent ids too!

        mush::fs::path cleave_target = cleaved_root / make_cleave_dirname(cleavage_values[i]);

        path_record[cleave_target] = cleaved_state;

        cautious_create_directory(cleave_target);
        mush::cu::xtal::write_poscar(cleaved_structures[i], cleave_target / "POSCAR");
        std::cout << "Wrote cleaved structure to " << cleave_target << "\n";
    }

    return path_record;
}

std::unordered_map<std::string, MultiRecord> write_shifter_structures(const std::vector<Structure>& shifted_structures,
                                                                      const std::vector<mush::ShiftRecord>& shift_records,
                                                                      const std::vector<std::vector<std::size_t>>& equivalence_map,
                                                                      const mush::fs::path& shifted_root,
                                                                      std::ostream& log,
                                                                      const MultiRecord& starting_state)
{
    std::unordered_map<std::string, MultiRecord> path_record;
    log << "Write shifted structures to..." << shifted_root << std::endl;

    for (int i = 0; i < shifted_structures.size(); ++i)
    {
        MultiRecord shifted_state = starting_state;
        shifted_state.a_index = shift_records[i].a;
        shifted_state.b_index = shift_records[i].b;

        for(auto ix : equivalence_map[i])
        {
            MultiRecord equivalent_state=starting_state;
            equivalent_state.a_index=shift_records[ix].a;
            equivalent_state.b_index=shift_records[ix].b;
            shifted_state.equivalent_structures.push_back(equivalent_state.id());
        }

        mush::fs::path shift_target = shifted_root / make_shift_dirname(shifted_state.a_index,shifted_state.b_index);

        path_record[shift_target] = shifted_state;

        cautious_create_directory(shift_target);
        mush::cu::xtal::write_poscar(shifted_structures[i], shift_target / "POSCAR");
        std::cout << "Wrote shifted structure to " << shift_target << "\n";
    }

    return path_record;
}
