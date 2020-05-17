#include "./slice.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include <cassert>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>

void write_slicer_structures(const mush::Slicer& slicer, const mush::fs::path& slices_path, std::ostream& log)
{
    log << "Back up prim used...\n";
    mush::cu::xtal::write_poscar(slicer.prim, slices_path / "prim.vasp");

    log << "Write sliced prim...\n";
    mush::cu::xtal::write_poscar(slicer.sliced_prim, slices_path / "sliced_prim.vasp");

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
    return "cleave__"+cleavestream.str();
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

        mush::fs::path cleave_target = cleaved_root / make_cleave_dirname(cleavage_values[i]);

        path_record[cleave_target]=cleaved_state;

        cautious_create_directory(cleave_target);
        mush::cu::xtal::write_poscar(cleaved_structures[i],cleave_target/"POSCAR");
        std::cout<<"Wrote cleaved structure to "<<cleave_target<<"\n";
    }

    return path_record;
}
