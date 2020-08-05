#include "./slice.hpp"
#include "./misc.hpp"
#include "casmutils/xtal/structure.hpp"
#include "multishift/slice_settings.hpp"
#include <casmutils/mush/shift.hpp>
#include <casmutils/mush/twist.hpp>
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
    auto aligned_sliced_lat = mush::make_aligned_lattice(slicer.sliced_prim.lattice());
    mush::cu::xtal::Structure aligned_sliced_prim = slicer.sliced_prim;
    // TODO: Is this really how it works?
    aligned_sliced_prim.set_lattice(aligned_sliced_lat, mush::cu::xtal::FRAC);
    mush::cu::xtal::write_poscar(aligned_sliced_prim, slices_path / "aligned_sliced_prim.vasp");

    log << "Write all possible basis translations...";
    for (int i = 0; i < slicer.floored_sliced_prims.size(); ++i)
    {
        log << "..";
        const auto& floored_slice = slicer.floored_sliced_prims[i];
        mush::fs::path target = slices_path / ("sliced_prim_floor." + std::to_string(i) + ".vasp");
        mush::cu::xtal::write_poscar(floored_slice, target);

        auto aligned_floored_lat = mush::make_aligned_lattice(floored_slice.lattice());
        mush::cu::xtal::Structure aligned_floored_slice = floored_slice.set_lattice(aligned_floored_lat, mush::cu::xtal::FRAC);
        mush::fs::path aligned_target = slices_path / ("aligned_sliced_prim_floor." + std::to_string(i) + ".vasp");
        mush::cu::xtal::write_poscar(aligned_floored_slice, aligned_target);
    }
    log << std::endl;

    return;
}

void setup_subcommand_slice(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* slice_sub =
        app.add_subcommand("slice", "Slice unit cell to expose desired plane. Creates input structures for shift/cleave/twist.");
    slice_sub->add_option("-s,--settings", *settings_path_ptr, "Settings file slicing primitive structure.")->required();

    slice_sub->callback([settings_path_ptr]() { run_subcommand_slice(*settings_path_ptr, std::cout); });
}

void run_subcommand_slice(const mush::fs::path& settings_path, std::ostream& log)
{
    json settings = load_json(settings_path);
    std::string project_name = extract_name_from_settings(settings);

    log << "Project name: " << project_name << std::endl;
    auto slice_settings = mush::SliceSettings::from_json(settings);

    mush::Slicer slicer(Structure::from_poscar(slice_settings.prim_path), slice_settings.miller_indexes);

    mush::fs::path slices_path(project_name + ".slices");
    cautious_create_directory(slices_path);

    log << "Back up settings used...\n";
    write_json(settings, slices_path / "slice.json");

    write_slicer_structures(slicer, slices_path, log);
}
