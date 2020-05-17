#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>
#include <multishift/slab.hpp>
#include <multishift/shift.hpp>
#include <multishift/slicer.hpp>
#include <multishift/slice_settings.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <nlohmann/json.hpp>
#include <stdexcept>
#include <string>

#include "./misc.hpp"
#include "./slice.hpp"

using json = nlohmann::json;
using Structure=mush::cu::xtal::Structure;

int main(int argc, char** argv)
{
    mush::fs::path settings_path;

    CLI::App app{"Slice, stack, shift, and rotate crystal slabs"};
    app.add_option("-s,--settings", settings_path, "Settings file for all possible operations");

    CLI::App* slice_sub = app.add_subcommand("slice", "Slice unit cell to expose desired plane. Creates input structures for shift/twist.");
    slice_sub->add_option("-s,--settings", settings_path, "Settings file slicing primitive structure.")->required();

    CLI::App* cleave_sub = app.add_subcommand("cleave", "Cleave unit cell to insert/remove empty space at the interface layer.");
    cleave_sub->add_option("-s,--settings", settings_path, "Settings file with path to slab unit, slab thickness, and cleavage values.")->required();

    CLI::App* shift_sub = app.add_subcommand("shift", "Shift unit cell along interface to generate stacking faults and structures for gamma surface calculations.");
    shift_sub->add_option("-s,--settings", settings_path, "Settings file with path to slab unit, slab thickness, and grid density.")->required();

    CLI11_PARSE(app, argc, argv);
    auto& log=std::cout;

    if (app.count("--settings"))
    {
        log << "Asked for all (unimplemented) \n";
        return 0;
    }

    if (slice_sub->count("--settings"))
    {
        json settings = load_json(settings_path);
        std::string project_name=extract_name_from_settings(settings);

        log<<"Project name: "<<project_name<<std::endl;
        auto slice_settings=mush::SliceSettings::from_json(settings);
    
        mush::Slicer slicer(Structure::from_poscar(slice_settings.prim_path),slice_settings.miller_indexes);

        mush::fs::path slices_path(project_name+".slices");
        cautious_create_directory(slices_path);

        log<<"Back up settings used...\n";
        write_json(settings,slices_path/"slice.json");

        write_slicer_structures(slicer,slices_path,log);
    }

    if(cleave_sub->count("--settings"))
    {
        json settings = load_json(settings_path);
        std::string project_name=extract_name_from_settings(settings);

        log<<"Project name: "<<project_name<<std::endl;
        auto cleavage_settings=mush::CleavageSettings::from_json(settings);
        auto slaber_settings=mush::SlabSettings::from_json(settings);

        auto slab=mush::make_stacked_slab(Structure::from_poscar(slaber_settings.slab_unit_path),slaber_settings.stacks);
        auto cleaved_structures=mush::make_cleaved_structures(slab,cleavage_settings.cleavage_values);

        mush::fs::path cleaved_path(project_name+".cleave");
        cautious_create_directory(cleaved_path);

        log<<"Back up settings used...\n";
        write_json(settings,cleaved_path/"cleave.json");
        log<<"Back up slab...\n";
        mush::cu::xtal::write_poscar(slab,cleaved_path/"slab.vasp");

        auto path_record=write_cleaver_structures(cleaved_structures, cleavage_settings.cleavage_values, cleaved_path, log, MultiRecord());

        log<<"Save record of structures to "<<cleaved_path<<"\n";
        write_json(record_to_json(path_record),cleaved_path/"record.json");
    }


    return 0;
}
