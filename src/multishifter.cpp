#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>
#include <multishift/slab.hpp>
#include <multishift/slicer.hpp>
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
        json raw_slice_settings=extract_slice_subsettings(settings);

        log<<"Project name: "<<project_name<<std::endl;
        auto slice_settings=mush::SliceSettings::from_json(raw_slice_settings);
    
        mush::Slicer slicer(Structure::from_poscar(slice_settings.prim_path),slice_settings.miller_indexes);

        mush::fs::path slices_path(project_name+".slices");
        cautious_create_directory(slices_path);

        log<<"Back up settings used...\n";
        write_json(raw_slice_settings,slices_path/"slice.json");

        write_slicer_structures(slicer,slices_path,log);
    }


    return 0;
}
