#include <CLI/CLI.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <filesystem>
#include <memory>
#include <multishift/definitions.hpp>
#include <multishift/shift.hpp>
#include <multishift/slab.hpp>
#include <multishift/slice_settings.hpp>
#include <multishift/slicer.hpp>
#include <nlohmann/json.hpp>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "./misc.hpp"
#include "./slice.hpp"

using json = nlohmann::json;
using Structure = mush::cu::xtal::Structure;

enum class COMMAND
{
    SHIFT,
    CLEAVE,
    TWIST
};

template <COMMAND>
std::unordered_map<std::string, MultiRecord>
run(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log);

template <>
std::unordered_map<std::string, MultiRecord>
run<COMMAND::CLEAVE>(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log)
{
    auto cleavage_settings = mush::CleavageSettings::from_json(settings);
    auto slaber_settings = mush::SlabSettings::from_json(settings);

    auto slab = mush::make_stacked_slab(Structure::from_poscar(slaber_settings.slab_unit_path), slaber_settings.stacks);
    auto cleaved_structures = mush::make_cleaved_structures(slab, cleavage_settings.cleavage_values);

    cautious_create_directory(root);

    log << "Back up settings used...\n";
    write_json(settings, root / "cleave.json");
    log << "Back up slab...\n";
    mush::cu::xtal::write_poscar(slab, root / "slab.vasp");

    return write_cleaver_structures(cleaved_structures, cleavage_settings.cleavage_values, root, log, starting_state);
}

template <>
std::unordered_map<std::string, MultiRecord>
run<COMMAND::SHIFT>(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log)
{
    throw std::runtime_error("Not implemented");
std::unordered_map<std::string, MultiRecord> tmp;
return tmp;
}

template <>
std::unordered_map<std::string, MultiRecord>
run<COMMAND::TWIST>(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log)
{
    throw std::runtime_error("Not implemented");
std::unordered_map<std::string, MultiRecord> tmp;
return tmp;
}

int main(int argc, char** argv)
{
    mush::fs::path settings_path;

    CLI::App app{"Slice, stack, shift, and rotate crystal slabs"};
    app.add_option("-s,--settings", settings_path, "Settings file for all possible operations");

    CLI::App* slice_sub = app.add_subcommand("slice", "Slice unit cell to expose desired plane. Creates input structures for shift/twist.");
    slice_sub->add_option("-s,--settings", settings_path, "Settings file slicing primitive structure.")->required();

    CLI::App* cleave_sub = app.add_subcommand("cleave", "Cleave unit cell to insert/remove empty space at the interface layer.");
    cleave_sub->add_option("-s,--settings", settings_path, "Settings file with path to slab unit, slab thickness, and cleavage values.")
        ->required();

    CLI::App* shift_sub = app.add_subcommand(
        "shift", "Shift unit cell along interface to generate stacking faults and structures for gamma surface calculations.");
    shift_sub->add_option("-s,--settings", settings_path, "Settings file with path to slab unit, slab thickness, and grid density.")
        ->required();

    CLI11_PARSE(app, argc, argv);
    auto& log = std::cout;

    if (app.count("--settings"))
    {
        /* log << "Asked for all (unimplemented) \n"; */
        /* return 0; */

        constexpr auto* run_cleave = run<COMMAND::CLEAVE>;
        constexpr auto* run_shift = run<COMMAND::SHIFT>;
        constexpr auto* run_twist = run<COMMAND::TWIST>;

        std::unordered_map<std::string, decltype(run_cleave)> command_dispacher{
            {"cleave", run_cleave}, {"shift", run_shift}, {"twist", run_twist}};

        json settings = load_json(settings_path);
        std::string project_name = extract_name_from_settings(settings);

        std::vector<std::string> executions=settings["execute"];
        mush::fs::path true_root(project_name + ".chain");
        std::unordered_map<std::string, MultiRecord> root_dirs{{true_root, MultiRecord()}};

        log << "Project name: " << project_name << std::endl;

        std::unordered_map<std::string, MultiRecord> full_record_data;
        for (const std::string& command : executions)
        {
            decltype(root_dirs) next_dirs;
            for (const auto& [root, state] : root_dirs)
            {
                auto partial_record_data = command_dispacher[command](settings, root, state, log);

                full_record_data.insert(partial_record_data.begin(), partial_record_data.end());
                next_dirs.insert(partial_record_data.begin(), partial_record_data.end());
            }
            root_dirs = next_dirs;
        }

        log << "Save record of structures to " << true_root << "\n";
        write_json(record_to_json(full_record_data), true_root / "record.json");
    }

    if (slice_sub->count("--settings"))
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

    if (cleave_sub->count("--settings"))
    {
        json settings = load_json(settings_path);
        std::string project_name = extract_name_from_settings(settings);
        mush::fs::path cleaved_path(project_name + ".cleave");

        log << "Project name: " << project_name << std::endl;
        log << "Save record of structures to " << cleaved_path << "\n";
        auto path_record = run<COMMAND::CLEAVE>(settings, cleaved_path, MultiRecord(), log);

        write_json(record_to_json(path_record), cleaved_path / "record.json");
    }

    return 0;
}
