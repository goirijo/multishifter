#include <CLI/CLI.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <filesystem>
#include <memory>
#include <multishift/definitions.hpp>
#include <multishift/shift.hpp>
#include <multishift/shifter.hpp>
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
#include "multishift/shifter.hpp"

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
run(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists);

template <>
std::unordered_map<std::string, MultiRecord> run<COMMAND::CLEAVE>(
    const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists)
{
    auto cleavage_settings = mush::CleavageSettings::from_json(settings);
    auto slaber_settings = mush::SlabSettings::from_json(settings);

    auto slab = mush::make_stacked_slab(Structure::from_poscar(slaber_settings.slab_unit_path), slaber_settings.stacks);

    log << "Cleaving structures...\n";
    auto cleaved_structures = mush::make_cleaved_structures(slab, cleavage_settings.cleavage_values);

    if (root_exists)
    {
        mush::fs::create_directory(root);
    }
    else
    {
        cautious_create_directory(root);
    }

    log << "Back up settings used...\n";
    write_json(settings, root / "cleave.json");
    log << "Back up slab...\n";
    mush::cu::xtal::write_poscar(slab, root / "slab.vasp");

    return write_cleaver_structures(cleaved_structures, cleavage_settings.cleavage_values, root, log, starting_state);
}

template <>
std::unordered_map<std::string, MultiRecord> run<COMMAND::SHIFT>(
    const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists)
{
    auto shift_settings = mush::ShiftSettings::from_json(settings);
    auto slaber_settings = mush::SlabSettings::from_json(settings);

    auto slab = mush::make_stacked_slab(Structure::from_poscar(slaber_settings.slab_unit_path), slaber_settings.stacks);

    log << "Shifting structures, please be patient...\n";
    mush::Shifter shifter(slab, shift_settings.a_density, shift_settings.b_density);

    if (root_exists)
    {
        mush::fs::create_directory(root);
    }
    else
    {
        cautious_create_directory(root);
    }

    log << "Back up settings used...\n";
    write_json(settings, root / "cleave.json");
    log << "Back up slab...\n";
    mush::cu::xtal::write_poscar(slab, root / "slab.vasp");

    return write_shifter_structures(shifter.shifted_structures, shifter.shift_records, shifter.equivalence_map, root, log, starting_state);
}

template <>
std::unordered_map<std::string, MultiRecord> run<COMMAND::TWIST>(
    const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists)
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
        constexpr auto* run_cleave = run<COMMAND::CLEAVE>;
        constexpr auto* run_shift = run<COMMAND::SHIFT>;
        constexpr auto* run_twist = run<COMMAND::TWIST>;

        std::unordered_map<std::string, decltype(run_cleave)> command_dispacher{
            {"cleave", run_cleave}, {"shift", run_shift}, {"twist", run_twist}};

        json settings = load_json(settings_path);
        std::string project_name = extract_name_from_settings(settings);

        std::vector<std::string> executions = settings["execute"];
        mush::fs::path true_root(project_name + ".chain");
        bool root_already_exists = false;
        std::unordered_map<std::string, MultiRecord> root_dirs{{true_root, MultiRecord()}};

        log << "Project name: " << project_name << std::endl;

        std::vector<std::unordered_map<std::string, MultiRecord>> full_records_data;
        for (const std::string& command : executions)
        {
            decltype(root_dirs) next_dirs;
            decltype(full_records_data)::value_type single_record;
            for (const auto& [root, state] : root_dirs)
            {
                if (command != executions[0])
                {
                    settings["stacks"] = 1;
                    // TODO: Make the slab.json path a function, since it's used several places
                    settings["slab_unit"] = mush::fs::path(root) / "POSCAR";
                    root_already_exists = true;
                }
                auto partial_record_data = command_dispacher[command](settings, root, state, log, root_already_exists);

                single_record.insert(partial_record_data.begin(),
                                     partial_record_data.end());
                next_dirs.insert(partial_record_data.begin(), partial_record_data.end());
            }
            root_dirs = next_dirs;
            full_records_data.push_back(single_record);
        }

        log << "Save record of structures to " << true_root << "\n";

        json full_record_json;
        std::string chained_command;
        assert(full_records_data.size()==executions.size());
        for(int i=0; i<full_records_data.size(); ++i)
        {
            chained_command+=executions[i];
            full_record_json[chained_command]=record_to_json(full_records_data[i]);
            if(i!=full_records_data.size()-1)
            {
                chained_command+="-";
            }
        }
        write_json(full_record_json, true_root / "record.json");
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
        auto path_record = run<COMMAND::CLEAVE>(settings, cleaved_path, MultiRecord(), log, false);

        log << "Save record of structures to " << cleaved_path << "\n";
        write_json(record_to_json(path_record), cleaved_path / "record.json");
    }

    if (shift_sub->count("--settings"))
    {
        json settings = load_json(settings_path);
        std::string project_name = extract_name_from_settings(settings);
        mush::fs::path shifted_path(project_name + ".shift");

        log << "Project name: " << project_name << std::endl;
        auto path_record = run<COMMAND::SHIFT>(settings, shifted_path, MultiRecord(), log, false);

        log << "Save record of structures to " << shifted_path << "\n";
        write_json(record_to_json(path_record), shifted_path / "record.json");
    }

    return 0;
}
