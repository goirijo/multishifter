#include "./shift.hpp"
#include "multishift/shifter.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include <multishift/slice_settings.hpp>

template <COMMAND>
std::unordered_map<std::string, MultiRecord>
run(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists);

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

std::string make_shift_dirname(int a, int b) { return "shift__" + std::to_string(a) + "." + std::to_string(b); }

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

        for (auto ix : equivalence_map[i])
        {
            MultiRecord equivalent_state = starting_state;
            equivalent_state.a_index = shift_records[ix].a;
            equivalent_state.b_index = shift_records[ix].b;
            shifted_state.equivalent_structures.push_back(equivalent_state.id());
        }

        mush::fs::path shift_target = shifted_root / make_shift_dirname(shifted_state.a_index, shifted_state.b_index);

        path_record[shift_target] = shifted_state;

        cautious_create_directory(shift_target);
        mush::cu::xtal::write_poscar(shifted_structures[i], shift_target / "POSCAR");
        std::cout << "Wrote shifted structure to " << shift_target << "\n";
    }

    return path_record;
}

void setup_subcommand_shift(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* shift_sub = app.add_subcommand(
        "shift", "Shift unit cell along interface to generate stacking faults and structures for gamma surface calculations.");
    shift_sub->add_option("-s,--settings", *settings_path_ptr, "Settings file with path to slab unit, slab thickness, and grid density.")
        ->required();

    shift_sub->callback([settings_path_ptr]() { run_subcommand_shift(*settings_path_ptr, std::cout); });
}

void run_subcommand_shift(const mush::fs::path& settings_path, std::ostream& log)
{
    json settings = load_json(settings_path);
    std::string project_name = extract_name_from_settings(settings);
    mush::fs::path shifted_path(project_name + ".shift");

    log << "Project name: " << project_name << std::endl;
    auto path_record = run<COMMAND::SHIFT>(settings, shifted_path, MultiRecord(), log, false);

    log << "Save record of structures to " << shifted_path << "\n";
    write_json(record_to_json(path_record), shifted_path / "record.json");
}
