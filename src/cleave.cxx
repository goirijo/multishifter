#include "./cleave.hpp"
#include "./misc.hpp"
#include <casmutils/xtal/structure_tools.hpp>
#include <multishift/shift.hpp>
#include <multishift/slice_settings.hpp>

std::string make_cleave_dirname(double cleavage)
{
    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << cleavage;
    return "cleave__" + cleavestream.str();
}

template <COMMAND>
std::unordered_map<std::string, MultiRecord>
run(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists);

template <>
std::unordered_map<std::string, MultiRecord> run<COMMAND::CLEAVE>(
    const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists)
{
    auto cleavage_settings = mush::CleavageSettings::from_json(settings);
    auto slaber_settings = mush::SlabSettings::from_json(settings);

    log << "Stacking slab...\n";
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

        // TODO: You have to update the cleavage values of the equivalent ids too!

        mush::fs::path cleave_target = cleaved_root / make_cleave_dirname(cleavage_values[i]);

        path_record[cleave_target] = cleaved_state;

        cautious_create_directory(cleave_target);
        mush::cu::xtal::write_poscar(cleaved_structures[i], cleave_target / "POSCAR");
        std::cout << "Wrote cleaved structure to " << cleave_target << "\n";
    }

    return path_record;
}

void setup_subcommand_cleave(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* cleave_sub = app.add_subcommand("cleave", "Cleave unit cell to insert/remove empty space at the interface layer.");
    cleave_sub
        ->add_option("-s,--settings", *settings_path_ptr, "Settings file with path to slab unit, slab thickness, and cleavage values.")
        ->required();

    cleave_sub->callback([settings_path_ptr]() { run_subcommand_cleave(*settings_path_ptr, std::cout); });
}

void run_subcommand_cleave(const mush::fs::path& settings_path, std::ostream& log)
{
    json settings = load_json(settings_path);
    std::string project_name = extract_name_from_settings(settings);
    mush::fs::path cleaved_path(project_name + ".cleave");

    log << "Project name: " << project_name << std::endl;
    auto path_record = run<COMMAND::CLEAVE>(settings, cleaved_path, MultiRecord(), log, false);

    log << "Save record of structures to " << cleaved_path << "\n";
    write_json(record_to_json(path_record), cleaved_path / "record.json");
}

