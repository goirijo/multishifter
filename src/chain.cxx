#include "./chain.hpp"
#include "./misc.hpp"
#include <nlohmann/json.hpp>

template <COMMAND>
std::unordered_map<std::string, MultiRecord>
run(const json& settings, const mush::fs::path& root, const MultiRecord& starting_state, std::ostream& log, bool root_exists);

void setup_subcommand_chain(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* chain_sub = app.add_subcommand("chain", "Combine slice, stack, cleave, and shift commands for gamma surface calculations.");
    chain_sub->add_option("-s,--settings", *settings_path_ptr, "Union of slice, cleave, and shift settings");

    chain_sub->callback([settings_path_ptr]() { run_subcommand_chain(*settings_path_ptr, std::cout); });
}

void run_subcommand_chain(const mush::fs::path& settings_path, std::ostream& log)
{
    constexpr auto* run_cleave = run<COMMAND::CLEAVE>;
    constexpr auto* run_shift = run<COMMAND::SHIFT>;

    std::unordered_map<std::string, decltype(run_cleave)> command_dispacher{
            {"cleave", run_cleave}, {"shift", run_shift}};

    json settings = load_json(settings_path);
    std::string project_name = extract_name_from_settings(settings);

    std::vector<std::string> executions{"shift","cleave"};
    mush::fs::path true_root(project_name + ".chain");
    bool root_already_exists = false;
    std::unordered_map<std::string, MultiRecord> root_dirs{{true_root, MultiRecord()}};

    log << "Project name: " << project_name << std::endl;

    std::vector<std::unordered_map<std::string, MultiRecord>> full_records_data;
    for (auto command : executions)
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

            single_record.insert(partial_record_data.begin(), partial_record_data.end());
            next_dirs.insert(partial_record_data.begin(), partial_record_data.end());
        }
        root_dirs = next_dirs;
        full_records_data.push_back(single_record);
    }

    log << "Save record of structures to " << true_root << "\n";

    json full_record_json;
    std::string chained_command;
    assert(full_records_data.size() == executions.size());
    for (int i = 0; i < full_records_data.size(); ++i)
    {
        chained_command += executions[i];
        full_record_json[chained_command] = record_to_json(full_records_data[i]);
        if (i != full_records_data.size() - 1)
        {
            chained_command += "-";
        }
    }
    write_json(full_record_json, true_root / "record.json");
}

