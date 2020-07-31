#ifndef CLEAVE_SUBCOMMAND_HH
#define CLEAVE_SUBCOMMAND_HH

#include "misc.hpp"
#include <multishift/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <CLI/CLI.hpp>

std::unordered_map<std::string, MultiRecord> write_cleaver_structures(const std::vector<Structure>& cleaved_structures,
                                                                         const std::vector<double>& cleavage_values,
                                                                         const mush::fs::path& cleaved_root,
                                                                         std::ostream& log,
                                                                         const MultiRecord& starting_state);


void setup_subcommand_cleave(CLI::App& app);
void run_subcommand_cleave(const mush::fs::path& settings_path, std::ostream& log);
#endif
