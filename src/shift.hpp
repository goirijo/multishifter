#ifndef SHIFT_SUBCOMMAND_HH
#define SHIFT_SUBCOMMAND_HH

#include "misc.hpp"
#include <casmutils/mush/shift.hpp>
#include <multishift/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <CLI/CLI.hpp>

std::unordered_map<std::string, MultiRecord> write_shifter_structures(const std::vector<Structure>& shifted_structures,
                                                                      const std::vector<mush::ShiftRecord>& shift_records,
                                                                      const std::vector<std::vector<std::size_t>>& equivalence_map,
                                                                      const mush::fs::path& shifted_root,
                                                                      std::ostream& log,
                                                                      const MultiRecord& starting_state);

void setup_subcommand_shift(CLI::App& app);
void run_subcommand_shift(const mush::fs::path& settings_path, std::ostream& log);


#endif
