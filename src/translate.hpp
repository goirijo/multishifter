#ifndef TRANSLATE_COMMAND_HH
#define TRANSLATE_COMMAND_HH

#include <CLI/CLI.hpp>
#include <filesystem>
#include <multishift/definitions.hpp>

void setup_subcommand_translate(CLI::App& app);
void run_subcommand_translate(const mush::fs::path& input_path, const mush::fs::path& output_path, const std::vector<double>& shift, int floor_ix, bool frac, std::ostream& log);

#endif
