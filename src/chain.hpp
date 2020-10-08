#ifndef CHAIN_SUBCOMMAND_HH
#define CHAIN_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>

void setup_subcommand_chain(CLI::App& app);
void run_subcommand_chain(const mush::fs::path& input_path, const mush::fs::path& output_dir, const std::vector<double>& cleavages, const std::vector<int>& grid_dims, std::ostream& log);

#endif
