#ifndef TWIST_SUBCOMMAND_HH
#define TWIST_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>

void setup_subcommand_twist(CLI::App& app);
void run_subcommand_twist(const mush::fs::path& input_path, const mush::fs::path& output_dir, const std::vector<double>& angles, int max_lattice_sites, double error_tol, std::string zone, std::string supercells, std::ostream& log);

#endif
