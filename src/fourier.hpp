#ifndef FOURIERCOMMAND_HH
#define FOURIERCOMMAND_HH

#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>

void setup_subcommand_fourier(CLI::App& app);
void run_subcommand_fourier(const mush::fs::path& data_path, double cleavage_slice, const std::string& value_key, double crush_value, std::ostream& log);

#endif
