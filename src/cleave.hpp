#ifndef CLEAVE_SUBCOMMAND_HH
#define CLEAVE_SUBCOMMAND_HH

#include "./chain.hpp"
#include <multishift/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <CLI/CLI.hpp>


void setup_subcommand_cleave(CLI::App& app);
//Callback is implemented in chain command

#endif
