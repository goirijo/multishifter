#ifndef SHIFT_SUBCOMMAND_HH
#define SHIFT_SUBCOMMAND_HH

#include "misc.hpp"
#include <casmutils/mush/shift.hpp>
#include <multishift/definitions.hpp>
#include <casmutils/xtal/structure.hpp>
#include <CLI/CLI.hpp>

void setup_subcommand_shift(CLI::App& app);
//Callback is implemented in chain command

#endif
