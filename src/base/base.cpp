#include "casm/casm_io/jsonParser.hh"
/* #include "casm/external/Eigen/Core" */
/* #include "casmutils/frankenstein.hpp" */
#include "cxxopts.hpp"
#include <iomanip>
#include <iostream>
#include <string>

#include "multishift/define.hpp"
#include "multishift/exceptions.hpp"
#include "multishift/io.hpp"
#include "multishift/settings.hpp"
#include "multishift/misc.hpp"

/**
 * multishift-base takes a primitive cell and some settings as input.
 * The primitive cell is sliced along the specified miller indexes, stacked a number of times
 * to create a slab, and then the slab has its basis translated so that the
 * gamma surface ends up between the correct atomic layer.
 */

int main(int argc, char* argv[])
{
    using CASM::fs::path;

    cxxopts::Options options("multishift-base", "Slice the primitive cell and create a slab.");

    // clang-format off
    options.add_options()
        /* ("p,prim", "Path to the primitive cell (VASP format)",cxxopts::value<path>()) */
        ("s,settings","Path to the settings file (json).",cxxopts::value<path>())
        /* ("n,nuke","Delete output directory and write everything out again.") */
        ("h,help","Print available options.");
    // clang-format on

    try
    {
        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0;
        }

        else
        {
            /* cxxopts::required_argument_notify(result, std::vector<std::string>{"prim", "settings"}); */
            cxxopts::required_argument_notify(result, std::vector<std::string>{"settings"});
        }

        const auto& settings_path = result["settings"].as<path>();
        /* const auto& prim_path=result["prim"].as<path>(); */

        auto all_settings = mush::FullSettings::from_path(settings_path);
        std::cout << "Using '" << all_settings.name() << "' as prefix." << std::endl;

        const auto& base_settings = all_settings.base_settings();
        std::cout << "Search in " << base_settings.prim_path() << " for primitive structure." << std::endl;
        std::cout << "Requested slicing along the " << base_settings.millers().transpose() << " direction."
                  << std::endl;
        if (!base_settings.floor_slab_atom_index())
        {
            std::cout << "No basis translation requested." << std::endl;
        }
        else
        {
            std::cout << "Will translate basis to bring atom " << base_settings.floor_slab_atom_index()
                      << " down (indexing starts at 1!)." << std::endl;
            std::cout << "Stack sliced unit " << base_settings.stacks() << " times." << std::endl;
        }

        loggy::divider();

        auto base = mush::MultiBase::from_settings(all_settings.base_settings());
        std::cout << "Structures constructed." << std::endl;
        std::cout << std::setw(20) << "Primitive cell: " << std::setw(20) << base.primitive().basis.size() << " atoms."
                  << std::endl;
        std::cout << std::setw(20) << "Shift unit: " << std::setw(20) << base.shift_unit().basis.size() << " atoms."
                  << std::endl;
        std::cout << std::setw(20) << "Slab: " << std::setw(20) << base.floored_slab().basis.size() << " atoms."
                  << std::endl;

        mush::MultiIO writer(all_settings.name());

        writer.drop_base(base);

        base_settings.to_json().write(writer.base_target() / (mush::BaseSettings::tag() + ".json"));

        loggy::divider();

        std::cout << "Structures written to " << writer.base_target() << " along with a backup of the settings used."
                  << std::endl;
    }

    catch (const cxxopts::missing_argument_exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Use '--help' to see option arguments" << std::endl;
        return 1;
    }

    catch (const mush::except::RequiredArgumentMissing& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Use '--help' to see options" << std::endl;
        return 2;
    }

    catch (const mush::except::UnspecifiedSettings& e)
    {
        // TODO: Be more helpful
        std::cerr << e.what() << std::endl;
        return 3;
    }

    catch (const mush::except::FileExists& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Rename it or delete it if you don't need it." << std::endl;
        return 4;
    }

    //TODO: Catch file doesn't exist (bad settings path)

    return 0;
}
