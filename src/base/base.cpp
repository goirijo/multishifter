#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
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

int _main()
{
    auto all_settings = mush::FullSettings::from_path("./mush.json");
    std::cout << all_settings.name() << std::endl;

    //***************

    const auto& base_settings = all_settings.base_settings();
    std::cout << base_settings.prim_path() << std::endl;
    std::cout << base_settings.millers().transpose() << std::endl;
    std::cout << base_settings.floor_slab_atom_index() << std::endl;
    std::cout << base_settings.stacks() << std::endl;

    auto base = mush::MultiBase::from_settings(all_settings.base_settings());

    mush::MultiIO writer(all_settings.name());
    writer.drop_base(base);

    base_settings.to_json().write(writer.base_target() / (mush::BaseSettings::tag() + ".json"));

    //***************

    const auto& shift_settings = all_settings.shift_settings();
    std::cout << shift_settings.slab_path() << std::endl;
    std::cout << shift_settings.a_points() << std::endl;
    std::cout << shift_settings.b_points() << std::endl;
    for (const auto& val : shift_settings.cleavage_values())
    {
        std::cout << val << " ";
    }
    std::cout << std::endl;

    auto shifter = mush::MultiShift::from_settings(shift_settings);
    writer.drop_shifts(shifter);

    shift_settings.to_json().write(writer.shift_target() / (mush::ShiftSettings::tag() + ".json"));
    Simplicity::write_poscar(shifter.reference_slab(), writer.shift_target() / "slab.json");

    //***************

    return 0;
}

namespace cxxopts
{
void required_argument_notify(const ParseResult& result, const std::vector<std::string>& required_arguments)
{

    for (const auto& arg : required_arguments)
    {
        if (!result.count(arg))
        {
            throw mush::except::RequiredArgumentMissing(arg);
        }
    }

    return;
}
} // namespace cxxopts

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
        std::cout << "Requested slcicing along the " << base_settings.millers().transpose() << " direction."
                  << std::endl;
        if (!base_settings.floor_slab_atom_index())
        {
            std::cout << "No basis translation requested." << std::endl;
        }
        else
        {
            std::cout << "Will translate basis to bring atom " << base_settings.floor_slab_atom_index()
                      << " down (indexing starts at 1!)." << std::endl;
            std::cout << "Stack sliced unit " << base_settings.stacks() << "times." << std::endl;
        }

        std::cout << "..........................................................." << std::endl;
        auto base = mush::MultiBase::from_settings(all_settings.base_settings());
        std::cout << "Structures constructed." << std::endl;
        std::cout << std::setw(15) << "Primitive cell: " << std::setw(15) << base.primitive().basis.size() << " atoms."
                  << std::endl;
        std::cout << std::setw(15) << "Shift unit: " << std::setw(15) << base.shift_unit().basis.size() << " atoms."
                  << std::endl;
        std::cout << std::setw(15) << "Slab: " << std::setw(15) << base.floored_slab().basis.size() << " atoms."
                  << std::endl;

        mush::MultiIO writer(all_settings.name());

        writer.drop_base(base);

        base_settings.to_json().write(writer.base_target() / (mush::BaseSettings::tag() + ".json"));

        std::cout << "..........................................................." << std::endl;
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
}
