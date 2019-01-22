#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/Lattice.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/external/Eigen/Core"
/* #include "casmutils/frankenstein.hpp" */
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

#include "./base.hpp"
#include "./io.hpp"

namespace mush
{

class ShiftSettings
{
public:
private:
};

class FullSettings
{
public:
    FullSettings(const std::string& init_name, const BaseSettings& init_base_settings);

    static FullSettings from_json(const CASM::jsonParser& init_json);
    static FullSettings from_path(const CASM::fs::path& init_path);

    const std::string& name() const { return m_name; }

    const BaseSettings& base_settings() const { return m_base_settings; }

private:
    /// Name for the particular set of shifts you're going to do
    std::string m_name;

    /// Specifies all the settings involved with creating a slab
    BaseSettings m_base_settings;
};

FullSettings::FullSettings(const std::string& init_name, const BaseSettings& init_base_settings)
    : m_name(init_name), m_base_settings(init_base_settings)
{
}

FullSettings FullSettings::from_json(const CASM::jsonParser& init_json)
{
    return FullSettings(init_json["name"].get<std::string>(), BaseSettings::from_json(init_json[BaseSettings::tag()]));
}

FullSettings FullSettings::from_path(const CASM::fs::path& init_path)
{
    CASM::jsonParser settings_dump(init_path);
    return FullSettings::from_json(settings_dump);
}

/// Given the raw shift unit, return the same structure with the basis translated such that
/// the given coordinate ends up on the shift plane
/* CASM::Structure cooked_shift_unit(const CASM::Structure& raw_shift_unit, const CASM::Coordinate& coord_to_floor) */
/* { */
/*     Rewrap::Structure cooked_shift_unit(raw_shift_unit); */
/*     Frankenstein::shift_coords_by(&cooked_shift_unit, -coord_to_floor.cart()); */
/*     return cooked_shift_unit; */
/* } */

} // namespace mush

/// Throws an exception if a bad combination of program options was given
void process_program_options(CASM::po::variables_map& vm)
{
    CASM::po::notify(vm);

    return;
}

int main_candidate(int argc, char* argv[])
{
    // clang-format off
    CASM::po::options_description desc("multishift options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("settings,s", CASM::po::value<CASM::fs::path>()->required(),"Path to file specifying how to perform shifts");
    // clang-format on

    CASM::po::variables_map vm;
    CASM::po::store(CASM::po::parse_command_line(argc, argv, desc), vm);

    if (vm.count("help"))
    {
        std::cout << desc << std::endl;
        return 0;
    }

    try
    {
        process_program_options(vm);
    }

    catch (CASM::po::error const& e)
    {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    std::cout << "lol" << std::endl;
}

int main()
{
    auto all_settings=mush::FullSettings::from_path("./mush.json");
    std::cout<<all_settings.name()<<std::endl;

    const auto& base_settings = all_settings.base_settings();
    std::cout << base_settings.prim_path() << std::endl;
    std::cout << base_settings.millers().transpose() << std::endl;
    std::cout << base_settings.floor_slab_atom_index() << std::endl;
    std::cout << base_settings.stacks() << std::endl;

    auto base=mush::MultiBase::from_settings(all_settings.base_settings());

    mush::MultiIO writer(all_settings.name());
    writer.drop_base(base_settings,base);

    return 0;
}
