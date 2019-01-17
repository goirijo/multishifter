#include "casm/casm_io/jsonParser.hh"
#include "casm/crystallography/Structure.hh"
#include "casm/external/Eigen/Core"
#include <boost/program_options.hpp>
#include <iostream>
#include <string>

namespace mush
{

/**
 * Simple structure to hold any settings required to generate the shifted
 * structures, starting from the primitive cell.
 * A bit on nomenclature:
 * "primitive" means the smallest unit that you begin your project with.
 * "shift unit" is the smallest structure possible that exposes the plane you want
 * to shift along as one of its unit cell facets, it is a supercell of the
 * primitive.
 * "stack" or "slab" is a 1 dimensional supercell of the slab, constructed by stacking
 * multiple slabs along the shift plane.
 */

class SlabSettings
{
public:
    SlabSettings(const CASM::jsonParser& init_settings);

    const CASM::fs::path& prim_path() const { return m_prim_path; }
    const Eigen::Vector3i& millers() const { return m_millers; }
    const int& floor_slab_atom() const { return m_floor_slab_atom; }
    const int& stacks() const { return m_stacks; }

private:
    /// The path to the primitive cell
    CASM::fs::path m_prim_path;

    /// The miller indexes for the desired shift plane relative to the primitive structure
    Eigen::Vector3i m_millers;

    /// Which atom of the slab to bring down to the base of the plane.
    /// You need to specify this if you have different possible planes the shifts can occur in.
    int m_floor_slab_atom;

    /// Specifies how thick the slabs your shifting should be, i.e. how many times you want to
    /// stacl the shift unit
    int m_stacks;
};

class ShiftSettings
{
public:
private:
};

class FullSettings
{
public:
    FullSettings(const CASM::jsonParser& init_settings);

    const std::string& name() const { return m_name; }

    const SlabSettings& slab_settings() const { return m_slab_settings; }

private:
    /// Name for the particular set of shifts you're going to do
    std::string m_name;

    /// Specifies all the settings involved with creating a slab
    SlabSettings m_slab_settings;
};

///Returns smallest possible structure with the shift plane exposed along the a-b vectors, but
///does not translate the basis in any way
CASM::Structure raw_shift_unit(const CASM::Structure& prim, const Eigen::Vector3i& miller_indexes)
{
    auto shift_lattice=prim.lattice().get_lattice_in_plane(miller_indexes);
    CASM::Structure raw_shift_unit(shift_lattice);
    raw_shift_unit.fill_supercell(prim);
    return raw_shift_unit;
}

CASM::Structure shift_unit(const CASM::Structure&, const CASM::Coordinate& coord_to_floor)
{
}

} // namespace mush

///Throws an exception if a bad combination of program options was given
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

    CASM::Structure test("./POS");
    return 0;
}

int main()
{
    return 0;
}
