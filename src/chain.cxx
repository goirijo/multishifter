#include "./chain.hpp"
#include "./common_options.hpp"
#include "./misc.hpp"
#include <multishift/shifter.hpp>
#include <nlohmann/json.hpp>

namespace cu = casmutils;

void setup_subcommand_chain(CLI::App& app)
{
    auto input_path_ptr = std::make_shared<mush::fs::path>();
    auto output_path_ptr = std::make_shared<mush::fs::path>();
    auto celavages_ptr = std::make_shared<std::vector<double>>();
    auto grid_dims_ptr = std::make_shared<std::vector<int>>();

    CLI::App* chain_sub = app.add_subcommand("chain", "Combine cleave and shift commands for gamma surface calculations.");

    populate_subcommand_input_option(chain_sub, input_path_ptr.get());
    populate_subcommand_output_option(chain_sub, output_path_ptr.get());

    chain_sub->add_option("-c,--cleave", *celavages_ptr, "List of cleavage values to insert between slabs.")->required();
    chain_sub
        ->add_option("-s,--shift",
                     *grid_dims_ptr,
                     "Grid dimensions to divide the ab-lane of the slab into. Each grid point will correspond to a particular shift "
                     "vector. Periodic image not counted in grid.")
        ->expected(2)
        ->required();

    chain_sub->callback([=]() { run_subcommand_chain<mush::SUBCOMMAND::CHAIN>(*input_path_ptr, *output_path_ptr, *celavages_ptr, *grid_dims_ptr, std::cout); });
}

mush::MultiRecord make_multirecord(const double cleave, const mush::Shifter& shifter, int ix)
{
    auto make_partial = [](double _cleave, const mush::Shifter& _shifter, int _ix) {
        mush::MultiRecord record;
        record.cleavage = _cleave;
        record.a_index = _shifter.shift_records[_ix].a;
        record.b_index = _shifter.shift_records[_ix].b;
        return record;
    };

    auto record = make_partial(cleave, shifter, ix);
    for (int eqx : shifter.equivalence_map[ix])
    {
        record.equivalent_structures.emplace_back(make_partial(cleave, shifter, eqx).id());
    }

    return record;
}

template<>
mush::fs::path mush::make_target_directory<mush::SUBCOMMAND::CHAIN>(const mush::MultiRecord& record)
{
    return mush::fs::path(mush::make_shift_dirname(record.a_index, record.b_index)) / mush::make_cleave_dirname(record.cleavage);
}

mush::json serialize(const mush::MultiRecord& record)
{
    mush::json j;

    j["grid_point"] = std::vector<int>{record.a_index, record.b_index};
    j["cleavage"] = record.cleavage;
    j["equivalent_structures"] = record.equivalent_structures;

    return j;
}

std::array<double,2> make_aligned_shift_vector(const mush::Shifter& shifter, int ix)
{
    auto aligned_lat=mush::make_aligned(shifter.shifted_structures[ix].lattice());
    Eigen::Vector3d a_shift=static_cast<double>(shifter.shift_records[ix].a)/shifter.grid_dims[0]*aligned_lat.a();
    Eigen::Vector3d b_shift=static_cast<double>(shifter.shift_records[ix].b)/shifter.grid_dims[1]*aligned_lat.b();
    Eigen::Vector3d aligned_shift=a_shift+b_shift;

    assert(cu::almost_equal(aligned_shift(2),0.0));

    return {aligned_shift(0),aligned_shift(1)};
}

std::array<std::array<double,2>,2> make_shift_units(const mush::Shifter& shifter)
{
    auto aligned_lat=mush::make_aligned(shifter.shifted_structures[0].lattice());
    Eigen::Vector3d a_shift=1.0/shifter.grid_dims[0]*aligned_lat.a();
    Eigen::Vector3d b_shift=1.0/shifter.grid_dims[1]*aligned_lat.b();

    std::array<double,2> a{a_shift(0),a_shift(1)};
    std::array<double,2> b{b_shift(0),b_shift(1)};

    return {a,b};
}


