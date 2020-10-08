#include "./chain.hpp"
#include "./common_options.hpp"
#include "./misc.hpp"
#include <casmutils/xtal/structure_tools.hpp>
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

    chain_sub->callback([=]() { run_subcommand_chain(*input_path_ptr, *output_path_ptr, *celavages_ptr, *grid_dims_ptr, std::cout); });
}

std::string make_cleave_dirname(double cleavage)
{
    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << cleavage;
    return "cleave__" + cleavestream.str();
}

std::string make_shift_dirname(int a, int b) { return "shift__" + std::to_string(a) + "." + std::to_string(b); }

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

mush::fs::path make_target_directory(const mush::MultiRecord& record)
{
    return mush::fs::path(make_shift_dirname(record.a_index, record.b_index)) / make_cleave_dirname(record.cleavage);
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

void run_subcommand_chain(const mush::fs::path& input_path,
                          const mush::fs::path& output_dir,
                          const std::vector<double>& cleavages,
                          const std::vector<int>& grid_dims,
                          std::ostream& log)
{
    mush::cautious_create_directory(output_dir);

    log << "Reading slab from " << input_path << "...\n";
    auto slab = cu::xtal::Structure::from_poscar(input_path);

    mush::json full_record;
    full_record["grid"] = grid_dims;
    full_record["cleavages"] = cleavages;

    log << "Shifting structures for " << grid_dims[0] << "x" << grid_dims[1] << " grid (please be patient)...\n";

    mush::Shifter shifter(slab, grid_dims[0], grid_dims[1]);
    assert(shifter.grid_dims[0] == grid_dims[0] && shifter.grid_dims[1] == grid_dims[1]);
    full_record["shift_units"]=make_shift_units(shifter);

    std::vector<std::vector<std::string>> unique_equivalent_groups;
    std::unordered_map<int,int> equivalence_map_ix_to_group_label;
    for (double cleave : cleavages)
    {
        std::unordered_set<int> recorded_equivalents;
        log << "Cleaving " << cleave << " angstroms...\n";
        for (int i = 0; i < shifter.size(); ++i)
        {
            auto cleaved_shifted_structures = mush::make_cleaved_structures(shifter.shifted_structures[i], cleavages);
            auto report = make_multirecord(cleave, shifter, i);

            if (recorded_equivalents.count(i) == 0)
            {
                for(int e : shifter.equivalence_map[i])
                {
                    equivalence_map_ix_to_group_label[e]=unique_equivalent_groups.size();
                }
                unique_equivalent_groups.push_back(report.equivalent_structures);
                recorded_equivalents.insert(shifter.equivalence_map[i].begin(), shifter.equivalence_map[i].end());
            }

            auto dir = make_target_directory(report);
            log << "Write structure to " << output_dir / dir << "...\n";
            mush::fs::create_directories(output_dir / dir);
            cu::xtal::write_poscar(shifter.shifted_structures[i], output_dir / dir / "POSCAR");

            auto chunk = serialize(report);
            chunk["directory"] = dir;
            chunk["shift"]=make_aligned_shift_vector(shifter, i);
            chunk["group"]=equivalence_map_ix_to_group_label[i];

            full_record["ids"][report.id()] = chunk;
        }
    }

    full_record["equivalents"] = unique_equivalent_groups;

    log << "Back up slab structure to " << output_dir / "slab.vasp"
        << "...\n";
    cu::xtal::write_poscar(slab, output_dir / "slab.vasp");
    log << "Save record to "<<output_dir/"record.json"<<"...\n";
    mush::write_json(full_record,output_dir/"record.json");
}

