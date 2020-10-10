#ifndef CHAIN_SUBCOMMAND_HH
#define CHAIN_SUBCOMMAND_HH

#include <CLI/CLI.hpp>
#include <multishift/definitions.hpp>
#include "./misc.hpp"
#include "multishift/shifter.hpp"
#include <casmutils/xtal/structure_tools.hpp>

void setup_subcommand_chain(CLI::App& app);

namespace cu = casmutils;

mush::MultiRecord make_multirecord(const double cleave, const mush::Shifter& shifter, int ix);
mush::json serialize(const mush::MultiRecord& record);
std::array<double,2> make_aligned_shift_vector(const mush::Shifter& shifter, int ix);
std::array<std::array<double,2>,2> make_shift_units(const mush::Shifter& shifter);

//Used for cleave, shift,and chain subcommands. The only difference between them
//is the output directory layout (single layer vs two layers)
template<mush::SUBCOMMAND subcommand>
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

    if(subcommand!=mush::SUBCOMMAND::CLEAVE)
    {
        log << "Shifting structures for " << grid_dims[0] << "x" << grid_dims[1] << " grid (please be patient)...\n";
    }

    mush::Shifter shifter(slab, grid_dims[0], grid_dims[1]);
    assert(shifter.grid_dims[0] == grid_dims[0] && shifter.grid_dims[1] == grid_dims[1]);
    full_record["shift_units"]=make_shift_units(shifter);

    std::vector<std::vector<std::string>> unique_equivalent_groups;
    std::unordered_map<int,int> equivalence_map_ix_to_group_label;
    for (double cleave : cleavages)
    {
        std::unordered_set<int> recorded_equivalents;
        if(subcommand!=mush::SUBCOMMAND::SHIFT)
        {
            log << "Cleaving " << cleave << " angstroms...\n";
        }
        for (int i = 0; i < shifter.size(); ++i)
        {
            auto cleaved_shifted_structure = mush::make_cleaved_structure(shifter.shifted_structures[i], cleave);
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

            auto dir = mush::make_target_directory<subcommand>(report);
            auto target_file=output_dir/dir/"POSCAR";
            log << "Write structure to " << target_file << "...\n";
            mush::fs::create_directories(output_dir / dir);
            cu::xtal::write_poscar(cleaved_shifted_structure, target_file);

            auto chunk = serialize(report);
            chunk["directory"] = dir;
            chunk["shift"]=make_aligned_shift_vector(shifter, i);
            chunk["orbit"]=equivalence_map_ix_to_group_label[i];

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

#endif
