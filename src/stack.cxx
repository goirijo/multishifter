#include "./stack.hpp"
#include "./common_options.hpp"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "casmutils/xtal/frankenstein.hpp"
#include "casmutils/mush/slab.hpp"
#include "multishift/slice_settings.hpp"
#include <filesystem>
#include <memory>

namespace cu=casmutils;

void setup_subcommand_stack(CLI::App& app)
{
    auto input_paths_ptr=std::make_shared<std::vector<mush::fs::path>>();
    auto output_path_ptr=std::make_shared<mush::fs::path>();

    CLI::App* stack_sub=app.add_subcommand("stack", "Stack multiple structures along the c direction.");
    stack_sub->add_option("-i,--inputs",*input_paths_ptr,"Paths to the structure files to stack together.")->required();
    populate_subcommand_output_option(stack_sub,output_path_ptr.get());
    
    stack_sub->callback([input_paths_ptr,output_path_ptr](){run_subcommand_stack(*input_paths_ptr,*output_path_ptr,std::cout);});
}

void run_subcommand_stack(const std::vector<mush::fs::path>& input_paths, const mush::fs::path& output_path, std::ostream& log)
{
    log << "Loading structures ...\n";
    std::vector<mush::cu::xtal::Structure> strucs;
    for(const auto& p : input_paths)
    {
        auto struc_in=cu::xtal::Structure::from_poscar(p);
        strucs.emplace_back(mush::make_aligned(struc_in));
    }

    log << "Stacking structures ...\n";
    auto stacked=mush::orthogonalize_c_vector(mush::cu::frankenstein::stack(strucs));
    stacked.within();

    log << "Write to "+output_path.string()<<std::endl;
    mush::cu::xtal::write_poscar(stacked,output_path);
}
