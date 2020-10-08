#include "./align.hpp"
#include "./common_options.hpp"
#include "casmutils/xtal/structure.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include <multishift/definitions.hpp>
#include <casmutils/mush/twist.hpp>
#include <casmutils/mush/slab.hpp>
#include <memory>

namespace cu=casmutils;

void setup_subcommand_align(CLI::App& app)
{
    auto input_path_ptr=std::make_shared<mush::fs::path>();
    auto output_path_ptr=std::make_shared<mush::fs::path>();
    auto prismatic_ptr=std::make_shared<bool>(false);

    CLI::App* align_sub = app.add_subcommand("align", "Reorient the structure so that the a and b vectors lie on the xy-plane.");

    populate_subcommand_input_option(align_sub,input_path_ptr.get());
    populate_subcommand_output_option(align_sub,output_path_ptr.get());
    align_sub->add_flag("-p,--prismatic",*prismatic_ptr,"Force the c vector to be perpendicular to the ab-plane (may break periodicity)");

    align_sub->callback([input_path_ptr,output_path_ptr,prismatic_ptr]() { run_subcommand_align(*input_path_ptr,*output_path_ptr,*prismatic_ptr,std::cout); });
}

void run_subcommand_align(const mush::fs::path& input_path, const mush::fs::path& output_path, bool prismatic, std::ostream& log)
{
    auto struc=mush::cu::xtal::Structure::from_poscar(input_path);
    mush::make_aligned(&struc);

    if(prismatic)
    {
        struc.set_lattice(mush::make_prismatic_lattice(struc.lattice()),mush::cu::xtal::CART);
    }

    mush::cu::xtal::write_poscar(struc,output_path);
    return;
}
