#include "./mutate.hpp"
#include "./common_options.hpp"
#include "casmutils/xtal/lattice.hpp"
#include "casmutils/xtal/site.hpp"
#include "casmutils/xtal/structure.hpp"
#include <casmutils/mush/shift.hpp>
#include <memory>
#include <vector>
#include <casmutils/xtal/structure_tools.hpp>
#include <casmutils/xtal/coordinate.hpp>

void setup_subcommand_mutate(CLI::App& app)
{
    auto input_path_ptr=std::make_shared<mush::fs::path>();
    auto output_path_ptr=std::make_shared<mush::fs::path>();
    auto mutation_ptr=std::make_shared<std::vector<double>>();
    auto frac_ptr = std::make_shared<bool>(false);

    CLI::App* mutate_sub=app.add_subcommand("mutate", "Break periodicity by altering the c vector to create a shiftor cleave.");
    populate_subcommand_input_option(mutate_sub,input_path_ptr.get(),true);
    populate_subcommand_output_option(mutate_sub,output_path_ptr.get(),true);
    auto opt_m=mutate_sub->add_option("-m,--mutation",*mutation_ptr,"Value to add to the c lattice vector.")->required()->expected(3);

    populate_subcommand_fractional(mutate_sub, frac_ptr.get(), opt_m);
    
    mutate_sub->callback([=](){run_subcommand_mutate(*input_path_ptr,*output_path_ptr,*mutation_ptr,*frac_ptr,std::cout);});
}

void run_subcommand_mutate(const mush::fs::path& input_path, const mush::fs::path& output_path, const std::vector<double>& mutation, bool frac, std::ostream& log)
{
    /* log << "Loading structure ...\n"; */
    /* auto struc=mush::cu::xtal::Structure::from_poscar(input_path); */

    /* Eigen::Vector3d vec(mutation[0],mutation[1],mutation[2]); */

    /* log << "Using coordinate mode "<<coord_mode<<" ...\n"; */
    /* if(coord_mode=="frac") */
    /* { */
    /*     vec=mush::cu::xtal::Coordinate::from_fractional(vec,struc.lattice()).cart(); */
    /* } */

    /* struc=mush::mutate(struc,vec); */

    /* /1* Eigen::Vector3d new_c_vector=struc.lattice().c()+vec; *1/ */
    /* /1* mush::cu::xtal::Lattice new_lattice(struc.lattice().a(),struc.lattice().b(),new_c_vector); *1/ */
    /* /1* struc.set_lattice(new_lattice,mush::cu::xtal::CART); *1/ */

    /* log << "Write to "+output_path.string()<<std::endl; */
    /* mush::cu::xtal::write_poscar(struc,output_path); */
}
