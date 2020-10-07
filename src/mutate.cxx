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
    log << "Loading structure...\n";
    auto struc=mush::cu::xtal::Structure::from_poscar(input_path);

    Eigen::Vector3d vec_cart(mutation[0],mutation[1],mutation[2]);
    Eigen::Vector3d vec_frac=mush::cu::xtal::Coordinate::from_fractional(vec_cart,struc.lattice()).cart();

    if(frac)
    {
        vec_frac=Eigen::Vector3d(mutation[0],mutation[1],mutation[2]);
        vec_cart=mush::cu::xtal::Coordinate::from_fractional(vec_frac,struc.lattice()).cart();
    }

    else
    {
        vec_cart=Eigen::Vector3d(mutation[0],mutation[1],mutation[2]);
        vec_frac=mush::cu::xtal::Coordinate(vec_cart).frac(struc.lattice());
    }

    log << "Mutate c-vector by "<<vec_cart.transpose()<<" (Cartesian) or "<<vec_frac.transpose()<<" (fractional)...\n";
    struc=mush::mutate(struc,vec_cart);

    log << "Write to "+output_path.string()<<"...\n";
    mush::cu::xtal::write_poscar(struc,output_path);
}
