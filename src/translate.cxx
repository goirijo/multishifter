#include "./translate.hpp"
#include "./common_options.hpp"
#include "casmutils/xtal/structure_tools.hpp"
#include "misc.hpp"
#include <casmutils/mush/shift.hpp>
#include <casmutils/mush/slab.hpp>
#include <casmutils/xtal/frankenstein.hpp>
#include <casmutils/xtal/coordinate.hpp>
#include <filesystem>
#include <memory>
#include <multishift/slice_settings.hpp>
#include <stdexcept>
#include <string>

void setup_subcommand_translate(CLI::App& app)
{
    auto input_path_ptr = std::make_shared<mush::fs::path>();
    auto output_path_ptr = std::make_shared<mush::fs::path>();
    auto shift_ptr = std::make_shared<std::vector<double>>(3, 0.0);
    auto floor_ix_ptr = std::make_shared<int>(0);
    auto frac_ptr = std::make_shared<bool>(false);

    CLI::App* translate_sub =
        app.add_subcommand("translate", "Rigidly translate the entire basis of a structure to change the atomic layer at the interface.");

    auto v_opt=translate_sub->add_option("-v,--value", *shift_ptr, "Vector specifying by which amount to translate the basis")->expected(3);
    
    translate_sub
        ->add_option("-f,--floor",
                     *floor_ix_ptr,
                     "Index of basis atom that should end up at the origin on the unit cell after translating (indexing begins at 1).")
        ->excludes(v_opt);

    populate_subcommand_input_option(translate_sub, input_path_ptr.get());
    populate_subcommand_output_option(translate_sub, output_path_ptr.get());
    populate_subcommand_fractional(translate_sub, frac_ptr.get(), v_opt);

    translate_sub->callback(
        [=]() { run_subcommand_translate(*input_path_ptr, *output_path_ptr, *shift_ptr, *floor_ix_ptr, *frac_ptr, std::cout); });
}

void run_subcommand_translate(const mush::fs::path& input_path,
                              const mush::fs::path& output_path,
                              const std::vector<double>& _shift,
                              int floor_ix,
                              bool frac,
                              std::ostream& log)
{
    log<<"Reading "<<input_path<<"...\n";
    auto struc = mush::cu::xtal::Structure::from_poscar(input_path);
    Eigen::Vector3d shift(_shift[0],_shift[1],_shift[2]);

    if (floor_ix != 0)
    {
        int num_sites = struc.basis_sites().size();
        if (floor_ix > num_sites || floor_ix < 0)
        {
            throw std::runtime_error("Floor index out of range! " + std::to_string(num_sites) + " in basis but received index " +
                                     std::to_string(floor_ix));
        }

        log<<"Determining shift value for atom "<<floor_ix<<"...\n";
        shift=struc.basis_sites()[floor_ix-1].cart();
    }

    Eigen::Vector3d frac_shift;

    if(frac)
    {
        frac_shift=shift;
        Eigen::Vector3d shift=mush::cu::xtal::fractional_to_cartesian(shift,struc.lattice());
    }
    else
    {
        frac_shift=mush::cu::xtal::cartesian_to_fractional(shift,struc.lattice());
    }


    log<<"Translating basis by "<<shift.transpose()<<" (Cartesian) or "<<frac_shift.transpose()<<" (fractional)...\n";
    auto translated_struc=mush::cu::frankenstein::translate_basis(struc,shift);
    translated_struc.within();

    log<<"Write final structure to "<<output_path<<"...\n";
    mush::cu::xtal::write_poscar(translated_struc,output_path);

    return;
}
