#include "./slice.hpp"
#include "./common_options.hpp"
#include "./misc.hpp"
#include <casmutils/mush/shift.hpp>
#include <casmutils/mush/twist.hpp>
#include <casmutils/xtal/structure.hpp>
#include <casmutils/xtal/structure_tools.hpp>
#include <filesystem>
#include <memory>
#include <multishift/slice_settings.hpp>
#include <stdexcept>
#include <string>
#include <vector>


void setup_subcommand_slice(CLI::App& app)
{
    auto input_path_ptr = std::make_shared<mush::fs::path>();
    auto output_path_ptr = std::make_shared<mush::fs::path>();
    auto miller_indexes_ptr = std::make_shared<std::vector<int>>();
    auto align_ptr = std::make_shared<bool>(false);

    CLI::App* slice_sub =
        app.add_subcommand("slice", "Slice unit cell to expose desired plane. Use output to construct slabs of a desired thickness.");
    slice_sub
        ->add_option("-m,--millers",
                     *miller_indexes_ptr,
                     "Miller indexes define the plane of your input structure that will be exposed on the facet of the output.")
        ->required();
    slice_sub->add_flag(
        "-x,--dont-align", *align_ptr, "Prevent rigidnly rotating output structure so that the exposed plane is in the xy Cartesian plane.");

    populate_subcommand_input_option(slice_sub, input_path_ptr.get(), true);
    populate_subcommand_output_option(slice_sub, output_path_ptr.get(), true);

    slice_sub->callback([=]() { run_subcommand_slice(*input_path_ptr, *output_path_ptr, *miller_indexes_ptr, *align_ptr, std::cout); });
}

void run_subcommand_slice(
    const mush::fs::path& input_path, const mush::fs::path& output_path, const std::vector<int>& millers, bool align, std::ostream& log)
{
    log << "Reading " << input_path << "...\n";
    auto prim = mush::cu::xtal::Structure::from_poscar(input_path);

    if (millers.size() != 3)
    {
        throw std::runtime_error("Exactly 3 miller indexes are required, but received " + std::to_string(millers.size()) + ".");
    }

    log << "Slice along (" << millers[0]<<", "<<millers[1]<<", "<<millers[2]<<")...\n";
    auto sliced_prim=mush::cu::xtal::slice_along_plane(prim, Eigen::Vector3i(millers[0],millers[1],millers[2]));

    if(align)
    {
        log<<"Align exposed plane to xy plane...\n";
        mush::make_aligned(&sliced_prim);
    }

    log << "Write final structure to "<<output_path<<"...\n";
    mush::cu::xtal::write_poscar(sliced_prim,output_path);
}
