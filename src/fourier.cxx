#include "./fourier.hpp"
#include "./misc.hpp"
#include <casmutils/mush/slab.hpp>
#include <filesystem>
#include <memory>
#include <multishift/fourier.hpp>
#include <multishift/slice_settings.hpp>
#include <ostream>
#include <tuple>
#include <vector>

namespace cu = casmutils;

namespace
{
std::vector<mush::InterPoint> read_unrolled_data(const mush::json record, double cleavage_slice, const std::string& value_key)
{
    int amax=record["grid"][0];
    int bmax=record["grid"][1];

    std::vector<mush::InterPoint> unrolled_data;
    for (const auto& id : record["ids"])
    {
        if(!cu::almost_equal(cleavage_slice,id["cleavage"].get<double>(),1e-8))
        {
            continue;
        }

        int a=id["grid_point"][0];
        int b=id["grid_point"][1];
        double v=id[value_key];
        unrolled_data.emplace_back(static_cast<double>(a)/amax,static_cast<double>(b)/bmax,v);
    }

    return unrolled_data;
}

cu::xtal::Lattice extract_pseudo_slab_lattice(const mush::json& record) {
    auto su=record["shift_units"];
    Eigen::Vector3d au(su[0][0],su[0][1],0.0);
    Eigen::Vector3d bu(su[1][0],su[1][1],0.0);
    Eigen::Vector3d c(0,0,1);

    int amax=record["grid"][0];
    int bmax=record["grid"][1];

    return cu::xtal::Lattice(au*amax,bu*bmax,c);
}

} // namespace

//*************************************************************************************//

void setup_subcommand_fourier(CLI::App& app)
{
    auto data_path_ptr = std::make_shared<mush::fs::path>();
    auto cleavage_slice_ptr = std::make_shared<double>();
    auto entry_key_ptr = std::make_shared<std::string>();
    auto crush_ptr = std::make_shared<double>();

    CLI::App* fourier_sub = app.add_subcommand("fourier", "Perform Fourier decomposition and get analytical expression for data set.");
    fourier_sub
        ->add_option("-d,--data",
                     *data_path_ptr,
                     "Amended 'record.json' like file, with an additional entry for values to interpolate for each structure id.")
        ->required();
    fourier_sub->add_option("-c,--cleavage-slice", *cleavage_slice_ptr, "Take data from these cleavage values.")->default_val(0.0);
    fourier_sub->add_option("-k,--key", *entry_key_ptr, "Key of the value that is being interpolated.")->required();
    fourier_sub->add_option("-x,--crush", *crush_ptr, "Basis functions that fall within this threshold will get added together, reducing the total number of basis functions.")->default_val(0.0)->default_val(1e-9);

    fourier_sub->callback([=]() { run_subcommand_fourier(*data_path_ptr, *cleavage_slice_ptr, *entry_key_ptr, *crush_ptr, std::cout); });
}

void run_subcommand_fourier(const mush::fs::path& data_path, double cleavage_slice, const std::string& value_key, double crush_value, std::ostream& log)
{
    log << "Load data from "<<data_path<<"...\n";
    auto record = mush::load_json(data_path);

    const auto& slab_lattice = ::extract_pseudo_slab_lattice(record);
    /* log << "Inferred surface vectors as:\n"; */
    /* log << "    a: "<<slab_lattice.a().transpose()<<"\n"; */
    /* log << "    b: "<<slab_lattice.b().transpose()<<"\n"; */

    log << "Extract values at cleavage "<<std::fixed << std::setprecision(6)<<cleavage_slice<<"...\n";
    auto unrolled_data = ::read_unrolled_data(record, cleavage_slice, value_key);
    mush::Interpolator ipolator(slab_lattice, unrolled_data);

    const auto& start_values = ipolator.sampled_values();
    auto [lat, ipolvalues] = ipolator.interpolate(start_values.size(), start_values[0].size());

    mush::Analytiker analyzer(ipolator);

    log << "Crushing functions smaller than " << std::fixed << std::setprecision(9)<<crush_value << "...\n";
    auto [real_functions, imag_functions] = analyzer.python_cart("x", "y", "np", crush_value);

    log << "Real functions:\n";
    log << real_functions;
    log << "\n\n\n";
    log << "Imaginary functions:\n";
    log << imag_functions;
    log << "\n";
    log << "Surface lattice vectors:\n";
    log << "    a: " << ipolator.real_lattice().a()(0) << ", " << ipolator.real_lattice().a()(1) << "\n";
    log << "    b: " << ipolator.real_lattice().b()(0) << ", " << ipolator.real_lattice().b()(1) << "\n";

    return;
}
