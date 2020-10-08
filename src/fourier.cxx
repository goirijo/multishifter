#include "./fourier.hpp"
#include "./misc.hpp"
#include <filesystem>
#include <memory>
#include <multishift/fourier.hpp>
#include <casmutils/mush/slab.hpp>
#include <multishift/slice_settings.hpp>
#include <ostream>
#include <tuple>
#include <vector>

namespace
{
std::vector<mush::InterPoint> read_unrolled_data(const mush::fs::path& data_path)
{
    auto data_dump = mush::load_json(data_path);

    std::vector<double> a_frac = data_dump["a_frac"];
    std::vector<double> b_frac = data_dump["b_frac"];
    std::vector<double> values = data_dump["values"];

    if (a_frac.size() != b_frac.size() || a_frac.size() != values.size())
    {
        throw std::runtime_error("Data was not properly formatted");
    }

    std::vector<mush::InterPoint> unrolled_data;
    for (int i = 0; i < a_frac.size(); ++i)
    {
        unrolled_data.emplace_back(a_frac[i], b_frac[i], values[i]);
    }

    return unrolled_data;
}
} // namespace

//*************************************************************************************//

void setup_subcommand_fourier(CLI::App& app)
{
    auto settings_path_ptr = std::make_shared<mush::fs::path>();

    CLI::App* fourier_sub = app.add_subcommand("fourier", "Perform Fourier decomposition and get analytical expression for data set.");
    fourier_sub->add_option("-s,--settings", *settings_path_ptr, "Settings file with slab, resolution, and lattice periodic data.");

    fourier_sub->callback([settings_path_ptr]() { run_subcommand_fourier(*settings_path_ptr, std::cout); });
}

void run_subcommand_fourier(const mush::fs::path& settings_path, std::ostream& log)
{
    auto settings = mush::load_json(settings_path);

    auto slaber_settings = mush::SlabSettings::from_json(settings);
    auto fourier_settings = mush::FourierSettings::from_json(settings);

    log << "Reading surface lattice vectors from " << slaber_settings.slab_unit_path << "...\n";
    auto slab_unit = mush::cu::xtal::Structure::from_poscar(slaber_settings.slab_unit_path);
    const auto& slab_lattice = slab_unit.lattice();

    log << "Reading data from " << fourier_settings.data_path << "...\n";
    std::vector<mush::InterPoint> unrolled_data = ::read_unrolled_data(fourier_settings.data_path);
    mush::Interpolator ipolator(slab_lattice, unrolled_data);

    const auto& start_values = ipolator.sampled_values();
    auto [lat, ipolvalues] = ipolator.interpolate(start_values.size(), start_values[0].size());

    mush::Analytiker analyzer(ipolator);

    log << "Crushing functions to a resolution of " << fourier_settings.resolution << "...\n";
    auto [real_functions, imag_functions] = analyzer.python_cart("x", "y", "np", fourier_settings.resolution);

    log << "Real functions:\n";
    log << real_functions;
    log << "\n\n\n";
    log << "Imaginary functions:\n";
    log << imag_functions;
    log << "\n";
    log << "Surface lattice vectors:\n";
    log << "a: " << ipolator.real_lattice().a()(0) << ", " << ipolator.real_lattice().a()(1) << "\n";
    log << "b: " << ipolator.real_lattice().b()(0) << ", " << ipolator.real_lattice().b()(1) << "\n";

    return;
}
