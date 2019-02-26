#include "multishift/fourier.hpp"
#include "multishift/io.hpp"
#include "casm/casm_io/json_io/container.hh"
#include "casmutils/structure.hpp"
#include "cxxopts.hpp"
#include "multishift/define.hpp"
#include "multishift/exceptions.hpp"
#include "multishift/misc.hpp"
#include <iostream>
#include <unordered_set>

/* mush::Interpolator::InterGrid best_reshaped(const std::vector<mush::InterPoint>& unrolled_data, double ka_dim, */
/*                                             double kb_dim) */
/* { */
/*     std::vector<std::pair<int, int>> kratio_candidates; */
/*     for (int i = 1; i < unrolled_data.size() + 1; ++i) */
/*     { */
/*         if (unrolled_data.size() % i == 0) */
/*         { */
/*             kratio_candidates.emplace_back(i, unrolled_data.size() / i); */
/*         } */
/*     } */

/*     auto best = kratio_candidates[0]; */
/*     for (const auto& candidate : kratio_candidates) */
/*     { */
/*     } */
/*     assert(false); */
/* } */


/**
 * multishift-fourier reads up a slab structure and a regular
 * grid of values on the ab-plane, and finds a fourier basis to
 * reproduce the data.
 * Coefficients for the Fourier interpolation
 * are then printed out.
 */

int main(int argc, char* argv[])
{
    using mush::fs::path;

    cxxopts::Options options("multishift-fourier", "Interpolate a periodic surface.");

    // clang-format off
    options.add_options()
        ("d,data","Path to the data file (json).",cxxopts::value<path>())
        ("l,lattice","Path to the structure whose ab-vectors the data is relative to (poscar).",cxxopts::value<path>())
        ("h,help","Print available options.");
    // clang-format on

    try
    {
        auto result = options.parse(argc, argv);

        if (result.count("help"))
        {
            std::cout << options.help() << std::endl;
            return 0;
        }

        else
        {
            cxxopts::required_argument_notify(result, std::vector<std::string>{"data"});
            cxxopts::required_argument_notify(result, std::vector<std::string>{"lattice"});
        }


        const auto& data_path=result["data"].as<path>();
        const auto& lat_path=result["lattice"].as<path>();
        //**************************************************************************//

        CASM::jsonParser test_json(data_path);

        //save me c++17
        std::vector<double> as,bs,vals;
        std::tie(as,bs,vals)=mush::MultiIO::unrolled_frac_grid_data_from_json(test_json);

        auto slab = Rewrap::Structure::from_poscar(lat_path);
        const CASM::Lattice& slab_lat = slab.lattice();
        mush::Interpolator ipolator(slab_lat, as, bs, vals);

        const auto& real_lat=ipolator.real_lattice();
        const auto& real_vals=ipolator.sampled_values();
        auto real_dump=mush::Interpolator::grid_to_json(real_lat,real_vals);
        real_dump.write(mush::fs::path("./real.json"));

        const auto& rcal_lat=ipolator.reciprocal_lattice();
        const auto& rcal_vals=ipolator.k_values();
        auto rcal_dump=mush::Interpolator::grid_to_json(rcal_lat,rcal_vals);
        rcal_dump.write(mush::fs::path("./rcal.json"));

        auto recover=ipolator.interpolate(50,50);
        auto ipol_dump=mush::Interpolator::grid_to_json(recover.first,recover.second);
        ipol_dump.write(mush::fs::path("./ipol.json"));

    }

    catch(const std::exception& e)
    {
        std::cout<<e.what()<<std::endl;
        throw e;
    }

    return 0;
}
