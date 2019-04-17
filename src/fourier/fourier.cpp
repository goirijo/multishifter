#include "multishift/fourier.hpp"
#include "multishift/io.hpp"
#include "multishift/settings.hpp"
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

std::pair<int,int> parse_grid_argument(const std::string& dimstring)
{
    std::string adim,bdim;
    std::string* target=&adim;

    for(int i=0; i<dimstring.size(); ++i)
    {
        if(std::isdigit(dimstring[i]))
        {
            target->push_back(dimstring[i]);
        }
        else
        {
            target=&bdim;
        }
    }

    return std::make_pair(std::stoi(adim),std::stoi(bdim));
}


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
        ("s,settings","Path to the settings file (json).",cxxopts::value<path>())
        ("i,interpolate","Give a grid density to calculate values for (e.g. '13x23'). Also requires --interpolator and --output.",cxxopts::value<std::string>())
        ("p,interpolator","Path to the file with the interpolation data that was generated when using --settings.",cxxopts::value<path>())
        ("o,output","Path to the file where interpolation should be written to.",cxxopts::value<path>())
        ("a,analytic","Print an analytic expression for the surface in terms of sin/cos in the specified format (python-cart,python-frac,latex).",cxxopts::value<std::string>())
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

        else if(result.count("interpolate"))
        {
            cxxopts::required_option_notify(result, std::vector<std::string>{"interpolator","output"});
            auto dims=parse_grid_argument(result["interpolate"].as<std::string>());
            auto ipol_path=result["interpolator"].as<path>();
            auto target=result["output"].as<path>();

            CASM::jsonParser ipol_dump(ipol_path);
            auto ipolator=mush::Interpolator::deserialize(ipol_dump);

            auto ipol_datas=ipolator.interpolate(dims.first,dims.second);
            auto ipolated_dump=mush::Interpolator::grid_to_json(ipol_datas.first,ipol_datas.second);

            mush::io::doesnt_exist_or_throw(target);
            ipolated_dump.write(target);
        }

        else if(result.count("analytic"))
        {
            cxxopts::required_option_notify(result, std::vector<std::string>{"interpolator"});
            cxxopts::invalid_parameter_notify(result, "analytic", std::vector<std::string>{"python-cart"});

            auto ipol_path=result["interpolator"].as<path>();
            CASM::jsonParser ipol_dump(ipol_path);
            auto ipolator=mush::Interpolator::deserialize(ipol_dump);

            mush::Analytiker anal(ipolator.k_values());
        }

        else
        {
            cxxopts::required_option_notify(result, std::vector<std::string>{"settings"});

            const auto& settings_path=result["settings"].as<path>();

            auto all_settings = mush::FullSettings::from_path(settings_path);
            std::cout << "Using '" << all_settings.name() << "' as prefix." << std::endl;

            const auto& fourier_settings = all_settings.fourier_settings();
            std::cout << "Search in " << fourier_settings.data_path() << " for data to fit to." << std::endl;

            std::cout<<"Using "<<fourier_settings.lattice_path()<<" as reference lattice."<<std::endl;
            auto slab = Rewrap::Structure::from_poscar(fourier_settings.lattice_path());
            const CASM::Lattice& slab_lat = slab.lattice();

            CASM::jsonParser fourier_data(fourier_settings.data_path());
            //save me c++17
            std::vector<double> as,bs;
            std::vector<std::vector<double>> vals;
            std::tie(as,bs,vals)=mush::MultiIO::unrolled_frac_grid_data_from_json(fourier_data,fourier_settings.values_to_interpolate());

            std::vector<mush::Interpolator> ipolators;
            ipolators.reserve(vals.size());
            assert(vals.size()==fourier_settings.values_to_interpolate().size());
            for(int i=0; i<vals.size(); ++i)
            {
                std::cout<<"Creating interpolator for '"<<fourier_settings.values_to_interpolate()[i]<<"' data..."<<std::endl;
                ipolators.emplace_back(slab_lat, as, bs, vals[i]);
            }

            loggy::divider();

            mush::MultiIO writer(all_settings.name());

            std::cout<<"Created "<<vals.size()<<" interpolators. Saving them to:"<<std::endl;
            for(int i=0; i<vals.size(); ++i)
            {
                auto target=writer.drop_interpolator(ipolators[i],fourier_settings.values_to_interpolate()[i]);
                std::cout<<"        "<<target<<" (also saving raw real and reciprocal gridpoint data)"<<std::endl;
            }

            std::cout<<"Making backups of your settings..."<<std::endl;
            auto json_settings=all_settings.fourier_settings().to_json();
            json_settings.write(writer.target<mush::FourierSettings>()/(mush::FourierSettings::docs.tag()+".json"));
            std::cout<<"Saving data and lattice structure..."<<std::endl;
            fourier_data.write(writer.target<mush::FourierSettings>()/"data.json");
            Simplicity::write_poscar(slab,writer.target<mush::FourierSettings>()/"lattice.vasp");



            loggy::divider();
            std::cout<<"All fourier coefficients have been saved."<<std::endl;
            std::cout<<"Use the files printed above as the --interpolator argument when using different options."<<std::endl;

        }

    }

    catch (const cxxopts::missing_argument_exception& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Use '--help' to see option arguments" << std::endl;
        return 1;
    }

    catch (const mush::except::RequiredOptionMissing& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Use '--help' to see options" << std::endl;
        return 2;
    }

    catch (const mush::except::UnspecifiedSettings& e)
    {
        // TODO: Be more helpful
        std::cerr << e.what() << std::endl;
        return 3;
    }

    catch (const mush::except::FileExists& e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Rename it or delete it if you don't need it." << std::endl;
        return 4;
    }

    catch (const mush::except::InvalidParameter& e)
    {
        std::cerr << e.what() << std::endl;
        return 5;
    }


    return 0;
}
