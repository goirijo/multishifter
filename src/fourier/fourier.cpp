#include "multishift/fourier.hpp"
#include "multishift/define.hpp"
#include "multishift/misc.hpp"
#include "casm/casm_io/json_io/container.hh"
#include "cxxopts.hpp"
#include "multishift/define.hpp"
#include "multishift/exceptions.hpp"
#include <iostream>
#include <unordered_set>
#include "casmutils/structure.hpp"

struct pair_hash {
    inline std::size_t operator()(const std::pair<long,long> & v) const {
        return v.first*31+v.second;
    }
};

void sanity_check(const std::vector<mush::InterPoint>& unrolled_data)
{
    long prec=1e8;
    std::unordered_set<std::pair<long,long>,pair_hash> unique_values;
    for(const auto& p : unrolled_data)
    {
        unique_values.insert(std::make_pair(static_cast<long>(p.a_frac*prec+0.5),static_cast<long>(p.b_frac*prec+0.5)));
    }

    if(unique_values.size()!=unrolled_data.size())
    {
        throw mush::except::BadData("There are multiple values defined per grid point!");
    }

    return;
}

mush::Interpolator::InterGrid direct_reshape(const std::vector<mush::InterPoint>& unrolled_data, int ka_dim, int kb_dim)
{
    if(ka_dim*kb_dim!=unrolled_data.size())
    {
        throw mush::except::DimensionalMismatch(ka_dim*kb_dim,unrolled_data.size(),"Cannot reshape vector to the specified grid.");
    }

    int i=0;
    mush::Interpolator::InterGrid final_grid;
    for(int ka=0; ka<ka_dim; ++ka)
    {
        mush::Interpolator::InterGrid::value_type kb_row;
        for(int kb=0; kb<kb_dim; ++kb, ++i)
        {
            kb_row.push_back(unrolled_data[i]);
        }
        final_grid.emplace_back(std::move(kb_row));
    }

    return final_grid;
}

mush::Interpolator::InterGrid best_reshaped(const std::vector<mush::InterPoint>& unrolled_data, double ka_dim, double kb_dim)
{
    std::vector<std::pair<int,int>> kratio_candidates;
    for(int i=1; i<unrolled_data.size()+1; ++i)
    {
        if(unrolled_data.size()%i==0)
        {
            kratio_candidates.emplace_back(i,unrolled_data.size()/i);
        }
    }

    auto best=kratio_candidates[0];
    for(const auto& candidate : kratio_candidates)
    {
    }
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
    mush::fs::path test_path("./dumb.json");
    CASM::jsonParser test_json(test_path);

    // I hate jsonParser
    std::vector<double> as;
    CASM::from_json(as, test_json["a_frac"]);
    /* auto as=test_json.get<std::vector<double>>(); */
    std::vector<double> bs;
    CASM::from_json(bs, test_json["b_frac"]);
    std::vector<double> vals;
    CASM::from_json(vals, test_json["value"]);

    if (as.size() != bs.size())
    {
        throw mush::except::DimensionalMismatch(
            as.size(), bs.size(), "Cannot create a mesh grid with the provided data points along a and b.");
    }
    if (as.size() != vals.size())
    {
        throw mush::except::DimensionalMismatch(
            as.size(), vals.size(), "The number of grid points does not match the number of values to fit to.");
    }
    
    // You have to ensure there's only a single data point per grid point
    std::vector<mush::InterPoint> unrolled_data;
    for(int i=0; i<vals.size(); ++i)
    {
        unrolled_data.emplace_back(as[i],bs[i],vals[i]);
    }

    sanity_check(unrolled_data);


    for (auto a : as)
    {
        std::cout << a << " ";
    }
    std::cout << std::endl << std::endl;
    ;
    for (auto b : bs)
    {
        std::cout << b << " ";
    }
    std::cout << std::endl;
    loggy::divider();

    auto slab=Rewrap::Structure::from_poscar("./slab.vasp");
    const CASM::Lattice& slab_lat=slab.lattice();
    auto reshaped=direct_reshape(unrolled_data,25,25);

    mush::Interpolator(slab_lat,reshaped);




    return 0;
}
