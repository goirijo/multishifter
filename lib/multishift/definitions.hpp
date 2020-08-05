#ifndef MUSH_DEFINITIONS_HH
#define MUSH_DEFINITIONS_HH

#include <casm/misc/CASM_math.hh>
#include <casm/misc/CASM_Eigen_math.hh>
#include <casmutils/definitions.hpp>
#include <nlohmann/json.hpp>

namespace casmutils
{
    namespace xtal
    {
        class Lattice;
    }
    namespace mush
    {
    }
}

namespace mush
{
    namespace cu=casmutils;
    using CASM::almost_equal;
    namespace fs=casmutils::fs;
    using nlohmann::json;

    using Lattice=cu::xtal::Lattice;

    using namespace cu::mush;

    //TODO
    /* double degrees_to_radians(double degrees); */
    /* double radians_to_degrees(double rad); */
}

#endif
