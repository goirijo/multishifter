#ifndef DEFINITIONS_HH
#define DEFINITIONS_HH

#include "casm/CASM_global_definitions.hh"

namespace Rewrap
{
    class Structure;
}

namespace CASM
{
    class Lattice;
}

namespace mush
{
    typedef Rewrap::Structure Structure;
    typedef CASM::Lattice Lattice;          //Maybe this will be Rewrap one day
    namespace fs=CASM::fs;
}

#endif
