#include "./exceptions.hpp"
#include "cxxopts.hpp"
#include "casm/casm_io/jsonParser.hh"

namespace cxxopts
{
void required_argument_notify(const ParseResult& result, const std::vector<std::string>& required_arguments)
{

    for (const auto& arg : required_arguments)
    {
        if (!result.count(arg))
        {
            throw mush::except::RequiredArgumentMissing(arg);
        }
    }

    return;
}
} // namespace cxxopts

namespace lazy
{
}

namespace loggy
{
void divider() { std::cout << "..........................................................." << std::endl; }
} // namespace loggy

