#ifndef MULTISHIFTMISC_HH
#define MULTISHIFTMISC_HH

#include "./exceptions.hpp"
#include "cxxopts.hpp"

namespace cxxopts
{
// TODO: split into cxx someday
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

namespace loggy
{
void divider() { std::cout << "..........................................................." << std::endl; }
} // namespace loggy

#endif
