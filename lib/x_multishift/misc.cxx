#include "./exceptions.hpp"
#include "cxxopts.hpp"
#include "casm/casm_io/jsonParser.hh"

namespace cxxopts
{
void required_option_notify(const ParseResult& result, const std::vector<std::string>& required_options)
{

    for (const auto& opt : required_options)
    {
        if (!result.count(opt))
        {
            throw mush::except::RequiredOptionMissing(opt);
        }
    }

    return;
}

void invalid_parameter_notify(const ParseResult& result, const std::string& option, const std::vector<std::string>& available_parameters)
{
    auto parameter=result[option].as<std::string>();
    for(const auto& par : available_parameters)
    {
        if(parameter==par)
        {
            return;
        }
    }

    throw mush::except::InvalidParameter(option,parameter,available_parameters);
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

