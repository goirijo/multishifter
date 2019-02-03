#ifndef MULTISHIFTMISC_HH
#define MULTISHIFTMISC_HH

#include "./exceptions.hpp"
#include "cxxopts.hpp"
#include "casm/casm_io/jsonParser.hh"

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

namespace lazy
{
    template<typename T>
    T get_or_value(const jsonParser& json, const std::string& key, const T& default_value)
    {
        if(json.contains(key))
        {
            return json[key].get<T>();
        }

        return default_value;
    }
}

namespace loggy
{
void divider() { std::cout << "..........................................................." << std::endl; }
} // namespace loggy

#endif
