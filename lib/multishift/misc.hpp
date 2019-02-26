#ifndef MULTISHIFTMISC_HH
#define MULTISHIFTMISC_HH

#include "./exceptions.hpp"
#include "casm/casm_io/jsonParser.hh"
#include "cxxopts.hpp"

namespace cxxopts
{
// TODO: split into cxx someday
void required_argument_notify(const ParseResult& result, const std::vector<std::string>& required_arguments);
} // namespace cxxopts

namespace lazy
{
template <typename T>
T get_or_value(const CASM::jsonParser& json, const std::string& key, const T& default_value)
{
    if (json.contains(key))
    {
        return json[key].get<T>();
    }

    return default_value;
}
} // namespace lazy

namespace loggy
{
void divider();
} // namespace loggy

namespace hashing
{
template <typename T1, typename T2>
struct pair_hash
{
    std::size_t operator()(const std::pair<T1, T2>& v) const { return v.first * 31 + v.second; }
};

} // namespace hashing

#endif
