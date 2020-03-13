#ifndef MULTISHIFTMISC_HH
#define MULTISHIFTMISC_HH

#include "./exceptions.hpp"
#include "casm/casm_io/jsonParser.hh"
#include "cxxopts.hpp"

namespace cxxopts
{
void required_option_notify(const ParseResult& result, const std::vector<std::string>& required_options);
void invalid_parameter_notify(const ParseResult& result, const std::string& option, const std::vector<std::string>& available_parameters);
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

template <typename T>
bool almost_equal(T lhs, T rhs, T tol=1e-12)
{
    return std::abs(lhs-rhs)<tol;
}

template <typename T>
bool almost_zero(T val, T tol=1e-12)
{
    return almost_equal(val, 0.0, tol);
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
