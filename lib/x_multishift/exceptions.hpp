#ifndef MULTIEXCEPT_HH
#define MULTIEXCEPT_HH

#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    std::string commatize(const std::vector<std::string>& words)
    {
        std::string commatized;
        for(const auto& w : words)
        {
            commatized+="'";
            commatized+=w;
            commatized+="', ";
        }

        commatized.pop_back();
        commatized.pop_back();

        return commatized;
    }
}

namespace mush
{
namespace except
{
class UnspecifiedSettings : public std::runtime_error
{
public:
    UnspecifiedSettings(const std::string& entry)
        : std::runtime_error("The settings entry '" + entry + "' was not specified.")
    {
    }
};
class RequiredOptionMissing : public std::runtime_error
{
public:
    RequiredOptionMissing(const std::string& option)
        : std::runtime_error("Expected option '--" + option + "' from command line.")
    {
    }
};

class InvalidParameter : public std::runtime_error
{
public:
    InvalidParameter(const std::string& option, const std::string& parameter, const std::vector<std::string>& available_parameters)
        : std::runtime_error("Bad parameter '"+parameter+"' for option '--"+option+"'. Available parameters are: "+commatize(available_parameters)+".")
    {
    }
};

class FileExists : public std::runtime_error
{
public:
    FileExists(const std::string& target_path) : std::runtime_error("File '" + target_path + "' already exists!") {}
};

class LeftHandedLattice : public std::runtime_error
{
public:
    LeftHandedLattice()
        : std::runtime_error("The lattice of a structure must be right handed. The lattice vectors must be redefined. ")
    {
    }
};

class DimensionalMismatch : public std::runtime_error
{
public:
    DimensionalMismatch(int expected, int got, const std::string& elaborate)
        : std::runtime_error("Expected dimension to be " + std::to_string(expected) + " but got " +
                             std::to_string(got) + " instead.\n" + elaborate)
    {
    }
};

class BadData : public std::runtime_error
{
public:
    BadData(const std::string& elaborate)
        : std::runtime_error("Bad input data.\n" + elaborate)
    {
    }
};

class BadSetting : public std::runtime_error
{
public:
    BadSetting(const std::string& key, const std::string& elaborate)
        : std::runtime_error("Bad settigs parameter for '" + key + "'\n" + elaborate)
    {
    }
};

class SettingMustBeGreaterThanZero : public BadSetting
{
public:
    SettingMustBeGreaterThanZero(const std::string& key) : BadSetting(key, "Value must be greater than zero") {}
};

class SettingMustBeGreaterEqualThanZero : public BadSetting
{
public:
    SettingMustBeGreaterEqualThanZero(const std::string& key)
        : BadSetting(key, "Value must be greater or equal to zero")
    {
    }
};

} // namespace except
} // namespace mush

#endif
