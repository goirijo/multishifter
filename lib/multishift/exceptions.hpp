#ifndef MULTIEXCEPT_HH
#define MULTIEXCEPT_HH

#include <stdexcept>
#include <string>

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
class RequiredArgumentMissing : public std::runtime_error
{
public:
    RequiredArgumentMissing(const std::string& argument)
        : std::runtime_error("Expected argument '" + argument + "' from command line.")
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
