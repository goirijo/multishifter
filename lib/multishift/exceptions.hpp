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
        : std::runtime_error("The settings entry '" + entry+"' was not specified.")
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
} // namespace except
} // namespace mush

#endif
