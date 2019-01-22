#ifndef MULTIEXCEPT_HH
#define MULTIEXCEPT_HH

#include <stdexcept>
#include <string>

namespace mush
{
namespace except
{
class FileExists : public std::runtime_error
{
public:
    FileExists(std::string target_path) : std::runtime_error("File '" + target_path + "' already exists!") {}
};
} // namespace except
} // namespace mush

#endif
