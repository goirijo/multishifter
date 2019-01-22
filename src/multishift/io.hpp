#ifndef MULTIIO_HH
#define MULTIIO_HH

#include <string>

namespace mush
{

class MultiBase;
class BaseSettings;

/**
 * Given classes defined in the mush namespace, this class
 * can write out the necessary files to create gamma surfaces,
 * collect UBER data, etc.
 */

class MultiIO
{
public:
    MultiIO(const std::string& init_name);

    std::string name() const {return m_name;}

    /// Writes out preliminary structures needed before the shifts are done:
    /// primitive, shift unit, raw slab (no translation), slab
    /// Pass the associated settings along with the structures
    void drop_base(const BaseSettings& settings, const MultiBase& preshift_structures);

private:
    /// Name of the calculations directories
    std::string m_name;
};

} // namespace mush

#endif
