#ifndef MULTIIO_HH
#define MULTIIO_HH

#include "./define.hpp"
#include <string>

namespace mush
{


class MultiBase;
class MultiShift;
struct SurfacePoint;

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

    /// Target directory for the base files
    fs::path base_target() const;

    /// Writes out preliminary structures needed before the shifts are done:
    /// primitive, shift unit, raw slab (no translation), slab
    /// Pass the associated settings along with the structures
    void drop_base(const MultiBase& preshift_structures);

    /// Target directory for the shift files
    fs::path shift_target() const;

    /// Writes out all the shifted slabs in nested directories
    void drop_shifts(const MultiShift& shifted_structures);

private:
    /// Name of the calculations directories
    std::string m_name;

    /// Directory path for a particular Surface coordinate
    fs::path _surface_point_path(const SurfacePoint& point) const;

};

} // namespace mush

#endif
