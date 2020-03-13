#ifndef MULTISHIFTSHIFT_HH
#define MULTISHIFTSHIFT_HH

#include "./define.hpp"
#include "./autodocs.hpp"
#include "./surf.hpp"
#include "casmutils/structure.hpp"
#include <tuple>
#include <vector>

namespace mush
{

/**
 * Holds values that specify how to grid up the shift plane and how much
 * cleavage to put in between the slabs.
 */

class ShiftSettings
{
public:
    ShiftSettings(const fs::path& init_slab_path, int init_a_points, int init_b_points,
                  const std::vector<double>& init_cleavage);
    static ShiftSettings from_json(const CASM::jsonParser& init_json);
    CASM::jsonParser to_json() const;

    fs::path slab_path() const { return m_slab_path; };
    int a_points() const { return m_a_points; };
    int b_points() const { return m_b_points; };
    const std::vector<double>& cleavage_values() const { return m_cleavage_values; };

    ///Information about each of the settings entries required to construct *this
    static const docs::SettingsInfo docs;

private:
    /// Where to find the slab that's supposed to get shifted
    fs::path m_slab_path;

    /// How many points to sample along the a-vector
    int m_a_points;

    /// How many points to sample along the b-vector
    int m_b_points;

    /// cleavage values. These are the values for how much vacuum you want to insert between the slabs in Angstrom
    std::vector<double> m_cleavage_values;

    /// Generate the documentation for the settings, used to construct static member
    static docs::SettingsInfo _initialized_documentation();

};

/**
 * Initialzed with a slab, shift densities, and cleavage values, creates a set of structures
 * corresponding to shifting the slab along the a-b-plane and inserts different amounts of vaccuum
 * between the slabs.
 */

class MultiShift
{
public:
    MultiShift(const Structure& init_slab, int a_density, int b_density, const std::vector<double>& cleavage_values);
    static MultiShift from_settings(const ShiftSettings& init_settings);

    const std::vector<std::pair<SurfacePoint, Structure>>& shifted_slabs() const { return this->m_shifted_slabs; };

    const Structure reference_slab() const {return this->m_reference_slab;};

private:
    /// The starting slab structure, this corresponds to no shift and no cleavage
    Structure m_reference_slab;

    /// All the shifted structures, paired up with values indicating their grid point and cleavage value
    std::vector<std::pair<SurfacePoint, Structure>> m_shifted_slabs;

    /// Given the number of shift points and cleavage values, run through all combinations
    std::vector<SurfacePoint> _multishift_coordinates(const Structure& ref_slab, int a_density, int b_density,
                                                      const std::vector<double>& cleavage_values);
};
} // namespace mush

#endif
