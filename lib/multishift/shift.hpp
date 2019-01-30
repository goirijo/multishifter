#ifndef MULTISHIFTSHIFT_HH
#define MULTISHIFTSHIFT_HH

#include "./define.hpp"
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

    static std::string tag() { return "shift"; };

private:
    /// Where to find the slab that's supposed to get shifted
    fs::path m_slab_path;

    /// How many points to sample along the a-vector
    int m_a_points;

    /// How many points to sample along the b-vector
    int m_b_points;

    /// cleavage values. These are the values for how much vacuum you want to insert between the slabs in Angstrom
    std::vector<double> m_cleavage_values;
};

/**
 * Just a collection of variables:
 * the integer shift values for a and b, along with the associated lattice
 * vector fraction of a and b. Also the cleavage value.
 */

class SurfacePoint
{
public:
    static std::vector<SurfacePoint> multishift_coordinates(const Structure& init_slab, int a_density, int b_density,
                                                            const std::vector<double>& cleavage_values);

     const int a;
     const int b;
     const double a_frac;
     const double b_frac;
     const double x_cart;
     const double y_cart;
     const double cleavage;

private:
    SurfacePoint(int a, int b, double a_frac, double b_frac, double x_cart, double y_cart, double cleavage)
        : a(a), b(b), a_frac(a_frac), b_frac(b_frac), x_cart(x_cart), y_cart(y_cart), cleavage(cleavage){};

    ///Given a lattice column matrix, return a lattice column matrix that has
    ///rotated the values such that the a vector points along the x-axis,
    ///axb points along the z axis, and information about the c vector is lost.
    static Eigen::Matrix3d _phony_lat_mat_reorient(const Eigen::Matrix3d& real_lat_mat);
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
