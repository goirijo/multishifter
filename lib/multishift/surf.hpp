#ifndef MULTISHIFTSURF_HH
#define MULTISHIFTSURF_HH

#include "casm/external/Eigen/Dense"
#include <complex>
#include <vector>

namespace CASM
{
class Lattice;
class jsonParser;
}

namespace mush
{

/// Given a lattice column matrix, return a lattice column matrix that has
/// rotated the values such that the a vector points along the x-axis,
/// axb points along the z axis, and information about the c vector is lost.
Eigen::Matrix3d phony_lat_mat_reorient(const Eigen::Matrix3d& real_lat_mat);

/**
 * Just a collection of variables:
 * the integer shift values for a and b, along with the associated lattice
 * vector fraction of a and b. Also the cleavage value.
 */

class SurfacePoint
{
public:
    static std::vector<SurfacePoint> multishift_coordinates(const CASM::Lattice& surf_lattice, int a_density,
                                                            int b_density, const std::vector<double>& cleavage_values);

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
};

/**
 * One point on the surface that should be used for the interpolation.
 * Only fractional coordinates needed.
 */

class InterPoint
{
public:
    InterPoint(double init_a_frac, double init_b_frac, std::complex<double> init_value)
        : a_frac(init_a_frac), b_frac(init_b_frac), value(init_value)
    {
    }

    double a_frac;
    double b_frac;
    std::complex<double> value;

    /// Return the cartesian coordinate relative to the given lattice. Position will land on the ab-plane, since there is no c-vector component.
    /// This could probably just return a 2d-vector, since your lattice is probably oriented such that ab-vectors are on xy-plane, but
    /// just in case it isn't, it'll stay as a 3d-vector;
    Eigen::Vector3d cart(const CASM::Lattice& ref_lat) const;
    
    /// Serialize into a jsonParser for immediate reconstruction
    CASM::jsonParser serialize() const;

    /// Reconstruct from a serialized jsonParser
    static InterPoint deserialize(const CASM::jsonParser& serialized);

private:
};

} // namespace mush

#endif
