#include "./surf.hpp"
#include "casm/crystallography/Lattice.hh"

namespace mush
{

Eigen::Matrix3d phony_lat_mat_reorient(const Eigen::Matrix3d& real_lat_mat)
{
    // We'll start with a zeroed matrix and fill in values
    Eigen::Matrix3d phony_lat_mat = Eigen::Matrix3d::Zero();
    auto dotprod = real_lat_mat.col(0).dot(real_lat_mat.col(1));
    auto bb = real_lat_mat.col(1).dot(real_lat_mat.col(1));

    // This is pretty dumb, but here it goes:
    // a-vector needs to get oriented along x, so the values are just
    // its length for the x component, and zeros for y and z
    phony_lat_mat(0, 0) = real_lat_mat.col(0).norm();

    // We know that the dot product between a and b should remain the same
    // a1b1+a2b2+a3b3=dotprod
    // but also a2=a3=b3=0, so
    // b1=dotprod/a1
    phony_lat_mat(0, 1) = dotprod / phony_lat_mat(0, 0);

    // The dot product of the b vector with itself should remain the same too
    // b1b1+b2b2+b3b3=bb
    // b2b2=bb-b1b1
    phony_lat_mat(1, 1) = std::sqrt(bb - phony_lat_mat(0, 1) * phony_lat_mat(0, 1));

    // Set the c-vector to just be a unit perpendicular to the ab-vectors
    phony_lat_mat.col(2) = phony_lat_mat.col(0).cross(phony_lat_mat.col(1)).normalized();
    return phony_lat_mat;
}

std::vector<SurfacePoint> SurfacePoint::multishift_coordinates(const CASM::Lattice& surf_lattice, int a_density,
                                                               int b_density,
                                                               const std::vector<double>& cleavage_values)
{
    /* auto phony_lat_mat = SurfacePoint::_phony_lat_mat_reorient(init_slab.lattice().lat_column_mat()); */
    auto phony_lat_mat = phony_lat_mat_reorient(surf_lattice.lat_column_mat());
    auto oriented_a_vec = phony_lat_mat.col(0);
    auto oriented_b_vec = phony_lat_mat.col(1);

    std::vector<SurfacePoint> all_coords;
    for (int a = 0; a < a_density; ++a)
    {
        for (int b = 0; b < b_density; ++b)
        {
            for (const auto cleave : cleavage_values)
            {
                double a_frac = static_cast<double>(a) / static_cast<double>(a_density);
                double b_frac = static_cast<double>(b) / static_cast<double>(b_density);

                auto cart_2d = oriented_a_vec * a_frac + oriented_b_vec * b_frac;
                double x_cart = cart_2d(0);
                double y_cart = cart_2d(1);

                // can't emplace because private constructor
                all_coords.push_back(SurfacePoint(a, b, a_frac, b_frac, x_cart, y_cart, cleave));
            }
        }
    }

    return all_coords;
}

//***********************************************************************

Eigen::Vector3d InterPoint::cart(const CASM::Lattice& ref_lat) const
{
    auto vec = a_frac * ref_lat[0] + b_frac * ref_lat[1] + 0.0 * ref_lat[2];
    return vec;
}

} // namespace mush
