#include "./misc.hpp"
#include "./shift.hpp"
#include "./exceptions.hpp"

namespace mush
{
ShiftSettings::ShiftSettings(const CASM::fs::path& init_slab_path, int init_a_points, int init_b_points,
                             const std::vector<double>& init_cleavage)
    : m_slab_path(init_slab_path),
      m_a_points(init_a_points),
      m_b_points(init_b_points),
      m_cleavage_values(init_cleavage)
{
    if(init_a_points<1)
    {
        throw except::SettingMustBeGreaterThanZero("a");
    }

    if(init_b_points<1)
    {
        throw except::SettingMustBeGreaterThanZero("b");
    }
}

ShiftSettings ShiftSettings::from_json(const CASM::jsonParser& init_json)
{
    auto aval=lazy::get_or_value<int>(init_json,"a",1);
    auto bval=lazy::get_or_value<int>(init_json,"b",1);
    auto cleave=lazy::get_or_value<std::vector<double>>(init_json,"cleavage",std::vector<double>{0.0});

    return ShiftSettings(init_json["slab"].get<CASM::fs::path>(), aval, bval, cleave);
}

CASM::jsonParser ShiftSettings::to_json() const
{
    CASM::jsonParser serialized;
    serialized["slab"] = m_slab_path.string();
    /* serialized["millers"]=m_millers; */
    serialized["a"] = m_a_points;
    serialized["b"] = m_b_points;
    serialized["cleavage"] = m_cleavage_values;
    return serialized;
}

//***********************************************************************

Eigen::Matrix3d SurfacePoint::_phony_lat_mat_reorient(const Eigen::Matrix3d& real_lat_mat)
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

    return phony_lat_mat;
}

std::vector<SurfacePoint> SurfacePoint::multishift_coordinates(const Structure& init_slab, int a_density, int b_density,
                                                               const std::vector<double>& cleavage_values)
{
    auto phony_lat_mat = SurfacePoint::_phony_lat_mat_reorient(init_slab.lattice().lat_column_mat());
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

MultiShift MultiShift::from_settings(const ShiftSettings& init_settings)
{
    return MultiShift(Structure::from_poscar(init_settings.slab_path()), init_settings.a_points(),
                      init_settings.b_points(), init_settings.cleavage_values());
}

MultiShift::MultiShift(const Structure& init_slab, int a_density, int b_density,
                       const std::vector<double>& cleavage_values)
    : m_reference_slab(init_slab)
{
    if (!m_reference_slab.lattice().is_right_handed())
    {
        throw except::LeftHandedLattice();
    }

    auto mush_coords = SurfacePoint::multishift_coordinates(init_slab, a_density, b_density, cleavage_values);

    const auto a_vector = std::get<0>(this->m_reference_slab.lattice().vectors());
    const auto b_vector = std::get<1>(this->m_reference_slab.lattice().vectors());
    const auto c_vector = std::get<2>(this->m_reference_slab.lattice().vectors());
    // This is the unit vector for the direction perpendicular to the shift plane
    const auto normal_unit = a_vector.cross(b_vector).normalized();

    for (const auto& mush_coord : mush_coords)
    {
        const auto new_c_vector =
            c_vector + mush_coord.a_frac * a_vector + mush_coord.b_frac * b_vector + mush_coord.cleavage * normal_unit;
        CASM::Lattice new_lattice(a_vector, b_vector, new_c_vector);

        this->m_shifted_slabs.push_back(std::make_pair(mush_coord, m_reference_slab));
        this->m_shifted_slabs.back().second.set_lattice(new_lattice, CASM::CART);
    }
}

/* std::vector<SurfacePoint> MultiShift::_multishift_coordinates(const Structure& ref_slab, int a_density, int
 * b_density, */
/*                                                               const std::vector<double>& cleavage_values) */
/* { */

/*     std::vector<SurfacePoint> all_coords; */
/*     for (int a = 0; a < a_density; ++a) */
/*     { */
/*         for (int b = 0; b < b_density; ++b) */
/*         { */
/*             for (const auto cleave : cleavage_values) */
/*             { */
/*                 double a_frac = static_cast<double>(a) / static_cast<double>(a_density); */
/*                 double b_frac = static_cast<double>(b) / static_cast<double>(b_density); */

/*                 all_coords.emplace_back(a, b, a_frac, b_frac, cleave); */
/*             } */
/*         } */
/*     } */

/*     return all_coords; */
/* } */

} // namespace mush
