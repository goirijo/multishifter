#include "./shift.hpp"
#include "./exceptions.hpp"
#include "./misc.hpp"

namespace mush
{
ShiftSettings::ShiftSettings(const CASM::fs::path& init_slab_path, int init_a_points, int init_b_points,
                             const std::vector<double>& init_cleavage)
    : m_slab_path(init_slab_path),
      m_a_points(init_a_points),
      m_b_points(init_b_points),
      m_cleavage_values(init_cleavage)
{
    if (init_a_points < 1)
    {
        throw except::SettingMustBeGreaterThanZero("a");
    }

    if (init_b_points < 1)
    {
        throw except::SettingMustBeGreaterThanZero("b");
    }
}

docs::SettingsInfo ShiftSettings::_initialized_documentation()
{
    docs::SettingsInfo docs("shift");
    // chain things here
    return docs;
}

const docs::SettingsInfo ShiftSettings::docs(ShiftSettings::_initialized_documentation());

ShiftSettings ShiftSettings::from_json(const CASM::jsonParser& init_json)
{
    auto aval = lazy::get_or_value<int>(init_json, "a", 1);
    auto bval = lazy::get_or_value<int>(init_json, "b", 1);
    auto cleave = lazy::get_or_value<std::vector<double>>(init_json, "cleavage", std::vector<double>{0.0});

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

    auto mush_coords = SurfacePoint::multishift_coordinates(init_slab.lattice(), a_density, b_density, cleavage_values);

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
