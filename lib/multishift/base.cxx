#include "./base.hpp"
#include "./exceptions.hpp"
#include "./misc.hpp"
#include "casmutils/frankenstein.hpp"

namespace mush
{
BaseSettings::BaseSettings(const CASM::fs::path& init_prim_path, const Eigen::Vector3i& init_millers,
                           int init_floor_slab_ix, int init_stacks)
    : m_prim_path(init_prim_path),
      m_millers(init_millers),
      m_floor_slab_atom_ix(init_floor_slab_ix),
      m_stacks(init_stacks)
{
    if (m_floor_slab_atom_ix < 0)
    {
        throw except::SettingMustBeGreaterEqualThanZero("floor_slab_index");
    }

    if (m_stacks < 1)
    {
        throw except::SettingMustBeGreaterThanZero("stacks");
    }
}

docs::SettingsInfo BaseSettings::_initialized_documentation()
{
    docs::SettingsInfo docs("base");
    // chain things here
    return docs;
}

const docs::SettingsInfo BaseSettings::docs(BaseSettings::_initialized_documentation());

BaseSettings BaseSettings::from_json(const CASM::jsonParser& init_settings)
{
    auto floor_index = lazy::get_or_value<int>(init_settings, "floor_slab_index", 0);
    auto stacks = lazy::get_or_value<int>(init_settings, "stacks", 1);

    // TODO: Make exception safe. The CASM exceptions are way to generic to forward though.

    return BaseSettings(init_settings["prim"].get<CASM::fs::path>(), init_settings["millers"].get<Eigen::Vector3i>(),
                        floor_index, stacks);
}

CASM::jsonParser BaseSettings::to_json() const
{
    CASM::jsonParser serialized;
    serialized["prim"] = m_prim_path.string();
    /* serialized["millers"]=m_millers; */
    serialized["millers"] =
        std::vector<int>{m_millers(0), m_millers(1),
                         m_millers(2)}; // This is to avoid a really dumb bug in casm that writes Eigen stuff stupid
    serialized["floor_slab_index"] = m_floor_slab_atom_ix;
    serialized["stacks"] = m_stacks;
    return serialized;
}

/* BaseSettings BaseSettings::from_path(const CASM::fs::path& init_path) */
/* { */
/*     CASM::jsonParser settings_dump(init_path); */
/*     return BaseSettings::from_json(settings_dump); */
/* } */

//*********************************************************************************************************

MultiBase::MultiBase(const Structure& init_prim, const Eigen::Vector3i& init_millers, int floor_atom_ix, int stacks)
    : m_prim(init_prim),
      m_shift_unit(this->_shift_unit_from_primitive(init_prim, init_millers)),
      /* m_floored_shift_unit(this->_floored_shift_unit(m_shift_unit, floor_atom_ix)), */
      m_slab(this->_stacked_units(m_shift_unit, stacks)),
      m_floored_slab(this->_floored_structure(m_slab, floor_atom_ix))
{
    if (!m_prim.lattice().is_right_handed())
    {
        throw except::LeftHandedLattice();
    }
}

MultiBase MultiBase::from_settings(const BaseSettings& init_settings)
{
    return MultiBase(Structure::from_poscar(init_settings.prim_path()), init_settings.millers(),
                     init_settings.floor_slab_atom_index(), init_settings.stacks());
}

Structure MultiBase::_shift_unit_from_primitive(const Structure& init_prim, const Eigen::Vector3i& init_millers)
{
    // The 0 means "get the smallest cell possible", and I wish it was the default
    auto shift_lattice = init_prim.lattice().get_lattice_in_plane(init_millers, 0);
    CASM::Structure raw_shift_unit(shift_lattice);
    raw_shift_unit.fill_supercell(init_prim);
    return Structure(raw_shift_unit);
}

Structure MultiBase::_floored_structure(const Structure& shiftable_struc, int floor_atom_ix)
{
    auto shifted_structure = shiftable_struc;

    // Index 0 means no translation
    if (floor_atom_ix != 0)
    {
        const auto& frac_coord = shiftable_struc.basis[floor_atom_ix].frac();
        Frankenstein::shift_coords_by(&shifted_structure, frac_coord);
    }

    return shifted_structure;
}

Structure MultiBase::_stacked_units(const Structure& stackable_struc, int stacks)
{
    // Always stack along c-direction
    Eigen::Matrix3i stack_mat;
    stack_mat << 1, 0, 0, 0, 1, 0, 0, 0, stacks;

    return Simplicity::make_super_structure(stackable_struc, stack_mat);
}

} // namespace mush
