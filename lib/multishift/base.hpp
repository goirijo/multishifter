#ifndef MULTIBUNDLE_HH
#define MULTIBUNDLE_HH

#include "./define.hpp"
#include "casmutils/structure.hpp"

namespace mush
{
/**
 * Simple structure to hold any settings required to generate the shifted
 * structures, starting from the primitive cell.
 * A bit on nomenclature:
 * "primitive" means the smallest unit that you begin your project with.
 * "shift unit" is the smallest structure possible that exposes the plane you want
 * to shift along as one of its unit cell facets, it is a supercell of the
 * primitive.
 * "stack" or "slab" is a 1 dimensional supercell of the slab, constructed by stacking
 * multiple slabs along the shift plane.
 */

class BaseSettings
{
public:
    BaseSettings(const fs::path& init_prim_path, const Eigen::Vector3i& init_millers, int init_floor_slab_ix, int init_stacks);
    static BaseSettings from_json(const CASM::jsonParser& init_settings);
    CASM::jsonParser to_json() const;
    /* static BaseSettings from_path(const fs::path& init_path); */

    const fs::path& prim_path() const { return m_prim_path; }
    const Eigen::Vector3i& millers() const { return m_millers; }
    const int& floor_slab_atom_index() const { return m_floor_slab_atom_ix; }
    const int& stacks() const { return m_stacks; }

    static std::string tag() {return "base";};

private:
    /// The path to the primitive cell
    fs::path m_prim_path;

    /// The miller indexes for the desired shift plane relative to the primitive structure
    Eigen::Vector3i m_millers;

    /// Which atom of the slab to bring down to the base of the plane. Index is relative to the jumbo slab structure.
    /// You need to specify this if you have different possible planes the shifts can occur in.
    int m_floor_slab_atom_ix;

    /// Specifies how thick the slabs your shifting should be, i.e. how many times you want to
    /// stacl the shift unit
    int m_stacks;
};


/**
 * Handles constructing preliminary structures given the primitive
 * cell and some settings parameters.
 * Creates the smallest cell with the shift plane exposed, and stacks
 * it into slabs.
 * The coordinates of the slab atoms can be translated so that the shift
 * plane appears between different layers of atoms. The translation is
 * specified by giving the index (starting from 1) of the atom you want to translate down
 * to the a-b plane. Index 0 indicates no translation is necessary.
 */

class MultiBase
{
    public:

        ///Explicit construction using all necessary parameters
        MultiBase(const Structure& init_prim, const Eigen::Vector3i& init_millers, int floor_atom_ix, int stacks);

        ///Construct via settigs object
        static MultiBase from_settings(const BaseSettings& init_settings);

        Structure primitive() const {return m_prim;};
        Structure shift_unit() const {return m_shift_unit;};
        Structure raw_slab() const {return m_slab;};
        Structure floored_slab() const {return m_floored_slab;};

    private:

        ///The true primitive cell
        Structure m_prim;

        ///The smallest cell whose a-b vectors define the shift plane
        Structure m_shift_unit;

        ///Same as the shift unit, but after translating the basis such that the specified floor atom ends up with c=0 (at the base of the a-b plane)
        /* Structure m_floored_shift_unit; */

        ///Supercell of the shift unit along the c-direction, stacked the specified number of times
        Structure m_slab;

        ///Supercell of the floored shift unit along the c-direction, stacked the specified number of times.
        ///This is the slab that will get shifted around after you start specifiyng the shift grid.
        Structure m_floored_slab;

        ///Returns smallest possible structure with the shift plane exposed along the a-b vectors, but
        ///does not translate the basis in any way
        static Structure _shift_unit_from_primitive(const Structure& init_prim, const Eigen::Vector3i& init_millers);

        ///Translate the atoms of the Structure so that the coordinates of the specified atom index end up on the a-b plane
        static Structure _floored_structure(const Structure& shiftable_struc, int floor_atom_ix);

        ///Stack the given Structure along the c direction the specified number of times
        static Structure _stacked_units(const Structure& stackable_struc, int stacks);
};

}

#endif
