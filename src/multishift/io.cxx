#include "./io.hpp"
#include "./base.hpp"
#include "./exceptions.hpp"
#include "casm/CASM_global_definitions.hh"

namespace mush
{
namespace io
{
std::string prefixed_dirname(const std::string& prefix, const std::string& value) { return prefix+"."+value; }
} // namespace io

MultiIO::MultiIO(const std::string& init_name) : m_name(init_name) {}

    void MultiIO::drop_base(const BaseSettings& settings, const MultiBase& preshift_structures)
    {
        CASM::fs::path target=io::prefixed_dirname(this->name(),BaseSettings::tag());
        auto exists=!CASM::fs::create_directory(target);

        if(exists)
        {
            throw except::FileExists(target.string());
        }

        auto settings_dump=settings.to_json();
        settings_dump.write(target/(BaseSettings::tag()+".json"));

        Simplicity::write_poscar(preshift_structures.primitive(),target/"prim.vasp");
        Simplicity::write_poscar(preshift_structures.shift_unit(),target/"shift_unit.vasp");
        Simplicity::write_poscar(preshift_structures.raw_slab(),target/"raw_slab.vasp");
        Simplicity::write_poscar(preshift_structures.floored_slab(),target/"final_slab.vasp");

        return;
    }

} // namespace mush
