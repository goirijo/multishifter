#include "casmutils/structure.hpp"

#include "./io.hpp"
#include "./base.hpp"
#include "./exceptions.hpp"
#include "./shift.hpp"

#include <sstream>

namespace mush
{
namespace io
{
std::string prefixed_dirname(const std::string& prefix, const std::string& value) { return prefix + "." + value; }

void mkdirs_or_throw(const fs::path& target)
{
    auto exists = !fs::create_directories(target);

    if (exists)
    {
        throw except::FileExists(target.string());
    }
}

} // namespace io

MultiIO::MultiIO(const std::string& init_name) : m_name(init_name) {}

fs::path MultiIO::base_target() const
{
    auto target = io::prefixed_dirname(this->name(), BaseSettings::tag());
    return target;
}

void MultiIO::drop_base(const MultiBase& preshift_structures)
{
    auto target = this->base_target();
    io::mkdirs_or_throw(target);

    /* auto settings_dump=settings.to_json(); */
    /* settings_dump.write(target/(BaseSettings::tag()+".json")); */

    Simplicity::write_poscar(preshift_structures.primitive(), target / "prim.vasp");
    Simplicity::write_poscar(preshift_structures.shift_unit(), target / "shift_unit.vasp");
    Simplicity::write_poscar(preshift_structures.raw_slab(), target / "raw_slab.vasp");
    Simplicity::write_poscar(preshift_structures.floored_slab(), target / "final_slab.vasp");

    return;
}

fs::path MultiIO::shift_target() const
{
    auto target = io::prefixed_dirname(this->name(), ShiftSettings::tag());
    return target;
}

fs::path MultiIO::_surface_point_path(const SurfacePoint& point) const
{
    fs::path shift_dir("shift_" + std::to_string(point.a) + "." + std::to_string(point.b));

    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << point.cleaveage;
    fs::path cleave_dir("cleave_" + cleavestream.str());

    return shift_dir / cleave_dir;
}

void MultiIO::drop_shifts(const MultiShift& shifted_structures) {
    std::vector<int> as;
    std::vector<int> bs;
    std::vector<double> a_fracs;
    std::vector<double> b_fracs;
    std::vector<double> x_carts;
    std::vector<double> y_carts;
    std::vector<double> cleavages;
    std::vector<std::string> pathnames;


    for(const auto point_struc_pair : shifted_structures.shifted_slabs())
    {
        const auto& shift_point=point_struc_pair.first;
        const auto& shifted_slab=point_struc_pair.second;

        const auto target=this->shift_target()/this->_surface_point_path(shift_point);
        io::mkdirs_or_throw(target);
    
        Simplicity::write_poscar(shifted_slab,target/"POSCAR");

        as.push_back(shift_point.a);
        bs.push_back(shift_point.b);
        a_fracs.push_back(shift_point.a_frac);
        b_fracs.push_back(shift_point.b_frac);
        x_carts.push_back(shift_point.x_cart);
        y_carts.push_back(shift_point.y_cart);
        cleavages.push_back(shift_point.cleaveage);
        pathnames.push_back(target.string());
    }

    //Keep a record of everything
    CASM::jsonParser shift_record;
    shift_record["a"]=as;
    shift_record["b"]=bs;
    shift_record["a_frac"]=a_fracs;
    shift_record["b_frac"]=b_fracs;
    shift_record["x_cart"]=x_carts;
    shift_record["y_cart"]=y_carts;
    shift_record["cleaveage"]=cleavages;
    shift_record["path"]=pathnames;
    shift_record.write(this->shift_target()/"record.json");
    

    return;
}

} // namespace mush
