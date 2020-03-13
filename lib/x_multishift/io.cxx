#include "casmutils/structure.hpp"

#include "./base.hpp"
#include "./exceptions.hpp"
#include "./fourier.hpp"
#include "./io.hpp"
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

void doesnt_exist_or_throw(const fs::path& target)
{
    auto exists = fs::exists(target);

    if (exists)
    {
        throw except::FileExists(target.string());
    }
}

} // namespace io

MultiIO::MultiIO(const std::string& init_name) : m_name(init_name) {}

/* fs::path MultiIO::base_target() const */
/* { */
/*     auto target = io::prefixed_dirname(this->name(), BaseSettings::docs.tag()); */
/*     return target; */
/* } */

void MultiIO::drop_base(const MultiBase& preshift_structures)
{
    auto target = this->target<BaseSettings>();
    io::mkdirs_or_throw(target);

    /* auto settings_dump=settings.to_json(); */
    /* settings_dump.write(target/(BaseSettings::tag()+".json")); */

    Simplicity::write_poscar(preshift_structures.primitive(), target / "prim.vasp");
    Simplicity::write_poscar(preshift_structures.shift_unit(), target / "shift_unit.vasp");
    Simplicity::write_poscar(preshift_structures.raw_slab(), target / "raw_slab.vasp");
    Simplicity::write_poscar(preshift_structures.floored_slab(), target / "final_slab.vasp");

    return;
}

fs::path MultiIO::drop_interpolator(const Interpolator& interpolated_data, const std::string& value_tag) const
{
    std::string filename=value_tag+ ".interpolator.json";
    auto target = this->target<FourierSettings>() / filename;
    auto serialized=interpolated_data.serialize();

    fs::create_directories(this->target<FourierSettings>());
    io::doesnt_exist_or_throw(target);
    serialized.write(target);

    std::string r_filename=value_tag+ ".real.json";
    auto r_target = this->target<FourierSettings>() / r_filename;
    io::doesnt_exist_or_throw(r_target);
    const auto& r_lat=interpolated_data.real_lattice();
    const auto& r_val=interpolated_data.sampled_values();
    Interpolator::grid_to_json(r_lat,r_val).write(r_target);

    std::string k_filename=value_tag+ ".rcal.json";
    auto k_target = this->target<FourierSettings>() / k_filename;
    io::doesnt_exist_or_throw(k_target);
    const auto& k_lat=interpolated_data.reciprocal_lattice();
    const auto& k_val=interpolated_data.k_values();
    Interpolator::grid_to_json(k_lat,k_val).write(k_target);

    return target;
}

/* fs::path MultiIO::shift_target() const */
/* { */
/*     auto target = io::prefixed_dirname(this->name(), ShiftSettings::docs.tag()); */
/*     return target; */
/* } */

fs::path MultiIO::_surface_point_path(const SurfacePoint& point) const
{
    fs::path shift_dir("shift_" + std::to_string(point.a) + "." + std::to_string(point.b));

    std::stringstream cleavestream;
    cleavestream << std::fixed << std::setprecision(6) << point.cleavage;
    fs::path cleave_dir("cleave_" + cleavestream.str());

    return shift_dir / cleave_dir;
}

void MultiIO::drop_shifts(const MultiShift& shifted_structures) const
{
    std::vector<int> as;
    std::vector<int> bs;
    std::vector<double> a_fracs;
    std::vector<double> b_fracs;
    std::vector<double> x_carts;
    std::vector<double> y_carts;
    std::vector<double> cleavages;
    std::vector<std::string> pathnames;

    for (const auto point_struc_pair : shifted_structures.shifted_slabs())
    {
        const auto& shift_point = point_struc_pair.first;
        const auto& shifted_slab = point_struc_pair.second;

        const auto target = this->target<ShiftSettings>() / this->_surface_point_path(shift_point);
        io::mkdirs_or_throw(target);

        Simplicity::write_poscar(shifted_slab, target / "POSCAR");

        as.push_back(shift_point.a);
        bs.push_back(shift_point.b);
        a_fracs.push_back(shift_point.a_frac);
        b_fracs.push_back(shift_point.b_frac);
        x_carts.push_back(shift_point.x_cart);
        y_carts.push_back(shift_point.y_cart);
        cleavages.push_back(shift_point.cleavage);
        pathnames.push_back(target.string());
    }

    // Keep a record of everything
    CASM::jsonParser shift_record;
    shift_record["a"] = as;
    shift_record["b"] = bs;
    shift_record["a_frac"] = a_fracs;
    shift_record["b_frac"] = b_fracs;
    shift_record["x_cart"] = x_carts;
    shift_record["y_cart"] = y_carts;
    shift_record["cleavage"] = cleavages;
    shift_record["path"] = pathnames;
    shift_record.write(this->target<ShiftSettings>() / "record.json");

    return;
}

std::tuple<std::vector<double>, std::vector<double>, std::vector<std::vector<double>>>
MultiIO::unrolled_frac_grid_data_from_json(const CASM::jsonParser& json_data,
                                           const std::vector<std::string>& value_tags)
{
    // I hate jsonParser
    std::vector<double> as;
    CASM::from_json(as, json_data["a_frac"]);
    /* auto as=json_data.get<std::vector<double>>(); */
    std::vector<double> bs;
    CASM::from_json(bs, json_data["b_frac"]);

    std::vector<std::vector<double>> vals;
    for (const auto& tag : value_tags)
    {
        std::vector<double> tag_data;
        CASM::from_json(tag_data, json_data[tag]);
        vals.emplace_back(std::move(tag_data));
    }

    return std::make_tuple(as, bs, vals);
}

} // namespace mush
