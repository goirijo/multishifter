#include "./fourier.hpp"
#include "./exceptions.hpp"
#include "./misc.hpp"
#include "casm/casm_io/jsonParser.hh"
#include "casm/casm_io/json_io/container.hh"
#include "casm/crystallography/Coordinate.hh"
#include <math.h>
#include <unordered_set>

namespace
{
int num_unique(std::vector<double> vals, long prec = 1e8)
{
    std::unordered_set<long> unique_values;
    for (auto v : vals)
    {
        unique_values.insert(static_cast<long>(v * prec + 0.5));
    }
    return unique_values.size();
}

void sanity_check(const std::vector<mush::InterPoint>& unrolled_data)
{
    long prec = 1e8;
    std::unordered_set<std::pair<long, long>, hashing::pair_hash<long, long>> unique_values;
    for (const auto& p : unrolled_data)
    {
        unique_values.insert(
            std::make_pair(static_cast<long>(p.a_frac * prec + 0.5), static_cast<long>(p.b_frac * prec + 0.5)));
    }

    if (unique_values.size() != unrolled_data.size())
    {
        throw mush::except::BadData("There are multiple values defined per grid point!");
    }

    return;
}

/* CASM::jsonParser serialize(const mush::Lattice& lat) */
/* { */
/*     CASM::jsonParser serialized; */
/*     std::vector<double> a{lat[0](0),lat[0](1),lat[0](2)}; */
/*     std::vector<double> b{lat[1](0),lat[1](1),lat[1](2)}; */
/*     std::vector<double> c{lat[2](0),lat[2](1),lat[2](2)}; */

/*     serialized["a"]=a; */
/*     serialized["b"]=b; */
/*     serialized["c"]=c; */

/*     return serialized; */
/* } */

} // namespace

namespace mush
{

const docs::SettingsInfo FourierSettings::docs(FourierSettings::_initialized_documentation());

docs::SettingsInfo FourierSettings::_initialized_documentation()
{
    docs::SettingsInfo docs("fourier");
    // chain things here
    return docs;
}

FourierSettings::FourierSettings(const fs::path& init_data_path, const fs::path& init_lattice_path,
                                 const std::vector<std::string>& init_value_tags)
    : m_data_path(init_data_path), m_lattice_path(init_lattice_path), m_value_tags(init_value_tags)
{
}

FourierSettings FourierSettings::from_json(const CASM::jsonParser& init_json)
{
    auto data_path = init_json["data"].get<fs::path>();
    auto lattice_path = init_json["lattice"].get<fs::path>();
    auto value_tags = init_json["values"].get<std::vector<std::string>>();

    return FourierSettings(data_path, lattice_path, value_tags);
}

CASM::jsonParser FourierSettings::to_json() const
{
    CASM::jsonParser serialized;
    serialized["data"] = m_data_path.string();
    serialized["lattice"] = m_lattice_path.string();
    serialized["values"] = m_value_tags;
    return serialized;
}

//********************************************************************************************

Interpolator::InterGrid Interpolator::_r_weigner_seitz(const InterGrid& init_values, const Lattice& real_lattice)
{
    InterGrid r_values;

    for (const auto& a_row : init_values)
    {
        InterGrid::value_type r_row;
        for (const auto& val : a_row)
        {
            CASM::Coordinate r_coord(val.a_frac, val.b_frac, 0.0, real_lattice, CASM::FRAC);
            r_coord.voronoi_within();
            r_row.emplace_back(r_coord.const_frac()(0), r_coord.const_frac()(1), val.value);
        }
        r_values.emplace_back(std::move(r_row));
    }

    return r_values;
}

Interpolator::InterGrid Interpolator::_k_grid(const InterGrid& init_values, const Lattice& reciprocal_lattice)
{
    InterGrid k_values;

    auto num_as = init_values.size();
    auto num_bs = init_values[0].size();

    // Yes, you want integers. Anything that doesn't fall on the reciprical lattice points is
    // gonna mess your interpolation up.
    // At this point grids are always even, because when you made the grid from the data, you
    // enforced periodicity by "repeating" the values at the edges.
    assert(num_as % 2 == 1);
    assert(num_bs % 2 == 1);
    int ka_centrize = num_as / 2;
    int kb_centrize = num_bs / 2;
    CASM::Coordinate center_shift(ka_centrize, kb_centrize, 0.0, reciprocal_lattice, CASM::FRAC);

    for (int a = 0; a < init_values.size(); ++a)
    {
        InterGrid::value_type kb_row;
        double ka_frac = static_cast<double>(a) / num_as;
        ka_frac = a;

        const auto& a_row = init_values[a];
        if (a_row.size() != num_bs)
        {
            throw except::DimensionalMismatch(a_row.size(), num_bs,
                                              "The grid of values to interpolate is not uniform.");
        }

        for (int b = 0; b < a_row.size(); ++b)
        {
            double kb_frac = static_cast<double>(b) / num_bs;
            kb_frac = b;
            CASM::Coordinate k_coord(ka_frac, kb_frac, 0.0, reciprocal_lattice, CASM::FRAC);
            k_coord -= center_shift;
            /* k_coord.voronoi_within(); */

            kb_row.emplace_back(k_coord.const_frac()(0), k_coord.const_frac()(1), 0.0);
        }
        k_values.emplace_back(std::move(kb_row));
    }
    return k_values;
}

std::vector<InterPoint*> Interpolator::_unrolled_values(InterGrid* grid_values_ptr)
{
    std::vector<InterPoint*> unrolled_values;
    auto& grid_values = *grid_values_ptr;
    for (auto& row : grid_values)
    {
        for (auto& val : row)
        {
            unrolled_values.push_back(&val);
        }
    }
    return unrolled_values;
}

std::vector<const InterPoint*> Interpolator::_unrolled_values(const InterGrid& grid_values)
{
    // TODO: There's some clever way to avoid the code duplication with some template or const_cast magic, but
    // I can't be bothered right now.
    std::vector<const InterPoint*> unrolled_values;
    for (auto& row : grid_values)
    {
        for (auto& val : row)
        {
            unrolled_values.push_back(&val);
        }
    }
    return unrolled_values;
}

std::pair<Lattice, Interpolator::InterGrid> Interpolator::interpolate(int a_dim, int b_dim) const
{
    std::complex<double> im(0, 1);
    InterGrid interpolated_values;
    const auto k_vals = _unrolled_values(this->m_k_values);
    for (int a = 0; a < a_dim; ++a)
    {
        InterGrid::value_type inter_row;
        for (int b = 0; b < b_dim; ++b)
        {
            InterPoint ipoint(static_cast<double>(a) / a_dim, static_cast<double>(b) / b_dim, 0.0);
            auto r_vec = ipoint.cart(this->m_real_lat);
            for (const auto& k_val : k_vals)
            {
                auto k_vec = k_val->cart(m_recip_lat);
                ipoint.value += k_val->value * k_val->weight * std::exp(im * r_vec.dot(k_vec));
            }
            inter_row.emplace_back(std::move(ipoint));
        }
        interpolated_values.emplace_back(std::move(inter_row));
    }

    return std::make_pair(m_real_lat, interpolated_values);
}

void Interpolator::_take_fourier_transform()
{
    std::complex<double> im(0, 1);
    std::complex<double> normalization = 1.0 / this->size();

    auto real_vals = _unrolled_values(m_real_ipoints);
    auto k_vals = _unrolled_values(&m_k_values);

    assert(real_vals.size() == k_vals.size() && real_vals.size() == this->size());

    /* for(auto k_val : k_vals) */
    /* { */
    /*     k_val->value=0.0; */
    /*     auto k_vec=k_val->cart(m_recip_lat); */
    /*     for(auto r_val : real_vals) */
    /*     { */
    /*         auto r_vec=r_val->cart(m_real_lat); */
    /*         k_val->value+=normalization*r_val->value*std::exp(-im*k_vec.dot(r_vec)); */
    /*     } */
    /* } */

    for (auto& k_row : m_k_values)
    {
        for (auto& k_val : k_row)
        {
            k_val.value = 0.0;
            auto k_vec = k_val.cart(this->m_recip_lat);

            for (const auto& r_row : m_real_ipoints)
            {
                for (const auto& r_val : r_row)
                {
                    auto r_vec = r_val.cart(this->m_real_lat);
                    k_val.value += normalization * r_val.value * r_val.weight * std::exp(-im * k_vec.dot(r_vec));
                }
            }
            /* std::cout << k_val.value << "    " << k_vec.transpose() << std::endl; */
        }
    }

    return;
}

Interpolator::Interpolator(const Lattice& init_lat, const std::vector<double>& a_fracs,
                           const std::vector<double>& b_fracs, const std::vector<double>& vals)
    : Interpolator(init_lat, Interpolator::_grid_from_unrolled_data(a_fracs, b_fracs, vals))
{
}

Interpolator::Interpolator(const Lattice& init_lat, const InterGrid& init_values)
    : m_real_lat(Lattice(phony_lat_mat_reorient(init_lat.lat_column_mat()))),
      m_recip_lat(this->m_real_lat.get_reciprocal()),
      m_real_ipoints(init_values),
      /* m_real_ipoints(this->_r_weigner_seitz(init_values,this->m_real_lat)), */
      m_k_values(this->_k_grid(this->m_real_ipoints, this->m_recip_lat))
{
    this->_take_fourier_transform();
}

std::pair<int, int> Interpolator::dims() const
{
    for (const auto& k_row : m_k_values)
    {
        assert(k_row.size() == m_k_values[0].size());
    }

    // TODO: Probaby overkill but you can also assert that the real
    // values have the right dimensions too

    assert(m_k_values.size() == m_real_ipoints.size());

    return std::make_pair(m_k_values.size(), m_k_values[0].size());
}

int Interpolator::size() const
{
    auto dims = this->dims();
    return dims.first * dims.second;
}

CASM::jsonParser Interpolator::grid_to_json(const Lattice& ref_lat, const InterGrid& grid_values)
{
    Lattice reorient_lat = Lattice(phony_lat_mat_reorient(ref_lat.lat_column_mat()));
    CASM::jsonParser data_dump;
    // TODO: Should this include integer values a and b?

    std::vector<double> a_frac;
    std::vector<double> b_frac;
    std::vector<double> x_cart;
    std::vector<double> y_cart;
    std::vector<double> re_val;
    std::vector<double> im_val;

    for (const auto& row : grid_values)
    {
        for (const auto& val : row)
        {
            auto vec = val.cart(ref_lat);

            a_frac.push_back(val.a_frac);
            b_frac.push_back(val.b_frac);
            x_cart.push_back(vec(0));
            y_cart.push_back(vec(1));
            assert(std::abs(vec(2)) < 0.000001);
            re_val.push_back(val.value.real());
            im_val.push_back(val.value.imag());
        }
    }

    data_dump["a_frac"] = a_frac;
    data_dump["b_frac"] = b_frac;
    data_dump["x_cart"] = x_cart;
    data_dump["y_cart"] = y_cart;
    data_dump["re(value)"] = re_val;
    data_dump["im(value)"] = im_val;

    return data_dump;
}

Interpolator::InterGrid Interpolator::_direct_reshape(const std::vector<mush::InterPoint>& unrolled_data, int ka_dim,
                                                      int kb_dim)
{
    if (ka_dim * kb_dim != unrolled_data.size())
    {
        throw mush::except::DimensionalMismatch(ka_dim * kb_dim, unrolled_data.size(),
                                                "Cannot reshape vector to the specified grid.");
    }

    int i = 0;
    mush::Interpolator::InterGrid final_grid;
    for (int ka = 0; ka < ka_dim; ++ka)
    {
        mush::Interpolator::InterGrid::value_type kb_row;
        for (int kb = 0; kb < kb_dim; ++kb, ++i)
        {
            kb_row.push_back(unrolled_data[i]);
        }
        final_grid.emplace_back(std::move(kb_row));
    }

    return final_grid;
}

void Interpolator::_make_unrolled_data_odd(std::vector<mush::InterPoint>* unrolled_data, int* final_adim,
                                           int* final_bdim)
{
    if (*final_adim % 2 == 0)
    {
        auto size = unrolled_data->size();
        for (int i = 0; i < size; ++i)
        {
            if (lazy::almost_equal(unrolled_data->operator[](i).a_frac, 0.0))
            {
                unrolled_data->operator[](i).weight /= 2.0;
                auto new_point = unrolled_data->operator[](i);
                new_point.a_frac = 1.0;
                unrolled_data->emplace_back(std::move(new_point));
            }
        }

        ++(*final_adim);
    }

    if (*final_bdim % 2 == 0)
    {
        auto size = unrolled_data->size();
        for (int i = 0; i < size; ++i)
        {
            if (lazy::almost_equal(unrolled_data->operator[](i).b_frac, 0.0))
            {
                unrolled_data->operator[](i).weight /= 2.0;
                auto new_point = unrolled_data->operator[](i);
                new_point.b_frac = 1.0;
                unrolled_data->emplace_back(std::move(new_point));
            }
        }

        ++(*final_bdim);
    }
    return;
}

Interpolator::InterGrid Interpolator::_grid_from_unrolled_data(const std::vector<double>& as,
                                                               const std::vector<double>& bs,
                                                               const std::vector<double>& vals)
{
    if (as.size() != bs.size())
    {
        throw mush::except::DimensionalMismatch(
            as.size(), bs.size(), "Cannot create a mesh grid with the provided data points along a and b.");
    }
    if (as.size() != vals.size())
    {
        throw mush::except::DimensionalMismatch(
            as.size(), vals.size(), "The number of grid points does not match the number of values to fit to.");
    }

    // Start by storing the measured grid point values
    std::vector<mush::InterPoint> unrolled_data;
    for (int i = 0; i < vals.size(); ++i)
    {
        unrolled_data.emplace_back(as[i], bs[i], vals[i]);
    }

    int final_adim = num_unique(as);
    int final_bdim = num_unique(bs);

    // In the case of even grids, you need to repeat values at the edges, so that you
    // can properly center the k-point grid later on
    Interpolator::_make_unrolled_data_odd(&unrolled_data, &final_adim, &final_bdim);

    sanity_check(unrolled_data);
    std::sort(unrolled_data.begin(), unrolled_data.end());
    return Interpolator::_direct_reshape(unrolled_data, final_adim, final_bdim);
}

CASM::jsonParser Interpolator::serialize() const
{
    CASM::jsonParser serialized;

    CASM::jsonParser rlat_serialized, klat_serialized;
    CASM::to_json(this->m_real_lat, rlat_serialized);
    CASM::to_json(this->m_recip_lat, klat_serialized);

    serialized["r_lattice"] = rlat_serialized;
    serialized["k_lattice"] = klat_serialized;

    auto dims = this->dims();
    serialized["a_dim"] = dims.first;
    serialized["b_dim"] = dims.second;

    std::vector<CASM::jsonParser> r_unroll, k_unroll;

    for (int a = 0; a < dims.first; ++a)
    {
        for (int b = 0; b < dims.second; ++b)
        {
            r_unroll.push_back(this->sampled_values()[a][b].serialize());
            k_unroll.push_back(this->k_values()[a][b].serialize());
        }
    }

    serialized["r_unrolled"] = r_unroll;
    serialized["k_unrolled"] = k_unroll;

    return serialized;
}

Interpolator Interpolator::deserialize(const CASM::jsonParser& serialized)
{
    Lattice r_lat, k_lat;
    CASM::from_json(r_lat, serialized["r_lattice"]);
    CASM::from_json(k_lat, serialized["k_lattice"]);

    int adim = serialized["a_dim"].get<int>();
    int bdim = serialized["b_dim"].get<int>();

    InterGrid r_grid, k_grid;
    int i = 0;
    for (int a = 0; a < adim; ++a)
    {
        r_grid.push_back(InterGrid::value_type());
        k_grid.push_back(InterGrid::value_type());
        for (int b = 0; b < bdim; ++b, ++i)
        {
            r_grid.back().emplace_back(InterPoint::deserialize(serialized["r_unrolled"][i]));
            k_grid.back().emplace_back(InterPoint::deserialize(serialized["k_unrolled"][i]));
        }
    }

    return Interpolator(std::move(r_lat), std::move(k_lat), std::move(r_grid), std::move(k_grid));
}

Interpolator::Interpolator(Lattice&& init_real, Lattice&& init_recip, InterGrid&& init_rpoints,
                           InterGrid&& init_kpoints)
    : m_real_lat(std::move(init_real)),
      m_recip_lat(std::move(init_recip)),
      m_real_ipoints(std::move(init_rpoints)),
      m_k_values(std::move(init_kpoints))
{
}

//********************************************************************************************

Analytiker::Analytiker(const Interpolator& init_ipolator)
    : m_formula_bits(this->_formula_bits(init_ipolator.k_values())), m_recip_lat(init_ipolator.reciprocal_lattice())
{
    //Default tolerance too strict
    /* assert(lazy::almost_zero(m_recip_lat[0](1))); */
    /* assert(lazy::almost_zero(m_recip_lat[0](2))); */
    /* assert(lazy::almost_zero(m_recip_lat[1](1))); */
    /* assert(lazy::almost_zero(m_recip_lat[2](0))); */
    /* assert(lazy::almost_zero(m_recip_lat[2](1))); */
    /* assert(lazy::almost_equal(m_recip_lat[2](2),1.0)); */
}

std::vector<Analytiker::FormulaBit> Analytiker::_formula_bits(const Interpolator::InterGrid& k_values)
{
    std::vector<FormulaBit> formula_bits;

    int adim = k_values.size();
    int bdim = k_values[0].size();

    // Keep track of which points you visited via inversion
    std::vector<std::vector<bool>> visited(bdim, std::vector<bool>(adim, false));

    assert(adim % 2 == 1);
    assert(bdim % 2 == 1);

    int acex = adim / 2;
    int bcex = bdim / 2;

    // First deal with the gamma point, which is the center of the grid, and has no inversion "twin"
    auto gamma_ptr = std::make_shared<const InterPoint>(k_values[acex][bcex]);
    assert(lazy::almost_zero(gamma_ptr->a_frac));
    assert(lazy::almost_zero(gamma_ptr->b_frac));

    formula_bits.emplace_back(gamma_ptr->value.real() * gamma_ptr->weight, FormulaBitBasis::RECOS, gamma_ptr);
    formula_bits.emplace_back(0.0, FormulaBitBasis::IMSIN, gamma_ptr);
    formula_bits.emplace_back(gamma_ptr->value.imag() * gamma_ptr->weight, FormulaBitBasis::IMCOS, gamma_ptr);
    formula_bits.emplace_back(0.0, FormulaBitBasis::RESIN, gamma_ptr);

    visited[acex][bcex] = true;

    // You only need to loop over half
    // We're pairing up k-points together because we can reduce the number
    // of functions if we remember that sin(-k)=-sin(k) and
    // cos(-k)=cos(k)
    for (int a = 0; a <= acex; ++a)
    {
        for (int b = -bcex; b <= bcex; ++b)
        {
            int cura = acex + a;
            int curb = bcex + b;
            int inva = acex - a;
            int invb = bcex - b;

            if (visited[cura][curb])
            {
                continue;
            }

            const auto& curk = k_values[cura][curb];
            const auto& invk = k_values[inva][invb];
            auto curk_ptr = std::make_shared<const InterPoint>(curk);

            // clang-format off
            formula_bits.emplace_back( curk.value.real() * curk.weight + invk.value.real() * invk.weight,FormulaBitBasis::RECOS, curk_ptr);
            formula_bits.emplace_back( curk.value.real() * curk.weight - invk.value.real() * invk.weight,FormulaBitBasis::IMSIN, curk_ptr);
            formula_bits.emplace_back( curk.value.imag() * curk.weight + invk.value.imag() * invk.weight,FormulaBitBasis::IMCOS, curk_ptr);
            formula_bits.emplace_back(-curk.value.imag() * curk.weight + invk.value.imag() * invk.weight,FormulaBitBasis::RESIN, curk_ptr);
            // clang-format on

            visited[cura][curb] = true;
            visited[inva][invb] = true;
        }
    }

    assert(formula_bits.size() == 4 * (adim * bdim + 1) / 2);
    return formula_bits;
}

bool Analytiker::_can_ignore_imaginary() const
{
    for (const auto& bit : m_formula_bits)
    {
        switch (std::get<1>(bit))
        {
        case FormulaBitBasis::IMSIN:
            if (!lazy::almost_zero(std::get<0>(bit)))
            {
                return false;
            }
            break;

        case FormulaBitBasis::IMCOS:
            if (!lazy::almost_zero(std::get<0>(bit)))
            {
                return false;
            }
            break;
        }
    }

    return true;
}

std::pair<std::string, std::string> Analytiker::python_cart(std::string x_var, std::string y_var,
                                                            std::string numpy) const
{
    std::string real_formula("0"), imag_formula("0");

    for (const auto& bit : m_formula_bits)
    {
        auto value = std::get<0>(bit);
        const auto& kpoint = *std::get<2>(bit).get();
        auto kcart = kpoint.cart(m_recip_lat);
        assert(lazy::almost_zero(kcart(2)));

        if (lazy::almost_zero(value))
        {
            continue;
        }

        // There's three parts to each entry :
        // VALUE*FUNCTION(X1*K1+X2*K2)
        // First convert the value to a string
        std::string value_str;
        if (value > 0)
        {
            value_str.push_back('+');
        }
        value_str += std::to_string(value);

        // Next work on the first entry within the function
        std::string dots = "(" + std::to_string(kcart(0)) + "*" + x_var;

        // Then work on the second entry within the function
        if (kcart(1) >= 0)
        {
            dots.push_back('+');
        }
        dots += std::to_string(kcart(1)) + "*" + y_var + ")";

        switch (std::get<1>(bit))
        {
        case FormulaBitBasis::RECOS:
            real_formula += value_str + "*" + numpy + ".cos" + dots;
            break;
        case FormulaBitBasis::IMSIN:
            imag_formula += value_str + "*" + numpy + ".sin" + dots;
            break;
        case FormulaBitBasis::IMCOS:
            imag_formula += value_str + "*" + numpy + ".cos" + dots;
            break;
        case FormulaBitBasis::RESIN:
            real_formula += value_str + "*" + numpy + ".sin" + dots;
            break;
        }
    }

    return std::make_pair(std::move(real_formula), std::move(imag_formula));
}

} // namespace mush
