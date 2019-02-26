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

} // namespace

namespace mush
{
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

    //Yes, you want integers. Anything that doesn't fall on the reciprical lattice points is
    //gonna mess your interpolation up.
    //TODO: Should you worry about even grids? It could potentially cause trouble that the k-points
    //on the edges don't have a "twin"
    int ka_centrize = num_as / 2;
    int kb_centrize = num_bs / 2;
    CASM::Coordinate center_shift(ka_centrize, kb_centrize, 0.0, reciprocal_lattice, CASM::FRAC);
    std::cout<<center_shift.const_frac().transpose()<<std::endl;

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
                ipoint.value += k_val->value * std::exp(im * r_vec.dot(k_vec));
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
                    k_val.value += normalization * r_val.value * std::exp(-im * k_vec.dot(r_vec));
                }
            }
            /* std::cout << k_val.value << "    " << k_vec.transpose() << std::endl; */
        }
    }

    return;
}

Interpolator::Interpolator(const Lattice& init_lat, const std::vector<double>& a_fracs,
                           const std::vector<double>& b_fracs, const std::vector<double>& vals):
    Interpolator(init_lat, Interpolator::_grid_from_unrolled_data(a_fracs,b_fracs,vals))
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

    // You have to ensure there's only a single data point per grid point
    std::vector<mush::InterPoint> unrolled_data;
    for (int i = 0; i < vals.size(); ++i)
    {
        unrolled_data.emplace_back(as[i], bs[i], vals[i]);
    }

    sanity_check(unrolled_data);
    return Interpolator::_direct_reshape(unrolled_data, num_unique(as), num_unique(bs));
}

} // namespace mush
