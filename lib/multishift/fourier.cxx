#include "./fourier.hpp"
#include <casmutils/mush/twist.hpp>
#include <casmutils/mush/slab.hpp>
#include "casmutils/xtal/site.hpp"
#include <casmutils/xtal/coordinate.hpp>
#include <stdexcept>
#include <unordered_set>
#include <utility>

namespace
{

template <typename T>
bool almost_equal(T lhs, T rhs, T tol = 1e-12)
{
    return std::abs(lhs - rhs) < tol;
}

std::pair<int, int> num_unique(const std::vector<mush::InterPoint>& unrolled_data, long prec = 1e8)
{
    std::unordered_set<long> unique_a;
    std::unordered_set<long> unique_b;
    for (const auto& v : unrolled_data)
    {
        unique_a.insert(static_cast<long>(v.a_frac * prec + 0.5));
        unique_b.insert(static_cast<long>(v.b_frac * prec + 0.5));
    }
    return std::make_pair(unique_a.size(), unique_b.size());
}
} // namespace
namespace mush
{
Eigen::Vector3d InterPoint::cart(const cu::xtal::Lattice& ref_lat) const
{
    Eigen::Vector3d vec = a_frac * ref_lat.a() + b_frac * ref_lat.b() + 0.0 * ref_lat.c();
    return vec;
}

bool InterPoint::operator<(const InterPoint& rhs) const
{
    // If a-fraction is the same, consider b-fraction
    if (almost_equal(this->a_frac, rhs.a_frac))
    {
        return this->b_frac < rhs.b_frac;
    }

    return this->a_frac < rhs.a_frac;
}

//*********************************************************************************//

Interpolator::Interpolator(const Lattice& init_lat, const std::vector<InterPoint>& real_data)
    : Interpolator(init_lat, _grid_from_unrolled_data(real_data))
{
}

Interpolator::Lattice make_phony_aligned_lattice(const Interpolator::Lattice& real_lat)
{
    Interpolator::Lattice real_aligned=make_aligned(real_lat);
    return Interpolator::Lattice(real_aligned.a(),real_aligned.b(),Eigen::Vector3d(0,0,1));
}

Interpolator::Interpolator(const Lattice& init_lat, const InterGrid& init_values)
    : m_real_lat(make_phony_aligned_lattice(init_lat)),
      m_recip_lat(cu::xtal::make_reciprocal(this->m_real_lat)),
      m_real_ipoints(init_values),
      /* m_real_ipoints(this->_r_weigner_seitz(init_values,this->m_real_lat)), */
      m_k_values(this->_k_grid(this->m_real_ipoints, this->m_recip_lat))
{
    this->_take_fourier_transform();
}

void Interpolator::_take_fourier_transform()
{
    std::complex<double> im(0, 1);
    double weight_sum=0.0;
    for (const auto& r_row : m_real_ipoints)
    {
        for (const auto& r_val : r_row)
        {
            weight_sum+=r_val.weight;
        }
    }
    std::complex<double> normalization = 1.0 / weight_sum;

    auto real_vals = _unrolled_values(m_real_ipoints);
    auto k_vals = _unrolled_values(&m_k_values);

    assert(real_vals.size() == k_vals.size() && real_vals.size() == this->size());

    for (auto& k_row : m_k_values)
    {
        for (auto& k_val : k_row)
        {
            k_val.value = 0.0;
            Eigen::Vector3d k_vec = k_val.cart(this->m_recip_lat);

            for (const auto& r_row : m_real_ipoints)
            {
                for (const auto& r_val : r_row)
                {
                    Eigen::Vector3d r_vec = r_val.cart(this->m_real_lat);
                    k_val.value += normalization * r_val.value * r_val.weight * std::exp(-im * k_vec.dot(r_vec));
                }
            }
            std::cout<<"DEBUGGING: k_val.value is "<<k_val.value<<std::endl;
            
        }
    }

    return;
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


Interpolator::InterGrid Interpolator::_grid_from_unrolled_data(std::vector<InterPoint> unrolled_data)
{
    auto [final_adim, final_bdim] = ::num_unique(unrolled_data);

    // In the case of even grids, you need to repeat values at the edges, so that you
    // can properly center the k-point grid later on
    Interpolator::_make_unrolled_data_odd(&unrolled_data, &final_adim, &final_bdim);

    std::sort(unrolled_data.begin(), unrolled_data.end());
    return Interpolator::_direct_reshape(unrolled_data, final_adim, final_bdim);
}

void Interpolator::_make_unrolled_data_odd(std::vector<mush::InterPoint>* unrolled_data, int* final_adim, int* final_bdim)
{
    if (*final_adim % 2 == 0)
    {
        auto size = unrolled_data->size();
        for (int i = 0; i < size; ++i)
        {
            if (::almost_equal(unrolled_data->operator[](i).a_frac, 0.0))
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
            if (::almost_equal(unrolled_data->operator[](i).b_frac, 0.0))
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

Interpolator::InterGrid Interpolator::_direct_reshape(const std::vector<mush::InterPoint>& unrolled_data, int ka_dim, int kb_dim)
{
    if (ka_dim * kb_dim != unrolled_data.size())
    {
        throw std::runtime_error("Cannot reshape vector to the specified grid.");
    }

    int i = 0;
    InterGrid final_grid;
    for (int ka = 0; ka < ka_dim; ++ka)
    {
        InterGrid::value_type kb_row;
        for (int kb = 0; kb < kb_dim; ++kb, ++i)
        {
            kb_row.push_back(unrolled_data[i]);
        }
        final_grid.emplace_back(std::move(kb_row));
    }

    return final_grid;
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
    Eigen::Vector3d center_shift(ka_centrize,kb_centrize,0.0);

    for (int a = 0; a < init_values.size(); ++a)
    {
        InterGrid::value_type kb_row;
        double ka_frac = static_cast<double>(a) / num_as;
        ka_frac = a;

        const auto& a_row = init_values[a];
        if (a_row.size() != num_bs)
        {
            throw std::runtime_error("The grid of values to interpolate is not uniform.");
        }

        for (int b = 0; b < a_row.size(); ++b)
        {
            double kb_frac = static_cast<double>(b) / num_bs;
            kb_frac = b;
            Eigen::Vector3d k_coord(ka_frac,kb_frac,0.0);
            k_coord -= center_shift;
            /* k_coord.voronoi_within(); */

            /* kb_row.emplace_back(std::move(k_coord)); */
            kb_row.emplace_back(k_coord(0), k_coord(1), 0.0);
        }
        k_values.emplace_back(std::move(kb_row));
    }

    assert(k_values.size()==init_values.size());
    assert(k_values[0].size()==init_values[0].size());

    /* for(int i=0; i<k_values.size(); ++i) */
    /* { */
    /*     for(int j=0; j<k_values[i].size(); ++j) */
    /*     { */
    /*         std::cout<<init_values[i][j].weight<<"    "; */
    /*         k_values[i][j].weight=init_values[i][j].weight; */
    /*     } */
    /*     std::cout<<"\n"; */
    /* } */

    return k_values;
}

std::pair<cu::xtal::Lattice, Interpolator::InterGrid> Interpolator::interpolate(int a_dim, int b_dim) const
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

//********************************************************************************************

Analytiker::Analytiker(const Interpolator& init_ipolator)
    : m_formula_bits(this->_formula_bits(init_ipolator.k_values())), m_recip_lat(init_ipolator.reciprocal_lattice())
{
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
    assert(almost_equal(gamma_ptr->a_frac,0.0));
    assert(almost_equal(gamma_ptr->b_frac,0.0));

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

    /* assert(formula_bits.size() == 4 * (adim * bdim + 1) / 2); */
    // TODO: Sort by magnitude
    std::sort(formula_bits.begin(), formula_bits.end(), [](const FormulaBit& lhs, const FormulaBit& rhs) -> bool {
        return std::abs(std::get<0>(lhs)) < std::abs(std::get<0>(rhs));
    });
    return formula_bits;
}

bool Analytiker::_can_ignore_imaginary() const
{
    for (const auto& bit : m_formula_bits)
    {
        switch (std::get<1>(bit))
        {
        case FormulaBitBasis::IMSIN:
            if (!almost_equal(std::get<0>(bit),0.0))
            {
                return false;
            }
            break;

        case FormulaBitBasis::IMCOS:
            if (!almost_equal(std::get<0>(bit),0.0))
            {
                return false;
            }
            break;

        default:
            break;
        }
    }

    return true;
}

std::string Analytiker::_to_string_formatted(double value, double precision)
{
    std::stringstream sstr;
    if (precision < 1e-8)
    {
        sstr << std::scientific << value;
    }

    else
    {
        sstr << std::fixed << std::setprecision(static_cast<int>(-std::log10(precision))) << value;
    }

    /* sstr<<value; */
    return sstr.str();
}

std::pair<std::string, std::string> Analytiker::python_cart(std::string x_var, std::string y_var, std::string numpy,
                                                            double precision) const
{
    std::string real_formula("0"), imag_formula("0");

    for (const auto& bit : m_formula_bits)
    {
        auto value = std::get<0>(bit);
        const auto& kpoint = *std::get<2>(bit).get();
        auto kcart = kpoint.cart(m_recip_lat);

        assert(almost_equal(kcart(2),0.0,precision));

        if (almost_equal(value, 0.0, precision))
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
        value_str += _to_string_formatted(value, precision);

        // Next work on the first entry within the function
        std::string dots = "(" + _to_string_formatted(kcart(0), precision) + "*" + x_var;

        // Then work on the second entry within the function
        if (kcart(1) >= 0)
        {
            dots.push_back('+');
        }
        dots += _to_string_formatted(kcart(1), precision) + "*" + y_var + ")";

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
