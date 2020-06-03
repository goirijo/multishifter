#ifndef FOURIER_HH
#define FOURIER_HH

#include "./definitions.hpp"
#include <complex>
#include <casmutils/xtal/lattice.hpp>

namespace mush
{
/**
 * One point on the surface that should be used for the interpolation.
 * Only fractional coordinates needed.
 */

class InterPoint
{
public:
    InterPoint(double init_a_frac, double init_b_frac, std::complex<double> init_value)
        : a_frac(init_a_frac), b_frac(init_b_frac), value(init_value), weight(1.0)
    {
    }

    double a_frac;
    double b_frac;
    std::complex<double> value;
    double weight;

    /// Return the cartesian coordinate relative to the given lattice. Position will land on the ab-plane, since there is no c-vector component.
    /// This could probably just return a 2d-vector, since your lattice is probably oriented such that ab-vectors are on xy-plane, but
    /// just in case it isn't, it'll stay as a 3d-vector;
    Eigen::Vector3d cart(const cu::xtal::Lattice& ref_lat) const;
    
    /// Priority is a fraction, then b fraction. Value is not considered.
    bool operator<(const InterPoint& rhs) const;

private:
};

/**
 * Given a lattice (only ab-vectors matter) and a list of values
 * and grid points, creates a reciprocal space from which to sample
 * plane waves uniformly, and calculates coefficients to
 * reproduce the given data.
 *
 * The reciprocal grid to sample in is determined by the shape of
 * the 2d arrangement of the interpolated values.
 */

class Interpolator
{
public:
    typedef std::vector<std::vector<InterPoint>> InterGrid;
    typedef mush::cu::xtal::Lattice Lattice;

    Interpolator(const Lattice& init_lat, const std::vector<InterPoint>& real_data);

    const InterGrid& sampled_values() const { return m_real_ipoints; }

    const InterGrid& k_values() const { return m_k_values; }

    const Lattice& real_lattice() const { return m_real_lat; }

    const Lattice& reciprocal_lattice() const { return m_recip_lat; }

    /// Number of k-points/data points
    int size() const;

    /// Number of k-points/data points along each direction
    std::pair<int, int> dims() const;

    /// Use the Fourier basis to reconstruct the signal at an arbitrary resolution
    std::pair<Lattice, InterGrid> interpolate(int a_dim, int b_dim) const;

private:
    /// The InterGrid must have specific dimensions, which should not be determined by outside forces
    Interpolator(const Lattice& init_lat, const InterGrid& init_values);

    /// This one is for when you call deserialize, don't use it for other stuff
    /* Interpolator(Lattice&& init_real, Lattice&& init_recip, InterGrid&& init_rpoints, InterGrid&& init_kpoints); */

    /// The real lattice where the gamma surface happened on the ab-plane
    /// This is not the true lattice of the structure, just one with the
    /// correct ab-vectors oriented along the xy-plane
    Lattice m_real_lat;

    /// The reciprocal lattice of the reduced real lattice
    Lattice m_recip_lat;

    /// Holds initial values at the real interpolation points
    InterGrid m_real_ipoints;

    /// Holds coefficients at each of the reciprocal k points
    InterGrid m_k_values;

    /// Constructs grid of k points, setting all coefficients to zero
    /// These k-points are integer multiples of the reciprocal lattice, because all
    /// the information of a gamma surface falls within a single real space unit
    /// cell (think about it, it's like the opposite of phonons, where you want 1st
    /// Brillouin zone)
    static InterGrid _k_grid(const InterGrid& init_values, const Lattice& reciprocal_lattice);

    /// Sets the coefficients for each of the k points
    void _take_fourier_transform();

    /// Helper to avoid so many loops. Unrolls 2d vector into a 1d vector.
    static std::vector<InterPoint*> _unrolled_values(InterGrid* grid_values);

    /// Helper to avoid so many loops. Unrolls 2d vector into a 1d vector, but pointers are const.
    static std::vector<const InterPoint*> _unrolled_values(const InterGrid& grid_values);

    /// Given an unrolled vector of grid data, reshape it to the specified dimensions
    static InterGrid _direct_reshape(const std::vector<mush::InterPoint>& unrolled_data, int ka_dim, int kb_dim);

    /// Given unrolled data for ab-grid fractions and values, reshape the data into a 2d-grid,
    /// making a few checks along the way to make sure it's even possible, given the data.
    /// The resulting grid always has odd dimensions. In the case of an even input, values at the
    /// periodic boundary are duplicated on the opposite edge in order for the Fourier transform
    /// to properly work.
    static InterGrid _grid_from_unrolled_data(std::vector<InterPoint> unrolled_data);

    /// Given unrolled data for ab-grid fractions and values, append additional values such that
    /// reshaping the grid results in odd dimensions along both a and b.
    /// Even dimensions are made odd by repeating the boundary values on the opposite side.
    static void _make_unrolled_data_odd(std::vector<mush::InterPoint>* unrolled_data, int* final_adim, int* final_bdim);
};

/**
 * Given an interpolator object, this class will express the plane waves of the
 * Fourier transform as an analytical formula
 */

class Analytiker
{
public:

    ///Used in the FormulaBit to assiciate a coefficient with the appropriate basis function
    /// Individual sin/cos basis in the following order:
    /// re(cos), im(sin), im(cos), re(sin), InterPoint
    enum class FormulaBitBasis
    {
        RECOS,
        IMSIN,
        IMCOS,
        RESIN,
    };

    /// Couples an InterPoint (k-point coefficient) with decomposed coefficients
    typedef std::tuple<double, FormulaBitBasis, std::shared_ptr<const InterPoint>> FormulaBit;

    /// Initialize with an interpolator
    Analytiker(const Interpolator& init_ipolator);

    /// Print the analytical expression in a Python compatible format, where the
    /// expected input are numpy like arrays for the x and y coordinates
    /// The first entry in the pair is made up of the real basis functions, while the
    /// second one has the imaginary ones (which are probably zero)
    std::pair<std::string,std::string> python_cart(std::string x_var, std::string y_var, std::string numpy, double precision) const;

private:
    /// Checks that the imaginary coefficients cancelled out when creating
    /// the formula bits
    bool _can_ignore_imaginary() const;

    /// Run through an entire k-point grid, and express the complex coefficients
    /// into individual coefficients for sin/cos.
    /// All the imaginary sin/cos terms will presumably be zero, but this routine does't
    /// check for that. The number of formula bits will be half of the k-point grid,
    /// since the terms can be halved by taking inversion symmetry of sin/cos into account.
    static std::vector<FormulaBit> _formula_bits(const Interpolator::InterGrid& k_values);

    /// Contains coefficients for each fo the basis functions
    std::vector<FormulaBit> m_formula_bits;

    /// Reciprocal lattice of the interpolator that *this was constructed with
    Interpolator::Lattice m_recip_lat;

    /// Turn a double into a string. If the precision is small enough, use scientific notation
    static std::string _to_string_formatted(double value, double precision);
};
}

#endif
