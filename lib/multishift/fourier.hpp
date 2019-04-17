#ifndef MULTISHIFTFOURIER_HH
#define MULTISHIFTFOURIER_HH

#include "./autodocs.hpp"
#include "./define.hpp"
#include "./surf.hpp"
#include "casm/crystallography/Lattice.hh"
#include <complex>
#include <vector>

namespace CASM
{
class jsonParser;
}

namespace mush
{

/**
 * Settings for interpolating values on a grid. Holds a path
 * to the data that should be interpolated.
 */

class FourierSettings
{
public:
    FourierSettings(const fs::path& init_data_path, const fs::path& init_lattice_path,
                    const std::vector<std::string>& init_value_tags);
    static FourierSettings from_json(const CASM::jsonParser& init_json);
    CASM::jsonParser to_json() const;

    const fs::path& data_path() const { return m_data_path; };
    const fs::path& lattice_path() const { return m_lattice_path; };
    const std::vector<std::string>& values_to_interpolate() const { return m_value_tags; };

    /// Information about each of the settings entries required to construct *this
    static const docs::SettingsInfo docs;

private:
    /// Path to json file that contains the data to interpolate. Each grid point
    /// is defined by its fractional coordinates, and must appear exactly once.
    fs::path m_data_path;

    /// Path to the structure whose lattice is associated with the fractional coordinates
    /// in the data file
    fs::path m_lattice_path;

    /// Fields in the data that should be interpolated
    std::vector<std::string> m_value_tags;

    /// Generate the documentation for the settings, used to construct static member
    static docs::SettingsInfo _initialized_documentation();
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

    Interpolator(const Lattice& init_lat, const std::vector<double>& a_fracs, const std::vector<double>& b_fracs,
                 const std::vector<double>& vals);

    /// Reconstruct from a serialized jsonParser
    static Interpolator deserialize(const CASM::jsonParser& serialized);

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

    /// Store values of an InterGrid for a_frac, b_frac, x_cart, y_cart, and the complex value as
    /// a json object
    static CASM::jsonParser grid_to_json(const Lattice& ref_lat, const InterGrid& grid_values);

    /// Serialize into a jsonParser for immediate reconstruction
    CASM::jsonParser serialize() const;

private:
    /// The InterGrid must have specific dimensions, which should not be determined by outside forces
    Interpolator(const Lattice& init_lat, const InterGrid& init_values);

    /// This one is for when you call deserialize, don't use it for other stuff
    Interpolator(Lattice&& init_real, Lattice&& init_recip, InterGrid&& init_rpoints, InterGrid&& init_kpoints);

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

    /// TODO: Places all the gridded up points within the weigner seitz cell
    static InterGrid _r_weigner_seitz(const InterGrid& init_values, const Lattice& real_lattice);

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
    static InterGrid _grid_from_unrolled_data(const std::vector<double>& as, const std::vector<double>& bs,
                                              const std::vector<double>& vals);

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
    std::pair<std::string,std::string> python_cart(std::string x_var, std::string y_var, std::string numpy) const;

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
    Lattice m_recip_lat;
};

} // namespace mush

#endif
