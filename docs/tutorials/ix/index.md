---
title: "Tutorial IX: What's in a record.json"
---
{% include head.html %}

<style type="text/css">
{% include warning.css %}
</style>

# Tutorial IX: What's in a `record.json?`
Certain commands available in `multishifter` generate many structures at once.
For these commands, a file named `record.json` is saved to the output directory.
This records contain information about each of the generated structures, which you will find useful.
There are two types of records:
records generated for `cleave`, `shift` or `chain`, and records generated for `twist`.

## `cleave` | `shift` | `chain`
There are 5 entries at the top level:
* cleavages
* grid
* shift_units
* ids
* equivalents

### "cleavages" and "grid"
These are simply the parameters when you called the relevant `multishift` command.
For `cleave` and `shift`, you will notice that default values have been inserted.
Indeed, these two commands are actually just `multishift chain` in disguise.
The directory layout will change slightly, but using `cleave` is equivalent to using `chain` with a $$1\times 1$$ grid, and using `shift` is like using `chain` with a single null cleavage value.

The cleavage and grid values used to generate all the slab structures are saved in these entries.
You should recognize their values from the values you passed to the command line.

### "shift_units"
When shifting structures, the $$ab$$ glide plane is divided into a regular grid.
The two vectors in this entry are integer divisions of the $$a$$ and $$b$$ vectors, determined by the grid values.
These vectors are 2-dimensional, and can be used to get the Cartesian representation of each grid point.

### "ids"
Each generated slab is assigned a unique ID in the form of "a:b:cleave".
"a" and "b" are always integers, and describe the position of the structure in the grid.
You can get the Cartesian position on the grid by multiplying "a" and "b" with the shift units.
The "cleave" values is the spacing added (or removed, if negative) between the slab images.
For example, the ID "0:4:5.000000" indicates that the slab has been shifted purely along the $$b$$ vector by 4 grid points, and has had $$5\\A$$ of vacuum inserted between the periodic images.

Eeach entry in the "ids" section uses one of these IDs as a key, and provides information about the slab structure:
* cleavage: the amount of spacing inserted between the periodic imagies of this particular structure.
* grid_point: integer representation of the point on the $$ab$$ plane grid divison.
* shift: the in plane shift that was applied to the input slab to get this structure.
* directory: the relative path where you'll find a structrue file in VASP format.
* equivalent_structures: a list of IDs that are symmetrically equivalent to this particular structure. You can expect degenerate energies if you calculate all of these.
* orbit: an index given to the group of structures in "equivalent_structures"

### "equivalents"
When shifting structures, groups of shifts often result in symmetrically equivalent structures getting generated.
The orbits of the structures are listed here, with each entry containing all the structures that are symmetrically equivalent to each other.
The more you can shrink the number of orbits by proper grid selection, the DFT fewer calculations you'll have to do.

## `twist`
In the twist reports, there are 4 entries at the top level:
* angles
* error_tolerance
* max_lattice_sites
* ids

### "angles", "error_tolerance" and "max_lattice_sites"
These are simply the parameters you passed through the command line that dictated how to find the Moir&#233; lattice.
Note that if you didn't pass values, the defaults will be written out.

### "ids"
Much like in the previous section, twisted bilayers are identifed by a unique ID.
The ID takes the form "ttwist:M" (always beginning with a t).
"twist" is the rotation angle, while "M" denotes the number of moirons present in a single unit cell of stacked layers.
For example, "t15.178179:3" indicates that the the slab was rotated by 15.178179 degrees, and that the final structures will contain 3 moirons (Moir&#233; lattice sites) when stacked together.

Every entry will contain information for two structures: "top" and "bottom".
The final twisted bilayer is constructed by stacking these two layers together.
These structures are supercells of the slab which, depending on the twist angle, may have undergone a small amount of deformation to achieve periodicity.

We'll start with the entries that are not related to the deformation:
* brillouin_zone: the brillouin zone used to determine the Moir&#233; lattice vectors, as explained in the [`twist` tutorial](../viii). This is determined at the CLI.
* layer: path to the VASP formatted structure of the twisted layer
* tile: path the VASP fromatted structure of the adjusted slab. A multiple of the tile is used to construct the final layer. For special angles, the tiles are either identical to the slab, or rigidly rotated by the twist angle.

#### Deformation adjustment metrics
When the Moir&#233; supercell is not perfectly coincident with the input (aligned) slab and the rotated slab, a deformation is applied to make commensurate with each other.
This deformation is referred to as $$\mathbf{F}$$.

The deformation gradient $$\mathbf{F}$$ can be factored into a product of a symmetric stretch matrix $$\mathbf{U}$$ and a rotation matrix $$\mathbf{R}$$ according to $$\mathbf{F}=\mathbf{R}\mathbf{U}$$
The stretch matrix $$\mathbf{U}$$ describes the deformation of the crystal, while the rotation matrix $$\mathbf{R}$$ in this situation corresponds to a rotation around lattice vectors $$\mathbf{A}\times \mathbf{B}$$

A single rotation can be easily extracted from $$\mathbf{R}$$, because slabs will always get [realigned](../i/), and the rotation will be purely about the $$z$$-axis.
To analyze the strain, we use the Biot strain defined as $$\mathbf{E}=\mathbf{U}-\mathbf{I}$$ where $$\mathbf{I}$$ is the identity matrix.
The strain is restricted to the two-dimensional space parallel to the $$xy$$ twist plane.
Convenient metrics of two-dimensional strain are

$$
\begin{equation}
    \eta_1=\frac{1}{\sqrt{2}}\left(E_{11}+E_{22}\right)
\end{equation}
$$

$$
\begin{equation}
    \eta_2=\frac{1}{\sqrt{2}}\left(E_{11}-E_{22}\right)
\end{equation}
$$

$$
\begin{equation}
    \eta_3=\sqrt{2}E_{12}
\end{equation}
$$

The first strain order parameter $$\eta_1$$ is a measure of symmetry preserving in-plane dilation.
The two remaining strain order parameters, $$\eta_2$$ and $$\eta_3$$, measure deviatoric deformations within the plane.
The value of these strain order parameters is what is used to determine what the "best" Moir&#233; supercell is.
Specifically, the supercell with the lowest value of $$\sqrt{\eta_1^{2}+\eta_2^{2}+\eta_3^{2}}$$, is considered the most commensurate.

The definitions above give meaning to the remaining entries:

* deformation_matrix: the deformation matrix $$\mathbf{F}$$
* rotation: the rotation matrix $$\mathbf{R}$$
* strain: the strain matrix $$\mathbf{U}$$
* angle_adjustment: the rotation angle extracted from $$\mathbf{R}$$.
* adjustment_strain_metrics: values for $$\eta_1$$, $$\eta_2$$, and $$\eta_1$$,
* dilation_adjustment:  dilation metric defined as $$\eta_1$$
* deviatoric_adjustment: deviatoric metric defined as $$\sqrt{\eta_2^{2}+\eta_3^{2}}$$

<div class="warning">
<b>WARNING:</b>
<br>
All matrix values were calculated using the a convention of <b>column</b> lattice vectors.
<br>
</div>
<div>
<br>
</div>
