---
title: "multishift guide"
---
{% include head.html %}

<style type="text/css">
{% include warning.css %}
</style>

<div class="error">
<b>THIS SITE IS UNDER CONSTRUCTION:</b>
<br>
Documentation is being uploaded before features are available.
<br>
</div>
<div>
<br>
</div>

# What is multishifter for?
A common technique for investigating dislocations, stacking faults, and slip planes in crystals involves the calculation of $$\gamma$$-surface energies.
The software provided in the `multishifter` toolbox facilitate the creation of crystal structures for calculations of UBER curves and surface energies along arbitrary slip planes.
The `multishifter` libraries also have tools that aid in the construction of twisted bilayers to create periodic Moir&#233; structures.

# Basic usage
The `multishift` executable is a suite of command line tools that can manipulate and create crystal structures for slab model calculations.
Each tool is specialized to complete a particular task, and accepts arguments directly from the command line.
Input and output files of crystal structures use the [VASP](https://www.vasp.at/wiki/index.php/POSCAR) format.

Starting from a primitive structure, several steps must be taken to create the slab models, each accomplished using a `multishifter` tool.
A typical workflow looks like
1. Slice the primitive structure to expose the desired surface plane (`slice`).
2. Specify the thickness of the slabs (`stack`).
3. Rigidly translate the basis of the slab, to expose the desired atomic layer/center the rotation axis (`translate`)
4. Generate a set of shifted slab structures with `shift` ($$\gamma$$-surface), cleave the structure at different intervals with `cleave` (UBER), or create commensurate structures that have been rotated with `twist` (Moir&#233;/twisted layers).

Each of these commands takes the form
```bash
multishift <command> --args
```
All available tools and a list on their input parameters can be found with
```bash
multishift --help
multhishift <command> --help
```

# Feature overview
`multishift` comes as a collection of tools, each meant to manipulate crystal structures a particular way.

## [slice](./tutorials/i)
`multishift slice` will create a supercell of a primitive cell that exposes a particular slip or surface plane.

### Parameters
- input: path to starting structure.
- output: output file path.
- millers: defines slip plane.

## [stack](./tutorials/ii)
`multishift stack` will concatenate a list of structures together.
The stack is created by fusing unit cells together along the $$ab$$ facet.

### Parameters
- input: list of structures to stack together.
- output: output file path.

## [translate](./tutorials/ii)
`multishift translate` will rigidly translate all the basis atoms of the given structure.

### Parameters
- input: path to starting structure.
- output: output file path.
- value: translation vector to apply to the basis.
- floor: index of atom that should end up at the origin after translating the basis.
- fractional: if given, "value" will be interpreted as a fractional value relative to the lattice vectors.

## [align](./tutorials/i)
`multishift align` will apply a rigid rotation to the structure, aligning the $a$ vector along the $x$ direction, and restricting the $b$ vector to the $xy$-plane.

### Parameters
- input: path to starting structure.
- output: output file path.
- prismatic: if given, the $$c$$ vector will also be rectified to be perpendicular to the $$ab$$ plane, creating a prismatic slab. This may break the periodicity of your crystal.

## [mutate](./tutorials/iii)
`multishift mutate` will alter the $$c$$ vector of the lattice, while retaining all Cartesian values of the basis.
Useful for individual stacking fault and surface energy calculations.

### Parameters
- input: path to starting structure.
- output: output file path.
- mutation: value to add to the third lattice vector.
- fractional: if given, "mutation" will be interpreted as a fractional value relative to the lattice vectors.

## [cleave](./tutorials/iv)
`multishift cleave` will generate a series of slab structures with a specified amount of vacuum between their periodic images.
Useful for Universal Binding Energy Relation (UBER) calculations.

### Parameters
- input: path to starting structure.
- output: output directory.
- values: vacuum spacings to insert between the periodic images of the slabs in $$\AA$$.
Negative values will bring periodic slabs closer to each other.

## [shift](./tutorials/v)
`multishift shift` will generate a series of slab structures that have been shifted relative to each other along the $$ab$$-plane.
Useful for $$\gamma$$-surface calculations.

### Parameters
- input: path to starting structure.
- output: output directory.
- grid: divisions to make along the $$ab$$-plane.
Each grid point will correspond to a particular shift applied to the slabs.
Periodic images are not counted in the grid.

## [chain](./tutorials/vi)
`multishift chain` combines the `shift` and `chain` commands.
For each perturbed slab generated by `shift`, all the values of `cleave` are applied to it.
Useful for $$\gamma$$-surface calculations.

### Parameters
- input: path to starting structure.
- output: output directory.
- grid: divisions to make along the $$ab$$-plane.
Each grid point will correspond to a particular shift applied to the slabs.
Periodic images are not counted in the grid.

## [fourier](./tutorials/vii)
`multishift fourier` reads DFT or UBER parameters that you've assigned to each gridpoint from a `chain` or `shift` command, and returns an analytical expression for your surface.
The expression is printed as python code.

### Parameters
- data: path to `json` file that holds the data to interpolate.
- cleavage-slice: specifies which cleavage values to read from the data set. Use this when your data set has multiple values per grid point.
- key: string used to access the values to interpolate in your `json` file.
- crush: minimum magnitude required for basis function to be included.

## [twist](./tutorials/viii)
`multishift twist` creates twisted commensurate supercells that can be combined to create Moir&#233; patterns when stacked together.

### Parameters
- input: path to starting structure.
- output: output directory.
- angles: rotation angles to apply to the slab. Rotation axis is always perpedicular to the $$ab$$-plane, and goes through the origin.
- max-lattice-sites: determines how large the search space for highly commensurate supercells should be. Larger values will result in less deformation of the final layers.
- error-tol: minimum improvement necessary to consider a larger supcercell [better](./tutorials/ix).
- brillouin-zone: select whether the rotated or the original (aligned) Brillouin zone should be used to map reciprocal vectors back into the first zone.
- supercells: determines whether only the best supercell, or supercells of different sizes should be outputted. 


# Tutorials
Follow this series of tutorials to acquaint yourself with all of the available features.
Each tutorial focuses on one or two of the available tools, and build on each other as they progress.

[Tutorial I: Slicing a structure](./tutorials/i/)<br/>
[Tutorial II: Stacking structures](./tutorials/ii/)<br/>
[Tutorial III: Stacking faults and surfaces](./tutorials/iii)<br/>
[Tutorial IV: Cleaving for UBER](./tutorials/iv)<br/>
[Tutorial V: Shifting for $$\gamma$$-surfaces](./tutorials/v)<br/>
[Tutorial VI: Chaining shifting and cleaving](./tutorials/vi)<br/>
[Tutorial VII: Fourier interpolation](./tutorials/vii)<br/>
[Tutorial VIII: Do the twist](./tutorials/viii)<br/>
[Tutorial IX: What's in a record.json?](./tutorials/ix)<br/>

