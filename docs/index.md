---
title: "multishift guide"
---
{% include head.html %}


# What is multishift for?
A common technique for investigating dislocations, stacking faults, and slip planes in crystals involves the calculation of $$\gamma$$-surface energies.
The software provided in the `multishifter` toolbox facilitate the creation of crystal structures for calculations of UBER curves and surface energies along arbitrary slip planes.
The `multishifter` libraries also have tools aid in the construction of twisted bilayers to create periodic Moir&#233; structures.

# Basic usage
The `multishift` executable is a suite of command line tools that can manipulate and create crystal structures for slab model calculations.
Each tool is specialized to complete a particular task, and accepts arguments either directly from the command line.
Input and output files of crystal structures use the [VASP](https://www.vasp.at/wiki/index.php/POSCAR) format.

A typical use case example is the calculation of a $$\gamma$$-surface, which requires generating slab structures that have been shifted along a particular crystallographic plane. 
Starting from a primitive structure, several steps must be taken to create the slab models, each accomplished using a `multishifter` tool:
1. Slice the primitive structure to expose the desired surface plane (`slice`).
2. Specify the thickness of the slabs (`stack`).
3. Rigidly translate the basis of the slab, to expose the desired atomic layer (`translate`).
4. Generate set of shifted slab structures, one for each point on the $$\gamma$$-surface (`shift`).

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
`multishift slice` will create a supercell of a primitive cell that exposes a particular slip plane.

### Parameters
- input: path to starting structure.
- output: output file path.
- millers: defines slip plane.

## [stack](./tutorials/ii)
`multishift stack` will concatenate a list of structures together.
The stack is created by fusing unit cells together along the $$ab$$ facet.

### Parameters
- inputs: list of structures to stack together.
- output: output file path.

## [translate](./tutorials/ii)
`multishift translate` will rigidly translate the basis atoms of the given structure.

### Parameters
- input: path to starting structure.
- output: output file path.
- value: translation vector to apply to the basis.
- floor: index of atom that should end up at the origin after translating the basis.
- fractional: if given, "value" will be interpreted as a fractional value relative to the lattice vectors.

## [align](./tutorials/i)
`multishift slice` will create a supercell of a primitive cell that exposes a particular slip plane.

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
- output: output directory (must not exist).
- values: vacuum spacings to insert between the periodic images of the slabs in $$\AA$$.
Negative values will bring periodic slabs closer to each other.


# Tutorials
[Tutorial I: Slicing a structure](./tutorials/i/)<br/>
[Tutorial II: Stacking structures](./tutorials/ii/)<br/>
[Tutorial III: Stacking faults and surfaces](./tutorials/iii)<br/>
[Tutorial IV: Cleaving for UBER](./tutorials/iv)<br/>
[Tutorial V: Shifting for $$\gamma$$-surfaces](./tutorials/v)<br/>
[Tutorial VI: Chaining shifting and cleaving](./tutorials/vi)<br/>
[Tutorial VII: Fourier interpolation](./tutorials/vii)<br/>
[Tutorial VIII: Do the twist](./tutorials/viii)<br/>
[Tutorial IX: What's in a record.json?](./tutorials/ix)<br/>

