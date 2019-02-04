---
title: multishift guide
---
{% include head.html %}


# What is multishift for?
A common technique for investigating dislocations, stacking faults, and slip planes in crystals involves the calculation of $$\gamma$$-surface energies.
The software provided in in the `multishifter` toolbox facilitate the creation of crystal structures for calculations of UBER curves and surface energies along arbitrary slip planes.

# Basic usage
Most `multishift` commands will take a settings file as input in a `json` format, where you will specify how to construct slabs, shift, interpolate, etc.
In order to keep data structured, a single settings file can be used for any of the `multishift` tools that requires one.
The settings file must contain a `name`, which will be used as a prefix for generated directories, and entries for each of the individual `multishift` tools you'll be calling.

For example, a settings file that will be used for `multishift-base` and `multishift-shift` might look like:

    {
        "name" : "Mg-pyramidal1"
        "base" :
        {
            "prim" : "Mg.vasp",
            "millers" : [1,0,1],
            "stacks" : 6
        },
        "shift" :
        {
            "slab" : "Mg-pyramidal1.vasp",
            "a" : 19,
            "b" : 21
        }
    }

# Feature overview
`multishift` comes as a collection of tools, each meant to handle different stages of $$\gamma$$-surface calculations.

## multishift-base
The `multishift-base` tool will create a supercell of a primitive cell that exposes a particular slip plane.
Energy surface calculations involve creating large slabs along a particular slip plane, and shifting the slabs past each other.
These slab structures are supercells of the primitive structure that have been constructed such that the ab-vectors of the slab lattice span the slip plane, while the c-vector is arbitrarily large and determines the thickness of the slab.

### Settings values
Use **base** as the settings key. Options are:
- prim: Path to the primitive structure (VASP format).
- millers: Array of 3 integers that specify the slip plane.
- slab_floor_index: The entire basis of the slab structure will be translated such that the atom with this index ends up on the very bottom (0 means no translation).
- stacks: Integer specifying how thick the slab should be, which is constructed by stacking the smallest possible unit that exposes the slip plane along the c-direction.

## multishift-shift
`multishift-shift` can create structures for both UBER curves and $$\gamma$$-surfaces.
Given a structure for a slab, this tool will modify by either changing amount of space perpendicular to the slip plane (for UBER calculations), shifting the periodic images of the slab along the slip plane ($$\gamma$$-surface), or both.

### Settings values
Use **shift** as the settings key. Options are:
- slab: Path to the slab structure that you want to create UBER or $$\gamma$$-surface calculations for.
- a: Number of grid points along the a-vector for the $$\gamma$$-surface
- b: Number of grid points along the b-vector for the $$\gamma$$-surface
- cleavage: List of values in $$\AA$$ of space that should be inserted between the slabs.

# Tutorials
[Tutorial I: Creating a slab structure](./tutorials/i/)<br/>
[Tutorial II: UBER of a layered structure](./tutorials/ii/)<br/>
[Tutorial III: A $$\gamma$$-surface for Mg]()<br/>
[Tutorial IV: A better $$\gamma$$-surface for Mg]()<br/>
