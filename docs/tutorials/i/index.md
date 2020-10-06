---
title: "Tutorial I: Slicing a Structure"
---
{% include head.html %}

# Tutorial I: Slicing a structure
In this tutorial we'll use `multishift slice` to create a structure of $$\mathrm{Mg}$$ with ab-vectors in the $$(1,0,1)$$ plane (pyramidal slip plane).

## Mg primitive cell
Current implementations of `multishifter` work exclusively with the [VASP](https://cms.mpi.univie.ac.at/vasp/guide/node59.html) structure format.
Download the primitive cell for $$\mathrm{Mg}$$ [here](./mg.vasp), or create a file called `mg.vasp` with the following crystallographic data:

    multishifter tutorial i
    1.00000000
                           1.59609453       2.76451683       0.00000000
                          -1.59609453       2.76451683       0.00000000
                           0.00000000       0.00000000       5.18401957
    Mg
    2
    Direct
                   0.6666667    0.6666667    0.7500000 Mg
                   0.3333333    0.3333333    0.2500000 Mg

The figure below on the left shows the primitive unit cell of $$\mathrm{Mg}$$, while the one on the right highlights the prismatic plane.
This is the slip plane we're interested in, and has miller indices $$(0,\bar{1},1)$$.

| ![Primitive structure of Mg](./mg.png){:width="70%"} | ![Pyramidal slip plane in HCP](./pyraslip.svg){: width="70%"} |

## Run "multishift slice"
From your working directory, run

```
multishift slice --input mg.vasp --millers 0 -1 1 --output mg_pyramidal.vasp
```

This will read `mg.vasp` as the starting structure, and slice it allong the specified miller indexes, saving the generated structure to `mg_pyramidal.vasp`.

## Output
You will now have a file named `mg_pyramidal.vasp`.
The plane defined by the `millers` argument is spanned by the $$ab$$ vectors of this structure.
By default, the `slice` command also performs a rigid rotation on the final structure, aligning the $$a$$ vector along the $$x$$ direction, and the $$b$$ vector in the $$xy$$ plane.
This final rotation can be supressed by using the flag `--dont-align`.
The rotation can also be applied to any structure using `multishift align`.

| ![Sliced Mg along 0-11](./mg_pyramidal.png){:width="100%"} | ![Pyramidal slip plane in HCP](./pyraslip_unit.svg){:width="70%"} |
