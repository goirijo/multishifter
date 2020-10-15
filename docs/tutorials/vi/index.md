---
title: "Tutorial VI: Chaining Shifting and Cleaving"
---
{% include head.html %}

<style type="text/css">
{% include warning.css %}
</style>

# Tutorial VI: Chaining shifting and cleaving
In this tutorial we'll use `multishift chain`, to multiply the output of `shift` and `cleave`.
We've already used [`multishift cleave`](../iv) to insert empty space between slabs, and we've used [`multishift shift`](../v) to slide the slabs parallel to each other.
What if we want to apply the `cleave` command to each structure generated by `shift` individually?


## Mg slab cell
We'll start with a slab of $$\mathrm{Mg}$$ with 8 atomic layers (4 unit cells tall).
This is the same slab from the [`mutate` tutorial](../iii).
You can download it [here]("./mg_stack4.vasp), or create a file called `mg_stack4.vasp` yourself with the follwing data:

    Mg stack
    1.00000000
       1.59609453    2.76451683    0.00000000
      -1.59609453    2.76451683    0.00000000
       0.00000000    0.00000000   20.73607828
    Mg
    8
    Direct
       0.66666670    0.66666670    0.18750000 Mg
       0.33333330    0.33333330    0.06250000 Mg
       0.66666670    0.66666670    0.43750000 Mg
       0.33333330    0.33333330    0.31250000 Mg
       0.66666670    0.66666670    0.68750000 Mg
       0.33333330    0.33333330    0.56250000 Mg
       0.66666670    0.66666670    0.93750000 Mg
       0.33333330    0.33333330    0.81250000 Mg

You can also create your own slab from a primitive cell using `multishift stack`, as explained in a [previous tutorial](../ii).

## Shift AND cleave your slab
The chain command is very similar to `shift` and `cleave`.
Suppose we want to shift the $$\mathrm{Mg}$$ slab the same way we did in the [`shift tutorial`](../v), but we also want to apply the spacings from the [`cleave tutorial`](../iv).
We call:

```bash
multishift chain --input mg_stack4.vasp --shift 3 3 --cleave -1.5 0.0 2.0 --output mg_chain
```

Once you've run this command, a similar direcotry structure appears, this time 2 levels deep:

```bash
mg_chain/
├── record.json
├── slab.vasp
├── shift__0.0
│   ├── cleave__0.000000
│   │   └── POSCAR
│   ├── cleave__-1.500000
│   │   └── POSCAR
│   └── cleave__2.000000
│       └── POSCAR
├── shift__0.1
│   ├── cleave__0.000000
│   │   └── POSCAR
│   ├── cleave__-1.500000
│   │   └── POSCAR
│   └── cleave__2.000000
│       └── POSCAR
├── ...
└── shift__2.2
    ├── cleave__0.000000
    │   └── POSCAR
    ├── cleave__-1.500000
    │   └── POSCAR
    └── cleave__2.000000
        └── POSCAR
```

It's as if the `shift` command had been run first, then for each generated structure, the `cleave` command was applied.