# About
For a tutorial on the features and workflow see the official [pages](https://goirijo.github.io/multishifter).

`multishift` is a [casm](https://github.com/prisms-center/CASMcode) and [casm-utilities](https://github.com/goirijo/casm-utilities) powered utility to create gamma surfaces and twisted bilayers along arbitrary crystallographic planes.

This work is a product of the [Van der Ven research group](https://labs.materials.ucsb.edu/vanderven/anton/) at the UCSB Materials department.
The repository is written and maintained by me, John Goiri.

# Citing
If this software helped your research please cite this repository and the following paper in addition to the papers listed in the [casm](https://github.com/prisms-center/CASMcode) repository.

```
Goiri, Jon Gabriel, and Anton Van der Ven. "MultiShifter: software to generate structural models of extended two-dimensional defects in 3D and 2D crystals." arXiv preprint arXiv:2011.10734 (2020).
```

```bibtex
@misc{goiri2020multishifter,
      title={MultiShifter: software to generate structural models of extended two-dimensional defects in 3D and 2D crystals}, 
      author={Jon Gabriel Goiri and Anton Van der Ven},
      year={2020},
      eprint={2011.10734},
      archivePrefix={arXiv},
      primaryClass={cond-mat.mtrl-sci}
}
```

The [paper](link/to/paper) outlines how the algorithms of `multishifter` operate, and applies them to metallic systems.
If you like this repository, consider giving me a star!

# Prerequisites
`multishift` is written on top of [casm-utilities](https://github.com/goirijo/casm-utilities), which requires a c++ compiler with c++17 support (`g++-9` recommended).

# Installation
You can install `multishift` by compiling from a tarball or zipball, or you can install it as a [casm-utilities](https://github.com/goirijo/casm-utilities) plugin.

## Installing from tarball or zipball
Head over to the [releases](https://github.com/goirijo/multishifter/releases) page and download the `tar.gz` or `zip` file of the latest version (*not* the "Source code" files, you want the attachment with the extension in the name).

Once you've downloaded it, unpack it and run

```bash
./configure
make
make install
```
    
The build involves compiling some relatively heavy libraries, so it may take a few minutes to finish. You can speed up the build running `make -j4`, which will parallize the compilation on 4 processors.
Don't forget to specify a compiler with c++17 support if it's no the default on your system!
You can specify any compiler flags you want during the `configure` step as well.

## Installing as a casm-utilities plugin
Documentation on the [`casm-utilities`](https://github.com/goirijo/casm-utilities) repository explains how to install plugins.
You'll have to clone the `casm-utilities` repository, and then this repository under the `plugins` directory.
Once this is done, compiling `casm-utilities` will also insall `multishifter`.

# Acknowledgements
Initial development of this software package was possible thanks to the financial support from NSF DMREF program under grant DMR-1534264.
Subsequent development was made possible by additional grants.
The scientific work was supported by the National Science Foundation DMREF grant: DMR-1729166 “DMREF/GOALI: Integrated ComputationalFramework for Designing Dynamically Controlled Alloy-Oxide Heterostructures”.
Software development was supported by the National Science Foundation, Grant No. OAC-1642433.
Computational resources provided by the National Energy Research Scientific Computing Center (NERSC), supported by the Office of Science and US Department of Energy under Contract No. DE-AC02-05CH11231, are gratefully acknowledged, in addition to support from the Center for Scientific Computing from the CNSI, MRL, and NSF MRSEC (No. DMR-1720256).

