# About
`mulstishift` is a [casm](https://github.com/prisms-center/CASMcode) powered utility to create gamma surfaces along arbitrary crystallographic planes.
This repository remains under heavy development.
Documentation and features are being continually added.

# Prerequisites
`multishift` requires a c++ compiler with c++11 support. Because it relies heavily on [casm](https://github.com/prisms-center/CASMcode), you will also need to have [boost](https://www.boost.org/doc/libs/1_66_0/more/getting_started/unix-variants.html) installed on your system.
The easiest way to do this is via your package manager, such as `apt` or `brew`, though you may choose to install it yourself.
Specifically, the boost dependencies are:
- boost system
- boost filesystem
- boost program options
- boost chrono
- boost regex

Future versions of `multishift` will reduce the dependencies as much as possible.
Ubuntu users can install all the dependencies with
```
apt-get install libboost-all-dev
```

# Installation
You can install `multishift` by compiling from a tarball or zipball, or by cloning the git repository.
## Installing from tarrball or zipball (recommended)
## Installing from git repository

# Usage
For a more or less complete tutorial on the features and workflow see the official [pages](https://goirijo.github.io/multishifter).
