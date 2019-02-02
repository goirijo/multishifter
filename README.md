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
## Installing from tarball or zipball (recommended)
Head over to the [releases](https://github.com/goirijo/multishifter/releases) page and download the `tar.gz` or `zip` file of the latest version (*not* the "Source code" files, you want the attachment with the extension in the name).

Once you've downloaded it, unpack it and run

    ./configure
    make
    make install
    
The build involves compiling some relatively heavy libraries, so it may take a few minutes to finish. You can speed up the build running `make -j4`, which will parallize the compilation on 4 processors.

## Installing from git repository
If you need to use features which are on the repository but have not been released yet, or want to compile from the git repository for whatever reason, you can clone whichever branch you're interested in, go inside, and run

    bash bootstrap.sh
    ./configure
    make
    make install
    
The `bootstrap.sh` script introduces additional dependencies, as it uses autotools to create the `configure` file. In order for this to work you'll need to have the following installed:

- autoconf
- automake
- libtool
- autoconf-archive

# Usage
For a more or less complete tutorial on the features and workflow see the official [pages](https://goirijo.github.io/multishifter).

Acknowledgements
