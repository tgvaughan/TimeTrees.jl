# TimeTrees

A tiny package that implements the TimeTree type for representing
fully-resolved phylogenetic time trees in Julia.  A constructor which generate
`TimeTree`s from Newick strings is provided, as are methods for manipulating
existing trees. In addition, a plot method is implemented which generates an
ASCII art depiction of a tree.

## Installation

TimeTrees is not yet a registered Julia package, so you'll need to install
it directly from the Github repository:

    Pkg.clone("http://github.com/tgvaughan/TimeTrees.jl")

## Documentation

Documentation is available through Julia's built-in help system.  To get started,
enter the following once the package is installed:

    using TimeTrees
    ?TimeTrees

For license information, see the LICENSE.md file in this directory.
