# TimeTrees

[![Build Status](https://travis-ci.org/tgvaughan/TimeTrees.jl.svg?branch=master)](https://travis-ci.org/tgvaughan/TimeTrees.jl)

A tiny package that implements the TimeTree type for representing
fully-resolved phylogenetic time trees in Julia.  A constructor which generate
`TimeTree`s from Newick strings is provided, as are methods for manipulating
existing trees. In addition, a plot method is implemented which generates an
ASCII depiction of a tree.

Here's an example of a TimeTree being generated from its Newick representation
and an ASCII visualization displayed in an interactive Julia session:

```
julia> using TimeTrees
julia> t = TimeTree(newick)
julia> plot(t)
                                              /----------------------* 1
              /-------------------------------+----------------------* 2
              |                         /----------------------------* 3
/-------------+                        /+        /-------------------* 4
|             |                        |\--------+                 /-* 5
|             \------------------------+         \-----------------+-* 6
|                                      |                   /---------* 7
|                                      |               /---+---------* 8
|                                      \---------------+/------------* 9
+                                                      \+           /* 10
|                                                       \-----------+* 11
|                                                                   \* 12
|                                     /------------------------------* 13
|      /------------------------------+   /--------------+-----------* 14
|      |                              |   |              \-----------* 15
|      |                              \---+          /--------------+* 16
\------+                                  \----------+              \* 17
       |                                             \---------------* 18
       |                                                           /-* 19
       \-----------------------------------------------------------+-* 20
```
(Assuming that the variable `newick` holds a string containing the Newick
representation.)

Trees with non-contemporaneous leaf ages are also supported:
```
julia> plot(t, labelLeaves = false)
                    /-------*
/-------------------+-----------------------------------------------*
|   /-----------------------------------------*
+   |                          /--------------------*
|   |      /-------------------+    /---------*
\---+      |                   \----+---------------+---*
    |      |                                        \---------*
    |      |                       /-------------*
    \------+-----------------------+        *
           |                       \--------+    /-+----*
           |                                |    | \------*
           |                                \----+-----------+-*
           +                                                 \------*
           |               /-----------------------------------------*
           |     /---------+     /--------------*
           |     |         \-----+--------------------------------*
           \-----+                       /--------------------------*
                 \------------+----------+----------*
                              |                            /-------*
                              \----------------------------+---------*
```

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
