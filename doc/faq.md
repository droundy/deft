---
title: Deft frequently asked questions.
layout: page
---


## How do I add a new C++ file to the deft core?

First create a new file in the `src/` directory.  Then
you edit the `Makefile.am` file, and search for a similar
file and make similar additions everywhere it shows up.


## How do I create a new C++ program to generate data?

We keep programs that generate data in `papers/PAPERNAME/figs/`,
so you'll want to create a `.cpp` file in that directory.  Then,
as in the case of a deft core file above, you edit `Makefile.am`,
and copy a similar file, making similar additions where ever it
shows up.  If you want your program to be automatically run when
building papers (which you should, if it runs quickly), then you
need to be sure to add a `.dat` dependency, which will tell make
to run your program.  Programs that generate data are built with
a `.mkdat` extension, which when run should generate a file with
a `.dat` extension, and the same base name.

## How do I speed up the make on srun?

You can run with multiple cores, designated by -j4

    srun -p debian -c4 make papers -j4

The `-c4` tells srun that make will use 4 cores.  The `-j4` tells make to use 4 cores.
If the srun job does not start immediately you will want to reduce the number of cores used.
The idea being that all the processes will run on a single node.

## How do I create a new density functional?

This requires a few steps, and is most easily done by copying and
modifying an existing functional.  You can begin by creating a new
haskell module for your functional by copying `src/haskell/SFMT.hs` to
a new `NewFunctional.hs` file for your new functional.  Edit that file
to rename SFMT at the top to your new filename (which should not
actually be NewFunctional, but should *start* with a capital letter).

Now you need to edit `src/haskell/functionals.hs` to make that program
use your functional to generate a new `src/NewFunctionalFast.cpp` file.

Next, you need to teach `Makefile.am` to generate and use the new
`src/NewFunctionalFast.cpp` file.  Edit this file, and add your new
filename, analogous to how `src/SoftFluidFast.cpp` shows up currently.
And this is all! You just need to go ahead and write code that uses
your functional.  Of course, making the haskell code implement the
functional you actually want requires more work.

I am currently transitioning to a new system for C++ code generation,
which is slightly different, and shows up in
`src/haskell/newfunctionals.hs`.
