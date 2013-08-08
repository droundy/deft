# Deft frequently asked questions.


## How do I add a new C++ file to the deft core?

First create a new file in the `src/` directory.  Then
you edit the `SConstruct` file, and search for a similar
file and add the new file at the same place.


## How do I create a new C++ program to generate data?

We keep programs that generate data in `papers/PAPERNAME/figs/`, so
you'll want to create a `.cpp` file in that directory.  Then, as in
the case of a deft core file above, you edit `SConstruct`, and copy a
similar file.  The primary distinction you want to make is: should
your program be automatically run, or is it very slow or a huge memory
hog? Also, does it just produce a single output file, or does it
producce a few, or many?  Programs that generate data are built with a
`.mkdat` extension, which when run should generate a file with a
`.dat` extension, and the same base name.

## How do I speed up the build on srun?

You can run with multiple cores, designated by -j4

    srun -c4 scons -j4

The `-c4` tells srun that make will use 4 cores.  The `-j4` tells scons to use 4 cores.
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

Next, you need to teach `SConstruct` to generate and use the new
`src/NewFunctionalFast.cpp` file.  Edit this file, and add your new
filename, analogous to how `src/SoftFluidFast.cpp` shows up currently.
And this is all! You just need to go ahead and write code that uses
your functional.  Of course, making the haskell code implement the
functional you actually want requires more work.

I am currently transitioning to a new system for C++ code generation,
which is slightly different, and shows up in
`src/haskell/newfunctionals.hs`.
