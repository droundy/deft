# Deft frequently asked questions.

[TOC]

## How do I add a new DFT C++ file to the deft core?

First create a new file in the `src/` directory.  Then you edit the
`fac/dft.py` file, and search for a similar file and add the new file
at the same place.

## How do I create a new Monte Carlo program?

To create a new Monte Carlo program, you want to copy an existing
program that is located in `src/MonteCarlo/`, but give it a new name.
Then you will add your new program to the list at the top of
`fac/monte-carlo.py`.  You should remember to `git add` your new
program, but if you forget then `fac` should remind you.


## How do I create a new C++ program to generate data?

We keep programs that generate data in `papers/PAPERNAME/figs/`, so
you'll want to create a `.cpp` file in that directory.  Then, as in
the case of a deft core file above, you edit `fac/dft.py`, and copy a
similar file.  The primary distinction you want to make is which files
it needs to be linked with.  Also you want to ask yourself whether the
program will create a single output file (and accept no arguments) or
whether it is something that should be given arguments to tell it what
to do.  Programs that generate data are always built (by convention)
with a `.mkdat` extension.  If they generate a single output (and have
no arguments), then when run they should generate a file with a `.dat`
extension, and the same base name.

If your program accepts arguments, then getting it run by fac is more
tricky, and requires putting some python code into a special comment
that fac reads.  You'll want to copy an existing example.

## How do I create a new density functional?

This requires a few steps, and is most easily done by copying and
modifying an existing functional.  You can begin by creating a new
haskell module for your functional by copying `src/haskell/SFMT.hs` to
a new `NewFunctional.hs` file for your new functional.  Edit that file
to rename SFMT at the top to your new filename (which should not
actually be NewFunctional, but should *start* with a capital letter).

Now you need to edit `src/haskell/create_generators.py` to make that
script generate a program using your functional to generate a new
`src/new/NewFunctionalFast.cpp` file.

Next, you need to teach `fac` to generate and use the new
`src/new/NewFunctionalFast.cpp` file.  Edit `fac/haskell.py` and add
your new functional to the list.  And this is all! You just need to go
ahead and write code that uses your functional.  Of course, making the
haskell code implement the functional you actually want requires more
work.
