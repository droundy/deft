# Code organization

Here is the beginning of a map for the source code.  At the top level
(in your `deft/` directory), there are a number of files, as well as a
few directories.

Of the files that live in the top-level directory, only a couple are
very important:

1. `autogen.sh` is run to initialize everything, and generate the
   `Makefile` that is used to build deft.  You only need to run
   `autogen.sh` once in each copy of deft that you have, so it's easy
   to overlook.

2. `Makefile.am` is the "makefile" that is used by `autogen.sh` (which
   itself runs [automake][] to generate the actual
   `Makefile`.  You will probably need to edit this file, if you add
   new source code files to deft.  Its syntax is a bit weird, so
   usually you just want to find something similar that is already
   there, and copy and paste, remembering to search for **every**
   occurance in `Makefile.am`

[automake]: http://www.gnu.org/software/automake/manual/automake.html

There are several directories that hold most of the interesting files:

1. [src](../src/src.html) is the directory where most of our core code is located.

2. [tests](../tests/tests.html) is where we keep test programs, that verify that our code
   behaves as we desire.

   You can run the tests by running:

       srun -p debian -c4 make -j4 check

3. [papers](../papers/papers.html) is where we keep LaTeX source for papers that use Deft.
   These directory also contain the specific code needed to generate
   the data for our papers, and the scripts to generate figures, etc.
   You can build the papers by running:

       srun -p debian -c4 make -j4 papers
