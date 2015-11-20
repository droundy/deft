# Code organization

Here is the beginning of a map for the source code.  At the top level
(in your `deft/` directory), there are a number of files, as well as a
few directories.
There are several directories that hold most of the interesting files:

1. [src](src.html) is the directory where most of our core code is located.
2. [tests](tests.html) is where we keep test programs, that verify that our code
   behaves as we desire.

   Currently the tests are not set up to work with fac, so they are
   not in use.

3. [papers](papers.html) is where we keep $\LaTeX$ source
   for papers that use Deft.  These directory also contain the
   specific code needed to generate the data for our papers, and the
   scripts to generate figures, etc.  You can build the papers (as
   well as everything else) by running:

         fac

4. The [fac](fac.html) directory holds most of the files needed to
   configure our build.
