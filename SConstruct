import os, time, string, glob

CacheDir(os.environ['HOME'] + '/.cache/scons')

# First, we want to set up the flags
env = Environment(CPPPATH=['src', 'include', 'tests'], LIBS=['fftw3'])
env.MergeFlags('-Wall -Werror -ansi')
env.MergeFlags('-Wno-unused-variable -Wno-unused-parameter -Wno-return-type -Wno-unused-local-typedefs')
env.MergeFlags('-O3')

# Configure git to run the test suite:
Default(env.Command(target = '.git/hooks/commit-msg',
                    source = 'git/commit-msg',
                    action = Copy("$TARGET", "$SOURCE")))
Default(env.Command(target = '.git/hooks/pre-commit',
                    source = 'git/pre-commit',
                    action = Copy("$TARGET", "$SOURCE")))

# Now we define utility functions for the tests.
passed_tests = 0
failed_tests = 0
total_tests = 0
OKBLUE = '\033[94m'
OKGREEN = '\033[92m'
WARNING = '\033[93m'
FAIL = '\033[91m'
ENDC = '\033[0m'
def test_function(target, source, env):
    global passed_tests, failed_tests
    testfile = str(source[0])
    logfile = str(target[0])
    name = testfile[6:len(testfile)-5]
    failfile = "tests/%s.failed" % name
    command = "%s > %s 2>&1" % (testfile, failfile)
    #print command
    start = time.clock()
    retval = os.system(command)
    end = time.clock()
    if retval == 0:
        os.system("mv %s %s" % (failfile, logfile))
        passed_tests += 1
        print "%sPASS: %s (%g seconds)%s" % (OKGREEN, name, end - start, ENDC)
    else:
        failed_tests += 1
        print "%sFAIL: %s (%g seconds)%s" % (FAIL, name, end - start, ENDC)
def check_function(target, source, env):
    global passed_tests, failed_tests
    run_tests = passed_tests + failed_tests
    already_passed = total_tests - run_tests
    already_string = ", %d tests already passed" % already_passed
    passed_string = "%d tests passed, " % passed_tests
    if passed_tests == 0: passed_string = ""
    if already_passed == 0: already_string = ""
    if run_tests == 0:
        print "%sAll %d tests already passed%s" % (OKGREEN, total_tests, ENDC)
        return 0
    if failed_tests == 1:
        print "%s%sone test failed%s%s" % ( passed_string, FAIL, ENDC, already_string )
        return 1
    if failed_tests != 0:
        print "%s%s%d tests failed%s%s" % ( passed_string, FAIL, failed_tests, ENDC, already_string )
        return 1
    print "%sAll %d tests passed%s%s" % ( OKGREEN, passed_tests, already_string, ENDC )
    return 0
check_runner = Builder(action = check_function)
env.Append(BUILDERS = {'Check': check_runner})
env.Check(target = 'check', source = 'monte-carlo')
env.Default('check')

def BuildTest(env, test, depends):
    global total_tests
    env.Program(target = 'tests/' + test + '.test',
                source = ['tests/' + test + '.cpp']+depends)
    testenv = env.Clone()
    #testenv.CacheDir(None) # we probably should not cache test results (which could be semi-random), but do
    test = testenv.Command(target = 'tests/' + test + '.log',
                           source = 'tests/' + test + '.test',
                           action = Action(test_function, '$SOURCE'))
    total_tests += 1
    Depends('check', test)
    return test
AddMethod(Environment, BuildTest)

def SVG(env, filename):
    filename = str(filename)
    if len(filename)>4 and filename[len(filename)-4:] == ".svg":
        filename = filename[:len(filename)-4]
    env.Command(target = filename+'.eps',
                source = filename+'.svg',
                action = 'inkscape --export-eps $TARGET $SOURCE')
    env.Command(target = filename+'.pdf',
                source = filename+'.eps',
                action = 'epstopdf $SOURCE')
AddMethod(Environment, SVG)

for name in Split(""" monte-carlo soft-monte-carlo pair-monte-carlo
                      triplet-monte-carlo polyhedra-monte-carlo
                      radial-distribution-monte-carlo """):
    env.Program(
        target=name,
        source=["src/Monte-Carlo/" + name + ".cpp", 'src/utilities.cpp'])
    Default(name)

# Generate source code (and pdfs) from haskell:

functionals = Builder(action = '$SOURCE $TARGET')
generate = Environment(BUILDERS = {'Functional' : functionals})
generated_sources = []
for source in Split(""" HardSpheresNoTensor2Fast TensorWhiteBearFast WhiteBearMarkIIFast
                        gSigmaA_by_handFast gSigmaA_automagicFast
                        TensorDensityXXFast n2DensityFast VectorDensityXFast
                        YuWuCorrelationFast SaftFluid2Fast EntropySaftFluid2Fast CorrelationGrossCorrectFast
                        gSigmaSm2Fast gSigmaAm2Fast gSigmaSFast gSigmaAFast
                        SoftFluidFast HardFluidFast WaterXFast HughesXFast """):
    filename = 'src/' + source + '.cpp'
    generated_sources.append(filename)
    generate.Functional(target = filename, source = 'src/haskell/functionals.exe')

newgenerated_sources = []
for source in Split(""" WhiteBearFast """):
    filename = 'src/new/' + source + '.cpp'
    newgenerated_sources.append(filename)
    generate.Functional(target = filename, source = 'src/haskell/newfunctionals.exe')

for pdf in Split(""" Association WhiteBear TensorWhiteBear WhiteBearMarkII Dispersion SaftFluid
                     SimpDispersion EntropySaftFluid GradDispersion JoinedGradDispersion
                     SimpGradDispersion """):
    generate.Functional(target = 'doc/' + pdf + '.pdf', source = 'src/haskell/latex-functionals.exe')
    Alias('pdf', 'doc/' + pdf + '.pdf')

# Here we have ordinary source code:

generic_sources = Split("""
  src/lattice.cpp src/utilities.cpp
  src/GridDescription.cpp src/Grid.cpp src/ReciprocalGrid.cpp
  src/IdealGas.cpp src/ChemicalPotential.cpp
  src/HardSpheres.cpp src/ExternalPotential.cpp
  src/Functional.cpp src/ContactDensity.cpp
  src/Gaussian.cpp src/Pow.cpp src/WaterSaftFast.cpp src/WaterSaft_by_handFast.cpp
  src/EffectivePotentialToDensity.cpp
  src/equation-of-state.cpp src/water-constants.cpp
  src/compute-surface-tension.cpp
  src/Minimizer.cpp src/Downhill.cpp
  src/Precision.cpp src/ConjugateGradient.cpp
  src/QuadraticLineMinimizer.cpp src/SteepestDescent.cpp

 """)
all_sources = generic_sources + generated_sources

# Here we have generic rules for our papers
for paper in Split(""" hughes-saft contact fuzzy-fmt pair-correlation water-saft
                       polyhedra """):
    p = env.PDF(target = 'papers/' + paper + '/paper.pdf',
                source = ['papers/' + paper + '/paper.tex'])
    Default(p)

Depends('index.html', 'papers/pair-correlation/figs/pretty-4.svg')

paper = Environment(tools=['gnuplot', 'matplotlib', 'mkdown'])
for paperfile in Glob('papers/*/paper.tex'):
    paperdir = str(paperfile)[:len(str(paperfile))-len('/paper.tex')]
    # first let's handle all gnuplot files
    for gpfile in Glob(paperdir + '/figs/*.gp'):
        gpfile = str(gpfile)[:len(str(gpfile))-3]
        x = paper.GplotGraph(gpfile, gp_chdir=paperdir)
        for node in x:
            outfile = str(node)
            basefile = outfile[:len(outfile)-4]
            if outfile[len(outfile)-4:] == ".eps":
                env.Command(target = basefile + ".pdf",
                            source = outfile,
                            action = 'epstopdf $SOURCE')
    # and now let's handle all python files
    for pyfile in Glob(paperdir + '/figs/*.py'):
        pyfile = str(pyfile)[:len(str(pyfile))-3]
        paper.Matplotlib(pyfile, py_chdir=paperdir)
    # and now we'll handle all svg files
    for svgfile in Glob(paperdir + '/figs/*.svg'):
        env.SVG(svgfile)

# #################### papers/hughes-saft ##################################################
env.Command(target = 'papers/hughes-saft/figs/single-rods-calculated-density.dat',
            source = ['papers/hughes-saft/figs/density_calc.py',
                      'papers/hughes-saft/figs/single-rod-in-water.dat'],
            action = 'cd papers/hughes-saft && python figs/density_calc.py')
env.Command(target = 'papers/hughes-saft/figs/single-rod-in-water.dat',
            source = Glob('papers/hughes-saft/figs/single-rod-*nm-energy.dat'),
            action = string.join(['cat '] +
                                 sorted(glob.glob('papers/hughes-saft/figs/single-rod-*nm-energy.dat')) +
                                 [' > $TARGET']))


# #################### papers/contact #######################################################

# #################### papers/water-saft ####################################################
env.Command(target = 'papers/water-saft/figs/single-rod-in-water.dat',
            source = Glob('papers/water-saft/figs/single-rod-*nm-energy.dat'),
            action = string.join(['cat '] +
                                 sorted(glob.glob('papers/water-saft/figs/single-rod-*nm-energy.dat')) +
                                 [' > $TARGET']))
env.Command(target = 'papers/water-saft/figs/hughes-single-rod-in-water.dat',
            source = Glob('papers/water-saft/figs/hughes-single-rod-*nm-energy.dat'),
            action = string.join(['cat '] +
                                 sorted(glob.glob('papers/water-saft/figs/hughes-single-rod-*nm-energy.dat')) +
                                 [' > $TARGET']))

# #################### papers/pair-correlation ##############################################

# #################### papers/fuzzy-fmt ##################################################


# The following programs generate a single .dat file that may be cached.
for mkdat in Split("""
	papers/water-saft/figs/surface-tension
	papers/water-saft/figs/pressure-with-isotherms
	papers/hughes-saft/figs/surface-tension
	papers/hughes-saft/figs/pressure-with-isotherms
	papers/contact/figs/gHS-vs-n
	papers/contact/figs/free-energy
      """):
    env.Program(target = mkdat + '.mkdat',
                source = [mkdat + '.cpp'] + all_sources)
    env.Command(mkdat + '.dat', mkdat + '.mkdat', './$SOURCE')

# The following programs generate several .dat files with different
# names, and therefore are unsafe to cache, since we do not list all
# these files here.
for mkdat in Split("""
	papers/hughes-saft/figs/single-rod-in-water-low-res
	papers/contact/figs/walls
	papers/pair-correlation/figs/walls
      """):
    env.Program(target = mkdat + '.mkdat',
                source = [mkdat + '.cpp'] + all_sources)
    # Do not cache output of mkdat, because other files also need to
    # be produced.
    NoCache(env.Command(mkdat + '.dat', mkdat + '.mkdat', './$SOURCE'))

for mkdat in Split("""
	papers/water-saft/figs/equation-of-state
	papers/hughes-saft/figs/equation-of-state
	papers/water-saft/figs/rods-in-water
	papers/water-saft/figs/four-rods-in-water
	papers/water-saft/figs/single-rod
	papers/water-saft/figs/hughes-single-rod
	papers/water-saft/figs/sphere
	papers/hughes-saft/figs/rods-in-water
	papers/hughes-saft/figs/four-rods-in-water
	papers/hughes-saft/figs/single-rod
	papers/hughes-saft/figs/single-rod-in-water-high-res
	papers/hughes-saft/figs/sphere
	papers/contact/figs/sphere
	papers/contact/figs/inner-sphere
	papers/contact/figs/test-particle-wall
	papers/pair-correlation/figs/sphere-with-wall
	papers/fuzzy-fmt/figs/walls
      """):
    env.Program(target = mkdat + '.mkdat',
                source = [mkdat + '.cpp'] + all_sources)

env.Command(target = ['papers/water-saft/figs/equation-of-state.dat',
                      'papers/water-saft/figs/experimental-equation-of-state.dat'],
            source = ['papers/water-saft/figs/equation-of-state.mkdat'],
            action = './$SOURCE')
env.Command(target = ['papers/hughes-saft/figs/equation-of-state.dat',
                      'papers/hughes-saft/figs/experimental-equation-of-state.dat'],
            source = ['papers/hughes-saft/figs/equation-of-state.mkdat'],
            action = './$SOURCE')

########################### Now let's build the website! ##################################

Default(paper.Markdown('index.md'))


# Here we have rules to build the haskell code

haskell = Environment(tools=['haskell'],
                      HSSEARCHPATH = ["src/haskell"],
                      HSPACKAGES = ["containers", "process", "HUnit"],
                      HSCFLAGS = ['-O2'])
haskell_source = []
for m in Split(""" LatexDouble Latex Expression CodeGen Statement HughesSaft WaterSaft Optimize
                   FMT WhiteBear IdealGas SFMT NewCode """):
    haskell_source.append("src/haskell/" + m + ".hs")

for program in Split("functionals newfunctionals test latex-functionals"):
    haskell.HaskellProgram(target = 'src/haskell/' + program + '.exe',
                           source = ['src/haskell/' + program  + '.hs'] + haskell_source)

haskell.Command(target = ['tests/generated-haskell/nice-sum.h',
                          'tests/generated-haskell/nice-sqrt.h',
                          'tests/generated-haskell/nice-sqrtandmore.h',
                          'tests/generated-haskell/nice-quad.h',
                          'tests/generated-haskell/nice-phi1.h',
                          'tests/generated-haskell/nice-phi2.h',
                          'tests/generated-haskell/nice-phi3.h',
                          'tests/generated-haskell/nice-n2xsqr.h',
                          'tests/generated-haskell/whitebear.tex',
                          'tests/generated-haskell/math.tex'],
                source = 'src/haskell/test.exe', action = '$SOURCE codegen')

haskell.Command(target = ['tests/new-generated-haskell/WhiteBear.h',
                          'tests/new-generated-haskell/integrate_sqr.h',
                          'tests/new-generated-haskell/volume_minus_one_sqr.h'],
                source = 'src/haskell/newfunctionals.exe', action = '$SOURCE tests')

################# Now do the test suite ##################################################
for test in Split(""" memory saft eos print-iter convolve-finite-difference precision
                      compare-gsigmas
                      convolve functional-of-double ideal-gas eps fftinverse generated-code  """):
    env.BuildTest(test, all_sources)

for test in Split(""" functional-arithmetic surface-tension new-generated-code """):
    env.BuildTest(test, generic_sources)

for test in Split(""" newcode """):
    env.BuildTest(test, ['src/new/Minimize.cpp'])

for test in Split(""" new-hard-spheres """):
    env.BuildTest(test, all_sources + newgenerated_sources)
