import os, time, string, glob, numpy

CacheDir(os.environ['HOME'] + '/.cache/scons')

# First, we want to set up the flags
env = Environment(CPPPATH=['src', 'include', 'tests'], LIBS=['fftw3', 'popt'])
env.MergeFlags('-Werror -ansi')

# The following is approximately equal to -Wall, but we use a additive
# set of statements, since the set of warnings encompased in -Wall
# sometimes changes in ways that cause incompatibilities.
env.MergeFlags('-Waddress -Warray-bounds -Wc++11-compat -Wchar-subscripts ' +
               '-Wenum-compare -Wcomment -Wformat -Wmaybe-uninitialized -Wmissing-braces ' +
               '-Wnonnull -Wparentheses -Wreorder -Wreturn-type ' +
               '-Wsequence-point -Wsign-compare -Wstrict-aliasing -Wstrict-overflow=1 ' +
               '-Wswitch -Wtrigraphs -Wuninitialized -Wunknown-pragmas -Wunused-function ' +
               '-Wunused-label -Wunused-value -Wunused-variable')

env.MergeFlags('-O3') # add -g and possibly -fno-inline here to enable debugging

# The following flags enable gcc to eliminate unused code from the
# final executable.  This reduces the size of the executable, and I
# hope it also means that executables are less likely to change when
# we add new code to deft that they do not use.
env.AppendUnique(LINKFLAGS=['-Wl,-gc-sections'],
                 CXXFLAGS=['-fdata-sections','-ffunction-sections'])

# Configure git to run the test suite:
Alias('git configuration',
      env.Command(target = '.git/hooks/commit-msg',
                  source = 'git/commit-msg',
                  action = Copy("$TARGET", "$SOURCE")))
Alias('git configuration',
      env.Command(target = '.git/hooks/pre-commit',
                  source = 'git/pre-commit',
                  action = Copy("$TARGET", "$SOURCE")))
Default('git configuration')

haskell = Environment(tools=['haskell'],
                      HSSEARCHPATH = ["src/haskell"],
                      HSPACKAGES = ["containers", "process", "HUnit"],
                      HSCFLAGS = ['-O2', '-Wall', '-Werror'])

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
    time.sleep(1) # work around sporadic race condition with "Text file busy" error.
    command = "%s 2>&1 > %s" % (testfile, failfile)
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
    env.Command(target = filename+'.pdf',
                source = filename+'.eps',
                action = 'epstopdf --outfile $TARGET $SOURCE')
    env.Command(target = filename+'.eps',
                source = filename+'.svg',
                action = 'inkscape --export-eps $TARGET $SOURCE')
    # env.Command(target = filename+'.pdf',
    #             source = filename+'.svg',
    #             action = 'inkscape --export-pdf $TARGET $SOURCE')
AddMethod(Environment, SVG)

for name in Split(""" monte-carlo soft-monte-carlo pair-monte-carlo
                      triplet-monte-carlo polyhedra-monte-carlo polyhedra-talk
                      square-well-monte-carlo
                      radial-distribution-monte-carlo """):
    env.Program(
        target=name,
        source=["src/Monte-Carlo/" + name + ".cpp", 'src/utilities.cpp', 'src/Monte-Carlo/polyhedra.cpp', 'src/Monte-Carlo/square-well.cpp', 'src/vector3d.cpp'])
    Alias('executables', name)
Default('executables')

# Generate source code (and pdfs) from haskell:

functionals = Builder(action = '$SOURCE $TARGET')
generate = Environment(BUILDERS = {'Functional' : functionals})
generated_sources = []
for source in Split(""" HardSpheresNoTensor2Fast TensorWhiteBearFast WhiteBearMarkIIFast
                        gSigmaA_by_handFast gSigmaA_automagicFast
                        TensorDensityXXFast n2DensityFast VectorDensityXFast
                        YuWuCorrelationFast SaftFluid2Fast EntropySaftFluid2Fast CorrelationGrossCorrectFast
                        gSigmaSm2Fast gSigmaAm2Fast gSigmaSFast gSigmaAFast
                        HughesHBFast
                        SoftFluidFast HardFluidFast HardRosenfeldFluidFast WaterXFast HughesXFast """):
    filename = 'src/' + source + '.cpp'
    generated_sources.append(filename)
    generate.Functional(target = filename, source = 'src/haskell/functionals.exe')

newgeneric_sources = Split(""" src/new/Minimize.cpp src/new/NewFunctional.cpp """)

newgenerated_sources = []
for name, module, hsfunctional, inputs in [
    # The following are just for testing purposes
    ("ExternalPotentialTest", "ExternalPotentialTest", "external_potential", '[(ER $ r_var "n", ER 1)]'),
    ("Quadratic", "Quadratic", "quadratic", '[(ER $ r_var "x", ER 1)]'),
    ("QuadraticN0", "QuadraticN0", "quadratic_n0", '[(ER $ r_var "x", ER 1)]'),
    ("QuadraticGaussian", "QuadraticGaussian", "quadratic_gaussian", '[(ER $ r_var "x", ER 1)]'),
    ("LogN0", "LogN0", "log_n0", '[(ER $ r_var "x", ER 1)]'),
    ("Phi1", "WhiteBear", "kTphi1", '[(ER $ r_var "x", ER 1)]'),
    ("Phi2", "WhiteBear", "kTphi2", '[(ER $ r_var "x", ER 1)]'),
    ("Phi3", "WhiteBear", "kTphi3", '[(ER $ r_var "x", ER 1)]'),
    ("SPhi1", "SFMT", "phi1", '[(ER $ r_var "x", ER 1)]'),
    ("SPhi2", "SFMT", "phi2", '[(ER $ r_var "x", ER 1)]'),
    ("SPhi3", "SFMT", "phi3", '[(ER $ r_var "x", ER 1)]'),
    # The rest are "real" functionals of sorts
    ("HomogeneousWhiteBear", "WhiteBear", "homogeneous_whitebear", '[]'),
    ("WhiteBear", "WhiteBear", "whitebear_n", '[(ER $ r_var "n", ER 1)]'),
    ("WhiteBearFluid", "WhiteBear", "whitebear_fluid_n", '[(ER $ r_var "n", ER 1)]'),
    ("SFMTFluid", "SFMT", "sfmt_fluid_n", '[(ER $ r_var "n", ER 1)]'),
    ("SFMTFluidVeff", "SFMT", "sfmt_fluid_Veff",
           '[(ER $ r_var "Veff", ER (exp(-r_var "Veff"/s_var "kT")))]'),
    ("HomogeneousSFMTFluid", "SFMT", "homogeneous_sfmt_fluid", '[]'),
    ("WaterSaft", "WaterSaft", "water_saft_n", '[]'), # no gradients:  for debugging!
    ("WaterSaftByHand", "WaterSaft", "water_saft_by_hand_n", '[]'), # no gradients:  for debugging!
    ("HomogeneousWaterSaft", "WaterSaft", "homogeneous_water_saft_n", '[]'),
    ("HomogeneousWaterSaftByHand", "WaterSaft", "homogeneous_water_saft_by_hand_n", '[]')]:
    # I'm sloppy and just recreate the generate_%s.hs files every time
    f = open('src/haskell/generate_%s.hs' % name, "w")
    f.write("""import NewCode
import %s ( %s )

main :: IO ()
main = createHeaderAndCppFiles %s %s "%s"
""" % (module, hsfunctional, hsfunctional, inputs, name))
    f.close()
    haskell.HaskellMake(target = 'src/haskell/generate_%s.exe' % name,
                        source = 'src/haskell/generate_%s.hs' % name)
    env.Command(target = ['src/new/%sFast.cpp' % name, 'src/new/%sFast.h' % name],
                source = 'src/haskell/generate_%s.exe' % name,
                action = './$SOURCE')
    newgenerated_sources += ['src/new/%sFast.cpp' % name]

for pdf in Split(""" Association WhiteBear TensorWhiteBear WhiteBearMarkII Dispersion SaftFluid
                     SimpDispersion EntropySaftFluid GradDispersion JoinedGradDispersion
                     SimpGradDispersion SFMT """):
    generate.Functional(target = 'doc/' + pdf + '.pdf', source = 'src/haskell/latex-functionals.exe')
    Alias('pdf', 'doc/' + pdf + '.pdf')
Default('pdf')

# Here we have ordinary source code:

generic_sources = Split("""

  src/lattice.cpp src/utilities.cpp src/Faddeeva.cc
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

env.AppendUnique(TARFLAGS = ['-c','-z'])
# Here we have generic rules for our papers
for paper in Split(""" hughes-saft contact fuzzy-fmt pair-correlation water-saft
                       square-well-liquid polyhedra renormalization electrostatics """):
    p = env.PDF(target = 'papers/' + paper + '/paper.pdf',
                source = ['papers/' + paper + '/paper.tex'])
    NoCache(p)
    Alias('papers', p)
    tar_env = Environment(TARFLAGS = r"-c -z --transform 's/papers\///'")
    tar_env.Tar(target = 'papers/' + paper + '/arxiv.tar.gz',
            source = ['papers/' + paper + '/paper.tex',
                      'papers/' + paper + '/paper.bbl'] +
                     Glob('papers/' + paper + '/figs/*.pdf') +
                     Glob('papers/' + paper + '/figs/*.jpg') +
                     Glob('papers/' + paper + '/figs/*.png') +
                     Glob('papers/' + paper + '/figs/*.tex'))
    Depends('papers/' + paper + '/arxiv.tar.gz', 'papers/' + paper + '/paper.pdf')
    Depends('papers/' + paper + '/arxiv.tar.gz', 'papers/' + paper + '/paper.bbl')
    Alias('papers', 'papers/' + paper + '/arxiv.tar.gz')
Default('papers')

Alias('papers', env.PDF('papers/thesis-hughes/project.tex'))
Alias('papers', env.PDF('papers/thesis-roth/project.tex'))
Alias('papers', env.PDF('papers/pair-correlation/notes.tex'))
Alias('papers', env.PDF('papers/pair-correlation/figs/ghs-analytics.tex'))
Alias('papers', env.PDF('papers/polyhedra/harmonics.tex'))
Alias('papers', env.PDF('papers/polyhedra/wigner-properties.tex'))

Depends('index.html', 'papers/pair-correlation/figs/pretty-4.svg')

paper = Environment(tools=['gnuplot', 'matplotlib', 'mkdown'])
for paperfile in Glob('papers/*/paper.tex') + Glob('papers/*/project.tex'):
    paperdir = os.path.dirname(str(paperfile))
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
    for pyfile in Glob(paperdir + '/figs/*.py') + Glob(paperdir + '/*.py'):
        pyfile = str(pyfile)[:len(str(pyfile))-3]
        paper.Matplotlib(pyfile, py_chdir=paperdir)
    # and now we'll handle all svg files
    for svgfile in Glob(paperdir + '/figs/*.svg'):
        env.SVG(svgfile)

# #################### papers/hughes-saft ##################################################
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
for atom in ['Ne', 'Ar', 'Kr', 'Xe']:
    env.Command(target = 'papers/water-saft/figs/lj-%s.dat' % atom,
                source = Glob('papers/water-saft/figs/lj-%s-*K-energy.dat' % atom),
                action = string.join(['cat '] +
                                     sorted(glob.glob('papers/water-saft/figs/lj-%s-*K-energy.dat' % atom)) +
                                     [' > $TARGET']))
for atom in ['Ne', 'Ar', 'Kr', 'Xe']:
    env.Command(target = 'papers/water-saft/figs/hughes-lj-%s.dat' % atom,
                source = Glob('papers/water-saft/figs/hughes-lj-%s-*K-energy.dat' % atom),
                action = string.join(['cat '] +
                                     sorted(glob.glob('papers/water-saft/figs/hughes-lj-%s-*K-energy.dat' % atom)) +
                                     [' > $TARGET']))

# #################### papers/pair-correlation ##############################################

# #################### papers/fuzzy-fmt ##################################################

# #################### papers/square-well-liquid ##################################################

Alias('papers', env.PDF('papers/square-well-liquid/histogram-paper.tex'))

# The following enables automagic monte-carlo generation of
# low-quality data for simple plots
for ff in [0.1, 0.2, 0.3, 0.4]:
    datadir = "papers/square-well-liquid/data/"
    for ww in [1.3, 1.5, 2.0, 3.0]:
        for N in [200]:
            env.Command(target = [datadir+"periodic-ww%04.2f-ff%04.2f-N%i-nw-E.dat" % (ww, ff, N),
                                  datadir+"periodic-ww%04.2f-ff%04.2f-N%i-nw-dos.dat" % (ww, ff, N),
                                  datadir+"periodic-ww%04.2f-ff%04.2f-N%i-nw-g.dat" % (ww, ff, N)],
                        source = 'square-well-monte-carlo',
                        action = './square-well-monte-carlo --nw --N %d --initialize=4000 --ff %g --ww %g  --iterations 10000' % (N, ff, ww))
            env.Command(target = [datadir+"periodic-ww%04.2f-ff%04.2f-N%i-flat-E.dat" % (ww, ff, N),
                                  datadir+"periodic-ww%04.2f-ff%04.2f-N%i-flat-dos.dat" % (ww, ff, N),
                                  datadir+"periodic-ww%04.2f-ff%04.2f-N%i-flat-g.dat" % (ww, ff, N)],
                        source = 'square-well-monte-carlo',
                        action = './square-well-monte-carlo --flat --N %d --initialize=4000 --ff %g --ww %g  --iterations 10000' % (N, ff, ww))
            for kT in [i*.1 for i in range(1,10)] + range(1,10):
                env.Command(target = [datadir+"periodic-ww%04.2f-ff%04.2f-N%i-kT%g-E.dat" % (ww, ff, N, kT),
                                      datadir+"periodic-ww%04.2f-ff%04.2f-N%i-kT%g-dos.dat" % (ww, ff, N, kT),
                                      datadir+"periodic-ww%04.2f-ff%04.2f-N%i-kT%g-g.dat" % (ww, ff, N, kT)],
                            source = 'square-well-monte-carlo',
                            action = './square-well-monte-carlo --kT %g --N %d --initialize=10000 --ff %g --ww %g  --iterations 10000' % (kT, N, ff, ww))


# #################### talks ##################################################


# and now we'll handle all svg files
for svgfile in Glob('talks/*/*/*.svg'):
    env.SVG(svgfile)
# and all py files
for talkdir in glob.glob('talks/*'):
    for pyfile in Glob(talkdir + '/*/*.py'):
        pyfile = str(pyfile)[:len(str(pyfile))-3]
        paper.Matplotlib(pyfile, py_chdir=talkdir)
    env.PDF(talkdir+'/slides.tex')
    #Alias('talks', talkdir+'/slides.pdf')

Alias('talks', ['talks/colloquium/slides.pdf'])
#Alias('talks', ['talks/colloquium/slides.pdf', 'talks/polyhedra/slides.pdf'])

###################### dependencies for polyhedra talk ####################
# servernum = 100
# def xserver():
#   # xvfb-run needs to run each command on a unique server,
#   # so using this function for the server number ensures that
#   global servernum
#   servernum -= 1
#   return servernum

# Depends('talks/polyhedra/slides.pdf',
#         env.Command(target = ['talks/polyhedra/figs/background.png'] +
#                     ['talks/polyhedra/figs/tet-%i.png' %i for i in xrange(3)] +
#                     ['talks/polyhedra/figs/ice-structure-%i.png' %i for i in xrange(5)],
#                     source = ['talks/polyhedra/generate-figs.py',
#                               'talks/polyhedra/dat/background.dat'] +
#                     ['talks/polyhedra/dat/tet-%i.dat' %i for i in xrange(3)] +
#                     ['talks/polyhedra/dat/ice-structure.dat'],
#                     action = 'cd talks/polyhedra && xvfb-run -n %i --server-args="-screen 0 1024x768x24" ./generate-figs.py' %xserver()))
# Depends('talks/polyhedra/slides.pdf',
#         NoCache(env.Command(target = ['talks/polyhedra/dat/background.dat'],
#                     source = 'polyhedra-talk',
#                     action = './polyhedra-talk')))
# Depends('talks/polyhedra/slides.pdf',
#         env.Command(target = ['talks/polyhedra/anim/mc-slow-%03i.pdf' % i
#                               for i in xrange(30)],
#                     source = 'talks/polyhedra/mc-tri-slow.py',
#                     action = 'cd talks/polyhedra && python mc-tri-slow.py'))
# Depends('talks/polyhedra/slides.pdf',
#         env.Command(target = ['talks/polyhedra/anim/mc100-%4.2f-%03i.png' % (.5, i)
#                               for i in xrange(30)],
#                     source = 'talks/polyhedra/mc-tri.py',
#                     action = 'cd talks/polyhedra && python mc-tri.py'))
# for ff in [0.42, 0.71]:
#   Depends('talks/polyhedra/slides.pdf',
#           env.Command(
#             target = ['papers/polyhedra/figs/anim/periodic-%04.2f-truncated_tetrahedron-216-%i.png' % (ff, i)
#                       for i in xrange(10)],
#             source = ['papers/polyhedra/figs/animate_polyhedra.py'] +
#             ['papers/polyhedra/figs/mc/vertices/periodic-%04.2f-vertices-truncated_tetrahedron-216-%i.dat' %(ff, i)
#              for i in xrange(10)],
#             action = 'cd papers/polyhedra && xvfb-run -n %i --server-args="-screen 0 1024x768x24" ./figs/animate_polyhedra.py %04.2f -p -N 216 -f 10 --save --hide --notext' %(xserver(), ff)))
#   Depends('talks/polyhedra/slides.pdf',
#           env.Command(
#             target = ['papers/polyhedra/figs/mc/vertices/periodic-%04.2f-vertices-truncated_tetrahedron-216-%i.dat'
#                       % (ff, i) for i in xrange(10)],
#                       source = 'polyhedra-monte-carlo',
#                       action = './polyhedra-monte-carlo --ff %04.2f --periodx 20 --periody 20 --periodz 20 --N 216 --iterations 100 --initialize_iterations 1000 --save_vertices 10' %ff))

###################### dependencies for colloquium ########################

Depends('talks/colloquium/slides.pdf',
        Split(""" talks/colloquium/figs/energy-solid.pdf
                  talks/colloquium/figs/energy-solid-gas-liquid.pdf
                  talks/colloquium/figs/energy-solid-gas.pdf """))
Depends('talks/colloquium/slides.pdf',
        env.Command(target = ['talks/colloquium/anim/mc-slow-%03i.pdf' % i
                              for i in xrange(31)],
                    source = 'talks/colloquium/mc-circle-slow.py',
                    action = 'cd talks/colloquium && python mc-circle-slow.py'))
Depends('talks/colloquium/slides.pdf',
        env.Command(target = ['talks/colloquium/anim/mc-density-%03i.pdf' % i
                              for i in xrange(31)],
                    source = 'talks/colloquium/mc-circle-slow.py',
                    action = 'cd talks/colloquium && python mc-circle-slow.py density'))
Depends('talks/colloquium/slides.pdf',
        env.Command(target = ['talks/colloquium/anim/mc-pair-%03i.pdf' % i
                              for i in xrange(31)],
                    source = 'talks/colloquium/mc-circle-slow.py',
                    action = 'cd talks/colloquium && python mc-circle-slow.py pair'))
Depends('talks/colloquium/slides.pdf',
        env.Command(target = ['talks/colloquium/anim/mc-gsigma-%03i.pdf' % i
                              for i in xrange(31)],
                    source = 'talks/colloquium/mc-circle-slow.py',
                    action = 'cd talks/colloquium && python mc-circle-slow.py gsigma'))

Default('talks')


# The following programs generate a single .dat file that may be cached.
for mkdat in Split("""
	papers/water-saft/figs/surface-tension
	papers/water-saft/figs/pressure-with-isotherms
	papers/hughes-saft/figs/surface-tension
	papers/hughes-saft/figs/pressure-with-isotherms
	papers/contact/figs/gHS-vs-n
	papers/contact/figs/free-energy
        papers/fuzzy-fmt/figs/homogeneous
      """):
    Alias('executables',
          env.Program(target = mkdat + '.mkdat',
                      source = [mkdat + '.cpp'] + all_sources))
    env.Command(mkdat + '.dat', mkdat + '.mkdat', './$SOURCE')

# The following programs generate several .dat files with different
# names, and therefore are unsafe to cache, since we do not list all
# these files here.
for mkdat in Split("""
	papers/hughes-saft/figs/single-rod-in-water-low-res
      """):
    Alias('executables',
          env.Program(target = mkdat + '.mkdat',
                      source = [mkdat + '.cpp'] + all_sources))
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
	papers/water-saft/figs/lj-atom
	papers/water-saft/figs/hughes-lj-atom
	papers/water-saft/figs/hughes-lj-atom-hs-density
	papers/hughes-saft/figs/rods-in-water
	papers/hughes-saft/figs/four-rods-in-water
	papers/hughes-saft/figs/single-rod
	papers/hughes-saft/figs/single-rod-in-water-high-res
	papers/hughes-saft/figs/sphere
	papers/contact/figs/sphere
	papers/contact/figs/inner-sphere
	papers/contact/figs/test-particle-wall
	papers/contact/figs/walls
	papers/pair-correlation/figs/sphere-with-wall
	papers/pair-correlation/figs/triplet-dft
	papers/pair-correlation/figs/walls
	papers/fuzzy-fmt/figs/walls
	papers/fuzzy-fmt/figs/soft-sphere
      """):
    Alias('executables',
          env.Program(target = mkdat + '.mkdat',
                      source = [mkdat + '.cpp'] + all_sources))

for mkdat in Split("""
	papers/fuzzy-fmt/figs/new-walls
      """):
    Alias('executables',
          env.Program(target = mkdat + '.mkdat',
                      source = [mkdat + '.cpp'] + generic_sources + newgeneric_sources +
                      ['src/new/SFMTFluidFast.cpp',
                       'src/new/SFMTFluidVeffFast.cpp', 'src/new/HomogeneousSFMTFluidFast.cpp']))
# rules for how to run fuzzy-fmt/figs/new-walls.mkdat:
for kT in [0.00001, 0.0001, 0.001, 0.01, 0.02, 0.03, 0.1]:
    for ff in [0.1, 0.2, 0.3, 0.4]:
        env.Command(target = "papers/fuzzy-fmt/figs/new-data/wall-%04.2f-%08.5g.dat" % (ff, kT),
                    source = ['papers/fuzzy-fmt/figs/new-walls.mkdat'],
                    action = './$SOURCE %g %g' % (ff, kT))

env.Command(target = ['papers/fuzzy-fmt/figs/walls.dat',
                      'papers/fuzzy-fmt/figs/wallshard-0.0000-0.10.dat',
                      'papers/fuzzy-fmt/figs/wallshard-0.0000-0.20.dat',
                      'papers/fuzzy-fmt/figs/wallshard-0.0000-0.30.dat',
                      'papers/fuzzy-fmt/figs/wallshard-0.0000-0.40.dat',
                      #'papers/fuzzy-fmt/figs/wallshard-0.0000-0.50.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0010-0.10.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0010-0.20.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0010-0.30.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0010-0.40.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0100-0.10.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0100-0.20.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0100-0.30.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0100-0.40.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0300-0.10.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0300-0.20.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0300-0.30.dat',
                      'papers/fuzzy-fmt/figs/wallssoft-0.0300-0.40.dat'],
                      #'papers/fuzzy-fmt/figs/wallssoft-0.0100-0.50.dat',
                      #'papers/fuzzy-fmt/figs/wallssoft-0.0200-0.50.dat',
                      #'papers/fuzzy-fmt/figs/wallssoft-0.0300-0.50.dat'],
            source = ['papers/fuzzy-fmt/figs/walls.mkdat'],
            action = './$SOURCE')

env.Command(target = ['papers/contact/figs/walls.dat',
                      'papers/contact/figs/wallsWB-0.10.dat',
                      'papers/contact/figs/wallsWB-0.20.dat',
                      'papers/contact/figs/wallsWB-0.30.dat',
                      'papers/contact/figs/wallsWB-0.40.dat'],
            source = ['papers/contact/figs/walls.mkdat'],
            action = './$SOURCE')

# The following lists the astonishing number of outputs generated by
# pair-correlation walls.mkdat.
env.Command(target = ['papers/pair-correlation/figs/walls.dat',
                      'papers/pair-correlation/figs/walls/x.dat',
                      'papers/pair-correlation/figs/walls/z.dat',
                      'papers/pair-correlation/figs/walls/inverse-sixth-dadz-this-work-0.30-rmax-5.dat',
                      'papers/pair-correlation/figs/walls/inverse-sixth-dadz-this-work-mc-0.30-rmax-5.dat',
                      'papers/pair-correlation/figs/walls/inverse-sixth-dadz-sokolowski-0.30-rmax-5.dat',
                      'papers/pair-correlation/figs/walls/square-well-dadz-this-work-0.30-1.790.dat',
                      'papers/pair-correlation/figs/walls/square-well-dadz-this-work-mc-0.30-1.790.dat',
                      'papers/pair-correlation/figs/walls/square-well-dadz-sokolowski-0.30-1.790.dat']+
            ['papers/pair-correlation/figs/wallsWB-%04.2f.dat' % ff
             for ff in [0.1, 0.2, 0.3]] +
            ['papers/pair-correlation/figs/walls/walls_daWB-%s-%04.2f-%05.3f.dat' % (method,ff,r)
             for ff in [0.1, 0.2, 0.3]
             for method in ['this-work', 'this-work-mc', 'sokolowski', 'fischer']
             for r in [2.005, 3.005]] +
            ['papers/pair-correlation/figs/walls/wallsWB-%s-pair-%04.2f-%04.2f.dat' % (method,ff,z0)
             for ff in [0.1, 0.2, 0.3]
             for method in ['this-work', 'this-work-mc', 'sokolowski', 'fischer']
             for z0 in numpy.arange(0.05, 3.95, 0.1)] +
            ['papers/pair-correlation/figs/walls/wallsWB-path-%s-pair-%04.2f-0.005.dat' % (method,ff)
             for ff in [0.1, 0.2, 0.3]
             for method in ['this-work', 'this-work-mc', 'sokolowski', 'fischer']],
            source = ['papers/pair-correlation/figs/walls.mkdat'],
            action = './$SOURCE')

env.Command(target = ['papers/water-saft/figs/equation-of-state.dat',
                      'papers/water-saft/figs/experimental-equation-of-state.dat'],
            source = ['papers/water-saft/figs/equation-of-state.mkdat'],
            action = './$SOURCE')
env.Command(target = ['papers/hughes-saft/figs/equation-of-state.dat',
                      'papers/hughes-saft/figs/experimental-equation-of-state.dat'],
            source = ['papers/hughes-saft/figs/equation-of-state.mkdat'],
            action = './$SOURCE')

########################### Now let's build the website! ##################################

Alias('webpage', paper.Markdown('index.md'))
Default('webpage')

# Here we have rules to build the haskell code

for hs in Glob("src/haskell/[A-Z]*.hs"):
    haskell.HaskellObject(hs)

for program in Split("functionals test latex-functionals"):
    haskell.HaskellMake(target = 'src/haskell/' + program + '.exe',
                        source = ['src/haskell/' + program  + '.hs'])

NoCache(
    haskell.Command(target = ['tests/generated-haskell/nice-sum.h',
                              'tests/generated-haskell/whitebear.tex',
                              'tests/generated-haskell/math.tex'],
                    source = 'src/haskell/test.exe',
                    # work around sporadic race condition giving "Text file busy" error.
                    action = 'sleep 1 && $SOURCE codegen'))

################# Now do the test suite ##################################################
for test in Split(""" memory saft eos print-iter convolve-finite-difference precision
                      compare-gsigmas
                      convolve functional-of-double ideal-gas eps fftinverse generated-code  """):
    env.BuildTest(test, all_sources)

for test in Split(""" new-fftinverse functional-arithmetic surface-tension """):
    env.BuildTest(test, generic_sources)

for test in Split(""" sfmt """):
    env.BuildTest(test, generic_sources + ['src/SoftFluidFast.cpp'])

for test in Split(""" newcode """):
    env.BuildTest(test, newgeneric_sources)

for test in Split(""" new-hard-spheres new-water-saft new-generated """):
    env.BuildTest(test, generic_sources + newgeneric_sources + newgenerated_sources)
