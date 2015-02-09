#!/usr/bin/python2

import re, string, os, sys

pyfs = sys.argv[1:]

import_re = re.compile(r"^import\s+(\w+)", re.M)

fixed_open = re.compile(r"=\s*open\(['\"]([^'\"\*]*)['\"]\s*,\s*['\"]w['\"]\s*\)")
open_changing_output = re.compile(r"open\(['\"]([^'\"\*]*)['\"]\s*%\s*(\(.*\))(\s*,[\w\s=]+)*\s*,\s*['\"]w['\"]\s*\)")
fixed_output = re.compile(r"savefig\(['\"]([^'\"\*]*)['\"](\s*,[\w\s=]+)*\)")
changing_output = re.compile(r"savefig\(['\"]([^'\"\*]*)['\"]\s*%\s*(\(.*\))(\s*,[\w\s=]+)*\s*\)")
arguments = re.compile(r"^#arg\s+(\w+)\s*=\s*(.*)$", re.M)

fixed_open_input = re.compile(r"=\s*open\(['\"]([^'\"\*]*)['\"]\s*,\s*['\"]r['\"]\s*\)")
fixed_input = re.compile(r"^[^\n#]*loadtxt\(['\"]([^'\"\*]*)['\"]\)", re.M)
changing_loadtxt_noparens = re.compile(r"^[^\n#]*loadtxt\(['\"]([^'\"\*]*)['\"]\s*%\s*([^\(\)\n]*)(\s*,[\w\s=]+)*\s*\)", re.M)
changing_loadtxt = re.compile(r"^[^\n#]*loadtxt\(['\"]([^'\"\*]*)['\"]\s*%\s*(\([^\)]*\))(\s*,[\w\s=]+)*\s*\)", re.M)
changing_input = re.compile(r"input:\s*['\"]([^'\"\*]*)['\"]\s*%\s*(\(.*\))")
list_comprehension_input = re.compile(r"input:\s*(\[\s*['\"][^'\"\*]*['\"].+for.*in.*\])\n")

def friendly_eval(code, context, local = None):
    try:
        return eval(code, local)
    except:
        print ("\nError evaluating '%s' from file %s" % (code, context))
        raise


for fname in pyfs:
    if fname == 'figs/*.py':
        continue

    f = open(fname, 'r')
    contents = f.read()
    f.close()

    inputs = fixed_input.findall(contents)
    inputs += fixed_open_input.findall(contents)

    imports = import_re.findall(contents)
    for i in imports:
        ipath = os.path.join(os.path.dirname(fname), i+'.py')
        if os.path.exists(ipath):
            inputs.append(ipath)

    outputs = fixed_output.findall(contents)
    outputs = [x[0] for x in outputs]
    outputs += fixed_open.findall(contents) # add any files created with open(...)
    argvals = arguments.findall(contents)
    if len(argvals) > 0:
        coutputs = changing_output.findall(contents) + open_changing_output.findall(contents)
        list_inputs = list_comprehension_input.findall(contents)
        cinputs = changing_input.findall(contents)
        cinputs += changing_loadtxt.findall(contents)
        cinputs += changing_loadtxt_noparens.findall(contents)
        allvalues = [{}]
        commandlineformat = ""
        commandlineargs = []
        for arg in argvals:
            commandlineformat += " \"%s\""
            commandlineargs.append(arg[0])
            newallvalues = []
            for v in friendly_eval(arg[1], fname):
                for av in allvalues:
                    newav = av.copy()
                    newav[arg[0]] = v
                    newallvalues.append(newav)
            allvalues = newallvalues
        for a in allvalues:
            aa = commandlineformat % friendly_eval('(' + string.join(commandlineargs, ",") + ')', 'file '+ fname, a)
            extrainputs = []
            extraoutputs = []
            for i in list_inputs:
                #print 'examining input:', i
                fnames = friendly_eval(i, fname, a)
                #print 'found', fnames
                extrainputs += fnames
            for i in cinputs:
                extrainputs.append(i[0] % friendly_eval(i[1], fname, a))
            for o in coutputs:
                extraoutputs.append(o[0] % friendly_eval(o[1], fname, a))
            print "? python2 %s %s" % (fname, aa)
            for i in inputs + extrainputs:
                print '<', i
            for o in outputs + extraoutputs:
                print '>', o
            print 'c .pyc'
            print 'C %s/.matplotlib' % os.getenv('HOME')
            print
    else:
        print '? python2', fname
        for i in inputs:
            print '<', i
        for o in outputs:
            print '>', o
        print 'c .pyc'
        print 'C %s/.matplotlib' % os.getenv('HOME')
        print
