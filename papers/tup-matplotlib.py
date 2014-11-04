#!/usr/bin/python2

import re, glob, string, os, sys

import_re = re.compile(r"^import\s+(\w+)", re.M)

fixed_open = re.compile(r"=\s*open\(['\"]([^'\"]*)['\"]\s*,\s*['\"]w['\"]\s*\)")
fixed_output = re.compile(r"savefig\(['\"]([^'\"]*)['\"](\s*,[\w\s=]+)*\)")
changing_output = re.compile(r"savefig\(['\"]([^'\"]*)['\"]\s*%\s*(\(.*\))(\s*,[\w\s=]+)*\s*\)")
arguments = re.compile(r"^#arg\s+(\w+)\s*=\s*(.*)$", re.M)

fixed_open_input = re.compile(r"=\s*open\(['\"]([^'\"]*)['\"]\s*,\s*['\"]r['\"]\s*\)")
fixed_input = re.compile(r"^[^\n#]*loadtxt\(['\"]([^'\"]*)['\"]\)", re.M)
changing_loadtxt_noparens = re.compile(r"^[^\n#]*loadtxt\(['\"]([^'\"]*)['\"]\s*%\s*([^\(\)\n]*)(\s*,[\w\s=]+)*\s*\)", re.M)
changing_loadtxt = re.compile(r"^[^\n#]*loadtxt\(['\"]([^'\"]*)['\"]\s*%\s*(\([^\)]*\))(\s*,[\w\s=]+)*\s*\)", re.M)
changing_input = re.compile(r"input:\s*['\"]([^'\"]*)['\"]\s*%\s*(\(.*\))")
list_comprehension_input = re.compile(r"input:\s*(\[.+for.*in.*\])\n")

def printrule(pyf, args, inputs, outputs):
    if len(outputs) == 0:
        return
    if len(inputs) == 0:
        print ': %s | *.pyc |> cd .. && python2 figs/%s %s |> %s' % (pyf, pyf, args, string.join(outputs))
        return
    print ': %s | *.pyc %s |> cd .. && python2 figs/%s %s |> %s' % (pyf, string.join(inputs), pyf, args, string.join(outputs))

def friendly_eval(code, context, local = None):
    try:
        return eval(code, local)
    except:
        print ("\nError evaluating '%s' from file %s" % (code, context))
        raise

pyfs = glob.glob('*.py')
if len(sys.argv) > 1:
    pyfs = sys.argv[1:]
for pyf in pyfs:
    f = open(pyf, 'r')
    contents = f.read()
    f.close()

    inputs = fixed_input.findall(contents)
    inputs += fixed_open_input.findall(contents)
    inputs = [i[5:] for i in inputs]

    imports = import_re.findall(contents)
    for i in imports:
        ipath = i+'.py'
        if os.path.exists(ipath):
            inputs.append(ipath)

    outputs = fixed_output.findall(contents)
    outputs = [x[0] for x in outputs]
    outputs += fixed_open.findall(contents) # add any files created with open(...)
    outputs = [o[5:] for o in outputs]
    argvals = arguments.findall(contents)
    if len(argvals) > 0:
        coutputs = changing_output.findall(contents)
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
            for v in friendly_eval(arg[1], pyf):
                for av in allvalues:
                    newav = av.copy()
                    newav[arg[0]] = v
                    newallvalues.append(newav)
            allvalues = newallvalues
        for a in allvalues:
            aa = commandlineformat % friendly_eval('(' + string.join(commandlineargs, ",") + ')', 'file '+ pyf, a)
            extrainputs = []
            extraoutputs = []
            for i in list_inputs:
                #print 'examining input:', i
                fnames = friendly_eval(i, pyf, a)
                #print 'found', fnames
                extrainputs += fnames
            for i in cinputs:
                #print 'examining input:', i[0], '%', '"%s"' % i[1]
                fname = i[0] % friendly_eval(i[1], pyf, a)
                extrainputs.append(fname)
            for o in coutputs:
                fname = o[0] % friendly_eval(o[1], pyf, a)
                #print 'generate', fname, 'from', aa
                extraoutputs.append(fname)
            extrainputs = [e[5:] for e in extrainputs]
            extraoutputs = [e[5:] for e in extraoutputs]
            printrule(pyf, aa, inputs + extrainputs, outputs + extraoutputs)
    else:
        printrule(pyf, "", inputs, outputs)
