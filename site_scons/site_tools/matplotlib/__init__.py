"""`matplotlib`

Tool specific initialization for python scripts using matplotlib.
"""

#
# Copyright (c) 2013 by David Roundy
# Copyright (c) 2013 by Pawel Tomulik
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE

__docformat__ = "restructuredText"

import re, glob, string
from SCons.Script import *

def runpython(env, pyfile, args, inputs, outputs, py_chdir):
    py_chdir = str(py_chdir)
    if py_chdir == "" or py_chdir == ".":
        act = "python %s %s" % (str(pyfile), args)
        #print "build", outputs, 'from', inputs, 'with', act
        env.Command(target = outputs,
                    source = [pyfile] + inputs,
                    action = act)
    else:
        pyrel = str(pyfile)[len(py_chdir)+1:]
        act = "cd %s && python %s %s" % (str(py_chdir), pyrel, args)
        goodoutputs = [py_chdir + '/' + o for o in outputs]
        goodinputs = [glob.glob(py_chdir + '/' + i) for i in inputs]
        goodinputs = []
        for i in inputs:
            globi = glob.glob(py_chdir + '/' + i)
            if len(globi) > 0:
                goodinputs.append(globi)
            elif not '*' in i and not '?' in i:
                goodinputs.append(py_chdir + '/' + i)
        #print "build", goodoutputs, 'from', goodinputs, 'with', act
        env.Command(target = goodoutputs,
                    source = [pyfile] + goodinputs,
                    action = act)

fixed_output = re.compile(r"savefig\(['\"](.*)['\"](\s*,[\w\s=]+)*\)")
changing_output = re.compile(r"savefig\(['\"](.*)['\"]\s*%\s*(\(.*\))(\s*,[\w\s=]+)*\s*\)")
arguments = re.compile(r"^#arg\s+(\w+)\s*=\s*(.*)$", re.M)

fixed_input = re.compile(r"^[^#]*loadtxt\(['\"](.*)['\"]\)", re.M)
changing_input = re.compile(r"input:\s*['\"](.*)['\"]\s*%\s*(\(.*\))")

def Matplotlib(env, source, py_chdir = ""):
    if len(source) < 4 or source[len(source)-3:] != ".py":
        source = source + '.py'
    node = File(source)
    contents = node.get_text_contents()

    outputs = fixed_output.findall(contents)
    outputs = [x[0] for x in outputs]
    inputs = fixed_input.findall(contents)
    argvals = arguments.findall(contents)
    if len(argvals) > 0:
        coutputs = changing_output.findall(contents)
        cinputs = changing_input.findall(contents)
        allvalues = [{}]
        commandlineformat = ""
        commandlineargs = []
        for arg in argvals:
            commandlineformat += " %s"
            commandlineargs.append(arg[0])
            newallvalues = []
            for v in eval(arg[1]):
                for av in allvalues:
                    newav = av.copy()
                    newav[arg[0]] = v
                    newallvalues.append(newav)
            allvalues = newallvalues
        for a in allvalues:
            aa = commandlineformat % eval('(' + string.join(commandlineargs, ",") + ')', a)
            extrainputs = []
            extraoutputs = []
            for i in cinputs:
                fname = i[0] % eval(i[1], a)
                extrainputs.append(fname)
            for o in coutputs:
                fname = o[0] % eval(o[1], a)
                extraoutputs.append(fname)
            runpython(env, source, aa, inputs + extrainputs, outputs + extraoutputs, py_chdir)
    else:
        runpython(env, source, "", inputs, outputs, py_chdir)

def generate(env):
    """Add Builders and construction variables to the Environment"""
    env.AddMethod(Matplotlib)

def exists(env):
    return True

# Local Variables:
# # tab-width:4
# # indent-tabs-mode:nil
# # End:
# vim: set syntax=python expandtab tabstop=4 shiftwidth=4 nospell:
