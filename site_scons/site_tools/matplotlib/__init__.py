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

import re, glob
from SCons.Script import *

def runpython(env, pyfile, args, inputs, outputs, py_chdir):
    py_chdir = str(py_chdir)
    if py_chdir == "" or py_chdir == ".":
        act = "python %s %s" % (str(pyfile), args)
        print "build", outputs, 'from', inputs, 'with', act
        env.Command(target = outputs,
                    source = [pyfile] + inputs,
                    action = act)
    else:
        pyrel = str(pyfile)[len(py_chdir)+1:]
        act = "cd %s && python %s %s" % (str(py_chdir), pyrel, args)
        goodoutputs = [py_chdir + '/' + o for o in outputs]
        okinputs = [py_chdir + '/' + i for i in inputs]
        goodinputs = [glob.glob(py_chdir + '/' + i) for i in inputs]
        print "build", goodoutputs, 'from', goodinputs, okinputs, 'with', act
        env.Command(target = goodoutputs,
                    source = [pyfile] + goodinputs,
                    action = act)


fixed_output = re.compile(r"savefig\(['\"](.*)['\"](\s*,[\w\s=]+)*\)")
changing1_output = re.compile(r"savefig\(['\"](.*)['\"]\s*%\s*\(\s*(\w+)\s*\)(\s*,[\w\s=]+)*\)")
arguments = re.compile(r"^#arg\s+(\w+)\s*=\s*(.*)$", re.M)

fixed_input = re.compile(r"^[^#]*loadtxt\(['\"](.*)['\"]\)", re.M)
changing1_input = re.compile(r"input:\s*['\"](.*)['\"]\s*%\s*\(\s*(\w+)\s*\)")

def Matplotlib(env, source, py_chdir = ""):
    if len(source) < 4 or source[len(source)-3:] != ".py":
        source = source + '.py'
    print 'examining', source
    node = File(source)
    contents = node.get_text_contents()

    outputs = fixed_output.findall(contents)
    outputs = [x[0] for x in outputs]
    inputs = fixed_input.findall(contents)
    if len(outputs) == 0:
        outputs1 = changing1_output.findall(contents)
        if outputs1 > 0:
            inputs1 = changing1_input.findall(contents)
            print 'inputs1 is', inputs1
            values = arguments.findall(contents)
            vmap = dict(values)
            if len(values) == 1:
                for v in eval(values[0][1]):
                    myinputs = []
                    myoutputs = []
                    print outputs1
                    for o in outputs1:
                        myoutputs.append(o[0] % v)
                    for i in inputs1:
                        myinputs.append(i[0] % v)
                    runpython(env, source, v, inputs + myinputs, outputs+myoutputs, py_chdir)
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
