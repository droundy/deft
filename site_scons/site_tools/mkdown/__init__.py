"""`markdowntool`

Tool specific initialization for creating website with markdown
"""

#
# Copyright (c) 2013 by David Roundy
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

import re, string, os
import markdown as mmdd
from SCons.Script import *


toupload = set()
def Upload(env, source):
    source = str(source)
    global toupload
    toupload.add('.layout/style.css')
    toupload.add('papers/pair-correlation/figs/pretty-4.svg') # background
    toupload.add(source)
    Depends('upload', source)

def doupload(target, source, env):
    host = 'science.oregonstate.edu'
    path = 'public_html/deft'
    os.system('ssh %s rm -rf %s' % (host, path))
    madedir = set()
    for f in toupload:
        dirname = os.path.join(path, os.path.dirname(f))
        if not dirname in madedir:
            print 'creating directory', dirname
            os.system('ssh %s mkdir -p %s' % (host, dirname))
            madedir.add(dirname)
        os.system('scp %s %s:%s/' % (f, host, os.path.join(path, os.path.dirname(f))))

def mkdown(target, source, env):
    base = ''
    sourcename = str(source[0])
    sourcedirname = os.path.dirname(sourcename)
    while sourcedirname != '':
        base = '../' + base
        sourcedirname = os.path.dirname(sourcedirname)

    mkstr = source[0].get_text_contents()
    titlere = re.compile(r"^\s*#\s*([^\n]*)(.*)", re.DOTALL)
    title = titlere.findall(mkstr)
    if len(title) == 0:
        title = "Deft"
    else:
        mkstr = title[0][1]
        title = mmdd.markdown(title[0][0])
        title = title[3:len(title)-4]
    mkstr = string.Template(mkstr).safe_substitute(base = base)
    sidebar = string.Template(source[1].get_text_contents()).safe_substitute(base = base)
    template = string.Template(source[2].get_text_contents())

    htmlfile = str(target[0])
    f = open(str(target[0]), 'w')
    f.write(template.safe_substitute(base = base,
                                     title = title,
                                     content = mmdd.markdown(mkstr, extensions=['mathjax']),
                                     sidebar = mmdd.markdown(sidebar, extensions=['mathjax'])))
    f.close()

linkre = re.compile(r"\[[^\]]*\]\(([^)]*)\)")
imgre = re.compile(r"<\s*[iI][mM][gG]\s+[sS][rR][cC]='([^']+)'")

def Markdown(env, source, sofar = set()):
    # sofar is the set of files we have already figured out how to
    # build.
    if len(source) < 4 or source[len(source)-3:] != ".md":
        source = source + '.md'
    if source in sofar:
        return [] # we already handled this!
    relto = os.path.dirname(source)
    node = File(source)
    contents = node.get_text_contents()

    links = linkre.findall(contents)

    sblinks = linkre.findall(string.Template(File('.layout/sidebar.md').get_text_contents()
                                             ).safe_substitute(base = ''))
    deps = []
    sofar.add(source)
    for link in sblinks + [os.path.normpath(os.path.join(relto,l)) for l in links]:
        origlink = link
        if len(link) > 5 and link[len(link)-5:] == '.html':
            link = link[:len(link)-5]
        if ':' in link:
            deps = deps # nothing to do
        elif (len(link) > 4 and (link[len(link)-4:] == '.pdf' or link[len(link)-4:] == '.svg')):
            env.Upload(origlink)
        else:
            deps += Markdown(env, link, sofar)

    for img in imgre.findall(contents):
        if ':' in img:
            deps = deps # nothing to do
        else:
            deps += [img]
            env.Upload(img)

    env.Upload(source[:len(source)-3]+'.html')
    return env.Command(target = source[:len(source)-3]+'.html',
                       source = [source,
                                 '.layout/sidebar.md',
                                 '.layout/page.html',
                                 'site_scons/mdx_mathjax.py'],
                       action = mkdown) + deps

def generate(env):
    """Add Builders and construction variables to the Environment"""
    env.AddMethod(Upload)
    env.AddMethod(Markdown)
    AlwaysBuild(env.Command(target = 'upload',
                            source = [],
                            action = doupload))

def exists(env):
    return True
