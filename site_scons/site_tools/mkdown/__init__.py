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

import re, string
import markdown as mmdd
from SCons.Script import *

def mkdown(target, source, env):
    base = 'doc/'
    if str(source[0])[:4] == base:
        base = ''

    mkstr = string.Template(source[0].get_text_contents()).safe_substitute(base = base)
    sidebar = string.Template(source[1].get_text_contents()).safe_substitute(base = base)
    template = string.Template(source[2].get_text_contents())

    htmlfile = str(target[0])
    f = open(str(target[0]), 'w')
    f.write(template.safe_substitute(base = base,
                                     content = mmdd.markdown(mkstr),
                                     sidebar = mmdd.markdown(sidebar)))
    f.close()

wikilinkre = re.compile(r"\[\[([^\]]*)\]\]")
linkre = re.compile(r"\[[^\]]*\]\(([^)]*)\)")

def Markdown(env, source):
    if len(source) < 4 or source[len(source)-3:] != ".md":
        source = source + '.md'
    node = File(source)
    contents = node.get_text_contents()

    wikilinks = wikilinkre.findall(contents)
    links = linkre.findall(contents)
    nodes = env.Command(target = source[:len(source)-3]+'.html',
                        source = [source, '.layout/sidebar.md', '.layout/page.html'],
                        action = mkdown)
    for link in wikilinks + links:
        if len(link) > 5 and link[len(link)-5:] == '.html':
            link = link[:len(link)-5]
        nodes += Markdown(env, 'doc/' + link)
    return nodes

def generate(env):
    """Add Builders and construction variables to the Environment"""
    env.AddMethod(Markdown)

def exists(env):
    return True
