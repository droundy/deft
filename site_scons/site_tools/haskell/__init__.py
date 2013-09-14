# -*- mode:python; coding:utf-8; -*-

#  A SCons tool to enable compilation of Haskell in SCons.
#
#  Copyright © 2009 Russel Winder
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY
# KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE
# WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
# LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
# OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
# WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#  This module is derived from the code to be found at http://www.scons.org/wiki/GhcBuilder.  That code has
#  no attribution and no licence statement, so the standard SCons licence has been applied.

from os.path import abspath, join
from SCons.Builder import Builder
from SCons.Scanner import Scanner
from SCons.Script  import *
import re

impdecl_regexp = re.compile(r"[\s]*import[\s]*(qualified)?[\s]*([\w\.]*).*")
foreign_decl_regexp = re.compile(r"[\s]*foreign[\s]*(export|import)[\s]*ccall(.*\"[\s]*wrapper[\s]*\")?.*::.*")

def find_foreign_declarations_that_require_stubs(node):
        foreigns = foreign_decl_regexp.findall(node.get_contents())
        for foreign in foreigns:
                if foreign[0] == "export":
                        return True
                if foreign[1] != "":
                        return True
        return False

def ghc_emitter(target, source, env):
        src_name = str(source[0])
        hi_name = src_name.replace(source[0].suffix, ".hi")
        target.append(hi_name) # mark .hi file as a target, just like .o file!
        stubs = []
        if find_foreign_declarations_that_require_stubs(source[0]):
                stub_name = str(source[0]).replace(source[0].suffix, "_stub")
                stubs.append(stub_name + ".o")
                env.Clean(target, stub_name + ".c")
                env.SideEffect(stub_name + ".h", target)
                env.Clean(target, stub_name + ".h")
        return (target + stubs, source)

def ghc_make_emitter(target, source, env):
        src_name = str(source[0])
        hi_names = [ src_name.replace(source[0].suffix, extension) for extension in [ '.o', '.hi' ] ]
        env.SideEffect(hi_names, target)
        env.Clean(target, hi_names)
        return (target, source)

def ghc_path_function(env, dir, targets, sources):
        if dir.path != ".":
                env.AppendUnique(HSSEARCHPATH = [dir.path])
        return tuple([dir.path] + env["HSSEARCHPATH"])

def ghc_scanner_function(node, env, path):
        from os.path import dirname, exists
        imports = impdecl_regexp.findall(node.get_contents())
        modules = map(lambda (qualified, module) : module, imports)
        interfaces = []
        for module in modules:
                module = module.replace(".", "/")
                interface = module + ".hi"
                hs_file   = module + ".hs"
                lhs_file  = module + ".lhs"
                for dir in path:
                        if exists(join(dir, hs_file)) or exists(join(dir, lhs_file)):
                                interfaces.append(interface)
                                break
        return interfaces


def ghc_make_scanner_function(node, env, path):
        from os.path import dirname, exists
        imports = impdecl_regexp.findall(node.get_contents())
        modules = map(lambda (qualified, module) : module, imports)
        interfaces = []
        for module in modules:
                module = module.replace(".", "/")
                interface = module + ".hi"
                o = module + ".o"
                hs_file   = module + ".hs"
                lhs_file  = module + ".lhs"
                for dir in path:
                        if exists(join(dir, hs_file)) or exists(join(dir, lhs_file)):
                                interfaces.append(interface)
                                interfaces.append(o)
                                break
        return interfaces

def exists(env):
        return WhereIs("ghc")

def generate(env):
        def _ghc_searchpath_opts(paths):
                if paths:
                        return reduce(lambda list, path: " -i" + path, paths, "")
                else:
                        return ""
        env["_ghc_searchpath_opts"] = _ghc_searchpath_opts
        def _ghc_package_opts(packages):
                return reduce(lambda list, package: list + " -package " + package, packages, "")
        env["_ghc_package_opts"] = _ghc_package_opts

        env["HSC"] = "ghc"
        env["HSCFLAGS"] = []
        env["HSLINKFLAGS"] = []
        env["HSSEARCHPATH"] = []
        env["HSPACKAGES"] = []
        env["_HSPACKAGE_OPTS"] = "${_ghc_package_opts(HSPACKAGES)}"
        env["_HSSEARCHPATH_OPTS"] = "${_ghc_searchpath_opts(HSSEARCHPATH)}"

        ghc_scanner = Scanner(
                function = ghc_scanner_function,
                skeys = [".hs", ".lhs"],
                path_function = ghc_path_function
                )
        ghc_make_scanner = Scanner(
                function = ghc_make_scanner_function,
                skeys = [".hs", ".lhs"],
                path_function = ghc_path_function
                )

        ghc_c_compiler = Builder(
                action = "$HSC $HSCFLAGS -c -o $TARGET $SOURCE",
                src_suffix = [ ".c" ],
                suffix = ".o",
                single_source = True
                )

        ghc_compiler = Builder(
                action = "$HSC $HSCFLAGS $_HSSEARCHPATH_OPTS -c -o $TARGET $SOURCE",
                src_suffix = [ ".hs", ".lhs" ],
                suffix = ".o",
                single_source = True,
                emitter = ghc_emitter,
                source_scanner = ghc_make_scanner
                )

        ghc_linker = Builder(
                action = "$HSC $HSLINKFLAGS $_HSPACKAGE_OPTS -o $TARGET $SOURCES",
                src_suffix = ".o",
                suffix = "$PROGSUFFIX",
                src_builder = [ ghc_compiler, ghc_c_compiler ]
                )

        ghc_make = Builder(
                action = "$HSC --make $HSCFLAGS $_HSSEARCHPATH_OPTS -o $TARGET $SOURCE",
                src_suffix = [ ".hs", ".lhs" ],
                suffix = "",
                single_source = True,
                emitter = ghc_make_emitter,
                source_scanner = ghc_make_scanner
                )

        env.Append( BUILDERS = {
                "HaskellProgram" : ghc_linker,
                "HaskellObject" : ghc_compiler ,
                "HaskellMake" : ghc_make
                } )
