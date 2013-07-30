"""`gnuplot`

Tool specific initialization for gnuplot.
"""

#
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

import SCons.Builder

_null = SCons.Builder._null

class _GplotRelTo(object):
    """Given a sequence of ``nodes`` return their paths relative to predefined
    ``base``.
    """
    def __init__(self, base):
        """Initializes the functional object
        
        **Arguments**

            - *base* - scons filesystem node representing base file or dir,
        """
        import SCons.Util
        self.base = base

    def __call__(self, nodes, *args, **kw):
        """Given a sequence of ``nodes`` return list of their paths relative to
           ``self.base``."""
        return [ self.base.rel_path(node) for node in nodes ]
   
def _GplotFvars(fdict, base):
    """Prepare list of gnuplot variables contaning file names.

    **Arguments**
        
        - *fdict* - dictionary with files (nodes) as returned by ``_gplot_fdict()``,
        - *base* - base directory (node),

    **Returns**

        a list of ``'variable="path"'`` strings, the ``variable`` s are keys
        from ``fdict`` and ``path``s are file names relative to ``base``
    """
    import SCons.Util
    if not fdict: return []
    return [ "'%s=\"%s\"'" % (k, base.rel_path(v)) for k,v in fdict.items() ]


def _gplot_arg2nodes(env, args, *args2, **kw):
    """Helper function. Convert arguments to a list of nodes.
    
    This function works similarly to ``env.arg2nodes()`` except it handles
    also dictionaries.
    
    **Arguments**

        - *env*   - SCons Environment object,
        - *args*  - arguments representing one or more files,
        - *args2* - other positional arguments (passed to ``env.arg2nodes()``),
        - *kw*    - keyword arguments (passed to ``env.arg2nodes()``).

    **Return**
        
        returns list of nodes.
    """
    import SCons.Util
    if SCons.Util.is_Dict(args):
        return env.arg2nodes(args.values(), *args2, **kw)
    else:
        return env.arg2nodes(args, *args2, **kw)

def _gplot_arg2nodes_dict(env, args, name = None, *args2, **kw):
    """Helper function. Convert arguments to a dict with file nodes.

    **Arguments**

        - *env*  - SCons Environment object,
        - *args* - arguments representing one or more files,
        - *name* - default name used when `args` is not a dictionary,
        - *args2* - other positional arguments (passed to ``env.arg2nodes()``),
        - *kw*    - keyword arguments (passed to ``env.arg2nodes()``).

    **Returns**
        
        dictionary of type ``{ 'key' : node }``,
    """
    import SCons.Util
    if SCons.Util.is_Dict(args):
        keys = args.keys()
        vals = args.values()
        nodes = dict( zip( keys, env.arg2nodes(vals, *args2, **kw) ) )
    elif SCons.Util.is_Sequence(args):
        keys = [ '%s%d' % (name, i+1) for i in xrange(0,len(args)) ]
        nodes = dict( zip( keys, env.arg2nodes(args, *args2, **kw) ) )
    elif args:
        nodes = { '%s%d' % (name,1) : env.arg2nodes(args,*args2,**kw)[0] }
    else:
        nodes = {}

    return nodes


def _gplot_fdict(env):
    """Helper function. Make a dictionary containing gnuplot command-line
    variables with input/output file names.

        Constuction variables ``$gp_inputs``, ``$gp_outputs``,
        and ``$gp_extoutputs`` are processed to create the specific dictionary.
       
        **Arguments**

            - *env* - SCons Environment object,

        **Returns**
            
            returns a ``{ 'name' : 'value' }`` dict where 'name's are
            gnuplot variable names and 'value's are corresponding values,

    """
    def fdict2(env, triples, f = None):
        nodes = {}
        for t in triples:
            try: args = env[t[0]]
            except KeyError: pass
            else:
                name = env.subst('$%s' % t[1])
                if not name: name = t[2]
                nodes.update(_gplot_arg2nodes_dict(env, args, name))

        if f is not None:
            for key, node in nodes.items():
                nodes[key] = f(node)
        return nodes
    result = {}

    triples = [ ('gp_outputs', 'GPLOTOUTVAR', 'output'),
                ('gp_extoutputs', 'GPLOTEOUTVAR', 'eoutput') ]
    result.update(fdict2(env, triples))
    triples = [ ('gp_inputs', 'GPLOTINVAR', 'input') ]
    result.update(fdict2(env, triples, lambda n : n.srcnode() ))
    return result

def _gplot_scan_for_outputs(env, base, source):
    """Helper function. Scan source files for 'set output' gnuplot commands.

    **Arguments**
        
        - *env* - the scons Environment object,
        - *base* - base directory (node) for the file names obtained from source,
        - *source* - list of source nodes to be scanned.

    **Return**
        
        list of output files (as nodes) extracted from the source files
    """
    import re
    _re = r'^\s*set\s+output\s+[\'"]([^\n\r\'"]+)[\'"](?:\s*;\s*)*#?.*$'
    _re = re.compile(_re, re.M)

    nodes = []
    # extract file names
    for src in source:
        contents = src.get_text_contents()
        names = _re.findall(contents)
        # convert to nodes
        nodes.extend(env.File(names, base))
    nodes = list(set(nodes))
    return nodes

def _GplotScanner(node, env, path, arg):
    """Scan gnuplot script for implicit dependencies.
    
    This scanner also handles the ``gp_inputs`` parameter.
    """
    import re
    dat_regexp = re.compile(r"^[\s]*'figs/([\w\.\-]+dat)'.*", re.M)
    csv_regexp = re.compile(r"^[\s]*'figs/([\w\.\-]+csv)'.*", re.M)
    #depfile_regexp = re.compile("(.*)")

    deps = dat_regexp.findall(node.get_contents())
    deps += csv_regexp.findall(node.get_contents())

    # Add input files to implicit dependencies
    try: inputs = env['_gp_input_nodes']
    except KeyError: pass
    else: deps.extend(inputs)

    return deps


def _GplotEmitter(target, source, env):
    """Append gp_outputs and gp_eoutputs to target list. The emitter also
    prepares the list of implicit dependencies for further processing in the
    scanner (see `_GplotScanner()`).
    """

    # scan source files for outpus
    outnodes = _gplot_scan_for_outputs(env, env['_gp_chdir'], source)

    try: outputs2 = env['gp_outputs']
    except KeyError: pass
    else: 
        outnodes2 =  _gplot_arg2nodes(env, outputs2)
        outnodes.extend([n for n in outnodes2 if n not in outnodes])

    try: outputs2 = env['gp_extoutputs']
    except KeyError: pass
    else:
        outnodes2 =  _gplot_arg2nodes(env, outputs2)
        outnodes.extend([n for n in outnodes2 if n not in outnodes])

    return  target + outnodes, source 

class _GplotBuilderObject (SCons.Builder.BuilderBase):
    """Gnuplot builder object"""    

    def _execute(self, env, target, source, *args):

        # Prepare our environment override a little bit
        
        try: inputs = env['gp_inputs']
        except KeyError: pass
        else: env['_gp_input_nodes'] = _gplot_arg2nodes(env, inputs)

        env['_gp_fdict'] = _gplot_fdict(env)
        sup = super(_GplotBuilderObject, self)
        return sup._execute(env, target, source, *args)

    def __call__(self, env, target=None, source=None, chdir=_null, **kw):
        import SCons.Node.FS

        # - by default change dir to the directory of calling SCons script,
        # - if gp_chdir is a string, interpret it as a path relative to the
        #   calling SCons script,
        # - if gp_chdir is None/False, revert SCons default behavior: run the
        #   command from the directory of the top level SConstruct.
        try: gp_chdir = kw['gp_chdir']
        except KeyError: gp_chdir = True # do chdir by default

        # I'think here is the place to convert gp_chdir to node,
        if gp_chdir:
            if SCons.Util.is_String(gp_chdir):
                kw['_gp_chdir'] = env.Dir(gp_chdir)
            elif not isinstance(gp_chdir, SCons.Node.FS.Base):
                # this is our default behavior
                kw['_gp_chdir'] = env.fs.getcwd()
            else:
                kw['_gp_chdir'] = gp_chdir
        else:
            kw['_gp_chdir'] = env.Dir('#') 
       
        if target is None: target = []

        sup = super(_GplotBuilderObject, self)
        return sup.__call__(env, target, source, chdir, **kw)

def _GplotBuilder(**kw):
    """A factory for gnuplot builder objects"""
    if 'action' in kw:
        kw['action'] = SCons.Action.Action(kw['action']) 
    return _GplotBuilderObject(**kw)

def _detect_gnuplot(env):
    if env.has_key('GNUPLOT'):
        return env['GNUPLOT']
    return env.WhereIs('gnuplot')

def generate(env):
    """Add Builders and construction variables to the Environment"""
    import SCons.Builder, SCons.Script

    gnuplot = _detect_gnuplot(env)
    if not gnuplot: gnuplot = 'gnuplot'
    env['GNUPLOT'] = gnuplot

    fvars = '$( ${_concat( "%s " % GPLOTVARPREFIX, ' \
          + '_GplotFvars( _gp_fdict, _gp_chdir), ' \
          + 'GPLOTVARSUFFIX, __env__ )} $)'

    srcs  = '$( ${_concat( "", SOURCES, "", __env__, ' \
          + '_GplotRelTo(_gp_chdir))} $)'

    com   = 'cd $_gp_chdir && $GNUPLOT $GNUPLOTFLAGS %s %s' % (fvars, srcs)
    env.SetDefault( GPLOTSUFFIX     = '.gp',
                    GPLOTINVAR      = 'input',
                    GPLOTOUTVAR     = 'output',
                    GPLOTEOUTVAR    = 'eoutput',
                    GPLOTVARPREFIX  = "-e",
                    GPLOTVARSUFFIX  = "",
                    _GplotFvars     = _GplotFvars,
                    _GplotRelTo     = _GplotRelTo,
                    GNUPLOTCOM      = com,
                    GNUPLOTCOMSTR   = '')
    try:
        env['BUILDER']['GplotGraph']
    except KeyError:
        scanner = SCons.Script.Scanner( function = _GplotScanner, argument = None )
        builder = _GplotBuilder( action = '$GNUPLOTCOM',
                                 src_suffix = '$GPLOTSUFFIX',
                                 emitter = _GplotEmitter,
                                 source_scanner = scanner )
        env.Append( BUILDERS = { 'GplotGraph' : builder } )

def exists(env):
    return _detect_gnuplot(env)

# Local Variables:
# # tab-width:4
# # indent-tabs-mode:nil
# # End:
# vim: set syntax=python expandtab tabstop=4 shiftwidth=4 nospell:
