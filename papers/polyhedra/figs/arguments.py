import argparse

parser = argparse.ArgumentParser(add_help = False)

parser.add_argument(
  'ff', metavar='ff', type=float, help='filling fraction')

parser.add_argument(
  '-N', metavar='N', type=int,default=0,
  help="""number of polyhedra, if not supplied then the first
  file with the proper filling fraction will be used""")

parser.add_argument(
  '-s', '--shape', metavar='S', default='truncated_tetrahedron',
  choices=['cube', 'tetrahedron', 'truncated_tetrahedron', 'cuboid'],
  help='type of polyhedron, defaults to truncated_tetrahedron')

parser.add_argument(
  '-r', '--ratio', type=float, default=0,
  help='ratio of edge lengths, required for cuboids and not used otherwise')

parser.add_argument(
  '-p', '--periodic', action='store_true',
  help='will use periodic cell - defaults to walls otherwise')

parser.add_argument(
  '--hide', action='store_true',
  help='will just save the plot and won\'t display it')
