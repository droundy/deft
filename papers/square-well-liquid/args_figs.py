import argparse

parser = argparse.ArgumentParser(add_help = False)

parser.add_argument(
		'-ww', metavar='FLOAT', type=float, default=0,
    help='Well width relative to ball diameter')

parser.add_argument(
		'-ff', metavar='FLOAT', type=float, default=0,
    help='Filling fraction')

parser.add_argument(
		'-i', metavar='INT', type=int, default=0,
    help='Number of interactions')

parser.add_argument(
		'-walls', metavar='INT', type=int, default=0,
    help='Number of walls')

parser.add_argument(
		'-ktemax', metavar='FLOAT', type=float, default=150,
    help='Maximum kt in heat capacity plot')

parser.add_argument(
		'--show', action='store_true',
    help='Show generated RDF figure')

parser.add_argument(
		'--all', action='store_true',
    help='Generate all figures (except RDF)')
