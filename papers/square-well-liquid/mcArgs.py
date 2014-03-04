import argparse

parser = argparse.ArgumentParser(add_help = False)

parser.add_argument(
		'-ff', metavar='FLOAT', type=float, nargs='+', default=[],
    help='Filling fraction(s)')

parser.add_argument(
		'-ww', metavar='FLOAT', type=float, nargs='+', default=[],
    help='Well width(s) relative to ball diameter')

parser.add_argument(
		'-N', metavar='INT', type=int, default=1000,
    help='Number of balls')

parser.add_argument(
		'-initialize', metavar='INT', type=int, default=500000,
		help='Number of iterations to run for initialization')

parser.add_argument(
		'-iterations', metavar='INT', type=int, default=2500000,
		help='Number of simulation iterations')

parser.add_argument(
		'-walls', metavar='INT', type=int, default=0, choices=[0,1,2,3],
    help='Number of walls')

parser.add_argument(
		'--weights', action='store_true', help='Use energy weights?')
