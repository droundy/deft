import argparse

parser = argparse.ArgumentParser(add_help = False)

parser.add_argument(
		'-N', metavar='INT', type=int, nargs='+', default=[],
    help='Number(s) of balls')

parser.add_argument(
		'-ww', metavar='FLOAT', type=float, nargs='+', default=[],
    help='Well width(s) relative to ball diameter')

parser.add_argument(
		'-ff', metavar='FLOAT', type=float, nargs='+', default=[],
    help='Filling fraction(s)')

parser.add_argument(
		'-walls', metavar='INT', type=int, default=0,
    choices=[0,1,2,3],
    help='Number(s) of dimensions in which to have walls')

parser.add_argument(
		'--weights', action='store_true', help='Use energy weights')
