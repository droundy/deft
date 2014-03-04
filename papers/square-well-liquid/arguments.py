import argparse

parser = argparse.ArgumentParser(add_help = False)

parser.add_argument(
		'-ff', metavar='ff', type=float, nargs='+', default=[],
    help='Filling fraction(s)')

parser.add_argument(
		'-interaction_scale', metavar='interaction_scale',
        type=float, nargs='+', default=[], help='Interaction scale(s)')

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
		'-walls', metavar='INT', type=int, default=0,
    help='Number of walls')

parser.add_argument(
		'--hide', action='store_true',
    help='Don\'t display plots when saved')
