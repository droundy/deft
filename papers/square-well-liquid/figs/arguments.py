#!/usr/bin/python2
import argparse

parser = argparse.ArgumentParser(add_help = False)

parser.add_argument(
		'-ww', type=float, nargs='+', default=[1.3],
    help='Well width(s) relative to ball diameter')

parser.add_argument(
		'-ff', type=float, nargs='+', default=[0.3],
    help='Filling fraction(s)')

parser.add_argument(
		'-N', type=int, nargs='+', default=[1000],
    help='Number(s) of balls')

parser.add_argument(
		'-kT', type=int, nargs='+', default=[0],
    help='Number(s) of balls')

parser.add_argument(
		'-walls', type=int, default=0,
    choices=[0,1,2,3],
    help='Number(s) of dimensions in which to have walls')

parser.add_argument(
		'--nw', action='store_true', help="Don't use energy weights")

parser.add_argument(
		'--flat', action='store_true', help="Use flat histogram method")

parser.add_argument(
		'--show', action='store_true', help='Show figure(s)')
