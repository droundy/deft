#!/usr/bin/python
# this should run and plot the free energy simulation eventually

import numpy
import sys
from ast import literal_eval as make_tuple

FILENAME = "periodic-ww1.00-ff0.30-N10-sf{0:f}-g.dat"

def run(args):
	sf = 0.99
	data = read_data_file_to_dict(FILENAME.format(sf))
	print data


def read_data_file_to_dict(filepath):
	meta_data = {}
	with open(filepath, 'r') as f:
		for line in f:
			if line.startswith('#') and ':' in line:
				sanitized = map(lambda s: s.strip(' #'), line.split(':'))
				meta_data[sanitized[0]] = num(sanitized[1])
	return meta_data


def num(s):
	# todo the way I hadnle the tuple case is terrible and makes me sad
	# I'm going to be lying awake thinking about it
    try:
        return int(s)
    except ValueError:
    	try:
    		return float(s)
    	except ValueError:
        	return map(num, make_tuple(s))




if __name__ == "__main__":
    run(sys.argv[1:])
