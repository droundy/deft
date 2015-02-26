#!/usr/bin/python
# this should run and plot the free energy simulation eventually

import numpy as np
import matplotlib.pyplot as plt
import sys, os, math
from subprocess import call
from ast import literal_eval as make_tuple

FILENAME = "periodic-ww1.00-ff{0}-N10-sf{1:f}"


def main(args):
	# run iteration
	# parse meta data
	# start iteration at "shrunk" filling fraction
	# goto start

	# here until line 60 I wrote while half asleep, don't judge me.
	ffs = []
	success_ratios = []

	sf = 0.99
	ff = 0.30

	ffs.append(ff)
	success_ratios.append(1) # for now. should be the absolute ratio

	for i in xrange(20):
		filename = FILENAME.format(ff, sf)

		call(['../free-energy-monte-carlo',
			'--sf', str(sf),
			'--ff', str(ff),
			'--filename', filename,
			])

 		# at some point I need to chdir so that the files aren't down a directory on accident.
		data = read_data_file_to_dict(os.path.join(os.getcwd(), 'free-energy-data', filename+"-g.dat"))
		next_ff = data['ff_small']
		total_checks = data['total checks of small cell']
		valid_checks = data['total valid small checks']
		success_ratio = (valid_checks * 1.0)/total_checks

		ffs.append(next_ff)
		success_ratios.append(success_ratio)

		ff = next_ff

	print ffs
	print success_ratios

	#do plot
	energy = -np.cumsum(map(math.log, success_ratios))
	print energy
	plt.scatter(np.array(ffs), energy)
	plt.ylabel('F/kT')
	plt.xlabel('ff')
	plt.show()



# I was awake when I wrote this part
def read_data_file_to_dict(filepath):
	meta_data = {}
	with open(filepath, 'r') as f:
		for line in f:
			if line.startswith('#') and ':' in line:
				sanitized = map(lambda s: s.strip(' #'), line.split(':'))
				meta_data[sanitized[0]] = num(sanitized[1])
	return meta_data


def num(s):
	# todo the way I handle the tuple case is terrible and makes me sad
	# I'm going to be lying awake thinking about it
    try:
        return int(s)
    except ValueError:
    	try:
    		return float(s)
    	except ValueError:
        	return map(num, make_tuple(s))




if __name__ == "__main__":
    main(sys.argv[1:])
