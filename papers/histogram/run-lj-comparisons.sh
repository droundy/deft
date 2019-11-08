set -ev

# This is how to run and compute the lj reference
python3 figs/lj-reference.py /home/droundy/src/sad-monte-carlo/tiny-lj-benchmark-0.001 lj-31

# This is how to run lj-multi-comparison.py
python3 figs/lj-multi-comparison.py /home/jordan/sad-monte-carlo/ data/lj-31-reference-lndos.dat lj/ 31 1 /home/jordan/sad-monte-carlo/lj-sad-31-bin001
