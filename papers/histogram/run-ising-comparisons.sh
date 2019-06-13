set -ev

# This is the reference for 32 X 32 system
#python figs/yaml-reference.py /home/jordan/sad-monte-carlo/ising-samc-1e5-32-reference data/ 2048 48

# This is the reference for 128 X 128 system

# This is for getting out the update factor (Gamma) plots
#python figs/parse-yaml-out.py /home/jordan/weniger-cp-data/sad-50-slow-s1 50 sad-50-slow

# This is for generating the main figures
python figs/yaml-multi-comparison.py /home/jordan/ising-cp-data/ /home/jordan/deft/papers/histogram/data/ising-32-reference-lndos.dat ising/N32 32 2048 48 true 8 ising-sad-32 ising-samc-1e5-32 ising-samc-1e6-32 ising-samc-1e7-32 ising-wl-32 ising-wl-inv-t-32

# ising-cp-data not weniger-cp-data
