set -ev

# This is for getting out the update factor (Gamma) plots
# This is the reference for 32 X 32 system
python figs/parse-ising-out.py /home/jordan/ising-cp-data/ising-wl-32-s1 32 ising-wl-32
python figs/parse-ising-out.py /home/jordan/ising-cp-data/ising-sad-32-s1 32 ising-sad-32
python figs/parse-ising-out.py /home/jordan/ising-cp-data/ising-wl-inv-t-32-s1 32 ising-wl-inv-t-32

# This is the reference for 128 X 128 system


# This is the reference for 32 X 32 system
python figs/yaml-multi-comparison.py /home/jordan/ising-cp-data/ /home/jordan/deft/papers/histogram/data/ising-32-reference-lndos.dat ising/N32 32 2040 48 true 8 ising-sad-32 ising-samc-1e5-32 ising-samc-1e6-32 ising-samc-1e7-32 ising-wl-32 ising-wl-inv-t-32
# This is the reference for 128 X 128 system
python figs/yaml-multi-comparison.py /home/jordan/ising-cp-data/ /home/jordan/deft/papers/histogram/data/ising-sad-128-s1-reference-lndos.dat ising/N128 128 32760 100 true 8 ising-sad-128 ising-samc-1e5-128 ising-samc-1e6-128 ising-samc-1e7-128 ising-wl-128 ising-wl-inv-t-128

# ising-cp-data not weniger-cp-data
