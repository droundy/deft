set -ev

# rsync -av knightley:src/deft/papers/square-well-fluid/data/lv/*-movie data/lv/

python figs/animate-dos.py lv ww1.30-ff0.20-100x10 tmi toe
python figs/animate-histogram.py lv ww1.30-ff0.20-100x10 tmi toe

# python figs/animate-dos.py lv ww1.30-ff0.20-100x20 tmi toe
# python figs/animate-histogram.py lv ww1.30-ff0.20-100x20 tmi toe

# # python figs/animate-dos.py lv N0050 tmi toe
python figs/animate-dos.py lv N0500-0.1 tmi toe
python figs/animate-dos.py lv N0500-0.2 tmi toe
python figs/animate-dos.py lv N0500-0.3 tmi toe

python figs/animate-histogram.py lv N0500-0.1 tmi toe
python figs/animate-histogram.py lv N0500-0.2 tmi toe
python figs/animate-histogram.py lv N0500-0.3 tmi toe

totem figs/movies/lv/ww1.30-ff0.20-100x20-dos/movie.mp4 figs/movies/lv/ww1.30-ff0.20-100x10-dos/movie.mp4

totem figs/movies/lv/ww1.30-ff0.20-100x20-hist/movie.mp4 figs/movies/lv/ww1.30-ff0.20-100x10-hist/movie.mp4

totem figs/movies/lv/N0500-*-dos/movie.mp4

totem figs/movies/lv/N0500-*-hist/movie.mp4
