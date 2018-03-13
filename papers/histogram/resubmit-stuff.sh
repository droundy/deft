#!/bin/bash

set -ev

ls data/lv

echo I am where I want to be

if ls -d data/lv/ww1.30-ff0.22-100x10-sad3*; then
    echo Please delete data first
    exit 1
fi

python data/run-sticky-wall.py 1.3 0.22 100 10 0.6 sad3
python data/run-sticky-wall.py 1.3 0.22 100 10 0.6 sad3 1
python data/run-sticky-wall.py 1.3 0.22 100 10 0.6 sad3 2

if ls -d data/lv/ww1.30-ff0.22-100x10-wltmmc-0.8-1e-10-s1*; then
    echo Please delete data first
    exit 1
fi

python data/run-sticky-wall.py 1.3 0.22 100 10 0.6 wltmmc 0.8 1e-10 2900 816 1

python data/run-sticky-wall.py 1.3 0.22 100 10 0.6 tmi3 golden
