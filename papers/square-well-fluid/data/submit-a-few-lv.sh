#!/bin/sh


set -v

cd data

set -e

python run-mc-liquid-vapor.py 1.3 0.2 100 10 0.5 tmi
python run-mc-liquid-vapor.py 1.3 0.2  80 10 0.5 tmi

python run-mc-liquid-vapor.py 1.3 0.2 100 10 0.5 toe
python run-mc-liquid-vapor.py 1.3 0.2  80 10 0.5 toe

python run-mc-liquid-vapor.py 1.3 0.2 100 10 0.5 tmmc
python run-mc-liquid-vapor.py 1.3 0.2  80 10 0.5 tmmc
