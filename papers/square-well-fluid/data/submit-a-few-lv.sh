#!/bin/sh


set -v

cd data

set -e

python run-mc-liquid-vapor.py 1.3 0.1 100 10 0.8 pessimistic tmi
python run-mc-liquid-vapor.py 1.3 0.2 100 10 0.8 pessimistic tmi
python run-mc-liquid-vapor.py 1.3 0.3 100 10 0.8 pessimistic tmi
python run-mc-liquid-vapor.py 1.3 0.1  80 10 0.8 pessimistic tmi
python run-mc-liquid-vapor.py 1.3 0.2  80 10 0.8 pessimistic tmi
python run-mc-liquid-vapor.py 1.3 0.3  80 10 0.8 pessimistic tmi

python run-mc-liquid-vapor.py 1.3 0.1 100 10 0.8 pessimistic
python run-mc-liquid-vapor.py 1.3 0.2 100 10 0.8 pessimistic
python run-mc-liquid-vapor.py 1.3 0.3 100 10 0.8 pessimistic
python run-mc-liquid-vapor.py 1.3 0.1  80 10 0.8 pessimistic
python run-mc-liquid-vapor.py 1.3 0.2  80 10 0.8 pessimistic
python run-mc-liquid-vapor.py 1.3 0.3  80 10 0.8 pessimistic
