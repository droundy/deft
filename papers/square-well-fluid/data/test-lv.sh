
set -ev

MINT=0.8

python run-mc-liquid-vapor.py 1.3 0.2 100 10 $MINT tmi
python run-mc-liquid-vapor.py 1.3 0.2 100 10 $MINT toe
# python run-mc-liquid-vapor.py 1.3 0.2 100 20 $MINT tmi
# python run-mc-liquid-vapor.py 1.3 0.2 100 20 $MINT toe

srun -J N0500-0.1-tmi ../../../liquid-vapor-monte-carlo --movies --N 500 --ff 0.1 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0500-0.1-tmi > lv/N0500-0.1-tmi.out &
srun -J N0500-0.2-tmi ../../../liquid-vapor-monte-carlo --movies --N 500 --ff 0.2 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0500-0.2-tmi > lv/N0500-0.2-tmi.out &
srun -J N0500-0.3-tmi ../../../liquid-vapor-monte-carlo --movies --N 500 --ff 0.3 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0500-0.3-tmi > lv/N0500-0.3-tmi.out &

srun -J N0500-0.1-toe ../../../liquid-vapor-monte-carlo --movies --N 500 --ff 0.1 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0500-0.1-toe > lv/N0500-0.1-toe.out &
srun -J N0500-0.2-toe ../../../liquid-vapor-monte-carlo --movies --N 500 --ff 0.2 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0500-0.2-toe > lv/N0500-0.2-toe.out &
srun -J N0500-0.3-toe ../../../liquid-vapor-monte-carlo --movies --N 500 --ff 0.3 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0500-0.3-toe > lv/N0500-0.3-toe.out &

srun -J N0050-0.1-tmi ../../../liquid-vapor-monte-carlo --movies --N 50 --ff 0.1 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0050-0.1-tmi > lv/N0050-0.1-tmi.out &
srun -J N0050-0.2-tmi ../../../liquid-vapor-monte-carlo --movies --N 50 --ff 0.2 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0050-0.2-tmi > lv/N0050-0.2-tmi.out &
srun -J N0050-0.3-tmi ../../../liquid-vapor-monte-carlo --movies --N 50 --ff 0.3 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0050-0.3-tmi > lv/N0050-0.3-tmi.out &

srun -J N0050-0.1-toe ../../../liquid-vapor-monte-carlo --movies --N 50 --ff 0.1 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0050-0.1-toe > lv/N0050-0.1-toe.out &
srun -J N0050-0.2-toe ../../../liquid-vapor-monte-carlo --movies --N 50 --ff 0.2 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0050-0.2-toe > lv/N0050-0.2-toe.out &
srun -J N0050-0.3-toe ../../../liquid-vapor-monte-carlo --movies --N 50 --ff 0.3 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0050-0.3-toe > lv/N0050-0.3-toe.out &

srun -J N0005-0.1-tmi ../../../liquid-vapor-monte-carlo --movies --N 5 --ff 0.1 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0005-0.1-tmi > lv/N0005-0.1-tmi.out &
srun -J N0005-0.2-tmi ../../../liquid-vapor-monte-carlo --movies --N 5 --ff 0.2 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0005-0.2-tmi > lv/N0005-0.2-tmi.out &
srun -J N0005-0.3-tmi ../../../liquid-vapor-monte-carlo --movies --N 5 --ff 0.3 --tmi --min-samples 100 --min-T $MINT --dir lv --filename N0005-0.3-tmi > lv/N0005-0.3-tmi.out &

srun -J N0005-0.1-toe ../../../liquid-vapor-monte-carlo --movies --N 5 --ff 0.1 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0005-0.1-toe > lv/N0005-0.1-toe.out &
srun -J N0005-0.2-toe ../../../liquid-vapor-monte-carlo --movies --N 5 --ff 0.2 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0005-0.2-toe > lv/N0005-0.2-toe.out &
srun -J N0005-0.3-toe ../../../liquid-vapor-monte-carlo --movies --N 5 --ff 0.3 --toe --min-samples 100 --min-T $MINT --dir lv --filename N0005-0.3-toe > lv/N0005-0.3-toe.out &
