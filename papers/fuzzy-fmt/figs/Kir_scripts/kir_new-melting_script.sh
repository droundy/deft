#!/bin/sh
#../new-melting.mkdat .325 0 1.1 2
#Best so far gwidth=.595 gives Diff=-28.4653
rm newmeltdataout.dat
#range of gwidth that work 1 .9 .8 .7 .6 .599 .598 .597 .596 .595 .594 .593 .592 .591  .59 .5 .4 
#for gwidth in .6 .599 .598 .597 .596 .595 .594 .593 .592 .591  .59 
#do
#../new-melting.mkdat .325 0 $gwidth 2
#done

#range of gwidth that work
#for gwidth in .72 .7 .68 .65 .64 .635 .63 .625 .62 .615 .61 .605 .6 .58 .56
#do
#../new-melting.mkdat .325 .1 $gwidth 2
#done

#range of gwidth that work
#for gwidth in 2 1.9 1.8 1.7 1.6 1.5 1.4 1.3 1.2 1.0 .9 .8 .75 .72 .7 .68 .65 .64 .635 .63 .625 .62 .615 .61 .605 .6
#for gwidth in 6 4 2  1.85 1.848 1.846 1.844 1.842 1.84 1.83  1 .5
#do
#../new-melting.mkdat .325 .1 $gwidth 2
#done

#range of gwidth that work
#for gwidth=.4576 gives N_crystal of 3.60013 for reduced number of spheres=3.6  DIFF=-38.5!
for gwidth in   .46 .459 .458 .4578 .4576 .4574 .4572 .457 .456   
do
../new-melting.mkdat 1.3 .1 $gwidth 2
done

./newmeltdataout
gnuplot newmeltdataout.gnu 
#pause -1
# this is a comment
#echo "$(($gwidth*2))"
#for ((i=0; i<4; i++)); do ls; done 
#gwidth=1.2
#end=1.6
#for ((i=0; i<4; i++)); do gwidth=start+1; done 
