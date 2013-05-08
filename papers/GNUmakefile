PAPERS=electrostatics contact fuzzy-fmt thesis-hughes water-saft hughes-saft pair-correlation polyhedra

all:
	for paper in $(PAPERS); do cd $$paper && $(MAKE) || exit 1; cd ..; done

clean:	
	for paper in $(PAPERS); do cd $$paper && $(MAKE) clean; cd ..; done
