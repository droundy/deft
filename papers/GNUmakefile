PAPERS=electrostatics contact thesis-hughes

all:
	for paper in $(PAPERS); do cd $$paper && $(MAKE) || exit 1; cd ..; done

clean:	
	for paper in $(PAPERS); do cd $$paper && $(MAKE) clean; cd ..; done
