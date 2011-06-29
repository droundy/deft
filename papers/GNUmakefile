PAPERS=electrostatics contact

all:
	for paper in $(PAPERS); do cd $$paper && $(MAKE); cd ..; done

clean:	
	for paper in $(PAPERS); do cd $$paper && $(MAKE) clean; cd ..; done
