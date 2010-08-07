# Build static html docs suitable for being shipped in the software
# package. This depends on ikiwiki being installed to build the docs.

ifeq ($(shell which ikiwiki),)
IKIWIKI=echo "** ikiwiki not found" >&2 ; echo ikiwiki
else
IKIWIKI=ikiwiki
endif

all:
	$(IKIWIKI) `pwd` html -v --wikiname FooBar --plugin=goodstuff \
		--exclude=html --exclude=Makefile

clean:
	rm -rf .ikiwiki html
