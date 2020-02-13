all: detector

%: %.cc
	g++ -Wall -Wextra -o $@ $< `geant4-config --cflags --libs`
