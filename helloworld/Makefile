# change the following values to reflect your Ariadne installation
INCLUDE_DIR = /usr/local/lib
LIB_DIR = /usr/local/lib

helloworld: helloworld.o
	g++ -o helloworld helloworld.o -L$(LIB_DIR) -lariadne

helloworld.o: helloworld.cc
	g++ -c helloworld.cc -frounding-math -I$(INCLUDE_DIR)

clean:
	rm -f helloworld helloworld.o
