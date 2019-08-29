CXXFLAGS = -g -O3 -I /usr/include/eigen3/
CXX = g++

carve: seamcarving.cpp
	$(CXX) $(CXXFLAGS) seamcarving.cpp -o carve

clean:
	rm -f *~
	rm -f *.o

run:
	./carve dog.jpg dog_2.jpg 500 420
