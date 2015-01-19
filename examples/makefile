default: all

include ../Makefile.in

I = -I../include

all: array vector

array: bin/array
	
bin/array: obj/array.o
	$(CXX) obj/array.o -L$(DIR_LIB) -o $@ $(LIB)
obj/array.o: array.cpp ../include/Base/Array.hpp ../include/Base/ArrayImpl.hpp
	$(CXX) $(CXXFLAGS) $(I) -c $< -o $@

vector: bin/vector
	
bin/vector: obj/vector.o
	$(CXX) obj/vector.o -L$(DIR_LIB) -o $@ $(LIB)
obj/vector.o: vector.cpp 
	$(CXX) $(CXXFLAGS) $(I) -c $< -o $@
