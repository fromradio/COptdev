# makefile for optimization framework
# written by Ruimin Wang
# powered by 'MathU'
# copyright@MathU

include ./Makefile.in




#objects = test.o basicmath.o

# INCLUDE += -I/usr/local/Cellar/Eigen/3.2.1/include/eigen3 -Iinclude/Frame/

# vpath %.cpp include/Frame

all:$(BIN_TARGET)


$(BIN_TARGET):$(OBJ)
	$(CC) $(OBJ) $(DIR_LIB) -o $@ $(LIBS)

$(DIR_OBJ)/%.o:%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY:clean
clean:
	find $(DIR_OBJ) -name *.o -exec rm -rf {} \;


# all debug
# all:
# 	@echo $(BIN_TARGET)
# 	@echo $(SRC)
# 	@echo $(notdir $(SRC))
# 	@echo $(OBJ)
# 	@echo $(TT)
# 	@echo "end"




# test for matrix
matrix: bin/matrix
bin/matrix: obj/matrix.o
	$(CC) obj/matrix.o $(DIR_LIB) -o $@ $(LIBS)
obj/matrix.o: test/matrix.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# test for simplex method
simplex: bin/simplex
bin/simplex: obj/simplex.o
	$(CC) obj/simplex.o $(DIR_LIB) -o $@ $(LIBS)
obj/simplex.o: test/simplex.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# test for matrix vector
vecmat: bin/vecmat
bin/vecmat: obj/vecmat.o
	$(CC) obj/vecmat.o $(DIR_LIB) -o $@ $(LIBS)
obj/vecmat.o: test/matrix_vector.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# test for sparse matrix
spmat: bin/spmat
bin/spmat: obj/spmat.o
	$(CC) obj/spmat.o $(DIR_LIB) -o $@ $(LIBS)
obj/spmat.o:test/spmat.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# test for umfpack
umfpack: bin/umfpack
bin/umfpack: obj/umfpack.o
	$(CC) obj/umfpack.o $(DIR_LIB) -o $@ $(LIBS)
obj/umfpack.o:test/umfpack.cpp
	$(CC) $(CFLAGS) -c $< -o $@

umfpackwrapper: bin/umfpackwrapper
bin/umfpackwrapper: obj/umfpackwrapper.o
	$(CC) obj/umfpackwrapper.o $(DIR_LIB) -o $@ $(LIBS)
obj/umfpackwrapper.o: test/umfpack_wrapper.cpp
	$(CC) $(CFLAGS) -c $< -o $@

omp: bin/omp
bin/omp: obj/omp.o
	$(CC) obj/omp.o $(DIR_LIB) -o $@ $(LIBS)
obj/omp.o: test/omp.cpp
	$(CC) $(CFLAGS) -c $< -o $@

ls: bin/ls
bin/ls: obj/ls.o
	$(CC) obj/ls.o $(DIR_LIB) -o $@ $(LIBS)
obj/ls.o: test/ls.cpp
	$(CC) $(CFLAGS) -c $< -o $@

eigen: bin/eigen
bin/eigen: obj/eigen.o
	$(CC) obj/eigen.o $(DIR_LIB) -o $@ $(LIBS)
obj/eigen.o: test/eigen_umfpack.cpp
	$(CC) $(CFLAGS) -c $< -o $@

lapack: bin/lapack
bin/lapack: obj/lapack.o
	$(CC) obj/lapack.o $(DIR_LIB) -o $@ $(LIBS)
obj/lapack.o: test/lapack_test.cpp
	$(CC) $(CFLAGS) -c $< -o $@

proximal: bin/proximal
bin/proximal: obj/proximal.o
	$(CC) obj/proximal.o $(DIR_LIB) -o $@ $(LIBS)
obj/proximal.o: test/proximal.cpp
	$(CC) $(CFLAGS) -c $< -o $@

mpitest: bin/mpitest
bin/mpitest: obj/mpitest.o
	mpic++ obj/mpitest.o $(DIR_LIB) -o $@ $(LIBS)
obj/mpitest.o: test/openmpi_test.cpp
	mpic++ $(CFLAGS) -c $< -o $@

lasso: bin/lasso
bin/lasso: obj/lasso.o
	$(CC) obj/lasso.o $(DIR_LIB) -o $@ $(LIBS)
obj/lasso.o: test/lasso.cpp include/Algorithms/Lasso.hpp
	$(CC) $(CFLAGS) -c $< -o $@

eigensolver: bin/eigensolver
bin/eigensolver: obj/eigensolver.o
	$(CC) obj/eigensolver.o $(DIR_LIB) -o $@ $(LIBS)
obj/eigensolver.o: test/eigensolver.cpp include/Algorithms/Eigenvalue.hpp
	$(CC) $(CFLAGS) -c $< -o $@

iotest: bin/iotest
bin/iotest: obj/iotest.o
	$(CC) obj/iotest.o $(DIR_LIB) -o $@ $(LIBS)
obj/iotest.o: test/io_test.cpp include/IOs/MtxFile.hpp
	$(CC) $(CFLAGS) -c $< -o $@
#$(TEST_BIN): $(TEST_OBJ)
#	$(CC) $(TEST_OBJ) $(DIR_LIB) -o $@ -lcblas -lblas

#$(TEST_OBJ):$(TEST_SRC)
#	$(CC) $(CFLAGS) -c $< -o $@

help: $(TEST_BIN) 
	
$(TEST_BIN): $(TEST_OBJ)
	$(CC) $(TEST_OBJ) $(DIR_LIB) -o $@ -lcblas -lblas
$(TEST_OBJ):$(TEST_SRC)
	$(CC) $(CFLAGS) -c $< -o $@



#test : $(objects)
#	$(CC) -o test $(objects)

#test.o : include/Frame/test.cpp
#	$(CC) -c $(INCLUDE) include/Frame/test.cpp

#basicmath.o : include/Frame/basicmath.cpp include/Frame/basicmath.h
#	$(CC) -c  $(INCLUDE) include/Frame/basicmath.cpp

#.PHONY : clean
#clean :
#	-rm test $(objects)