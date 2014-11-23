# makefile for optimization framework
# written by Ruimin Wang
# powered by 'MathU'
# copyright@MathU

DIR_INC += -I/usr/local/Cellar/Eigen/3.2.1/include/eigen3 -I./include/Base -I./include/LeastSquares -I./include/FunctionRepository -I./include/Algorithms -I./include/ThirdParty -I./include
DIR_SRC += ./src/Frame ./src/LeastSquares ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin
# DIR_LIB = -LE:/libs

SRC = $(foreach n, $(DIR_SRC),$(wildcard $(n)/*.cpp))
# TT = $(patsubst %.cpp,%.o,$(SRC))
# TT = $(notdir $(SRC))
OBJ = $(patsubst %.cpp,$(DIR_OBJ)/%.o,$(notdir $(SRC)))

TEST_BIN = $(DIR_BIN)/test
TEST_OBJ = $(DIR_OBJ)/test.o
TEST_SRC = ./test/simplex.cpp

TARGET = all
BIN_TARGET = $(DIR_BIN)/$(TARGET)

CC = g++
CFLAGS = -g -Wall $(DIR_INC)

vpath %.cpp ./src/Frame
vpath %.cpp ./src/LeastSquares
vpath %.cpp ./src

#objects = test.o basicmath.o

# INCLUDE += -I/usr/local/Cellar/Eigen/3.2.1/include/eigen3 -Iinclude/Frame/

# vpath %.cpp include/Frame

all:$(BIN_TARGET)


$(BIN_TARGET):$(OBJ)
	$(CC) $(OBJ) $(DIR_LIB) -o $@ -lcblas -lblas

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

test: $(TEST_BIN)

matrix: bin/matrix

# test for matrix
bin/matrix: obj/matrix.o
	$(CC) obj/matrix.o $(DIR_LIB) -o $@ -lcblas -lblas
obj/matrix.o: test/matrix.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# test for simplex method
simplex: bin/simplex
bin/simplex: obj/simplex.o
	$(CC) obj/simplex.o $(DIR_LIB) -o $@ -lcblas -lblas
obj/simplex.o: test/simplex.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# test for matrix vector
vecmat: bin/vecmat
bin/vecmat: obj/vecmat.o
	$(CC) obj/vecmat.o $(DIR_LIB) -o $@ -lcblas -lblas
obj/vecmat.o: test/matrix_vector.cpp
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