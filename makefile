# makefile for optimization framework
# written by Ruimin Wang
# powered by 'MathU'
# copyright@MathU

DIR_INC += -I./libs/Eigen -I./include/Base -I./include/LeastSquares -I./include/FunctionRepository -I./include/Algorithms -I./include/ThirdParty -I./include
DIR_SRC += ./src/Frame ./src/LeastSquares ./src
DIR_OBJ = ./obj
DIR_BIN = ./bin

SRC = $(foreach n, $(DIR_SRC),$(wildcard $(n)/*.cpp))
# TT = $(patsubst %.cpp,%.o,$(SRC))
# TT = $(notdir $(SRC))
OBJ = $(patsubst %.cpp,$(DIR_OBJ)/%.o,$(notdir $(SRC)))

TARGET = main
BIN_TARGET = $(DIR_BIN)/$(TARGET)

CC = g++
CFLAGS = -g -Wall $(DIR_INC)

vpath %.cpp ./src/Frame
vpath %.cpp ./src/LeastSquares
vpath %.cpp ./src

#objects = test.o basicmath.o

# INCLUDE += -I/usr/local/Cellar/Eigen/3.2.1/include/eigen3 -Iinclude/Frame/

# vpath %.cpp include/Frame


$(BIN_TARGET):$(OBJ)
	$(CC) $(OBJ) -o $@ -lcblas

$(DIR_OBJ)/%.o:%.cpp
	$(CC) $(CFLAGS) -c $< -o $@

.PHONY:clean
clean:
	find $(DIR_OBJ) -name *.o -exec rm -rf {} \;


# all debug
all:
	@echo $(BIN_TARGET)
	@echo $(SRC)
	@echo $(notdir $(SRC))
	@echo $(OBJ)
	@echo $(TT)
	@echo "end"


#test : $(objects)
#	$(CC) -o test $(objects)

#test.o : include/Frame/test.cpp
#	$(CC) -c $(INCLUDE) include/Frame/test.cpp

#basicmath.o : include/Frame/basicmath.cpp include/Frame/basicmath.h
#	$(CC) -c  $(INCLUDE) include/Frame/basicmath.cpp

#.PHONY : clean
#clean :
#	-rm test $(objects)