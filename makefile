SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj
EIGEN_DIR = /usr/include/eigen3

CC = g++
CPPFLAGS = -Iinclude
CFLAGS = -Wall
EIGEN = -I $(EIGEN_DIR)

LINK_TARGET = test, test_rng

SRC = $(wildcard $(SRC_DIR)/*.cpp)
OBJ_ALL = $(patsubst %.cpp, %.o, $(addprefix $(OBJ_DIR)/, $(notdir $(SRC))))
OBJ_TAR = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(LINK_TARGET)))
OBJ = $(filter-out $(OBJ_TAR), $(OBJ_ALL))

REBUILDABLES = $(OBJ_ALL) $(LINK_TARGET) 


all : $(LINK_TARGET)

clean: 
	rm -f $(REBUILDABLES)
	rm -f log.txt


$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cpp
	$(CC) -g  $(CPPFLAGS) $(CFLAGS) $(EIGEN) -o $@ -c $< 

test : $(OBJ_DIR)/test.o 
	$(CC) -g $(CPPFLAGS) $(CFLAGS) $(EIGEN) -o test $(OBJ_DIR)/test.o

test_rng : $(OBJ_DIR)/test_rng.o $(OBJ_DIR)/random.o
	$(CC) -g $(CPPFLAGS) $(CFLAGS) $(EIGEN) -o $@ $(OBJ_DIR)/test_rng.o $(OBJ_DIR)/random.o