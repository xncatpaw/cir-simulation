SRC_DIR = src
INC_DIR = include
OBJ_DIR = obj

CC = gcc
CPPFLAGS = -Iinclude
CFLAGS = -Wall

LINK_TARGET = 

SRC = $(wildcard $(SRC_DIR)/*.c)
OBJ_ALL = $(patsubst %.c, %.o, $(addprefix $(OBJ_DIR)/, $(notdir $(SRC))))
OBJ_TAR = $(addprefix $(OBJ_DIR)/, $(addsuffix .o, $(LINK_TARGET)))
OBJ = $(filter-out $(OBJ_TAR), $(OBJ_ALL))

REBUILDABLES = $(OBJ_ALL) $(LINK_TARGET) 


all : $(LINK_TARGET)

clean: 
	rm -f $(REBUILDABLES)
	rm -f log.txt

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.c
	$(CC) -g  $(CPPFLAGS) $(CFLAGS) -o $@ -c $< 