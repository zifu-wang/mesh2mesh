# -*- makefile -*-
.PHONY: all clean

COMPILER = 
FLAGS = 

MODULE =

all: ${MODULE}.so

${MODULE}.so: ${MODULE}.f90
	${COMPILER} ${FLAGS} -c $< -m ${MODULE}

clean:
	rm *.pyc *~
