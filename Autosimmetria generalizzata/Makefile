#!/bin/bash

target = all_zeros_test.out 					 \
		 all_ones_test.out  					 \
		 generalized_test_Filippo.out 			 \
		 generalized_test_Filippo_migliorato.out \
		 test_reduction.out		

OBJ  = parser.o autosymmetry.o equations.o

ccflags = -std=c99 #-Wall

.PHONY: all clean cleanall runall

pla/%.pla: all
	./run.sh $@

%.out: %.c $(OBJ)
	gcc  $(ccflags) $< $(OBJ) -o $@ -lcudd -lm
	
%.o : %.c %.h
	gcc $(ccflags) -c -o $@ $<

all: $(target)

clean:
	rm -r -f S.pla output.*.txt eq/* fk/* *.err *.log
	
cleanall: clean
	rm -f *.o *.out