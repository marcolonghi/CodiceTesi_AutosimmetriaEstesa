#!/bin/bash

target = all_zeros_test.out 					 \
		 all_ones_test.out  					 \
		 all_ones_test_new.out 			 		 \
		 all_zeros_test_new.out					 \
		 generalized_new.out                     \
		 all_zeros_test_affini.out				 \
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