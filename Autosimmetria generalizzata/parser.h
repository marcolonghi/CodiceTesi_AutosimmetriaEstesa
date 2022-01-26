#ifndef _PLA_PARSER_
#define _PLA_PARSER_

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <cudd.h>

#define MAX_LINE_LENGTH 1024

typedef struct boolean_function_t {
    int inputs;
    int outputs;
    int dcs;
    DdNode** on_set;
    DdNode** dc_set;
} boolean_function_t;

boolean_function_t* parse_pla(DdManager* manager, char* file_path, int alfa);

int** parse_equations(char* eq_path, int* rows, int* columns);

void printPla(DdManager* manager, char* outputfile, DdNode* bdd, int n_var);

#endif
