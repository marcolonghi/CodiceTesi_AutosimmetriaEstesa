#include <stdio.h>
#include <stdlib.h>
#include <cudd.h>
#include <math.h>
#include <time.h>

#include "string.h"
#include "autosymmetry.h"
#include "parser.h"

#define DD_APA_BITS	((int) sizeof(DdApaDigit) * 8)
#define DD_APA_BASE	((DdApaDoubleDigit) 1 << DD_APA_BITS)

typedef uint64_t DdApaDoubleDigit;

int main(int argc, char ** argv) 
{
    init();

    boolean_function_t* f = parse_pla(manager, argv[1], 1);
    int i;

    char *str = malloc(strlen(argv[1]));

    strncpy(str, argv[1] + 10, strlen(argv[1]));

    char buf[100000];
    sprintf(buf, "mkdir eq/%s", str);
    system(buf);
    sprintf(buf, "mkdir fk/%s", str);
    system(buf);

    //fprintf(stdout, "numero variabili: %d\n", f->inputs);

    DdNode* lf,*lf2;
    DdNode* s;
    DdNode* ls;
    DdNode* x;
    DdNode* tr;
    DdNode* new_f;
    DdNode* fk;
    int bit; //Indica se si esegue il test generalizzato o meno

    DdNode** vars = calloc(2*f->inputs, sizeof(DdNode*));

    double elements=0, degree=0, avg_degree=0;

    if(vars==NULL)
    {
        perror("Errore chiamata calloc: ");
        return -1;
    }

    for(int i=0; i<f->inputs; i++)
    {
        DdNode* tmp = Cudd_bddIthVar(manager, i);
        vars[i] = Cudd_bddIthVar(manager, (2*f->inputs)-1-i);
        vars[(2*f->inputs)-1-i] = tmp;
    }

    int** cv = calloc(f->outputs, sizeof(int*));

    if(cv==NULL)
    {
        perror("Errore chiamata calloc: ");
        free(vars);
        return -1;
    }

    //array per i nomi dei file contenenti le equazioni di riduzione e le pla delle forme ridotte.
    char** eq_file_names = calloc(f->outputs, sizeof(char*));
    char** red_file_names = calloc(f->outputs, sizeof(char*));

    eq_file_names[0] = calloc(50, sizeof(char));
    red_file_names[0] = calloc(50, sizeof(char));
    sprintf(eq_file_names[0], "eq/%s/eq0.re", str);
    sprintf(red_file_names[0], "fk/%s/fk0.pla", str);

    //array per le variabili canoniche
    cv[0] = calloc(f->inputs, sizeof(int));  

    clock_t start0, end0, start1, end1;
    double total=0, max_degree=0;
    int max_ns=0, min_ns = INT32_MAX, ns=0, out_num=0;

    //Calcolo di 2^{n-1}
    int digits2 = Cudd_ApaNumberOfDigits(f->inputs+1);
    DdApaNumber potenzaDueMenoUno = Cudd_NewApaNumber(digits2);
    Cudd_ApaPowerOfTwo(digits2, potenzaDueMenoUno, f->inputs-1);

    //per ciascun output della funzione trattata effettua il test di autosimmetria ed eventualmente il calcolo della forma ridotta e delle equazioni di riduzione

    start0 = clock();

    bit = 0;

    //Conto il numero di elementi dell'onset di f
    int digitsOn, digitsDC, digitsSum, digitsTest;
    DdApaNumber countOn, countDC, countSum, countTest;
    int compare;
    countOn = Cudd_ApaCountMinterm(manager,f->on_set[0],f->inputs,&digitsOn);
    countDC = Cudd_ApaCountMinterm(manager,f->dc_set[0],f->inputs,&digitsDC);

    DdNode* sum; 
    sum = Cudd_bddOr(manager, f->on_set[0], f->dc_set[0]);
    Cudd_Ref(sum);
    countSum = Cudd_ApaCountMinterm(manager,sum,f->inputs,&digitsSum);

/*         int dimensione1 = (int) (digitsOn * log10((double) DD_APA_BASE)) + 1;
    int dimensione2 = (int) (digitsDC * log10((double) DD_APA_BASE)) + 1;
    int dimensione3 = (int) (digitsSum * log10((double) DD_APA_BASE)) + 1;

    char stringaOn[dimensione1];
    char stringaDC[dimensione2];
    char stringaSum[dimensione3];

    strcpy(stringaOn, Cudd_ApaStringDecimal(digitsOn, countOn));
    strcpy(stringaDC, Cudd_ApaStringDecimal(digitsDC, countDC));
    strcpy(stringaSum, Cudd_ApaStringDecimal(digitsSum, countSum)); */

    //fprintf(stdout, "i: %d | potenza2meno1: %s, onset: %s, dcset: %s, sum: %s\n", i, stringa, stringaOn, stringaDC, stringaSum);

    compare = Cudd_ApaCompare(digits2, potenzaDueMenoUno, digitsSum, countSum);

    if(compare == 1){
        //potenzaDueMenoUno > countSum | Classico test generalizzato    
        s = buildS(sum, f->on_set[0], f->inputs);
        ls = build_Ls_3(s, f->inputs, f->on_set[0], f->dc_set[0]);
        //calcola la dimensione di ls
        elements = Cudd_CountMinterm(manager, ls, f->inputs);
        degree = log(elements)/log(2);
        new_f = buildNewF(f->on_set[0], ls, f->inputs);
        bit = 1;
    }
    else if(compare == 0){
        //potenzaDueMenoUno = countSum 
        //Test di nuova autosimmetria su sum
        lf = buildLfNew(sum, f->inputs);
        //Test generalizzato su f->on_set[0]
        s = buildS(sum, f->on_set[0], f->inputs);
        ls = build_Ls_3(s, f->inputs, f->on_set[0], f->dc_set[0]);
        //Si tiene quello che ha cardinalità maggiore
        if (Cudd_CountMinterm(manager, lf, f->inputs) >= Cudd_CountMinterm(manager, ls, f->inputs)){
            elements = Cudd_CountMinterm(manager, lf, f->inputs);
            degree = log(elements)/log(2);
            bit = 0;
        }
        else{
            elements = Cudd_CountMinterm(manager, ls, f->inputs);
            degree = log(elements)/log(2);
            new_f = buildNewF(f->on_set[0], ls, f->inputs);
            bit = 1;
        }
        //fprintf(stdout, "potenzaDueMenoUno = countSum | bit: %d\n", bit);
    }
    else if(compare == -1){
        //potenzaDueMenoUno < countSum | Classico test generalizzato
        s = buildS(sum, f->on_set[0], f->inputs);
        ls = build_Ls_3(s, f->inputs, f->on_set[0], f->dc_set[0]);
        new_f = buildNewF(f->on_set[0], ls, f->inputs);
        countTest = Cudd_ApaCountMinterm(manager,new_f,f->inputs,&digitsTest);
        if (Cudd_ApaCompare(digitsTest, countTest, digits2, potenzaDueMenoUno) == 0){
            lf = buildLfNew(new_f, f->inputs);
            //calcola la dimensione di ls
            elements = Cudd_CountMinterm(manager, lf, f->inputs);
            degree = log(elements)/log(2);
            bit = 0;
        }
        else{
            //calcola la dimensione di ls
            elements = Cudd_CountMinterm(manager, ls, f->inputs);
            degree = log(elements)/log(2);
            bit = 1;
        }
    }

    end0 = clock();

    //se la funzione è autosimmetrica di grado k>=1
    //calcola la forma ridotta e le equazioni di riduzione
    if(degree>=1 && degree<=f->inputs)
    {
        start1 = clock();

        if (bit == 1){
            cv[0] = computeCanonicalVariables(ls, 0, f->inputs, cv[0]);
            fk = restrictionFunction(new_f, cv[0], f->inputs);
        }
        else {
            cv[0] = computeCanonicalVariables(lf, 0, f->inputs, cv[0]);
            fk = restrictionFunction(f->on_set[0], cv[0], f->inputs);
        }

        if(Cudd_SupportSize(manager, fk)!=(f->inputs-(int)degree))
        {
            fprintf(stderr, "Costruzione di fk fallita\n");
            return -1;
        }

        EQ_manager* eq_manager = EQ_manager_init();
        Eqns_t* eq = new_eqns(eq_manager);

        if (bit == 1){
            eq = reductionEquations(ls, 0, f->inputs, eq_manager);
        }
        else {
            eq = reductionEquations(lf, 0, f->inputs, eq_manager);
        }
        
        complements_all(eq);

        end1 = clock();

        FILE *eq_file = fopen(eq_file_names[0], "w");
        print_equations(eq_file, eq);
        fclose(eq_file);
        EQ_manager_quit(eq_manager);

        printPla(manager, red_file_names[0], fk, 2*f->inputs);

        Cudd_RecursiveDeref(manager, f->on_set[0]);
        Cudd_RecursiveDeref(manager, f->dc_set[0]);
        if (bit == 0){
            Cudd_RecursiveDeref(manager, lf);
        }
        if (bit == 1){
            Cudd_RecursiveDeref(manager, s);
            Cudd_RecursiveDeref(manager, ls);
            Cudd_RecursiveDeref(manager, new_f);
        }
        Cudd_RecursiveDeref(manager, fk);

        if(degree>max_degree) max_degree = degree;

        avg_degree += degree;
        total += (double)(end1-start1);
    }
    else
    {
        Cudd_RecursiveDeref(manager, f->on_set[0]);
        Cudd_RecursiveDeref(manager, f->dc_set[0]);
        if (bit == 0){
            Cudd_RecursiveDeref(manager, lf);
        }
        if (bit == 1){
            Cudd_RecursiveDeref(manager, s);
            Cudd_RecursiveDeref(manager, ls);
        }
    }

    ns = Cudd_SupportSize(manager, f->on_set[0]);

    total += (double)(end0-start0);

    printf("%s %i %i %f %.2f %.2f %f\n", argv[1], f->inputs, f->outputs, (double)(f->dcs/f->outputs), (double)(avg_degree/f->outputs), max_degree, (total/CLOCKS_PER_SEC));

    free(eq_file_names[0]);
    free(red_file_names[0]);
    free(cv[0]);

    free(eq_file_names);
    free(red_file_names);
    free(f->dc_set);
    free(f->on_set);
    free(vars);
    free(cv);
    free(f);

    quit();

    return 0;
}

