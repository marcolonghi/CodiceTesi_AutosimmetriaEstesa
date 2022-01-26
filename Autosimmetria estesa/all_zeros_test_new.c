#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "autosymmetry.h"
#include "parser.h"

#define DD_APA_BITS	((int) sizeof(DdApaDigit) * 8)
#define DD_APA_BASE	((DdApaDoubleDigit) 1 << DD_APA_BITS)

typedef uint64_t DdApaDoubleDigit;

int divisibilePerDue(char* vettoreChar){
    int posizione = strlen(vettoreChar)-1;
    if(vettoreChar[posizione]=='0' || vettoreChar[posizione]=='2' || vettoreChar[posizione]=='4' || vettoreChar[posizione]=='6' || vettoreChar[posizione]=='8'){
        return 1;
    }
    else return 0;
}

int main(int argc, char** argv)
{
    init();

    boolean_function_t* f = parse_pla(manager, argv[1], 1);
    DdNode* lf, *lf1;
    DdNode* fk;

    int** cv = calloc(f->outputs, sizeof(int*));
    double elements=0, degree=0, avg_degree=0, total=0, max_degree=0;
    clock_t start0, end0, start1, end1;
    int i, ns=0;

    char *str = malloc(strlen(argv[1]));

    strncpy(str, argv[1] + 10, strlen(argv[1]));

    char buf[100000];
    sprintf(buf, "mkdir eq/%s", str);
    system(buf);
    sprintf(buf, "mkdir fk/%s", str);
    system(buf);

    //Calcolo di 2^{n-1}, per il primo preprocessing
    //Se il numero elementi dell'onset è uguale a 2^{n-1} allora si calcola L_f come L0 \cup L1
    int digits2 = Cudd_ApaNumberOfDigits(f->inputs+1);
    DdApaNumber potenzaDueMenoUno = Cudd_NewApaNumber(digits2);
    Cudd_ApaPowerOfTwo(digits2, potenzaDueMenoUno, f->inputs-1);

    //Memorizzo il valore di 2^{n-1} in un array di char
    int dimensione = (int) (digits2 * log10((double) DD_APA_BASE)) + 1;
    char stringa[dimensione];
    strcpy(stringa, Cudd_ApaStringDecimal(digits2, potenzaDueMenoUno));

    //array per i nomi dei file contenenti le equazioni di riduzione
    //e le pla delle forme ridotte.
    char** eq_file_names = calloc(f->outputs, sizeof(char*));
    char** red_file_names = calloc(f->outputs, sizeof(char*));

    if(cv==NULL)
    {
        perror("Errore chiamata calloc: ");
        return -1;
    }

    eq_file_names[0] = calloc(50, sizeof(char));
    red_file_names[0] = calloc(50, sizeof(char));
    sprintf(eq_file_names[0], "eq/%s/eq0.re", str);
    sprintf(red_file_names[0], "fk/%s/fk0.pla", str);

    //array per le variabili canoniche
    cv[0] = calloc(f->inputs, sizeof(int));

    //Conto il numero di elementi dell'onset di f
    //Se il numero di elementi dell'onset è divisibile per 2 allora si esegue il test di autosimmetria classico (solo L0)
    int digits;
    DdApaNumber count;
    int compare;
    char countStringa[dimensione];

    start0 = clock();

    //Conto il numero di elementi dell'onset di f
    count = Cudd_ApaCountMinterm(manager,f->on_set[0],f->inputs,&digits);
    strcpy(countStringa, Cudd_ApaStringDecimal(digits, count));

    //Vedere se è uguale a potenzaDueMenoUno
    //Decidere se tenere Cudd_ApaCompare o fare confronto di stringhe
    compare = Cudd_ApaCompare(digits2, potenzaDueMenoUno, digits, count);

    free(count);

    if (compare==0){
        //Il numero elementi dell'onset è uguale a 2^{n-1}
        //Calcolo lf nuovo
        lf = buildLfNew(f->on_set[0], f->inputs);
    }
    //else if (countStringa[strlen(countStringa)-1]=='0' || countStringa[strlen(countStringa)-1]=='2' || countStringa[strlen(countStringa)-1]=='4' || countStringa[strlen(countStringa)-1]=='6' || countStringa[strlen(countStringa)-1]=='8'){
    else if (divisibilePerDue(countStringa)==1){
        //Il numero elementi dell'onset è divisibile per due
        //Calcolo lf autosimmetria
        lf = buildLf(f->on_set[0], f->inputs);
    }
    else{
        //Il calcolo è inutile, non si trovano regolarità
        lf = buildLfZero();
    }

    elements = Cudd_CountMinterm(manager, lf, f->inputs);
    degree = log(elements)/log(2);

    end0 = clock();

    if(degree>=1 && degree<=f->inputs)
    {
        start1 = clock();

        cv[0] = computeCanonicalVariables(lf, 0, f->inputs, cv[0]);

        fk = restrictionFunction(f->on_set[0], cv[0], f->inputs);

        if(Cudd_SupportSize(manager, fk)!=(f->inputs-(int)degree))
        {
            fprintf(stderr, "Costruzione di fk fallita\n");
            return -1;
        }

        EQ_manager* eq_manager = EQ_manager_init();
        Eqns_t* eq = new_eqns(eq_manager);

        eq = reductionEquations(lf, 0, f->inputs, eq_manager);
        complements_all(eq);

        end1 = clock();

        FILE *eq_file = fopen(eq_file_names[0], "w");
        print_equations(eq_file, eq);
        fclose(eq_file);
        EQ_manager_quit(eq_manager);

        printPla(manager, red_file_names[0], fk, 2*f->inputs);

        Cudd_RecursiveDeref(manager, f->on_set[0]);
        Cudd_RecursiveDeref(manager, f->dc_set[0]);
        Cudd_RecursiveDeref(manager, lf);
        Cudd_RecursiveDeref(manager, fk);

        if(degree>max_degree) max_degree = degree;

        avg_degree += degree;

        total += (double)(end1-start1);
    }
    else
    {
        Cudd_RecursiveDeref(manager, f->on_set[0]);
        Cudd_RecursiveDeref(manager, f->dc_set[0]);
        Cudd_RecursiveDeref(manager, lf);
    }

    ns = Cudd_SupportSize(manager, f->on_set[0]);

    total += (double)(end0-start0);

    printf("%s %i %i %f %.2f %.2f %f\n", argv[1], f->inputs, f->outputs, (double)(f->dcs/f->outputs), (double)(avg_degree/f->outputs), max_degree, (total/CLOCKS_PER_SEC));

    free(eq_file_names[0]);
    free(red_file_names[0]);
    free(cv[0]);

    free(potenzaDueMenoUno);
    free(eq_file_names);
    free(red_file_names);
    free(f->dc_set);
    free(f->on_set);
    free(cv);
    free(f);   

    quit();

    return 0;
}
