#include <math.h>
#include <time.h>

#include "autosymmetry.h"
#include "parser.h"

int main(int argc, char ** argv) 
{
    init();

    boolean_function_t* f = parse_pla(manager, argv[1], 1);
    int i;

    DdNode* lf;
    DdNode* s;
    DdNode* u;
    DdNode* ls;
    DdNode* x;
    DdNode* tr;
    DdNode* new_f;
    DdNode* fk;

    char *str = malloc(strlen(argv[1]));

    strncpy(str, argv[1] + 12, strlen(argv[1]));

    char buf[100000];
    sprintf(buf, "mkdir eq/%s", str);
    system(buf);
    sprintf(buf, "mkdir fk/%s", str);
    system(buf);

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

    //array per i nomi dei file contenenti le equazioni di riduzione
    //e le pla delle forme ridotte.
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

    //per ciascun output della funzione trattata
    //effettua il test di autosimmetria ed eventualmente
    //il calcolo della forma ridotta e delle equazioni di riduzione

    start0 = clock();

    //spazio vettoriale lf
    //lf = buildLf(f->on_set[0], f->inputs);

    //unione tra on-set e dc-set
    u = Cudd_bddOr(manager, f->on_set[0], f->dc_set[0]);
    Cudd_Ref(u);

    //insieme S
    s = buildS(u, f->on_set[0], f->inputs);

    ls = build_Ls_3(s, f->inputs, f->on_set[0], f->dc_set[0]);

    //calcola la dimensione di ls
    elements = Cudd_CountMinterm(manager, ls, f->inputs);
    degree = log(elements)/log(2);

    end0 = clock();

    //se la funzione è autosimmetrica di grado k>=1
    //calcola la forma ridotta e le equazioni di riduzione
    if(degree>=1 && degree<=f->inputs)
    {
        //printf("%.3lf\n", degree);

        if(degree>max_degree) max_degree = degree;

        avg_degree += degree;
        total += (double)(end0-start0);

        new_f = buildNewF(f->on_set[0], ls, f->inputs);

        cv[0] = computeCanonicalVariables(ls, 0, f->inputs, cv[0]);

        fk = restrictionFunction(new_f, cv[0], f->inputs);

        if(Cudd_SupportSize(manager, fk)!=(f->inputs-(int)degree))
        {
            fprintf(stderr, "Costruzione di fk fallita\n");
            return -1;
        }

        EQ_manager* eq_manager = EQ_manager_init();
        Eqns_t* eq = new_eqns(eq_manager);

        eq = reductionEquations(ls, 0, f->inputs, eq_manager);
        complements_all(eq);

        end1 = clock();
        
        FILE *eq_file = fopen(eq_file_names[0], "w");
        print_equations(eq_file, eq);
        fclose(eq_file);
        EQ_manager_quit(eq_manager);

        printPla(manager, red_file_names[0], fk, 2*f->inputs);

        Cudd_RecursiveDeref(manager, f->on_set[0]);
        Cudd_RecursiveDeref(manager, f->dc_set[0]);
        Cudd_RecursiveDeref(manager, ls);
        Cudd_RecursiveDeref(manager, u);
        Cudd_RecursiveDeref(manager, fk);

        if(degree>max_degree) max_degree = degree;

        avg_degree += degree;

        total += (double)(end1-start1);
    }
    else
    {
        Cudd_RecursiveDeref(manager, f->on_set[0]);
        Cudd_RecursiveDeref(manager, f->dc_set[0]);
        Cudd_RecursiveDeref(manager, u);
        Cudd_RecursiveDeref(manager, ls);
    }

    ns = Cudd_SupportSize(manager, f->on_set[0]);

    total += (double)(end0-start0);

    printf("%s %i %i %.3f %.3f %f\n", argv[1], f->inputs, f->outputs, (double)(avg_degree/f->outputs), max_degree, (total/CLOCKS_PER_SEC));

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