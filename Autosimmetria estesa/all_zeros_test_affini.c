#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "autosymmetry.h"
#include "parser.h"

float massimo (float a, float b){
    if (a >= b){
        return a;
    }
    else{
        return b;
    }
}

int main(int argc, char** argv)
{
    init();

    boolean_function_t* f = parse_pla(manager, argv[1], 1);
    DdNode* lf, *fk;

    double potenzaDue = pow(2, f->inputs);
    //fprintf(stdout, "Potenza di due: %f\n", potenzaDue);
    fprintf(stdout, "f->inputs: %d\n", f->inputs);

    int intpotenzaDue = (int) potenzaDue;
    //fprintf(stdout, "Potenza di due int: %d\n", intpotenzaDue);

    //int matrice[intpotenzaDue][f->inputs];

    int *matrice[intpotenzaDue];
    //int *matrice = malloc(sizeof(int)*intpotenzaDue);
    for (int i = 0; i < intpotenzaDue; i++){
        int *row = malloc(sizeof(int)*f->inputs);
        matrice[i] = row;
    }    

    int temp, iteratore;

    for(int j=0; j<potenzaDue; j++){
        int variabili[f->inputs];
        temp = j;
        for(int k=0; k<f->inputs; k++){
            variabili[f->inputs-k-1] = temp%2;
            temp = temp/2;
        }
        for(int i=0; i<f->inputs; i++){
            matrice[j][i]=variabili[i];
        }   
    }

    double elements=0, degree=0, avg_degree=0, total=0, max_degree=0;
    clock_t start0, end0, start1, end1;
    int i, ns=0;

    //array per i nomi dei file contenenti le equazioni di riduzione
    //e le pla delle forme ridotte.
    char** eq_file_names = calloc(f->outputs, sizeof(char*));
    char** red_file_names = calloc(f->outputs, sizeof(char*));

    char buf[100000];
    system(buf);
    system("rm -r eq");
    system("rm -r fk");
    system("rm -r eqN");
    system("rm -r fkN");
    system("mkdir eq");
    system("mkdir fk");
    system("mkdir eqN");
    system("mkdir fkN");
    sprintf(buf, "for i in {0..%d}; do mkdir fk/fk$i; done", intpotenzaDue-1);
    system(buf);
    sprintf(buf, "for i in {0..%d}; do mkdir fkN/fk$i; done", intpotenzaDue-1);
    system(buf);
    sprintf(buf, "for i in {0..%d}; do mkdir eq/eq$i; done", intpotenzaDue-1);
    system(buf);
    sprintf(buf, "for i in {0..%d}; do mkdir eqN/eq$i; done", intpotenzaDue-1);
    system(buf);

/*     int** cv = calloc(f->outputs, sizeof(int*));
    if(cv==NULL)
    {
        perror("Errore chiamata calloc: ");
        return -1;
    }
    //array per le variabili canoniche
    for(i=0; i<f->outputs; i++){
        cv[i] = calloc(f->inputs, sizeof(int));
    } */

/*     for(i=0; i<f->outputs; i++){
        eq_file_names[i] = calloc(10, sizeof(char));
        red_file_names[i] = calloc(10, sizeof(char));
        sprintf(eq_file_names[i], "eq/eq%i.re",i);
        sprintf(red_file_names[i], "fk/fk%i.pla",i);
    } */



/*     for(int i=0; i<potenzaDue; i++){
        for(int j=0; j<f->inputs; j++){
            fprintf(stdout, "| %d |", matrice[i][j]);
        }
        fprintf(stdout, "\n");
    } */

    //Preparazione variabili
    DdNode* variabile[f->inputs];
    for(int k=0; k<f->inputs; k++){
        variabile[k] = Cudd_bddIthVar(manager, (2*k)+1);
        Cudd_Ref(variabile[k]);
    }

    for(int i=0; i<f->outputs; i++){

        start0 = clock();

        float max = 0;
        DdNode *varXor;
        fprintf(stdout, "\n-----------------------------");
        fprintf(stdout, "\n-----------------------------\n");
        //fprintf(stdout, "f->on_set[%i]: %d\n", i, Cudd_PrintMinterm(manager, f->on_set[i])); 
        fprintf(stdout, "i: %d, NumMintermini: %f\n", i, Cudd_CountMinterm(manager, f->on_set[i], f->inputs));

        for(int j=0; j<potenzaDue; j++){    
            for(int k=0; k<f->outputs; k++){
                eq_file_names[k] = calloc(17, sizeof(char));
                red_file_names[k] = calloc(17, sizeof(char));
                sprintf(eq_file_names[k], "eq/eq%i/eq%i.re",j,k);
                sprintf(red_file_names[k], "fk/fk%i/fk%i.pla",j,k);
            }

            int** cv = calloc(f->outputs, sizeof(int*));
            if(cv==NULL)
            {
                perror("Errore chiamata calloc: ");
                return -1;
            }
            //array per le variabili canoniche
            for(int k=0; k<f->outputs; k++){
                cv[k] = calloc(f->inputs, sizeof(int));
            }

            varXor = f->on_set[i];
            Cudd_Ref(f->on_set[i]);
            fprintf(stdout, "counter: %d\n", j); 
            for(int k=0; k<f->inputs; k++){
                if(matrice[j][k]==1){
                    varXor = Cudd_bddXor(manager, varXor, variabile[k]);
                    Cudd_Ref(varXor);
                    fprintf(stdout, "xor con x%d\n", k+1);
                }
            }

            //fprintf(stdout, "Mintermini varXor: %d\n", Cudd_PrintMinterm(manager, varXor)); 
            fprintf(stdout, "Numero mintermini varXor: %f\n", Cudd_CountMinterm(manager, varXor, f->inputs)); 
            lf = buildLf(varXor, f->inputs);
            fprintf(stdout, "Mintermini LF: %d\n", Cudd_PrintMinterm(manager, lf));
            fprintf(stdout, "Numero mintermini LF: %f\n", Cudd_CountMinterm(manager, lf, f->inputs));
            //fprintf(stdout, "\n");
            max = massimo(max, Cudd_CountMinterm(manager, lf, f->inputs));
            elements = Cudd_CountMinterm(manager, lf, f->inputs);
            //fprintf(stdout, "%d\n", Cudd_PrintMinterm(manager, lf));
            degree = log(elements)/log(2);

            end0 = clock();

            //fprintf(stdout, "Numero di elementi in Lf: %.2f\n", elements);
            //fprintf(stdout, "Grado di autosimmetria di f_out_%i: %.2f\n", i, degree);

            if(degree>=1 && degree<=f->inputs)
            {
                start1 = clock();

                cv[i] = computeCanonicalVariables(lf, 0, f->inputs, cv[i]);

                fk = restrictionFunction(varXor, cv[i], f->inputs);
                //fprintf(stdout, "Mintermini fk: %d\n", Cudd_PrintMinterm(manager, fk));
                //fprintf(stdout, "Cudd_SupportSize(manager, fk): %d\n", Cudd_SupportSize(manager, fk));
                

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

                FILE *eq_file = fopen(eq_file_names[i], "w");
                print_equations(eq_file, eq);
                fclose(eq_file);
                EQ_manager_quit(eq_manager);

                printPla(manager, red_file_names[i], fk, 2*f->inputs);

                Cudd_RecursiveDeref(manager, lf);
                Cudd_RecursiveDeref(manager, fk);

                if(degree>max_degree) max_degree = degree;

                avg_degree += degree;

                total += (double)(end1-start1);
            }
            else
            {
                Cudd_RecursiveDeref(manager, lf);
            }

            ns = Cudd_SupportSize(manager, f->on_set[i]);

            //printf("%s(%i) %i %i %i %f\n", argv[1], i, f->inputs, ns, (int)degree, (double)(f->dcs/f->outputs));

            total += (double)(end0-start0);
            Cudd_RecursiveDeref(manager, varXor);  

            for(int k=0; k<f->outputs; k++){
                free(cv[k]);
            }
            free(cv);   
            fprintf(stdout, "\n-----------------------------\n");   
        }

/*         fprintf(stdout, "\nCudd_not (xor 1)\n");

        for(int j=0; j<potenzaDue; j++){    
            for(int k=0; k<f->outputs; k++){
                eq_file_names[k] = calloc(18, sizeof(char));
                red_file_names[k] = calloc(18, sizeof(char));
                sprintf(eq_file_names[k], "eqN/eq%i/eq%i.re",j,k);
                sprintf(red_file_names[k], "fkN/fk%i/fk%i.pla",j,k);
            }

            int** cv = calloc(f->outputs, sizeof(int*));
            if(cv==NULL)
            {
                perror("Errore chiamata calloc: ");
                return -1;
            }
            //array per le variabili canoniche
            for(int k=0; k<f->outputs; k++){
                cv[k] = calloc(f->inputs, sizeof(int));
            }

            varXor = Cudd_Not(f->on_set[i]);
            Cudd_Ref(f->on_set[i]);
            fprintf(stdout, "counter: %d\n", j); 
            for(int k=0; k<f->inputs; k++){
                if(matrice[j][k]==1){
                    varXor = Cudd_bddXor(manager, varXor, variabile[k]);
                    Cudd_Ref(varXor);
                    fprintf(stdout, "f xor 1 + xor con x%d\n", k+1);
                }
            }
            //fprintf(stdout, "varXor: %d | NumMintermini: %f\n", Cudd_PrintMinterm(manager, varXor), Cudd_CountMinterm(manager, varXor, f->inputs));
            fprintf(stdout, "Mintermini varXor: %f\n", Cudd_CountMinterm(manager, varXor, f->inputs)); 
            lf = buildLf(varXor, f->inputs);
            fprintf(stdout, "Mintermini LF: %d\n", Cudd_PrintMinterm(manager, lf));
            fprintf(stdout, "Numero mintermini LF: %f\n", Cudd_CountMinterm(manager, lf, f->inputs));
            //fprintf(stdout, "\n");
            max = massimo(max, Cudd_CountMinterm(manager, lf, f->inputs));
            elements = Cudd_CountMinterm(manager, lf, f->inputs);
            //fprintf(stdout, "%d\n", Cudd_PrintMinterm(manager, lf));
            degree = log(elements)/log(2);

            end0 = clock();

            //fprintf(stdout, "Numero di elementi in Lf: %.2f\n", elements);
            //fprintf(stdout, "Grado di autosimmetria di f_out_%i: %.2f\n", i, degree);

            if(degree>=1 && degree<=f->inputs)
            {
                start1 = clock();

                cv[i] = computeCanonicalVariables(lf, 0, f->inputs, cv[i]);

                fk = restrictionFunction(varXor, cv[i], f->inputs);
                //fprintf(stdout, "Mintermini fk: %d\n", Cudd_PrintMinterm(manager, fk));
                //fprintf(stdout, "Cudd_SupportSize(manager, fk): %d\n", Cudd_SupportSize(manager, fk));
                

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

                FILE *eq_file = fopen(eq_file_names[i], "w");
                print_equations(eq_file, eq);
                fclose(eq_file);
                EQ_manager_quit(eq_manager);

                printPla(manager, red_file_names[i], fk, 2*f->inputs);

                Cudd_RecursiveDeref(manager, lf);
                Cudd_RecursiveDeref(manager, fk);

                if(degree>max_degree) max_degree = degree;

                avg_degree += degree;

                total += (double)(end1-start1);
            }
            else
            {
                Cudd_RecursiveDeref(manager, lf);
            }

            ns = Cudd_SupportSize(manager, f->on_set[i]);

            //printf("%s(%i) %i %i %i %f\n", argv[1], i, f->inputs, ns, (int)degree, (double)(f->dcs/f->outputs));

            total += (double)(end0-start0);
            Cudd_RecursiveDeref(manager, varXor);  

            for(int k=0; k<f->outputs; k++){
                free(cv[k]);
            }
            free(cv); 
            fprintf(stdout, "\n-----------------------------\n");     
        }
 */        fprintf(stdout, "Massimo: %f\n", max);
    }

    printf("%s %i %i %f %.2f %.2f %f\n", argv[1], f->inputs, f->outputs, (double)(f->dcs/f->outputs), (double)(avg_degree/f->outputs), max_degree, (total/CLOCKS_PER_SEC));

    for(i=0; i<f->outputs; i++)
    {
        free(eq_file_names[i]);
        free(red_file_names[i]);
        //free(cv[i]);
    }

    free(eq_file_names);
    free(red_file_names);
    free(f->dc_set);
    free(f->on_set);
    //free(cv);
    free(f);   

    quit();

    return 0;
}
