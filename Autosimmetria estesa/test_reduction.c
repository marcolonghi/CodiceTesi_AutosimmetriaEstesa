#include <stdio.h>
#include <stdlib.h>

#include "parser.h"

#define MAX_OUT 3

static void inline free_matrix(int** m, int rows)
{
    for(int i=0; i<rows; i++)
        free(m[i]);
    
    free(m);
    m=NULL;
}

int main(int argc, char** argv)
{
    if(argc!=4)
    {
        fprintf(stderr, "Utilizzo: %s -funzione -funzione ridotta -numero output\n", argv[0]);
        return -1;
    }

    DdManager* manager = Cudd_Init(0,0,CUDD_UNIQUE_SLOTS,CUDD_CACHE_SLOTS,0);

    int rows, columns;

    int** a = parse_equations(argv[3], &rows, &columns);

    if(a==NULL)
    {
        perror("Errore chiamata parse_equations: ");
        return -1;
    }

    //effettua il parsing della funzione ridotta senza aggiungere le variabili alfa
    boolean_function_t* fk = parse_pla(manager, argv[2], 0);

    if(fk==NULL)
    {
        fprintf(stderr, "Funzione ridotta non trovata\n");
        return -1;
    }

    DdNode* xor;
    DdNode* tmp1;
    DdNode* tmp2;
    DdNode* tmp3;
    DdNode* tmp4;
    int i,j;

    DdNode** vars = malloc((columns)*sizeof(DdNode*));
    
    //creo n nuove variabili, il numero di colonne
    //della matrice A è uguale al numero di variabili
    //della funzione originale
    for(i=0; i<columns; i++)
        vars[i] = Cudd_bddNewVar(manager);

    //costruisco l'input per la funzione ridotta
    //prendendo in considerazione le nuove variabili
    for(i=0; i<rows; i++)
    {
        xor = Cudd_ReadLogicZero(manager);
        Cudd_Ref(xor);

        for(j=0; j<columns; j++)
            if(a[i][j])
            {
                tmp1 = Cudd_bddXor(manager, xor, vars[j]);
                Cudd_Ref(tmp1);

                Cudd_RecursiveDeref(manager, xor);
                xor = tmp1;
            }

        tmp3 = Cudd_bddCompose(manager, fk->on_set[0], xor, i);
        Cudd_Ref(tmp3);

        Cudd_RecursiveDeref(manager, fk->on_set[0]);
        fk->on_set[0] = tmp3;
    }

    //effettua il parsing della funzione originale senza considerare
    //le variabili alfa
    boolean_function_t* f = parse_pla(manager, argv[1], 0);

    //permuta le variabili di fk in modo da avere
    //lo stesso ordinamento di quelle in f
    int perm[columns+rows];

    for(i=0; i<columns+rows; i++)
    {
        if(i<rows) perm[i] = rows+i+1;
        else perm[i] = i-rows;
    }

    tmp4 = Cudd_bddPermute(manager, fk->on_set[0], perm);
    Cudd_Ref(tmp4);

    Cudd_RecursiveDeref(manager, fk->on_set[0]);
    fk->on_set[0] = tmp4;

    //unione tra on set e dc set di f
    DdNode* u = Cudd_bddOr(manager, f->on_set[0], f->dc_set[0]);
    Cudd_Ref(u);
    
    //on-set di f deve essere incluso in quello di fk
    int cond1 = Cudd_bddLeq(manager, f->on_set[0], fk->on_set[0]);

    //on-set di fk deve essere incluso nell'unione tra
    //dc-set e on-set di f, in questo modo si verifica
    //che fk ha lo stesso comportamento di f
    int cond2 = Cudd_bddLeq(manager, fk->on_set[0], u);

    Cudd_RecursiveDeref(manager, u);

    free_matrix(a, rows);
    Cudd_Quit(manager);

    fprintf(stdout, "%s | Condizione 1: %d | Condizione 2: %d\n", argv[1], cond1, cond2);

    return (cond1 && cond2);
}
