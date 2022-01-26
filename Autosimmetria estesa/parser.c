#include "parser.h"

static DdNode* parse_product(DdManager* manager ,char* input_line, int inputs, int alfa)
{
    // Ciascun prodotto inizialmente è costante uguale ad uno
    DdNode* f = Cudd_ReadOne(manager);
    DdNode* tmp;
    DdNode* x;

    Cudd_Ref(f);

    // Legge la specifica di un termine prodotto
    // se il valore in posizione i è 1 la variabile appare affermata
    // se è 0 negata, se è - non compare.
    for(int i=inputs-1; i>=0; i--)
    {
        if(input_line[i]=='-') continue;

        if(alfa) x = Cudd_bddIthVar(manager, (2*i)+1);
        else x = Cudd_bddIthVar(manager, i);

        if(input_line[i]=='0')
            tmp = Cudd_bddAnd(manager, f, Cudd_Not(x));
        else
            tmp = Cudd_bddAnd(manager, f, x);

        Cudd_Ref(tmp);

        Cudd_RecursiveDeref(manager, f);
        f = tmp;
    }
    
    return f;
}

boolean_function_t* parse_pla(DdManager* manager, char* file_path, int alfa)
{
    if(file_path==NULL) return NULL;

    char line[MAX_LINE_LENGTH];
    int i, read=0, inputs=0, outputs=0, dcs=0;
    FILE* pla_file;

    if((pla_file = fopen(file_path, "r"))==NULL) return NULL;

    //on-set e dc-set di ciascun output della funzione
    DdNode** on_set=NULL;
    DdNode** dc_set=NULL;

    while((read = fscanf(pla_file, "%s\n", line))!=EOF)
    {
        // dbg_printf("Linea letta: %s\n", line);

        if(line[0]=='.' && line[1]=='i')
        {
            fscanf(pla_file, "%d\n", &inputs);

            //alfa è un paramentro che permette di creare le variabili alfa
            //direttamente al momento del parsing. Si veda anche la costruzione
            //di un termine prodotto.
            for(i=0; i<inputs+(alfa*inputs); i++)
                Cudd_bddIthVar(manager, i);

            // dbg_printf("Numero di input della funzione: %i\n", inputs);
        }
        
        else if(line[0]=='.' && line[1]=='o')
        {
            fscanf(pla_file, "%d\n", &outputs);
            // dbg_printf("Numero di output della funzione: %i\n", outputs);

            on_set = calloc(outputs, sizeof(DdNode*));
            dc_set = calloc(outputs, sizeof(DdNode*));

            if(on_set==NULL || dc_set==NULL)
            {
                perror("Errore chiamata calloc ");
                return NULL;
            }

            //inizialmente tutti gli on-set e dc-set sono vuoti
            for(i=0; i<outputs; i++)
            {
                on_set[i] = Cudd_ReadLogicZero(manager);
                dc_set[i] = Cudd_ReadLogicZero(manager);
                Cudd_Ref(on_set[i]);
                Cudd_Ref(dc_set[i]);
            }
        }
        
        else if((line[0]!='.' && line[1]!='o') || (line[0]!='.' && line[1]!='i') || line[0]!='.' && line[1]!='e')
        {
            char out_lines[MAX_LINE_LENGTH], output;
            DdNode* conj;
            DdNode* tmp;

            fscanf(pla_file, "%s\n", out_lines);

            for(i=0; i<outputs; i++)
            {
                output = out_lines[i];

                // dbg_printf("Output letto: %c\n", output);

                //il termine presente sulla linea letta deve essere
                //inserito nell on-set.
                if(output=='1')
                {
                    conj = parse_product(manager, line, inputs, alfa);

                    tmp = Cudd_bddOr(manager, on_set[i], conj);
                    Cudd_Ref(tmp);

                    Cudd_RecursiveDeref(manager, conj);
                    Cudd_RecursiveDeref(manager, on_set[i]);

                    on_set[i] = tmp;
                }

                //il termine presente sulla linea letta deve essere
                //inserito nel dc-set
                if(output=='*' || output=='-')
                {
                    conj = parse_product(manager, line, inputs, alfa);

                    tmp = Cudd_bddOr(manager, dc_set[i], conj);
                    Cudd_Ref(tmp);

                    Cudd_RecursiveDeref(manager, conj);
                    Cudd_RecursiveDeref(manager, dc_set[i]);

                    dc_set[i] = tmp;

                    //conta quanti don't cares sono presenti su tutti
                    //gli output della funzione. Serve per capire se la
                    //funzione trattata è non completamente specificata.
                    dcs++;
                }
            }
        }
    }

    boolean_function_t* f = calloc(1, sizeof(boolean_function_t));

    if(f==NULL)
    {
        perror("Errore chiamata calloc: ");

        for(i=0; i<outputs; i++)
        {
            Cudd_RecursiveDeref(manager, on_set[i]);
            Cudd_RecursiveDeref(manager, dc_set[i]);
        }
        
        free(on_set);
        free(dc_set);
    }

    f->inputs = inputs;
    f->outputs = outputs;
    f->dcs = dcs;
    f->on_set = on_set;
    f->dc_set = dc_set;
    fclose(pla_file);

    return f;
}

int** parse_equations(char* eq_path, int* rows, int* columns)
{
    if(eq_path==NULL || rows==NULL || columns==NULL) return NULL;

    FILE* f = fopen(eq_path, "r");

    if(f==NULL) return NULL;

    fscanf(f, "%d\n", rows);
    fscanf(f, "%d\n", columns);

    int** A = calloc((*rows), sizeof(int*));
    int i, j;

    for(i=0; i<(*rows); i++)
        A[i] = calloc((*columns), sizeof(int));

    for(i=0; i<(*rows); i++)
    {
        for(j=0; j<(*columns); j++)
        {
            if(!fscanf(f, "%d", &A[i][j]))
                break;
        }
    }

    return A;
}

void printPla(DdManager* manager, char* outputfile, DdNode* bdd, int n_var) {
	
    DdGen* gen;
	int* cube;
	CUDD_VALUE_TYPE value;
	int ssize;

	FILE *f = fopen(outputfile, "w");

	if(f == NULL) {
        fprintf(stderr, "Errore durante l'apertura del file: %s\n", outputfile);
        exit(-1);
	}

	ssize = Cudd_SupportSize(manager, bdd);

	if (ssize == 0) { // se la funzione è costante...
        

        fprintf(f, ".i 1\n");

        fprintf(f, ".o 1\n");

        /*
         * dalla documentazione di CUDD, la funzione Cudd_bddPickOneCube
         * restituisce un cubo dell on-set. se la funzione è costante uguale a
         * zero, la chiamata fallisce e restituisce zero. Funziona sempre ma non
         * ho trovato un modo più semplice per farlo.
         */

        char* cube_str = (char* )malloc(n_var*sizeof(char));
        if (Cudd_bddPickOneCube(manager, bdd, cube_str) == 1)
            fprintf(f, "- 1\n");
        else
            fprintf(f, "- 0\n");

        free(cube_str);

	} else {
        
        fprintf(f, ".i %d\n", ssize);
        fprintf(f, ".o 1\n");

        //indici di variabili da cui il bdd dipende
        int* support = Cudd_SupportIndex(manager, bdd);

        Cudd_ForeachCube(manager, bdd, gen, cube, value) {
        
            for (int i = 0; i < n_var; i++) {
            
                if (support[i]) {
                
                    if (cube[i] == 2)
                        fprintf(f, "%c", '-');
                    else
                        fprintf(f, "%d", cube[i]);
                }
            }

            fprintf(f, " 1\n");
        }

        free(support);
	}

	fclose(f);
}