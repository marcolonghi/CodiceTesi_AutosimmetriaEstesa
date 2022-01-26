#!/usr/bin/env bash

function run_generalized_Filippo {
    timeout 10m ./generalized_test_Filippo.out $1 >> generalized_test_Filippo.log 2>> generalized_test_Filippo.err

    if [ $? -eq 0 ]; then
        echo "generalized_test_Filippo su pla $1 eseguito correttamente"
    else
        echo "generalized_test_Filippo su pla $1 fallito"
    fi
}

function run_generalized_Filippo_migliorato {
    timeout 10m ./generalized_test_Filippo_migliorato.out $1 >> generalized_test_Filippo_migliorato.log 2>> generalized_test_Filippo_migliorato.err

    if [ $? -eq 0 ]; then
        echo "generalized_test_Filippo_migliorato su pla $1 eseguito correttamente"
    else
        echo "generalized_test_Filippo_migliorato su pla $1 fallito"
    fi
}

function run_all_zeros {
    timeout 10m ./all_zeros_test.out $1 >> all_zeros_test.log 2>> all_zeros_test.err

    if [ $? -eq 0 ]; then
        echo "all_zeros_test su pla $1 eseguito correttamente"
    else
        echo "all_zeros_test su pla $1 fallito"
    fi
}

function run_all_ones {
    timeout 10m ./all_ones_test.out $1 >> all_ones_test.log 2>> all_ones_test.err

    if [ $? -eq 0 ]; then
        echo "all_ones_test su pla $1 eseguito correttamente"
    else
        echo "all_ones_test su pla $1 fallito"
    fi
}

function run_test_reduction {
    timeout 10m ./test_reduction.out $1 $2 $3 >> test_reduction.log 2>> test_reduction.err

    if [ $? -eq 1 ]; then
        echo "test_reduction su pla $1 eseguito correttamente"
    else
        echo "test_reduction su pla $1 fallito"
    fi
}

if [ $1 -eq 0 ] 
then
	for filename in pla_splitDC/*.pla
    		do
	    		run_generalized_Filippo $filename
    		done
fi

if [ $1 -eq 1 ] 
then
	for filename in pla_splitDC/*.pla
    		do
	    		run_generalized_Filippo_migliorato $filename
    		done
fi

if [ $1 -eq 2 ] 
then
	for filename in pla_splitDC/*.pla
    		do
	    		run_all_zeros $filename
    		done
fi

if [ $1 -eq 3 ] 
then
	for filename in pla_splitDC/*.pla
    		do
	    		run_all_ones $filename
    		done
fi

if [ $1 -eq 4 ] 
then
	for filename in pla_splitDC/*.pla
    		do
                SUBSTRING=$(echo $filename| cut -d'/' -f 2)
                run_test_reduction $filename fk/$SUBSTRING/fk0.pla eq/$SUBSTRING/eq0.re
    		done
fi
