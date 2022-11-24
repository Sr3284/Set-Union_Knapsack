#!/bin/bash

# Scrip para rodar todos os casos de teste

for input in data/setI/*.sukp data/setII/*.sukp
do
    for config in output/*.config
    do
	if [ ! -f output/${input##data/}-mochila-${config##output/}.out ]; then
	    xargs --a $config ./bin/mochila $input
	fi
    done
done

