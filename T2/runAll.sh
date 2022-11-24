#!/bin/bash

# Scrip para rodar todos os casos de teste

for input in data/setI/*.sukp data/setII/*.sukp
do
    for tipo in 3 4 5  ; do
	tempo=5
	perc=50
	tam=0.10
	for lns in 0 1  ; do
	    if [ ! -f "${input##*/}-$tipo-$tempo-$lns-$perc-($tam).out" ]; then
		./bin/mochila $input $tipo $tempo $lns $perc $tam > "${input##*/}-$tipo-$tempo-$lns-$perc-($tam).out"
	    fi
	done
    done
done

