#!/bin/bash

for file in output/*-aleatoria.config.out
do
	cat $file >> aleatoria.csv
done