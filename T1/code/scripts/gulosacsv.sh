#!/bin/bash

for file in output/*-gulosa.config.out
do
	cat $file >> gulosa.csv
done