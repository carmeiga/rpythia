#!/bin/bash

# rm *.txt
# rm main00

pro=$1

make main00
./main00 $pro > saida_$pro.txt

grep -n 'Event Listing' saida_$pro.txt | cut -f1 -d: > separadores_$pro.txt


