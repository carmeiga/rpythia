#!/bin/bash

# rm *.txt
# rm main00

pro=$1

make main00
#./main00 $pro > saida_$pro.txt
./main00 s > saida_s.txt
./main00 h > saida_h.txt
./main00 t > saida_t.txt
./main00 w > saida_w.txt

#grep -n 'Event Listing' saida_$pro.txt | cut -f1 -d: > separadores_$pro.txt
grep -n 'Event Listing' saida_s.txt | cut -f1 -d: > separadores_s.txt
grep -n 'Event Listing' saida_h.txt | cut -f1 -d: > separadores_h.txt
grep -n 'Event Listing' saida_w.txt | cut -f1 -d: > separadores_w.txt
grep -n 'Event Listing' saida_t.txt | cut -f1 -d: > separadores_t.txt


