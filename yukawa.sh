make main00
./main00 > saida.dat

grep -n 'Event Listing' saida.dat | cut -f1 -d: > separadores.dat



