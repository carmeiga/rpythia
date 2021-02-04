rm *.txt
rm main00


make main00
./main00 > saida.txt

grep -n 'Event Listing' saida.txt | cut -f1 -d: > separadores.txt


