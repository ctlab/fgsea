The "runFgsea.R" script calculates P-values using the FGSEA-multilevel method and writes the results to "multilevelResults.tsv". The script also creates the input files "roundRanks.txt" and "inpPathways.txt" for the exact algorithm. 

The "exact.cpp" script takes as input files "inpPathways.txt." and "roundRanks.txt" and calculates the exact P-values. Exact values are saved to the file "exactResults.tsv".

Examples of commands for running scripts:
1) Rscript runFgsea.R
2) g++ -O2 exact.cpp -o exact.out
3) ./exact.out
