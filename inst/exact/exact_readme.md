The "runFgsea.R" script calculates the following P-values:
\[
\mathrm{P} \left( s_{r}^{+} \left(q \right) \geqslant \gamma \right)
\]
 using the FGSEA-multilevel method and writes the results to "multilevelResults.tsv". The script also creates the input files "roundRanks.txt" and "inpPathways.txt" for the exact algorithm. 

The "exact.cpp" script takes as input files "inpPathways.txt." and "roundRanks.txt" and calculates the same probabilities. The results are saved to the file "exactResults.tsv" the form of pairs $(p, \delta)$, where the exact value is guaranteed to lay no further than $\delta$ from $p$.

Examples of commands for running scripts:
1. Rscript runFgsea.R
2. g++ -O2 exact.cpp -o exact.out
3. ./exact.out
