#include "geseca.h"
#include "ScoreRuler.h"


NumericVector gesecaCpp(const NumericMatrix & E, const NumericVector & inpScores,
                        unsigned genesetSize, unsigned sampleSize, int seed, double eps){
    std::vector<double> scores = as<std::vector<double> > (inpScores);
    std::vector<std::vector<double> > expressionMatrix(E.nrow());
    for (unsigned i = 0; i < E.nrow(); i++){
        NumericVector currentRow = E(i, _);
        expressionMatrix[i] = as<std::vector<double> >(currentRow);
    }

    ScoreRuler ruler(expressionMatrix, sampleSize, genesetSize);
    std::vector<double> pvals;

    double maxScore = *std::max_element(scores.begin(), scores.end());
    ruler.extend(maxScore, seed, eps);
    for (auto score : scores){
        pvals.push_back(ruler.getPvalue(score, eps));
    }
    return wrap(pvals);
}
