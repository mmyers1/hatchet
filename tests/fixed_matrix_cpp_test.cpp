#include "fixed_matrix.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 2 && argc != 3)
    {
        return 2;
    }

    try
    {
        const int clones = argc == 3 ? std::atoi(argv[2]) : 3;
        std::map<std::string, int> clusters;
        clusters["cluster1"] = 0;
        clusters["cluster2"] = 1;
        clusters["1"] = 0;
        clusters["2"] = 1;
        clusters["3"] = 2;
        std::map<std::string, int> samples;
        samples["sample1"] = 0;
        samples["Tumor1"] = 0;
        samples["Tumor2"] = 1;
        samples["Tumor3"] = 2;

        FixedMatrix fixed = FixedMatrix::read(argv[1], clones, clusters, samples);
        DoubleArray weights(3, 1.0);
        CNMap copyNumbers;
        fixed.validate(6, 0.03, 1, false, -1, copyNumbers, weights);

        IntMatrix A(3, IntArray(clones, 0));
        IntMatrix B(3, IntArray(clones, 0));
        DoubleMatrix U(clones, DoubleArray(3, 0.0));
        for(int segment = 0; segment < 3; ++segment)
        {
            A[segment][0] = 1;
            B[segment][0] = 1;
        }
        fixed.applyCopyNumbers(A, B);
        fixed.conditionProportions(U, 0.03);
        fixed.verify(A, B, U);

        if(argc == 2 &&
           (A[0][1] != 2 || B[1][2] != 1 || std::fabs(U[0][0] - 0.2) > TOL ||
            std::fabs(U[1][0] - 0.3) > TOL || std::fabs(U[2][0] - 0.5) > TOL))
        {
            return 3;
        }
    }
    catch(const std::string& error)
    {
        std::cerr << error << std::endl;
        return 1;
    }
    return 0;
}
