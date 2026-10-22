#include "fixed_matrix.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc < 2 || argc > 5)
    {
        return 2;
    }

    if(std::string(argv[1]) == "--canonicalize")
    {
        // Burdens: clone1 = 10, clone2 = 5, clone3 = 5; the clone2/clone3 tie is kept in input order
        IntMatrix A = {{1, 3, 1, 1}, {1, 2, 1, 1}};
        IntMatrix B = {{1, 1, 0, 0}, {1, 1, 1, 1}};
        DoubleMatrix U = {{0.1}, {0.5}, {0.3}, {0.1}};
        canonicalizeCloneOrder(A, B, U);
        const IntMatrix expectedA = {{1, 1, 1, 3}, {1, 1, 1, 2}};
        const IntMatrix expectedB = {{1, 0, 0, 1}, {1, 1, 1, 1}};
        const DoubleMatrix expectedU = {{0.1}, {0.3}, {0.1}, {0.5}};
        return A == expectedA && B == expectedB && U == expectedU ? 0 : 3;
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
        int base = 1;
        bool ampdel = false;
        IntSet ampdelExempt;
        int cnstates = -1;
        if(argc >= 4)
        {
            const std::string scenario(argv[3]);
            if(scenario == "clonal")
            {
                copyNumbers[0] = CNState(1, 1);
            }
            else if(scenario == "ampdel")
            {
                base = 2;
                ampdel = true;
            }
            else if(scenario == "ampdel-exempt")
            {
                base = 2;
                ampdel = true;
                ampdelExempt.insert(0);
            }
            else if(scenario == "cnstates")
            {
                cnstates = 1;
            }
            else if(scenario == "zero-weight")
            {
                weights[0] = 0.0;
            }
        }
        const double maxHomdelSize = argc == 5 ? std::stod(argv[4]) : 0.005;
        fixed.validate(6, 0.03, base, ampdel, ampdelExempt, cnstates, copyNumbers, weights, maxHomdelSize);

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
            std::fabs(U[1][0] - 0.3) > TOL || std::fabs(U[2][0] - 0.5) > TOL || !fixed.addressesTumorClones()))
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
