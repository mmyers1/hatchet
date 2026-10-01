#include "sample_clonal.h"

#include <iostream>

int main(int argc, char** argv)
{
    if(argc != 2 && argc != 4)
    {
        return 2;
    }
    try
    {
        std::map<std::string, int> clusters;
        clusters["c1"] = 0;
        clusters["c2"] = 1;
        clusters["c3"] = 2;
        std::map<std::string, int> samples;
        samples["Tumor"] = 0;
        const SampleClonal sampleClonal = SampleClonal::read(argv[1], clusters, samples);
        if(argc == 4)
        {
            DoubleArray weights;
            weights.push_back(1.0);
            weights.push_back(99.0);
            if(std::string(argv[2]) == "zero-weight")
            {
                weights[0] = 0.0;
                weights[1] = 100.0;
            }
            DoubleArray purities(1, 0.5);
            sampleClonal.validate(6, 1, false, -1, weights, purities, std::stod(argv[3]));
            return 0;
        }
        DoubleMatrix R(3, DoubleArray(1, 5.0));
        R[0][0] = 2.0 / 3.0;
        R[1][0] = 2.5 / 3.0;
        const SampleClonalCalibration calibration = sampleClonal.calibrate(R, 1);
        if(sampleClonal.entries().size() < 2 || sampleClonal.maxTotalCopyNumber() != 3 ||
           std::fabs(calibration.gamma[0] - 3.0) > TOL ||
           std::fabs(calibration.purity[0] - 0.5) > TOL)
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
