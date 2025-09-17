#ifndef NKGUTILS_H
#define NKGUTILS_H

#include <vector>
#include <string>
#include <map>

struct Event {
    std::vector<double> EhitR;
    std::vector<double> MhitR;
    std::vector<double> Eweights;
    std::vector<double> Mweights;
    bool hasWeights = false;
};

extern Event gEM;

typedef double (*NKGFunc)(double r, double s, double Ne);

double NKG_EM(double r, double s, double Ne);
double NKG_MU(double r, double s, double Ne);

void InitializeData(const char* filename, double& theta, double& phi);

bool PerformFit(std::vector<double>& hits,
                NKGFunc nkg,
                const char* label,
                double& out_s,
                double& out_Ne,
                double& out_chi2ndf,
                const std::vector<double>* weights = nullptr);

void CalculateStatsForRcuts(const std::vector<double>& data,
                             const std::vector<int>& rcuts,
                             std::map<std::string, double>& results,
                             const std::vector<double>* weights = nullptr);

#endif // NKGUTILS_H

