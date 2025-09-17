#include "NKGUtils.h"
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <map>
#include <string>
#include <vector>

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " input.root output.root" << std::endl;
        return 1;
    }

    const char* inputFileName = argv[1];
    const char* outputFileName = argv[2];

    double theta = 0.0, phi = 0.0;
    InitializeData(inputFileName, theta, phi);

    double s_e = 0, Ne_e = 0, chi2ndf_e = -1;
    double s_m = 0, Ne_m = 0, chi2ndf_m = -1;

    // 拟合
    PerformFit(gEM.EhitR, NKG_EM, "Electron", s_e, Ne_e, chi2ndf_e,
               gEM.hasWeights ? &gEM.Eweights : nullptr);
    PerformFit(gEM.MhitR, NKG_MU, "Muon", s_m, Ne_m, chi2ndf_m,
               gEM.hasWeights ? &gEM.Mweights : nullptr);

    // 统计
    std::vector<int> rcuts = {200, 400, 600, 800};
    std::map<std::string, double> stats_e, stats_m;
    CalculateStatsForRcuts(gEM.EhitR, rcuts, stats_e,
                           gEM.hasWeights ? &gEM.Eweights : nullptr);
    CalculateStatsForRcuts(gEM.MhitR, rcuts, stats_m,
                           gEM.hasWeights ? &gEM.Mweights : nullptr);

    // 写入 ROOT 文件
    TFile* outputFile = TFile::Open(outputFileName, "RECREATE");
    if (!outputFile || outputFile->IsZombie()) {
        std::cerr << "Error creating output file: " << outputFileName << std::endl;
        return 1;
    }

    TTree* outputTree = new TTree("outputTree", "Fit and statistics");
    outputTree->Branch("s_e", &s_e);
    outputTree->Branch("esize", &Ne_e);
    outputTree->Branch("chi2ndf_e", &chi2ndf_e);
    outputTree->Branch("s_m", &s_m);
    outputTree->Branch("msize", &Ne_m);
    outputTree->Branch("chi2ndf_m", &chi2ndf_m);

    for (int rc : rcuts) {
        std::string tag = std::to_string(rc);
        outputTree->Branch(("e_mean_" + tag).c_str(), &stats_e["mean_" + tag]);
        outputTree->Branch(("e_var_" + tag).c_str(), &stats_e["var_" + tag]);
        outputTree->Branch(("e_skew_" + tag).c_str(), &stats_e["skew_" + tag]);
        outputTree->Branch(("e_kurt_" + tag).c_str(), &stats_e["kurt_" + tag]);
        outputTree->Branch(("m_mean_" + tag).c_str(), &stats_m["mean_" + tag]);
        outputTree->Branch(("m_var_" + tag).c_str(), &stats_m["var_" + tag]);
        outputTree->Branch(("m_skew_" + tag).c_str(), &stats_m["skew_" + tag]);
        outputTree->Branch(("m_kurt_" + tag).c_str(), &stats_m["kurt_" + tag]);
    }

    outputTree->Fill();
    outputTree->Write();
    outputFile->Close();

    std::cout << "Done. Output saved to " << outputFileName << std::endl;
    return 0;
}

