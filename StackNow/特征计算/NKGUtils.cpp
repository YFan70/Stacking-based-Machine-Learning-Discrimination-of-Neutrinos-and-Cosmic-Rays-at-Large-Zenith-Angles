#include "NKGUtils.h"
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TMinuit.h>
#include <iostream>
#include <cmath>

Event gEM;

static std::vector<double>* sHitR = nullptr;
static NKGFunc gCurrentNKG = nullptr;
static const std::vector<double>* sWeights = nullptr;

double NKG_EM(double r, double s, double Ne) {
    if (r <= 0 || Ne <= 0) return 0.0;
    const double rm = 120.0;
    double a = 1.269 * s + 1.623;
    double term1 = Ne / (2.0 * TMath::Pi() * rm * rm);
    double gamma_factor = TMath::Gamma(a - s) / (TMath::Gamma(s) * TMath::Gamma(a - 2.0 * s));
    return term1 * gamma_factor * pow(r / rm, s - 2.0) * pow(1.0 + r / rm, s - a);
}

double NKG_MU(double r, double s, double Ne) {
    if (r <= 0 || Ne <= 0) return 0.0;
    const double rm = 725.0;
    double a = 3.754 * s - 0.385;
    double term1 = Ne / (2.0 * TMath::Pi() * rm * rm);
    double gamma_factor = TMath::Gamma(a - s) / (TMath::Gamma(s) * TMath::Gamma(a - 2.0 * s));
    return term1 * gamma_factor * pow(r / rm, s - 2.0) * pow(1.0 + r / rm, s - a);
}

void InitializeData(const char* filename, double& theta, double& phi) {
    TFile* file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file: " << filename << std::endl;
        return;
    }

    TTree* pri = (TTree*)file->Get("pri");
    TTree* sec = (TTree*)file->Get("sec");

    if (!pri || !sec) {
        std::cerr << "Missing 'pri' or 'sec' tree." << std::endl;
        file->Close();
        return;
    }

    double theta_deg, phi_deg;
    pri->SetBranchAddress("priTheta", &theta_deg);
    pri->SetBranchAddress("priPhi", &phi_deg);
    pri->GetEntry(0);
    theta = theta_deg * TMath::DegToRad();
    phi = phi_deg * TMath::DegToRad();

    int secId;
    double secX, secY, secWeight = 1.0;
    bool useWeight = sec->GetBranch("secWeight") != nullptr;
    sec->SetBranchAddress("secId", &secId);
    sec->SetBranchAddress("secX", &secX);
    sec->SetBranchAddress("secY", &secY);
    if (useWeight) sec->SetBranchAddress("secWeight", &secWeight);

    gEM.EhitR.clear();
    gEM.MhitR.clear();
    gEM.Eweights.clear();
    gEM.Mweights.clear();
    gEM.hasWeights = useWeight;

    for (Long64_t i = 0; i < sec->GetEntries(); ++i) {
        sec->GetEntry(i);
        double x_m = secX / 100.0;
        double y_m = secY / 100.0;

        double vx = sin(theta) * cos(phi);
        double vy = sin(theta) * sin(phi);
        double vz = cos(theta);
        double cx = y_m * vz;
        double cy = -x_m * vz;
        double cz = x_m * vy - y_m * vx;
        double R = sqrt(cx * cx + cy * cy + cz * cz);

        if (secId == 1 || secId == 2 || secId == 3) {
            gEM.EhitR.push_back(R);
            if (useWeight) gEM.Eweights.push_back(secWeight);
        } else if (secId == 5 || secId == 6) {
            gEM.MhitR.push_back(R);
            if (useWeight) gEM.Mweights.push_back(secWeight);
        }
    }
    file->Close();
}

// Minuit-compatible FCN
void FCN_Generic(int& npar, double* gin, double& f, double* par, int iflag) {
    double s = par[0];
    double Ne = par[1];
    const int nbin = 40;
    const double rmin = 10.0, rmax = 600.0;
    const double bin_width = (rmax - rmin) / nbin;

    std::vector<double> bin_count(nbin, 0.0);
    for (size_t i = 0; i < sHitR->size(); ++i) {
        double r = (*sHitR)[i];
        if (r < rmin || r > rmax) continue;
        int bin = static_cast<int>((r - rmin) / bin_width);
        if (bin >= 0 && bin < nbin)
            bin_count[bin] += (sWeights ? (*sWeights)[i] : 1.0);
    }

    double chi2 = 0.0;
    for (int i = 0; i < nbin; ++i) {
        double r1 = rmin + i * bin_width;
        double r2 = r1 + bin_width;
        double r_center = 0.5 * (r1 + r2);
        double area = TMath::Pi() * (r2 * r2 - r1 * r1);
        double pred = gCurrentNKG(r_center, s, Ne) * area;
        double obs = bin_count[i];
        if (pred > 0) chi2 += pow(obs - pred, 2) / pred;
    }
    f = chi2;
}

bool PerformFit(std::vector<double>& hits, NKGFunc nkg, const char* label,
                double& out_s, double& out_Ne, double& out_chi2ndf,
                const std::vector<double>* weights) {
    if (hits.size() < 10) return false;

    sHitR = &hits;
    gCurrentNKG = nkg;
    sWeights = weights;

    TMinuit minuit(2);
    minuit.SetFCN(FCN_Generic);
    minuit.DefineParameter(0, "s", 1.3, 0.05, 0.5, 2.5);
    minuit.DefineParameter(1, "Ne", hits.size(), 10.0, 0.7 * hits.size(), 4.0 * hits.size());

    if (minuit.Migrad() != 0) return false;

    double err;
    minuit.GetParameter(0, out_s, err);
    minuit.GetParameter(1, out_Ne, err);

    // chi2/ndf post-fit calculation
    const int nbin = 40;
    const double rmin = 10.0, rmax = 600.0;
    const double bin_width = (rmax - rmin) / nbin;
    std::vector<double> bin_count(nbin, 0.0);
    for (size_t i = 0; i < hits.size(); ++i) {
        double r = hits[i];
        if (r < rmin || r > rmax) continue;
        int bin = static_cast<int>((r - rmin) / bin_width);
        if (bin >= 0 && bin < nbin)
            bin_count[bin] += (weights ? (*weights)[i] : 1.0);
    }

    double chi2 = 0.0;
    int ndf = 0;
    for (int i = 0; i < nbin; ++i) {
        double r1 = rmin + i * bin_width;
        double r2 = r1 + bin_width;
        double r_center = 0.5 * (r1 + r2);
        double area = TMath::Pi() * (r2 * r2 - r1 * r1);
        double obs = bin_count[i];
        double pred_density = nkg(r_center, out_s, out_Ne);
        double pred = pred_density * area;
        if (obs > 0 && pred > 1e-12) {
            chi2 += pow(obs - pred, 2) / pred;
            ++ndf;
        }
    }

    out_chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;
    return true;
}

void CalculateStatsForRcuts(const std::vector<double>& data,
                             const std::vector<int>& rcuts,
                             std::map<std::string, double>& results,
                             const std::vector<double>* weights) {
    for (int rc : rcuts) {
        double sum = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, count = 0.0;
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i] > rc) continue;
            double w = weights ? (*weights)[i] : 1.0;
            double x = data[i];
            sum += w * x;
            sum2 += w * x * x;
            sum3 += w * x * x * x;
            sum4 += w * x * x * x * x;
            count += w;
        }

        std::string tag = std::to_string(rc);
        if (count > 0.0) {
            double mean = sum / count;
            double var = (sum2 - sum * sum / count) / (count - 1);
            double skew = 0, kurt = 0;
            if (var > 1e-12) {
                double sigma = sqrt(var);
                double mu3 = (sum3 - 3 * sum2 * sum / count + 2 * pow(sum, 3) / (count * count)) / count;
                double mu4 = (sum4 - 4 * sum3 * sum / count + 6 * sum2 * sum * sum / (count * count) - 3 * pow(sum, 4) / pow(count, 3)) / count;
                skew = mu3 / pow(sigma, 3);
                kurt = mu4 / (var * var) - 3.0;
            }
            results["mean_" + tag] = mean;
            results["var_" + tag] = var;
            results["skew_" + tag] = skew;
            results["kurt_" + tag] = kurt;
        } else {
            results["mean_" + tag] = 0;
            results["var_" + tag] = 0;
            results["skew_" + tag] = 0;
            results["kurt_" + tag] = 0;
        }
    }
}

