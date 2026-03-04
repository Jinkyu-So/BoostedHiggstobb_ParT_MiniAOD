#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <cmath>

void check_MeanStd() {

    TFile *file = TFile::Open("/home/sjk9523/htobb/higgsToBB_testfiles/ntuple_merged_0.root");
    TTree *tree = (TTree*)file->Get("deepntuplizer/tree");

    float_t fj_pt;
    float_t fj_e;
    float_t fj_mass;
    float_t fj_eta;
    std::vector<float> *pfcand_dxysig = nullptr;
    std::vector<float> *pfcand_dzsig  = nullptr;
    std::vector<float> *pfcand_ptrel = nullptr;
    std::vector<float> *pfcand_erel = nullptr;
    std::vector<float> *pfcand_deltaR = nullptr;

    tree->SetBranchAddress("fj_pt", &fj_pt);
    tree->SetBranchAddress("pfcand_ptrel", &pfcand_ptrel);
    tree->SetBranchAddress("pfcand_erel", &pfcand_erel);
    tree->SetBranchAddress("pfcand_deltaR", &pfcand_deltaR);
    tree->SetBranchAddress("fj_eta", &fj_eta);
    tree->SetBranchAddress("fj_mass", &fj_mass);
    tree->SetBranchAddress("pfcand_dxysig", &pfcand_dxysig);
    tree->SetBranchAddress("pfcand_dzsig", &pfcand_dzsig);

    ///five variables

    // ---- accumulators ----
    double sum_pt_log = 0, sum2_pt_log = 0;
    double sum_ptrel_log = 0, sum2_ptrel_log = 0;
    double sum_e_log = 0, sum2_e_log = 0;
    double sum_erel_log = 0, sum2_erel_log = 0;
    double sum_deltaR = 0, sum2_deltaR = 0;
    double sum_dxysig = 0, sum2_dxysig = 0;
    double sum_dzsig = 0, sum2_dzsig = 0;

    long long count = 0;

    auto clip = [](double x, double minv, double maxv) {
        if (x < minv) return minv;
        if (x > maxv) return maxv;
        return x;
    };

    // ---- event loop ----
    Long64_t nEntries = tree->GetEntries();

    for (Long64_t i = 0; i < nEntries; ++i) {

        tree->GetEntry(i);

        double p = fj_pt * std::cosh(fj_eta);
        double fj_e = std::sqrt(p * p + fj_mass * fj_mass);

        // loop over PF candidates
        for (size_t j = 0; j < pfcand_ptrel->size(); ++j) {

            double ptrel = pfcand_ptrel->at(j);
            double erel  = pfcand_erel->at(j);
            double deltaR = pfcand_deltaR->at(j);
            double dxysig = pfcand_dxysig->at(j);
            double dzsig  = pfcand_dzsig->at(j);
            
            // ---- clipping ----
            dxysig = clip(dxysig, -20.0, 20.0);
            dzsig  = clip(dzsig,  -20.0, 20.0);
            

            // skip invalid or padded entries
            if (ptrel <= 0 || erel <= 0) continue;

            double pt_log       = std::log(ptrel * fj_pt);
            double ptrel_log    = std::log(ptrel);
            double e_log        = std::log(erel * fj_e);
            double erel_log     = std::log(erel);
    

            // accumulate
            sum_pt_log      += pt_log;
            sum2_pt_log     += pt_log * pt_log;

            sum_ptrel_log   += ptrel_log;
            sum2_ptrel_log  += ptrel_log * ptrel_log;

            sum_e_log       += e_log;
            sum2_e_log      += e_log * e_log;

            sum_erel_log    += erel_log;
            sum2_erel_log   += erel_log * erel_log;

            sum_deltaR      += deltaR;
            sum2_deltaR     += deltaR * deltaR;

            sum_dxysig  += dxysig;
            sum2_dxysig += dxysig * dxysig;
            
            sum_dzsig   += dzsig;
            sum2_dzsig  += dzsig * dzsig;

            count++;
        }
    }

    // ---- compute mean & std ----
    auto compute = [&](double sum, double sum2) {
        double mean = sum / count;
        double var  = sum2 / count - mean * mean;
        double std  = std::sqrt(var);
        return std::make_pair(mean, std);
    };

    auto pt_stats      = compute(sum_pt_log, sum2_pt_log);
    auto ptrel_stats   = compute(sum_ptrel_log, sum2_ptrel_log);
    auto e_stats       = compute(sum_e_log, sum2_e_log);
    auto erel_stats    = compute(sum_erel_log, sum2_erel_log);
    auto deltaR_stats  = compute(sum_deltaR, sum2_deltaR);
    auto dxysig_stats  = compute(sum_dxysig, sum2_dxysig);
    auto dzsig_stats   = compute(sum_dzsig, sum2_dzsig);

    // ---- print results ----
    std::cout << "[pfcand_pt_log, "      << pt_stats.first      << ", " << pt_stats.second      << "]\n";
    std::cout << "[pfcand_logptrel, "    << ptrel_stats.first   << ", " << ptrel_stats.second   << "]\n";
    std::cout << "[pfcand_e_log, "       << e_stats.first       << ", " << e_stats.second       << "]\n";
    std::cout << "[pfcand_logerel, "     << erel_stats.first    << ", " << erel_stats.second    << "]\n";
    std::cout << "[pfcand_deltaR, "      << deltaR_stats.first  << ", " << deltaR_stats.second  << "]\n";
    std::cout << "[pfcand_dxysig, "      << dxysig_stats.first  << ", " << dxysig_stats.second  << "]\n";
    std::cout << "[pfcand_dzsig, "      << dzsig_stats.first  << ", " << dzsig_stats.second  << "]\n";

        ///calculate for the five varibles


    }

