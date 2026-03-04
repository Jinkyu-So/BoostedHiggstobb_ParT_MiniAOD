#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <vector>
#include <cmath>

void check_energy() {

    TFile *file = TFile::Open("/home/sjk9523/htobb/higgsToBB_testfiles/ntuple_merged_0.root");
    TTree *tree = (TTree*)file->Get("deepntuplizer/tree");

    std::vector<float> *ptrel = nullptr;
    std::vector<float> *etarel = nullptr;
    std::vector<float> *erel = nullptr;

    tree->SetBranchAddress("pfcand_ptrel", &ptrel);
    tree->SetBranchAddress("pfcand_etarel", &etarel);
    tree->SetBranchAddress("pfcand_erel", &erel);

    Long64_t nentries = tree->GetEntries();

    double max_rel_diff = 0.0;
    double avg_rel_diff = 0.0;
    long long count = 0;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        for (size_t j = 0; j < ptrel->size(); j++) {

            double E_calc = (*ptrel)[j] * std::cosh((*etarel)[j]);
            double E_true = (*erel)[j];

            if (E_true == 0) continue;

            double rel_diff = std::abs(E_calc - E_true) / std::abs(E_true);

            avg_rel_diff += rel_diff;
            if (rel_diff > max_rel_diff)
                max_rel_diff = rel_diff;

            count++;
        }
    }

    avg_rel_diff /= count;

    std::cout << "Checked particles: " << count << std::endl;
    std::cout << "Average relative difference: " << avg_rel_diff << std::endl;
    std::cout << "Max relative difference: " << max_rel_diff << std::endl;

    file->Close();
}
