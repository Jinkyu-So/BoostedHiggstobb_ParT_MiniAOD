#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <filesystem>

namespace fs = std::filesystem;

void count_all_events() {

    std::string path = "/home/sjk9523/htobb/higgsToBB_files/";
    long long total_events = 0;

    for (const auto &entry : fs::directory_iterator(path)) {
        std::string filename = entry.path().string();

        // Skip non-root files
        if (filename.substr(filename.find_last_of(".") + 1) != "root") continue;

        TFile *file = TFile::Open(filename.c_str());
        if (!file || file->IsZombie()) {
            std::cerr << "Failed to open: " << filename << std::endl;
            continue;
        }

        TTree *tree = (TTree*)file->Get("deepntuplizer/tree");
        if (!tree) {
            std::cerr << "Tree not found in: " << filename << std::endl;
            file->Close();
            continue;
        }

        Long64_t nentries = tree->GetEntries();
        total_events += nentries;

        file->Close();
    }

    std::cout << "Total events across all files: " << total_events << std::endl;
}