#include <iostream>

using namespace std;

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TSystem.h"

void MergeOutput(Int_t nJobs=9){
  TChain* fChain_Mode0 = new TChain("StpNTuple");
  TChain* fChain_Mode1 = new TChain("StpNTuple");
  TChain* fChain_Final = new TChain("MNTTree");
  TString BasePath = std::getenv("IonCatchersPath");
  TString baseDir = BasePath; baseDir += "/build/";
  TString currDir, fileRun1, fileRun2, fileRun3, sMerge;
  sMerge = "cat ";
  for(Int_t iJob=0; iJob<nJobs; iJob++) {
    currDir = baseDir; currDir+= "run"; currDir += iJob; currDir+= "/";
    cout<<"Processing directory "<<currDir<<endl;
    fileRun1 = currDir;  fileRun1 += "NTuple_Mode0.root";
    fChain_Mode0->Add(fileRun1);
    fileRun2 = currDir;  fileRun2 += "NTuple_Mode1.root";
    fChain_Mode1->Add(fileRun2);
    fileRun3 = currDir;  fileRun3 += "Tree.root";
    fChain_Final->Add(fileRun3);
    sMerge += currDir; sMerge += "Output_Mode0.txt ";
  }
  sMerge += ">Output_Mode0.txt";
  
  fChain_Mode0->Merge("Ntuple_Mode0.root");
  fChain_Mode1->Merge("Ntuple_Mode1.root");
  fChain_Final->Merge("Tree.root");

  gSystem->Exec(sMerge);
}
