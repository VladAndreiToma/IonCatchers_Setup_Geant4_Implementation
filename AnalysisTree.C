// Code that loads the NTuple_ModeX.root (X=1,2) output files from Geant4
// and creates one output file with a TTree of custom events and tracks

#include <TFile.h>
#include <TString.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TUnixSystem.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <ctime>

using namespace std;

Bool_t printOut = true;

class Volume {
public:
  Int_t vol;     // volume number
  // track properties at entry in volume:
  Float_t ein, tin, pin, qin, xin, yin, zin; // energy, theta, phi, charge, x/y/z position
  //  track properties integrated over volume:
  Int_t nhit;               // number of hits/steps
  Float_t edep, slen; // total deposited energy, total path length
  //  track properties at exit from volume:
  Float_t eout, tout, pout, qout, xout, yout, zout; // energy, theta, phi, charge, x/y/z position

  void FillVolume(Int_t v1, Float_t v2, Float_t v3, Float_t v4, Float_t v5, Float_t v6, Float_t v7, Float_t v8,
		  Int_t v9, Float_t v10, Float_t v11, Float_t v12, Float_t v13, Float_t v14, Float_t v15, Float_t v16, Float_t v17, Float_t v18) {
    vol=v1; ein=v2; tin=v3; pin=v4; qin=v5; xin=v6; yin=v7; zin=v8;
    nhit=v9; edep=v10; slen=v11;
    eout=v12; tout=v13; pout=v14; qout=v15; xout=v16; yout=v17; zout=v18;
  }
  void Reset() { vol=-1; ein=-1.; tin=-1.; pin=-1.; qin=-1.; xin=-1.; yin=-1.; zin=-1.;
    nhit=-1; edep=-1.; slen=-1.; eout=-1.; tout=-1.; pout=-1.; qout=-1.; xout=-1.; yout=-1.; zout=-1.; }
};

class Track {
public:
  Int_t Z;      // nuclear charge
  Int_t A;      // nuclear mass number
  Float_t M;  // nuclear mass
  Int_t par;        // parent number
  UInt_t nvol; // number of volumes in track
  std::vector<Volume> tV; // check size against nvol!
  
  void FillTrack(Int_t v1, Int_t v2, Float_t v3, Int_t v4) { Z=v1; A=v2; M=v3; par=v4; }
  void FillNvol(Int_t v) { nvol=v; }
  void Reset() { Z=-1; A=-1; M=-1.; par=-1; nvol=0; tV.clear(); }
  void PrintTrack(Int_t type) {
    if(type==0) std::cout<<"BEAM: ";
    if(type==1) std::cout<<"FRAG1: ";
    if(type==2) std::cout<<"FRAG2: ";
    std::cout<<type<<" Z="<<this->Z<<", A="<<this->A<<", par="<<this->par<<", nvol="<<this->tV.size()<<std::endl;
      for(UInt_t iv=0; iv<this->tV.size(); iv++)
	std::cout<<"   "<<this->tV[iv].vol<<", "<<this->tV[iv].ein<<", "<<this->tV[iv].nhit<<", "<<this->tV[iv].eout<<std::endl;
  }
  void Copy(Track *source) {
    this->Z = source->Z;
    this->A = source->A;
    this->M = source->M;
    this->par = source->par;
    this->nvol = source->nvol;
    this->tV = source->tV;
  }
};

class Event : public TObject {
public:
  Int_t time;      // days since June 1, 2022
  Int_t run;        // run numer
  Int_t evt;        // event number
  Int_t trk;         // current track number
  Track beam;
  Track frag1;
  Track frag2;
  void Reset() { time=-1; run=-1; evt=-1; trk=-1; beam.Reset(); frag1.Reset(); frag2.Reset(); }
  ClassDef(Event,1)
};
ClassImp(Event)

Bool_t BuildTree(Int_t inRun);
Bool_t MergeTrees();

// MAIN MACRO:
Bool_t AnalysisTree(Bool_t TwoStep = true) {
  if(TwoStep) {
    Bool_t haveT0 = BuildTree(0);
    if(!haveT0) { std::cout<<"AnalysisTree - First tree not generated! Exiting..."<<std::endl; return kFALSE; }
    Bool_t haveT1 = BuildTree(1);
    if(!haveT1) { std::cout<<"AnalysisTree - Second tree not generated! Exiting..."<<std::endl; return kFALSE; }
    Bool_t haveT = MergeTrees();
    if(!haveT) { std::cout<<"AnalysisTree - Trees not merged! Exiting..."<<std::endl; return kFALSE; }
  } else {
    Bool_t haveT2 = BuildTree(2);
    if(!haveT2) { std::cout<<"AnalysisTree - Unique tree not generated! Exiting..."<<std::endl; return kFALSE; }
  }
  std::cout<<"AnalysisTree - DONE!"<<std::endl;
  
  return kTRUE;
}

// FUNCTION THAT BUILDS INDIVIDUAL TREES:
Bool_t BuildTree(Int_t inRun = 1){
  // GET INPUT NTUPLE
  TString inFile = "NTuple_Mode"; inFile += inRun; inFile += ".root";
  TFile* fin = new TFile(inFile);
  if(!fin) { std::cout<<"BuildTree - Couldn't get the file "<<inFile<<"! Exiting..."<<std::endl; return kFALSE; }
  TNtuple* InNtu = (TNtuple*)fin->Get("StpNTuple");
  if(!InNtu) { std::cout<<"BuildTree - Couldn't get the ntuple from "<<inFile<<"! Exiting..."<<std::endl; return kFALSE; }
  Int_t nEntries = InNtu->GetEntries();
  if(nEntries<1) { std::cout<<"BuildTree - Ntuple from "<<inFile<<" is empty! Exiting..."<<std::endl; return kFALSE; }
  std::cout<<"BuildTree - Looping over "<<nEntries<<" ntuple rows from "<<inFile<<std::endl;
  Double_t fHitE, fHitT, fHitP, fm, fx, fy, fz, fedep, fslen, fsq;
  Int_t fevt, ftrk, fpare, fpart, fq, fa, fstp, fdet;
  InNtu->SetBranchAddress("evt",&fevt);
  InNtu->SetBranchAddress("trk",&ftrk);
  InNtu->SetBranchAddress("pare",&fpare);
  InNtu->SetBranchAddress("part",&fpart);
  InNtu->SetBranchAddress("Z",&fq);
  InNtu->SetBranchAddress("A",&fa);
  InNtu->SetBranchAddress("M",&fm);
  InNtu->SetBranchAddress("E",&fHitE);
  InNtu->SetBranchAddress("T",&fHitT);
  InNtu->SetBranchAddress("P",&fHitP);
  InNtu->SetBranchAddress("stp",&fstp);
  InNtu->SetBranchAddress("x",&fx);
  InNtu->SetBranchAddress("y",&fy);
  InNtu->SetBranchAddress("z",&fz);
  InNtu->SetBranchAddress("det",&fdet);
  InNtu->SetBranchAddress("edep",&fedep);
  InNtu->SetBranchAddress("slen",&fslen);
  InNtu->SetBranchAddress("Q",&fsq);
  
  // SETUP THE OUTPUT TREE
  TString outFile = "Tree_Mode"; outFile += inRun; outFile += ".root";
  TFile* fout=new TFile(outFile,"RECREATE");
  TTree* OutTree = new TTree("OutTree","tree of tracks");
  Event* mnt = new Event;
  mnt->Reset();
  OutTree->Branch("e",&mnt);

  // GET RUN NUMBER AND DATE:
  std::string runDir = gSystem->GetWorkingDirectory();
  if(runDir.empty()) { std::cout<<"BuildTree - Cannot get run directory! Exiting..."<<std::endl; return kFALSE; }
  std::string runStr = runDir.substr(runDir.find("run")+3);
  Int_t runNum = -1;
  if(runStr.length()<3) runNum = std::stoi(runStr);
  time_t timeNow = time(0); // in seconds since 00:00, Jan 1, 1970
  struct tm timeRef = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  timeRef.tm_year = 122; timeRef.tm_mon = 5; timeRef.tm_mday = 1; // 00:00, June 1, 2022
  Int_t runTime = Int_t(difftime(timeNow, mktime(&timeRef))/24./3600.);
  
  // LOOP OVER INPUT NTUPLE
  Bool_t FirstVolume = kTRUE, FirstTrack = kTRUE, FirstEvent = kTRUE;
  Bool_t newEvent = kFALSE, newTrack = kFALSE, newVolume = kFALSE;
  Int_t currEvt=-1, currTrk=-1, currVol=-1, iVol=0;
  Int_t  inEvt=0, inPar=-1, inTrk=0, inVol=-1., inZ=0, inA=0, NHit=0;
  Float_t inE=0., inT=0., inP=0., inQ=0., inx=0., iny=0., inz=0.;
  Float_t  inM=0., EDep=0., SLen=0.;
  Float_t outE=0., outT=0., outP=0., outQ=0., outx=0., outy=0., outz=0.;
  div_t DivByMil;
  for(Int_t iRow=0; iRow<nEntries; iRow++){
    DivByMil = std::div(iRow,1000000);
    if(DivByMil.rem == 0)  std::cout<<"Row "<<iRow<<std::endl;
    InNtu->GetEntry(iRow);

    newEvent = kFALSE; newTrack = kFALSE; newVolume = kFALSE;
    if(currEvt!=fevt) { newEvent = kTRUE; currEvt = fevt; currTrk = -1; currVol = -1; } // if newEvt force newTrk, newVol
    if(currTrk!=ftrk)  { newTrack = kTRUE; currTrk = ftrk; currVol = -1; } // if newTrk force newVol
    if(currVol!=fdet) { newVolume = kTRUE; currVol = fdet; }
    
    if(newVolume) {
      if(!FirstVolume) {
	if(printOut) std::cout<<"Evt/Trk/Par/Vol:"<<inEvt<<"/"<<inTrk<<"/"<<inPar<<"/"<<inVol<<" iVol="<<iVol
			      <<" Z/A="<<inZ<<"/"<<inA<<"  Ein="<<inE<<" NHit="<<NHit<<" Eout="<<outE<<std::endl;
	mnt->time = runTime;
	mnt->run = runNum;
	mnt->evt = inEvt;
	mnt->trk = inTrk;
	Volume tempVol;
	tempVol.FillVolume(inVol, inE, inT, inP, inQ, inx, iny, inz, NHit, EDep, SLen, outE, outT, outP, outQ, outx, outy, outz);
	if(mnt->trk==0) {
	  mnt->beam.FillTrack(inZ, inA, inM, inPar);
	  mnt->beam.tV.push_back(tempVol);
	} else if(mnt->trk==1) {
	  mnt->frag1.FillTrack(inZ, inA, inM, inPar);
	  mnt->frag1.tV.push_back(tempVol);
	} else if(mnt->trk==2) {
	  mnt->frag2.FillTrack(inZ, inA, inM, inPar);
	  mnt->frag2.tV.push_back(tempVol);
	} else { std::cout<<"Error - track index: "<<inTrk<<std::endl; return kFALSE;}
	tempVol.Reset();
      }
      // set track parameters:
      inEvt = fevt; inTrk =  ftrk-1; inPar = fpare; inZ = fq; inA = fa; inM = fm;
      // set volume entry parameters:
      inVol = fdet; inE = fHitE; inT = fHitT; inP = fHitP; inQ = fsq; inx = fx; iny = fy; inz = fz;
      // initialize volume integrated parameters:
      NHit = 0; EDep = 0.; SLen = 0.;
      FirstVolume = kFALSE;
      iVol++;
    }
    // set volume integrated parameters:
    NHit += 1.; EDep += fedep; SLen += fslen;
    // set volume exit parameters:
    outE = fHitE; outE = fHitE; outT = fHitT; outP = fHitP; outQ = fsq; outx = fx; outy = fy; outz = fz;
    
    if(newTrack) {
      if(!FirstTrack) {
	if(mnt->trk==0) {
	  mnt->beam.FillNvol(iVol);
	  if(printOut) mnt->beam.PrintTrack(mnt->trk);
	} else if(mnt->trk==1) {
	  mnt->frag1.FillNvol(iVol);
	  if(printOut) mnt->frag1.PrintTrack(mnt->trk);
	} else if(mnt->trk==2) {
	  mnt->frag2.FillNvol(iVol);
	  if(printOut) mnt->frag2.PrintTrack(mnt->trk);
	}
      }
      FirstTrack = kFALSE;
      iVol = 0;
    }

    if(newEvent) {
      if(!FirstEvent) {
	OutTree->Fill();
	mnt->Reset();
      }
      FirstEvent = kFALSE;
    }
    
    // ftxt<<"Evt/Trk/Det/Stp="<<fevt<<"/"<<ftrk<<"/"<<fdet<<"/"<<fstp<<" HitE/T/P="<<fHitE<<"/"<<fHitT<<"/"<<fHitP<<" a/z/q="<<fa<<"/"<<fq<<"/"<<fsq<<" edep/slen="<<fedep<<"/"<<fslen<<" x/y/z="<<fx<<"/"<<fy<<"/"<<fz<<endl;
  } // loop over Ntuple steps

  std::cout<<"BuildTree - last event addition:"<<std::endl<<std::endl;
  OutTree->Fill();
  fout->Write();
  fout->Close();

  return kTRUE;
}

// FUNCTION THAT MERGES THE TWO TREES:
Bool_t MergeTrees(){
  // GET THE FIRST TREE
  TFile* fin1=new TFile("Tree_Mode0.root","READ");
  if(!fin1) { std::cout<<"MergeTrees - Couldn't get the file Tree_Mode0.root! Exiting..."<<std::endl; return kFALSE; }
  TTree* InTree1 = (TTree*)fin1->Get("OutTree");
  if(!InTree1) { std::cout<<"MergeTrees - Couldn't get the tree from Tree_Mode0.root! Exiting..."<<std::endl; return kFALSE; }
  Int_t nEntries = InTree1->GetEntries();
  if(nEntries<1) { std::cout<<"MergeTrees - Tree from Tree_Mode0.root is empty! Exiting..."<<std::endl; return kFALSE; }
  std::cout<<"MergeTrees - Looping over "<<nEntries<<" tree rows from Tree_Mode0.root"<<std::endl;
  Event* mnt1 = new Event;
  mnt1->Reset();
  InTree1->SetBranchAddress("e",&mnt1);

  // GET THE SECOND TREE AND BUILD AN EVENT INDEX
  TFile* fin2=new TFile("Tree_Mode1.root","READ");
  if(!fin2) { std::cout<<"MergeTrees - Couldn't get the file Tree_Mode1.root! Exiting..."<<std::endl; return kFALSE; }
  TTree* InTree2 = (TTree*)fin2->Get("OutTree");
  if(!InTree2) { std::cout<<"MergeTrees - Couldn't get the tree from Tree_Mode1.root! Exiting..."<<std::endl; return kFALSE; }
  Int_t nIdx = InTree2->BuildIndex("evt","run");
  if(nIdx<1) { std::cout<<"MergeTrees - Tree index from Tree_Mode1.root is empty! Exiting..."<<std::endl; return kFALSE; }
  std::cout<<"MergeTrees - Built tree index with "<<nIdx<<" entries from Tree_Mode1.root."<<std::endl;
  if(nEntries!=nIdx) { std::cout<<"MergeTrees - number of entries mismatch: "<<nEntries<<"/"<<nIdx<<std::endl; return kFALSE; }
  Event *mnt2 = new Event;
  mnt2->Reset();
  InTree2->SetBranchAddress("e",&mnt2);
  
  // SETUP THE OUTPUT TREE
  TFile* fout=new TFile("Tree.root","RECREATE");
  TTree* OutTree = new TTree("MNTTree","tree of tracks");
  Event* mnt = new Event;
  mnt->Reset();
  OutTree->Branch("e",&mnt);

  // LOOP OVER FIRST TREE AND INSERT EVENT FROM SECOND TREE
  div_t DivByMil;
  //  for(Int_t iRow=0; iRow<nEntries; iRow++){
  for(Int_t iRow=0; iRow<nEntries-1; iRow++){   // do not process last event!
    DivByMil = std::div(iRow,1000000);
    if(DivByMil.rem == 0)  std::cout<<"Row "<<iRow<<std::endl;
    InTree1->GetEntry(iRow);

    Int_t nRun = mnt1->run;
    Int_t nEvt = mnt1->evt;
    if(printOut)
      std::cout<<"File1: Evt="<<mnt1->evt<<" Run="<<mnt1->run
	       <<" Beam:Z/A/nvol="<<mnt1->beam.Z<<"/"<<mnt1->beam.A<<"/"<<mnt1->beam.nvol<<"-"<<mnt1->beam.tV.size()
	       <<" Frag1:Z/A/nvol="<<mnt1->frag1.Z<<"/"<<mnt1->frag1.A<<"/"<<mnt1->frag1.nvol<<"-"<<mnt1->frag1.tV.size()
	       <<" Frag2:Z/A/nvol="<<mnt1->frag2.Z<<"/"<<mnt1->frag2.A<<"/"<<mnt1->frag2.nvol<<"-"<<mnt1->frag2.tV.size()<<std::endl;

    Int_t nBytes = InTree2->GetEntryWithIndex(nEvt,nRun);
    if(nBytes<1) { std::cout<<"Event #"<<nEvt<<" from Run #"<<nRun<<" from second tree is empty! Exiting..."<<std::endl; return kFALSE; }
    if(printOut) 
      std::cout<<"File2: Evt="<<mnt2->evt<<" Run="<<mnt2->run
	       <<" Frag1:Z/A/nvol="<<mnt2->frag1.Z<<"/"<<mnt2->frag1.A<<"/"<<mnt2->frag1.nvol<<"-"<<mnt2->frag1.tV.size()
	       <<" Frag2:Z/A/nvol="<<mnt2->frag2.Z<<"/"<<mnt2->frag2.A<<"/"<<mnt2->frag2.nvol<<"-"<<mnt2->frag2.tV.size()<<std::endl;

    // check array sizes:
    if(mnt1->beam.nvol!=mnt1->beam.tV.size()) {std::cout<<"MergeTrees - Beam array size error!"<<std::endl; return kFALSE;}
    if(mnt1->frag1.nvol>0 && mnt1->frag1.nvol!=mnt1->frag1.tV.size()) {std::cout<<"MergeTrees - Frag11 array size error!"<<std::endl; return kFALSE;}
    if(mnt1->frag2.nvol>0 && mnt1->frag2.nvol!=mnt1->frag2.tV.size()) {std::cout<<"MergeTrees - Frag12 array size error!"<<std::endl; return kFALSE;}
    if(mnt2->frag1.nvol>0 && mnt2->frag1.nvol!=mnt2->frag1.tV.size()) {std::cout<<"MergeTrees - Frag21 array size error!"<<std::endl; return kFALSE;}
    if(mnt2->frag2.nvol>0 && mnt2->frag2.nvol!=mnt2->frag2.tV.size()) {std::cout<<"MergeTrees - Frag22 array size error!"<<std::endl; return kFALSE;}

    // copy event from 1st tree:
    mnt->time = mnt1->time;
    mnt->run = nRun;
    mnt->evt = nEvt;
    mnt->trk = mnt1->trk;
    mnt->beam.Copy(&mnt1->beam);
    mnt->frag1.Copy(&mnt1->frag1);
    mnt->frag2.Copy(&mnt1->frag2);
    
    // make fragment 1 association:
    if(mnt2->frag1.Z>0 && mnt2->frag1.A>0) { // fragment 1 exists (exit the target)
      if(mnt2->frag1.Z==mnt1->frag1.Z && mnt2->frag1.A==mnt1->frag1.A) { // insert MNT2.1 into MNT.1
	if(mnt1->frag1.nvol<1 || mnt2->frag1.nvol<1) {std::cout<<"MergeTrees - array size error 1!"<<std::endl; return kFALSE;}
	mnt->frag1.nvol += mnt2->frag1.nvol;
	for(UInt_t iv=0; iv<mnt2->frag1.tV.size(); iv++) mnt->frag1.tV.push_back(mnt2->frag1.tV[iv]);
      } else if(mnt2->frag1.Z==mnt1->frag2.Z && mnt2->frag1.A==mnt1->frag2.A) { // insert MNT2.1 into MNT.2
	if(mnt1->frag2.nvol<1 || mnt2->frag1.nvol<1) {std::cout<<"MergeTrees - array size error 2!"<<std::endl; return kFALSE;}
	mnt->frag2.nvol += mnt2->frag1.nvol;
	for(UInt_t iv=0; iv<mnt2->frag1.tV.size(); iv++) mnt->frag2.tV.push_back(mnt2->frag1.tV[iv]);
      } else  { std::cout<<"MergeTrees - Couldn't make fragment 1 association! Exiting..."<<std::endl; return kFALSE; }
    } else { if(printOut) std::cout<<"     Fragment 1 did not exit target..."<<std::endl; }
    
    // make fragment 2 association:
    if(mnt2->frag2.Z>0 && mnt2->frag2.A>0) { // fragment 2 exists (exit the target)
      if(mnt2->frag2.Z==mnt1->frag1.Z && mnt2->frag2.A==mnt1->frag1.A) { // insert MNT2.2 into MNT.1
	if(mnt1->frag1.nvol<1 || mnt2->frag2.nvol<1) {std::cout<<"MergeTrees - array size error 3!"<<std::endl; return kFALSE;}
	mnt->frag1.nvol += mnt2->frag2.nvol;
	for(UInt_t iv=0; iv<mnt2->frag2.tV.size(); iv++) mnt->frag1.tV.push_back(mnt2->frag2.tV[iv]);
      } else if(mnt2->frag2.Z==mnt1->frag2.Z && mnt2->frag2.A==mnt1->frag2.A) { // insert MNT2.2 into MNT.2
	if(mnt1->frag2.nvol<1 || mnt2->frag2.nvol<1) {std::cout<<"MergeTrees - array size error 4!"<<std::endl; return kFALSE;}
	mnt->frag2.nvol += mnt2->frag2.nvol;
	for(UInt_t iv=0; iv<mnt2->frag2.tV.size(); iv++) mnt->frag2.tV.push_back(mnt2->frag2.tV[iv]);
      } else  { std::cout<<"MergeTrees - Couldn't make fragment 2 association! Exiting..."<<std::endl; return kFALSE; }
    } else { if(printOut) std::cout<<"     Fragment 2 did not exit target..."<<std::endl; }

    OutTree->Fill();  
    mnt->Reset();
  }

  fout->Write();
  fout->Close();

  return kTRUE;
}
