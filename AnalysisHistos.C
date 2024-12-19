// Loops over the tree in Tree.root and makes plots

#include <TFile.h>
#include <TString.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH2F.h>
#include <TMath.h>
#include <TUnixSystem.h>

#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <ctime>

using namespace std;

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

// MAIN MACRO
// Config: 0=INCREASE+Cf252, 1=IGISOL+MNT  
Bool_t AnalysisHistos(TString fileName = "5MXePt_Tree.root", Int_t Config = 1){
  Int_t nBinsT=0, nBinsE=0, nBinsR=0, nBinsZ=0;
  Double_t maxT=0., maxE=0., maxR=0., maxZ=0.;
  Bool_t inGas = kFALSE;
  if(Config==0) { // Settings for INCREASE Cf252
  nBinsT=70; nBinsE=37; nBinsR=15; nBinsZ=50;
  maxT=140.; maxE=0.74; maxR=150.; maxZ=500.; // maxR&maxZ also define effective stopping volume in gas
  } else if(Config==1) { // Settings for IGISOL MNT
  nBinsT=90; nBinsE=30; nBinsR=50; nBinsZ=35;
  maxT=90.; maxE=1.5; maxR=100.; maxZ=70.;       // maxR&maxZ also define effective stopping volume in gas
  } else { std::cout<<"Configuration not implemented!"<<std::endl; return kFALSE; }
  
  // GET INPUT NTUPLE
  TFile* fin = new TFile(fileName,"READ");
  if(!fin) { std::cout<<"Couldn't get the input file"<<std::endl; return kFALSE; }
  TTree* InTree = (TTree*)fin->Get("MNTTree");
  if(!InTree) { std::cout<<"Couldn't get the tree! Exiting..."<<std::endl; return kFALSE; }
  Int_t nEvts = InTree->GetEntries();
  if(nEvts<1) { std::cout<<"Tree is empty! Exiting..."<<std::endl; return kFALSE; }
  std::cout<<"Looping over "<<nEvts<<" events"<<std::endl;
  Event* mnt = new Event;
  mnt->Reset();
  InTree->SetBranchAddress("e",&mnt);

  // SETUP OUTPUT HISTOGRAM FILE:
  TFile* fout = new TFile("Histos.root","RECREATE");
  if(!fout) { std::cout<<"Couldn't create the output file"<<std::endl; return kFALSE; }
  TH2F* hExTreleased = new TH2F("hExTreleased", "Energy vs. Angle released ions", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hExTstopped = new TH2F("hExTstopped", "Energy vs. Angle stopped ions", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hAxZstopped = new TH2F("hAxZstopped", "A vs Z stopped ions", 40,28.,68., 98,74.,172.);
  TH2F* hrxzstopped = new TH2F("hrxzstopped", "r vs z stopped ions", nBinsZ, 0., maxZ, nBinsR, 0., maxR);
  TH2F* hExTescaped = new TH2F("hExTescaped", "Energy vs. Angle hit escaped", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hExTAlpha = new TH2F("hExTAlpha", "Energy vs. Angle hit alpha sources", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hExTdumped = new TH2F("hExTdumped", "Energy vs. Angle hit beam dump", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hExTTiFrame = new TH2F("hExTTiFrame", "Energy vs. Angle hit Ti frame", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hExTMiniCage = new TH2F("hExTMiniCage", "Energy vs. Angle hit mini cages", nBinsT, 0., maxT, nBinsE, 0., maxE);
  TH2F* hExTTargHold = new TH2F("hExTTargHold", "Energy vs. Angle hit target wheel", nBinsT, 0., maxT, nBinsE, 0., maxE);
  
  // LOOP OVER INPUT TREE
  Float_t Rad2Deg = 180./TMath::Pi();
  div_t DivByMil;
  for(Int_t iRow=0; iRow<nEvts; iRow++){
    DivByMil = std::div(iRow,100);
    if(DivByMil.rem == 0)  std::cout<<"Processed "<<iRow<<" events"<<std::endl;
    InTree->GetEntry(iRow);
    
    // loop over fragment 1 volumes:
    Float_t ein = -1., tin = -1., Rout = -1., Zout = -1.;
    Bool_t inCell = kFALSE;
    Bool_t isReleased = kFALSE, isStopped = kFALSE, isEscaped = kFALSE;
    Bool_t isAlpha = kFALSE, isDumped = kFALSE, isTiFrame = kFALSE, isMiniCage = kFALSE, isTargHold = kFALSE;
    Float_t AF1 = mnt->frag1.A;
    Float_t ZF1 = mnt->frag1.Z;
    UInt_t nF1 = mnt->frag1.nvol;
    for(UInt_t iF1=0; iF1<nF1; iF1++) {
      Int_t Vol = mnt->frag1.tV[iF1].vol;
      if(Vol==5) { // layer after target
	isReleased = kTRUE;
	ein = mnt->frag1.tV[iF1].ein/AF1;
	tin = mnt->frag1.tV[iF1].tin*Rad2Deg;
      }
      if(Vol==1) { // gas cell
	Float_t xout = mnt->frag1.tV[iF1].xout;
	Float_t yout = mnt->frag1.tV[iF1].yout;
	Rout = sqrt(xout*xout+yout*yout);
	Zout = mnt->frag1.tV[iF1].zout;
      }
    }
    Int_t lastVol1 = mnt->frag1.tV[nF1-1].vol;
    if(Rout>0. && Zout>0.) inCell = (Rout<maxT && Zout<maxR); // defines effective stopping volume in gas
    if(lastVol1==1 && inCell) isStopped = kTRUE;
    if(lastVol1==0 || (lastVol1==1 && !inCell)) isEscaped = kTRUE;
    if(lastVol1>=11 && lastVol1<=12) isAlpha = kTRUE; // last volume is in the alpha sources
    if(lastVol1>=13 && lastVol1<=21) isDumped = kTRUE; // last volume is in the beam dump
    if(lastVol1>=24 && lastVol1<=25) isTiFrame = kTRUE; // last volume is in the Ti foil frame
    if(lastVol1>=26 && lastVol1<=28) isMiniCage = kTRUE; // last volume is in the Ti foil frame
    if(lastVol1>=29 && lastVol1<=31) isTargHold = kTRUE; // last volume is in the Ti foil frame
    if(ein>0 && tin>0 && isReleased) hExTreleased->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isStopped) {
      hExTstopped->Fill(tin,ein);
      hAxZstopped->Fill(ZF1,AF1);
      if(Rout>0. && Zout>0.) hrxzstopped->Fill(Zout, Rout);
    }
    if(ein>0 && tin>0 && isReleased && isEscaped) hExTescaped->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isAlpha) hExTAlpha->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isDumped) hExTdumped->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isTiFrame) hExTTiFrame->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isMiniCage) hExTMiniCage->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isTargHold) hExTTargHold->Fill(tin,ein);
    
    // loop over fragment 2 volumes:
    ein = -1.; tin = -1.; Rout = -1.; Zout = -1.;
    inCell = kFALSE;
    isReleased = kFALSE; isStopped = kFALSE; isEscaped = kFALSE;
    isAlpha = kFALSE; isDumped = kFALSE; isTiFrame = kFALSE; isMiniCage = kFALSE; isTargHold = kFALSE;
    Float_t AF2 = mnt->frag2.A;
    Float_t ZF2 = mnt->frag2.Z;
    UInt_t nF2 = mnt->frag2.nvol;
    for(UInt_t iF2=0; iF2<nF2; iF2++) {
      Int_t Vol = mnt->frag2.tV[iF2].vol;
      if(Vol==5) { // fragment 2 in layer after target
	isReleased = kTRUE;
	ein = mnt->frag2.tV[iF2].ein/AF2;
	tin = mnt->frag2.tV[iF2].tin*Rad2Deg;
      }
      if(Vol==1) { // gas cell
	Float_t xout = mnt->frag2.tV[iF2].xout;
	Float_t yout = mnt->frag2.tV[iF2].yout;
	Rout = sqrt(xout*xout+yout*yout);
	Zout = mnt->frag2.tV[iF2].zout;
      }
   }
    Int_t lastVol2 = mnt->frag2.tV[nF2-1].vol;
    if(Rout>0. && Zout>0.) inCell = (Rout<maxT && Zout<maxR); // defines effective stopping volume in gas
    if(lastVol2==1 && inCell) isStopped = kTRUE;
    if(lastVol2==0 || (lastVol2==1 && !inCell)) isEscaped = kTRUE;
    if(lastVol2>=11 && lastVol2<=12) isAlpha = kTRUE; // last volume is in the alpha sources
    if(lastVol2>=13 && lastVol2<=21) isDumped = kTRUE; // last volume is in the beam dump
    if(lastVol2>=24 && lastVol2<=25) isTiFrame = kTRUE; // last volume is in the Ti foil frame
    if(lastVol2>=26 && lastVol2<=28) isMiniCage = kTRUE; // last volume is in the Ti foil frame
    if(lastVol2>=29 && lastVol2<=31) isTargHold = kTRUE; // last volume is in the Ti foil frame
    if(ein>0 && tin>0 && isReleased) hExTreleased->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isStopped) {
      hExTstopped->Fill(tin,ein);
      hAxZstopped->Fill(ZF2,AF2);
      if(Rout>0. && Zout>0.) hrxzstopped->Fill(Zout, Rout);
    }
    if(ein>0 && tin>0 && isReleased && isEscaped) hExTescaped->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isAlpha) hExTAlpha->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isDumped) hExTdumped->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isTiFrame) hExTTiFrame->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isMiniCage) hExTMiniCage->Fill(tin,ein);
    if(ein>0 && tin>0 && isReleased && isTargHold) hExTTargHold->Fill(tin,ein);
  } // loop over tree events

  fout->Write();
  fout->Close();

  return kTRUE;
}
