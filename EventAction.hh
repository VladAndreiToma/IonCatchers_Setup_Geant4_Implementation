#ifndef EventAction_h
#define EventAction_h 1

#include "G4UserEventAction.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "globals.hh"
#include "Randomize.hh"

#include "RunInput.hh"
#include "RunAction.hh"

#include <vector>

// Author: Paul Constantin
// Event action class. Defines data members to hold the parameters to be stored.
// They are collected every step/track via the functions AddXYZ().
class EventAction : public G4UserEventAction {
public:
  EventAction(RunInput*, RunAction*);
  virtual ~EventAction();
  
  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);

  G4int GetEvt() { return eventID; }
  // step fillers:
  void AddPart(G4int v1, G4int v2, G4int v3) { sPrt.push_back(v1); sTrk.push_back(v2); sPre.push_back(v3); }
  void AddSZAM(G4int v1, G4int v2, G4double v3) { sZ.push_back(v1); sA.push_back(v2); sM.push_back(v3); }
  void AddSKin(G4double v1, G4double v2, G4double v3) { sE.push_back(v1); sT.push_back(v2); sP.push_back(v3); }
  void AddLabel(G4int v1, G4int v2) { sStp.push_back(v1); sDet.push_back(v2); }
  void AddPos(G4double v1, G4double v2, G4double v3) { sx.push_back(v1); sy.push_back(v2); sz.push_back(v3); }
  void AddELQ(G4double v1, G4double v2, G4double v3) { sEDep.push_back(v1); sSLen.push_back(v2); sQ.push_back(v3);}
  // track fillers:
  void AddTrk(G4int v1, G4int v2, G4int v3, G4int v4) { tTrk.push_back(v1); tPrt.push_back(v2); tPre.push_back(v3); tPrc.push_back(v4); }
  void AddTZAM(G4int v1, G4int v2, G4double v3) { tZ.push_back(v1); tA.push_back(v2); tM.push_back(v3); }
  void AddTKin(G4double v1, G4double v2, G4double v3) { tE.push_back(v1); tT.push_back(v2); tP.push_back(v3); }
  
private:
  RunInput* runInput;
  RunAction* fRunAction;
  G4int eventID;
  // step variables:
  std::vector<G4int>    sTrk;
  std::vector<G4int>    sPrt;
  std::vector<G4int>    sPre;
  std::vector<G4int>    sZ;
  std::vector<G4int>    sA;
  std::vector<G4double> sM;
  std::vector<G4double> sE;
  std::vector<G4double> sT;
  std::vector<G4double> sP;
  std::vector<G4int>    sStp;
  std::vector<G4double> sx;
  std::vector<G4double> sy;
  std::vector<G4double> sz;
  std::vector<G4int>    sDet;
  std::vector<G4double> sEDep;
  std::vector<G4double> sSLen;
  std::vector<G4double> sQ;
  // track variables:
  std::vector<G4int>    tTrk;
  std::vector<G4int>    tPrt;
  std::vector<G4int>    tPre;
  std::vector<G4int>    tPrc;
  std::vector<G4int>    tZ;
  std::vector<G4int>    tA;
  std::vector<G4double> tM;
  std::vector<G4double> tE;
  std::vector<G4double> tT;
  std::vector<G4double> tP;
 
  G4AnalysisManager* analysisManager;
  void FillStpNtuple(unsigned int, G4int, G4int, G4int);
  void FillTrkNtuple(unsigned int, G4int, G4int, G4int);
  void PrintStp();
  void PrintTrk();
};

#endif
    
