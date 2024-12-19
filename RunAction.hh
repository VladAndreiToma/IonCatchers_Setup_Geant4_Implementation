#ifndef RunAction_h
#define RunAction_h 1

#include <vector>
#include <deque>

#include "RunInput.hh"
#include "DetectorConstruction.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "G4Material.hh"
#include "G4HadronicProcessStore.hh"
#include "G4ProcessTable.hh"
#include "G4IonTable.hh"
#include "G4NistManager.hh"
#include "G4EmCalculator.hh"
#include "G4StableIsotopes.hh"

#include "Analysis.hh"

// Author: Paul Constantin
// Run action class. Histograms and ntuples are created in BeginOfRunAction().
// They are saved in the output file in a format selected in Analysis.hh.
class RunAction : public G4UserRunAction {
  public:
    RunAction(RunInput*, DetectorConstruction*);
    virtual ~RunAction();

    virtual void BeginOfRunAction(const G4Run*);
    virtual void EndOfRunAction(const G4Run*);

    G4int GetStpNtuple() {return stpNtu;}
    G4int GetTrkNtuple() {return trkNtu;}
    std::vector<G4String> GetVolumeList() { return fVolList; }

  private:
    RunInput* runInput;
    DetectorConstruction* setup;
    std::vector<G4String> fVolList; // list of geometry volumes
    G4int stpNtu, trkNtu;
    G4AnalysisManager* analysisManager;

    void PrintXSections();
    void PrintEMInfo();
};

#endif

