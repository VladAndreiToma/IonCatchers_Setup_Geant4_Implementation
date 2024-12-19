#ifndef SteppingAction_h
#define SteppingAction_h 1

#include <iostream>
#include <fstream>

#include "G4UserSteppingAction.hh"
#include "G4StepPoint.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"

#include "RunInput.hh"
#include "RunAction.hh"
#include "EventAction.hh"

// Author: Paul Constantin
// Collects the step parameters to be outputted by EventAction and RunAction
class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(RunAction*,  EventAction*, RunInput*);
  virtual ~SteppingAction();

  virtual void UserSteppingAction(const G4Step*);
    
private:
  RunAction*  fRunAction;
  EventAction*  fEventAction;
  RunInput* runInput;
  std::vector<G4String> fVolList;
  G4int iTest1p, iTest2p, iTest3p, iTest4p, iTest5p;
  G4int iTest1a, iTest2a, iTest3a, iTest4a, iTest5a;
  std::ofstream ionFile;

  void PrintSecondaries(const G4Step*);
};

#endif
