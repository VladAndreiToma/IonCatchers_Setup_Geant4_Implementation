#include "Randomize.hh"

#include <time.h>
#include <algorithm>
using namespace std;

#include "RunInput.hh"
#include "PhysicsList.hh"
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"

// #ifdef G4MULTITHREADED
// #include "G4MTRunManager.hh"
// #else
#include "G4RunManager.hh"
// #endif

#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "TROOT.h"
#include "TApplication.h"

// Author: Paul Constantin
int main(int argc, char** argv){
  TApplication *app = new TApplication("app", 0, 0); // needed for reading ROOT files
  RunInput* runInput = new RunInput();
  G4bool VisualOn = runInput->IsVisualOn();

  CLHEP::HepRandomEngine* randGen = new CLHEP::RanecuEngine;
  long int seed = static_cast<long int> (time(NULL));   // get seed from time
  G4Random::setTheEngine(randGen);   // choose the Random engine
  G4Random::setTheSeed(seed);

  // get number of events
  G4int numberOfEvent = runInput->GetNumEvents();
  if(runInput->RunMode1()) {
    std::ifstream ionFile;
    ionFile.open("./ReleasedIons.txt", ios::in);
    if(!ionFile) { G4cout<<"IonCatchers.cc - Error: Can't open the released ion file."<<G4endl; exit(-1); }
    numberOfEvent = std::count(std::istreambuf_iterator<char>(ionFile),  std::istreambuf_iterator<char>(), '\n');
    runInput->SetNumEvents(numberOfEvent);
  }
  
  // construct the default run manager (use G4MTRunManager for multi-threading)
// #ifdef G4MULTITHREADED
//   G4MTRunManager* runManager = new G4MTRunManager;
//   runManager->SetNumberOfThreads(8);
//   G4cout<<"Running in multi-threading mode!"<<G4endl;
// #else
  G4RunManager* runManager = new G4RunManager;
  G4cout<<"Running in sequential mode!"<<G4endl;
// #endif

  // Set mandatory initialization classes
  DetectorConstruction* setup = new DetectorConstruction(runInput);
  runManager->SetUserInitialization(setup);
  runManager->SetUserInitialization(new PhysicsList(runInput, randGen));
  runManager->SetUserInitialization(new ActionInitialization(runInput, randGen, setup));
  runManager->Initialize();                                    // initialize G4 kernel
  
  G4VisManager* visManager = new G4VisExecutive;
  if(VisualOn) {
    visManager->Initialize();
    // interactive mode : define UI session
    G4UIExecutive* UIexec = new G4UIExecutive(argc, argv);  // then "/control/execute vis.mac" and "/run/beamOn 1" 
    UIexec->SessionStart();
    delete UIexec;
  } else {
    // get the pointer to the UI manager and set verbosities
//     G4UImanager* UI = G4UImanager::GetUIpointer();
//     UI->ApplyCommand("/run/verbose 1");
//     UI->ApplyCommand("/event/verbose 1");
//     UI->ApplyCommand("/tracking/verbose 1");
//     UI->ApplyCommand("/physics/verbose 1");

    // start a run
    G4cout<<G4endl<<"***Running "<<numberOfEvent<<" events with random seed "<<seed<<G4endl<<G4endl;
    runManager->BeamOn(numberOfEvent);
  }
  
  // Free the store: user actions, physics list and detector description are owned and 
  // deleted by the run manager, so they should not be deleted in the main() program!
  delete app;
  delete runInput;
  delete visManager;
  delete runManager;  
  
  return 0;
}
