#ifndef ActionInitialization_h
#define ActionInitialization_h 1

#include "G4VUserActionInitialization.hh"
#include "RunInput.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"

// Author: Paul Constantin
// Action initialization class: the user defines here all the user action classes
class ActionInitialization : public G4VUserActionInitialization {
public:
  ActionInitialization(RunInput*, CLHEP::HepRandomEngine*, DetectorConstruction*);
  virtual ~ActionInitialization();
  
  //     virtual void BuildForMaster() const; // called by the master thread
  virtual void Build() const; // called by worker threads
  
private:
  RunInput* runInput;
  CLHEP::HepRandomEngine* randGen;
  DetectorConstruction* setup;
};

#endif

    
