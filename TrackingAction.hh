#ifndef TrackingAction_h
#define TrackingAction_h 1

#include "G4UserTrackingAction.hh"
#include "EventAction.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"
#include "G4UnitsTable.hh"

// Author: Paul Constantin
// Collects the track parameters to be outputted by EventAction and RunAction
// G4UserStackingAction controls track stacking: place in urgent or waiting stacks, postpone to next event or kill them.
class TrackingAction : public G4UserTrackingAction {
public:  
  TrackingAction(EventAction*, RunInput*);
  virtual ~TrackingAction() {};
  
  virtual void PreUserTrackingAction(const G4Track*);   
  virtual void PostUserTrackingAction(const G4Track*);
  
private:
  EventAction* fEventAction;
  RunInput* runInput;
};

#endif
