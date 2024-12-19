#include "ActionInitialization.hh"

ActionInitialization::ActionInitialization(RunInput* rinp, CLHEP::HepRandomEngine* rgen, DetectorConstruction* dc) : G4VUserActionInitialization(), runInput(rinp), randGen(rgen), setup(dc) {}

ActionInitialization::~ActionInitialization() {}

void ActionInitialization::Build() const {
  SetUserAction(new PrimaryGeneratorAction(runInput, randGen));
  RunAction* runAction = new RunAction(runInput, setup);
  SetUserAction(runAction);

  EventAction* eventAction = new EventAction(runInput, runAction);
  SetUserAction(eventAction);
  SetUserAction(new SteppingAction(runAction, eventAction, runInput));
  if(runInput->TrackNtuple())
    SetUserAction(new TrackingAction(eventAction, runInput));
}
