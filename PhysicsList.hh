#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "Randomize.hh"
#include "globals.hh"

#include "RunInput.hh"

// Author: Paul Constantin
// Defines the particles and the physics processes to be simulated.
class PhysicsList : public G4VModularPhysicsList {
public:
  PhysicsList(RunInput*, CLHEP::HepRandomEngine*);
  virtual ~PhysicsList();
  virtual void ConstructParticle();   // construction of particles
  virtual void ConstructProcess();    // construct processes and register them to particles
  virtual void SetCuts();
  
private:  
  RunInput* runInput;
  CLHEP::HepRandomEngine* randGen;
};

#endif

