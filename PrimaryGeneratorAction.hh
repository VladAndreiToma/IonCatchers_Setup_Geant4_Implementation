#ifndef PrimaryGeneratorAction_h
#define PrimaryGeneratorAction_h 1

#include <iostream>
#include <fstream>

#include "RunInput.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ChargedGeantino.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "Randomize.hh"
#include "globals.hh"

// Author: Paul Constantin
// Generates input (primary) beam (radiation)
class PrimaryGeneratorAction: public G4VUserPrimaryGeneratorAction {
public:
  PrimaryGeneratorAction(RunInput*, CLHEP::HepRandomEngine*);
  virtual ~PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event*);
  G4ParticleGun* GetParticleGun() {return fParticleGun;}

private:
  G4ParticleGun* fParticleGun;
  RunInput* runInput;
  CLHEP::HepRandomEngine* randGen;
  G4int bZ, bA, bEvt, bTrk, bPar;
  G4double bQ, bx, by, bz, bEne, bM, bR, bThe, bPhi;
  G4double dx, dy, dz;
  std::ifstream ionFile;
  G4RandGauss* enerRand;
  G4RandGauss* gausRand;
  G4RandFlat* flatRand;
};

#endif
