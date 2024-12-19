#ifndef MNTDataModel_h
#define MNTDataModel_h 1

#include "RunInput.hh"

#include "globals.hh"
#include "G4HadronicInteraction.hh"
#include "G4CascadeInterface.hh"
#include "G4HadProjectile.hh"
#include "G4Nucleus.hh"
#include "G4DynamicParticle.hh"
#include "G4PhysicsTable.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"
#include "Randomize.hh"

#include <cmath>
#include <cstdlib>
#include <vector>

#include "TFile.h"
#include "TTree.h"
#include "TH2D.h"

// Author: Paul Constantin
// Multi-nucleon transfer process with input data from Alexander Karpov theoretical model 
class MNTDataModel : public G4HadronicInteraction {
public:
  MNTDataModel(RunInput*, CLHEP::HepRandomEngine*);
  virtual ~MNTDataModel();
 
  virtual G4HadFinalState* ApplyYourself(const G4HadProjectile& aTrack, G4Nucleus& targetNucleus);

private:
  MNTDataModel & operator=(const MNTDataModel &right);
  MNTDataModel(const MNTDataModel&);
  
  G4bool IsValidProjectile(G4int Z, G4int A);
  G4bool IsValidTarget(G4int Z);
  void GetKinHist(G4int, G4int, G4int, G4int);
  G4int NearestEnergy(G4double);

  RunInput* runInput;
  CLHEP::HepRandomEngine* randGen;
  G4RandFlat*  flatRand;
  G4IonTable*    theTableOfIons;

  TFile* inFile;
  TTree* treeKin;
  TH2D*  histKin;

  G4int fZ, fA, fN;
  G4float fE, fS;
  G4float tKin[20000];   // adjust the size to that in the input tree
  G4float tThe[20000];
  G4float tSig[20000];
  G4double maxKin;
};

#endif
