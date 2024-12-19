#ifndef SpontaneousFission_h
#define SpontaneousFission_h 1

#include <vector>
#include <algorithm>
#include <CLHEP/Units/SystemOfUnits.h>

#include "G4ios.hh"
#include "globals.hh"
#include "G4VRestDiscreteProcess.hh"
#include "G4ParticleChangeForDecay.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4DynamicParticle.hh"
#include "G4DecayProducts.hh"
#include "G4Ions.hh"
#include "G4IonTable.hh"
#include "Randomize.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4HadronicProcessType.hh"

#include "TFile.h"
#include "TH2D.h"

// Author: Paul Constantin
// Generates spontaneous fission fragments
class SpontaneousFission : public G4VRestDiscreteProcess {
  public:
    SpontaneousFission(CLHEP::HepRandomEngine* randGen, const G4String& processName="SpontaneousFission");
    ~SpontaneousFission();

    virtual void ProcessDescription(std::ostream& outFile) const;
    G4bool IsApplicable(const G4ParticleDefinition&);
    void SelectAVolume(const G4String aVolume);
    void SelectAllVolumes();
    void BuildPhysicsTable(const G4ParticleDefinition &);
    G4VParticleChange* Apply(const G4Track& theTrack, const G4Step&  theStep);

  protected:
    G4double GetMeanFreePath(const G4Track& theTrack, G4double previousStepSize, G4ForceCondition* condition);
    G4double GetMeanLifeTime(const G4Track& theTrack, G4ForceCondition*);
  
  private:
    SpontaneousFission(const SpontaneousFission &right, CLHEP::HepRandomEngine*);
    SpontaneousFission & operator=(const SpontaneousFission &right);

    CLHEP::HepRandomEngine* randGen;
    G4bool isAllVolumesMode;
    G4bool isInitialised;
    std::vector<G4String> ValidVolumes;
    G4double fRemainderLifeTime;
    // ParticleChange for decay process
    G4ParticleChangeForDecay fParticleChangeForDecay;
    G4IonTable*  theTableOfIons;
    G4RandFlat*  flatRand;
    G4RandGauss* gaussRand;
    TFile* CfYieldFile;

    G4bool IsValid(const G4ParticleDefinition&);
    G4VParticleChange* Clear();
    G4double GetNeutMult(G4double);
    G4double GetMeanTKE(G4double);
    G4double GetSigTKE(G4double);

    // inline implementations 
    inline G4double AtRestGetPhysicalInteractionLength(const G4Track& track, G4ForceCondition* condition) {
      fRemainderLifeTime = G4VRestDiscreteProcess::AtRestGetPhysicalInteractionLength(track, condition);
      return fRemainderLifeTime;
    }
    inline G4VParticleChange* AtRestDoIt(const G4Track& theTrack, const G4Step& theStep) {
      return Apply(theTrack, theStep);
    }
    inline G4VParticleChange* PostStepDoIt(const G4Track& theTrack, const G4Step& theStep) {
      return Apply(theTrack, theStep);
    }
};

#endif

