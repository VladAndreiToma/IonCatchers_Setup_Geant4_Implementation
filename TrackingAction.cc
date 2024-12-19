#include "TrackingAction.hh"

using namespace std;

TrackingAction::TrackingAction(EventAction* eventAction, RunInput* rinp):G4UserTrackingAction(), fEventAction(eventAction), runInput(rinp) {}

void TrackingAction::PreUserTrackingAction(const G4Track* /*track*/) { }

void TrackingAction::PostUserTrackingAction(const G4Track* track) {
  if(!runInput->TrackNtuple()) return;
  
  G4double energy = track->GetVertexKineticEnergy();
  G4double theta = track->GetVertexMomentumDirection().theta();
  G4double phi = track->GetVertexMomentumDirection().phi();
  G4int trackID = track->GetTrackID();
  G4int parentID = track->GetParentID();
 
  const G4ParticleDefinition* particle = track->GetParticleDefinition();
  G4int        partZ = (G4int)particle->GetPDGCharge();
  G4int        partA = particle->GetAtomicMass();
  G4double partM = particle->GetPDGMass()/MeV;
  G4int    particleID = runInput->GetParticleID(track);
  G4int    processID = runInput->GetProcessID(track);
  
  fEventAction->AddTrk(trackID, particleID, parentID, processID);
  fEventAction->AddTZAM(partZ, partA, partM);
  fEventAction->AddTKin(energy, theta, phi);
}

/*
  G4Track: G4double -> GetLocalTime(), GetGlobalTime(), GetTotalEnergy(), GetKineticEnergy(), GetVelocity(), GetWeight(),
                       GetTrackLength(), GetStepLength(), GetVertexKineticEnergy().
           G4ThreeVector& -> GetPosition(), GetMomentumDirection(), GetMomentum(), GetPolarization(),
	                     GetVertexPosition(), GetVertexMomentumDirection().
	   G4int -> GetTrackID(), GetParentID(), GetCurrentStepNumber().
	   G4bool -> IsGoodForTracking(), IsBelowThreshold().
	   G4TouchableHandle& -> GetTouchableHandle(), GetNextTouchableHandle(), GetOriginTouchableHandle().
	   G4PhysicalVolume* -> GetVolume(), GetNextVolume().
	   G4LogicalVolume* -> GetLogicalVolumeAtvertex().
	   G4TrackStatus -> GetTrackStatus(). G4Material* -> GetMaterial().
	   G4VProcess* -> GetCreatorProcess(). G4ParticleDefinition* -> GetDefinition().
 */
