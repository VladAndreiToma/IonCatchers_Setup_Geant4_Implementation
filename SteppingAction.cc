#include "SteppingAction.hh"

using namespace std;

SteppingAction::SteppingAction(RunAction*  runAction, EventAction* eventAction, RunInput* rinp) : G4UserSteppingAction(), fRunAction(runAction), fEventAction(eventAction), runInput(rinp), fVolList(), iTest1p(0), iTest2p(0), iTest3p(0), iTest4p(0), iTest5p(0), iTest1a(0), iTest2a(0), iTest3a(0), iTest4a(0), iTest5a(0) {
  if(runInput->RunMode0()) {
    ionFile.open("./ReleasedIons.txt", ios::out);
    if(!ionFile) { G4cout<<"SteppingAction - Error: Can't open the output file."<<G4endl; exit(1); }
  }
}

SteppingAction::~SteppingAction() {
  for(unsigned int i=0; i<fVolList.size(); i++) G4cout<<"Volume #"<<i<<" is "<<fVolList[i]<<G4endl; G4cout<<G4endl;
  G4cout<<"IONS IN TEST LAYER 1: primary="<<iTest1p<<", secondary="<<iTest1a<<G4endl;
  G4cout<<"IONS IN TEST LAYER 2: primary="<<iTest2p<<", secondary="<<iTest2a<<G4endl;
  G4cout<<"IONS IN TEST LAYER 3: primary="<<iTest3p<<", secondary="<<iTest3a<<G4endl;
  G4cout<<"IONS IN TEST LAYER 4: primary="<<iTest4p<<", secondary="<<iTest4a<<G4endl;
  G4cout<<"IONS IN TEST LAYER 5: primary="<<iTest5p<<", secondary="<<iTest5a<<G4endl;
  if(ionFile) ionFile.close();
}

void SteppingAction::UserSteppingAction(const G4Step* step) {
  G4StepPoint* prePoint = dynamic_cast<G4StepPoint*>(step->GetPreStepPoint());
  if(!prePoint) { G4cout<<"No point associated to step - Exiting!"<<G4endl; exit(0); }
  G4ThreeVector posPt = prePoint->GetPosition();                // prePoint position
  // G4StepPoint* postPoint = dynamic_cast<G4StepPoint*>(step->GetPostStepPoint());
  // if(!postPoint) { G4cout<<"No point associated to step - Exiting!"<<G4endl; exit(0); }
  G4Track* track = dynamic_cast<G4Track*>(step->GetTrack());
  if(!track) { G4cout<<"No track associated to step - Exiting!"<<G4endl; exit(0); }
  //  G4ThreeVector momVtx = track->GetVertexMomentumDirection();   // track vertex direction
  G4ParticleDefinition* particle = track->GetDefinition();      // particle PDG data
  G4VPhysicalVolume* stepVol = dynamic_cast<G4VPhysicalVolume*>(prePoint->GetTouchableHandle()->GetVolume()); // current volume
  if(!stepVol) { G4cout<<"No volume associated to step - Exiting!"<<G4endl; exit(0); }
  G4String volName = stepVol->GetName();
  
  G4int iDet = -1;
  fVolList = fRunAction->GetVolumeList(); // cannot be in constructor
  for(unsigned int i=0; i<fVolList.size(); i++)
    if(volName==fVolList[i]) {
      iDet = i;
      break;
    }
  if(iDet<0) return;   // do not store steps in unregistered volumes (see RunAction::BeginOfRunAction)

  G4double eDep = step->GetTotalEnergyDeposit()/MeV;
  G4double stepL = step->GetStepLength()/mm;
  G4int        parentID = track->GetParentID();
  G4int        trackID = track->GetTrackID();
  G4int        stepID = track->GetCurrentStepNumber();
  //  G4double TrkE = track->GetVertexKineticEnergy()/MeV;
  //  G4double TrkT = momVtx.theta()/rad;
  G4double StpE = prePoint->GetKineticEnergy()/MeV;
  G4double StpT = prePoint->GetMomentumDirection().theta()/rad;
  G4double StpP = prePoint->GetMomentumDirection().phi()/rad;
  G4double StpX = posPt.getX()/mm; G4double StpY = posPt.getY()/mm; G4double StpZ = posPt.getZ()/mm;
  G4double preQ = prePoint->GetCharge();
  G4int        partZ = (G4int)particle->GetPDGCharge();
  G4int        partA = particle->GetAtomicMass();
  G4double partM = particle->GetPDGMass()/MeV;
  //  G4int    preProc = runInput->GetProcessID(prePoint);
  //  G4int    postProc = runInput->GetProcessID(postPoint);
  //  G4int    trackProc = runInput->GetProcessID(track);
  G4int        partID = runInput->GetParticleID(track);
  
  // ION STATS (see RunAction::BeginOfRunAction):
  if(iDet==3) { // 1st test layer
     if(parentID==0 && partID==6) iTest1p++;
     if(partID>6) iTest1a++;
  }
  if(iDet==4) { // 2nd test layer
     if(parentID==0 && partID==6) iTest2p++;
     if(partID>6) iTest2a++;
  }
  if(iDet==5) { // 3rd test layer
     if(parentID==0 && partID==6) iTest3p++;
     if(partID>6) iTest3a++;
  }
  if(iDet==6 && !runInput->IgisolMNT()) { // 4th test layer
     if(parentID==0 && partID==6) iTest4p++;
     if(partID>6) iTest4a++;
  }
  if(iDet==7 && !runInput->IgisolMNT()) { // 5th test layer
     if(parentID==0 && partID==6) iTest5p++;
     if(partID>6) iTest5a++;
  }

  if(partID>5) { // see EventAction::EndOfEventAction
    fEventAction->AddPart(partID, trackID, parentID);
    fEventAction->AddSZAM(partZ, partA, partM);
    fEventAction->AddSKin(StpE, StpT, StpP);
    fEventAction->AddLabel(stepID, iDet);
    fEventAction->AddPos(StpX, StpY, StpZ);
    fEventAction->AddELQ(eDep, stepL, preQ);
  }

  if(runInput->RunMode0() && iDet==5 && partID>6)   // store in a file the ions in test layer after target (MNT) or after Ti foil (Cf)
    ionFile<<fEventAction->GetEvt()<<" "<<trackID<<" "<<parentID<<" "<<partZ<<" "<<partA<<" "<<partM<<" "<<preQ<<" "<<StpE<<" "<<StpT<<" "<<StpP<<" "<<StpX<<" "<<StpY<<" "<<StpZ<<endl;
}
      
void SteppingAction::PrintSecondaries(const G4Step* step){
  const G4TrackVector* secList = step->GetSecondary();          // vector of secondaries
  G4int nSec = secList->size();
  G4int nGamma=0, nElect=0, nPosit=0, nNeutr=0, nProto=0, nAlpha=0;
  if(nSec>0) {
    for(G4int iSec=0; iSec<nSec; iSec++) {
      const G4Track* secTk = (*secList)[iSec];
      G4int secID = runInput->GetParticleID(secTk);
      if(secID==0) nGamma++;
      else if (secID==1) nElect++;
      else if (secID==2) nPosit++;
      else if (secID==3) nNeutr++;
      else if (secID==4) nProto++;
      else if (secID==5) nAlpha++;
      else if(secID>5) {
	G4ParticleDefinition* particle = const_cast<G4ParticleDefinition*>(secTk->GetDefinition());
	if(!particle) { G4cout<<"ERROR - SteppingAction::PrintSecondaries : no particle associated!"<<G4endl; exit(-1); }
	G4cout<<", ion="<<particle->GetParticleName();
      } else { G4cout<<"ERROR - SteppingAction::PrintSecondaries: unidentified secondary "<<secID<<G4endl; exit(-1); }
    }
  } else G4cout<<", no secondaries!";
  if(nGamma>0) G4cout<<", gamma="<<nGamma;
  if(nElect>0) G4cout<<", elect="<<nElect;
  if(nPosit>0) G4cout<<", posit="<<nPosit;
  if(nNeutr>0) G4cout<<", neutr="<<nNeutr;
  if(nProto>0) G4cout<<", proto="<<nProto;
  if(nAlpha>0) G4cout<<", alpha="<<nAlpha;
  G4cout<<G4endl;
}

/*
  Classes selected public methods:
  G4Step: G4StepPoint* -> GetPreStepPoint(), GetPostStepPoint().
          G4Track* -> GetTrack().
	  G4TrackVector* -> GetSecondary(), GetfSecondary().
	  G4double -> GetStepLength(), GetTotalEnergyDeposit(), GetDeltaTime(), GetDeltaEnergy().
	  G4ThreeVector -> GetDeltaPosition(), GetDeltaMomentum().
	  G4bool -> IsFirstStepInVolume(), IsLastStepInVolume().
  G4StepPoint: G4double -> GetLocalTime(), GetGlobalTime(), GetTotalEnergy(), GetKineticEnergy(), GetVelocity(), GetBeta(),
                           GetGamma(), GetMass(), GetCharge(), GetWeight(), GetMagneticMoment().
	       G4ThreeVector& -> GetPosition(), GetMomentumDirection(), GetMomentum(), GetPolarization().
	       G4VPhysicalVolume* -> GetPhysicalVolume().
	       G4VSensitiveDetector* -> GetSensitiveDetector().
	       G4TouchableHandle& -> GetTouchableHandle().
	       G4Material* -> GetMaterial().
	       G4VProcess* -> GetProcessDefinedStep().
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
  G4VTouchable: G4VPhysicalVolume* -> GetVolume(). G4VSolid* -> GetSolid(). G4int -> GetCopyNumber().
  G4VPhysicalVolume: G4String& -> GetName(). G4int -> GetCopyNo(), GetMultiplicity().
  G4TrackVector: typedef std::vector<G4Track*> G4TrackVector.
*/
