#include "PrimaryGeneratorAction.hh"

using namespace std;

PrimaryGeneratorAction::PrimaryGeneratorAction(RunInput* rinp, CLHEP::HepRandomEngine* rgen) : G4VUserPrimaryGeneratorAction(), fParticleGun(0), runInput(rinp), randGen(rgen), bZ(0), bA(0), bEvt(0), bTrk(0), bPar(0), bQ(0.*eplus), bx(0.*mm), by(0.*mm), bz(0.*mm), bEne(0.*MeV), bM(0.*MeV), bR(0.*mm), bThe(0.*rad), bPhi(0.*rad) {
  // default direction for beam ions:
  dx = 0.; dy = 0.; dz = 1.; //  positive z
  
  if(runInput->RunMode1()) {
    ionFile.open("./ReleasedIons.txt", ios::in);
    if(!ionFile) { G4cout<<"PrimaryGeneratorAction - Error: Can't open the output file."<<G4endl; exit(-1); }
    G4cout<<"PrimaryGeneratorAction: using ReleasedIons.txt as source..."<<G4endl;
  }

  enerRand = new G4RandGauss(randGen, runInput->GetBeamEMean(), runInput->GetBeamEWidth());
  gausRand = new G4RandGauss(randGen, 0., runInput->GetBeamPWidth());
  flatRand = new G4RandFlat(randGen);

  fParticleGun = new G4ParticleGun(1);    // with one particle
  G4ParticleDefinition* part = G4ChargedGeantino::ChargedGeantino();
  fParticleGun->SetParticleDefinition(part);  
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  // delete enerRand;
  delete fParticleGun;
  if(ionFile) ionFile.close();
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
  if(runInput->RunMode1()) {   // read reaction product ions from file
    ionFile>>bEvt>>bTrk>>bPar>>bZ>>bA>>bM>>bQ>>bEne>>bThe>>bPhi>>bx>>by>>bz;
    bM *= MeV; bQ *= eplus; bEne *= MeV; bThe *= rad; bPhi *= rad; bx *= mm;  by *= mm;  bz *= mm;
    runInput->SetBeamEvent(bEvt); runInput->SetBeamTrack(bTrk); runInput->SetBeamParent(bPar);
    dx = sin(bThe)*cos(bPhi);
    dy = sin(bThe)*sin(bPhi);
    dz = cos(bThe);
    //  G4cout<<"Generate released ion: Evt/Trk/Par="<<bEvt<<"/"<<bTrk<<"/"<<bPar<<", Z/A="<<bZ<<"/"<<bA<<G4endl;
  } else {   // get beam ions from distributions
    if(runInput->BeamName()=="Cf252") { // Cf252 spontaneous fission
      bZ = 98, bA = 252;
      bQ = 98.*eplus;
      bR = (runInput->GetSourceR())*sqrt(flatRand->fire()); bR *= mm;
      bPhi = CLHEP::twopi*flatRand->fire(); bPhi *= rad;
      bx = bR*cos(bPhi);
      by = bR*sin(bPhi);
      bz = runInput->GetSourceZ(); bz *= mm; // at interface between Pt foil and Au layer
      bEne = 0.*MeV; // at rest
    } else if(runInput->BeamName()=="U238") { // MNT inside INCREASE
      bZ = 92, bA = 238;
      bQ = 72.*eplus;
      bx = gausRand->fire(); bx *= mm;
      by = gausRand->fire(); by *= mm;
      bz = runInput->GetSourceZ(); bz *= mm;
      bEne = bA*enerRand->fire(); bEne*= MeV;
    } else if(runInput->BeamName()=="Xe136") { // MNT inside IGISOL
      bZ = 54, bA = 136;
      bQ = 31.*eplus;
      bR = (runInput->GetSourceR())*sqrt(flatRand->fire()); bR *= mm;
      bPhi = CLHEP::twopi*flatRand->fire(); bPhi *= rad;
      bx = bR*cos(bPhi);
      by = bR*sin(bPhi);
      bz = runInput->GetSourceZ(); bz *= mm; // 1mm behind target
      bEne = bA*enerRand->fire(); bEne*= MeV;
    } else { G4cout<<"PrimaryGeneratorAction - Error: undefined beam type."<<G4endl; exit(0); }
    // G4cout<<"Generate primary ion: E="<<bEne/bA/MeV<<"MeV/n, (x,y,z)=("<<bx/mm<<","<<by/mm<<","<<bz/mm<<")"<<G4endl;
  }
  
  if(bZ<1) {G4cout<<"PrimaryGeneratorAction - Error: wrong input!"<<G4endl; exit(-1);}
  G4ParticleDefinition* part = G4IonTable::GetIonTable()->GetIon(bZ, bA, 0.*keV); // ground state
  if(part) {
    fParticleGun->SetParticleDefinition(part);
    fParticleGun->SetParticleCharge(bQ);
  } else { G4cout<<"PrimaryGeneratorAction::GeneratePrimaries - ERROR: ion not defined!"<<G4endl; exit(0); }
  fParticleGun->SetParticlePosition(G4ThreeVector(bx,by,bz));
  fParticleGun->SetParticleEnergy(bEne);
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(dx,dy,dz));

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
