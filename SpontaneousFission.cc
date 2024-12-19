#include "SpontaneousFission.hh"

using namespace CLHEP;

SpontaneousFission::SpontaneousFission(CLHEP::HepRandomEngine* rgen, const G4String& processName) : G4VRestDiscreteProcess(processName, fDecay), randGen(rgen), isAllVolumesMode(true), isInitialised(false) {
  SetProcessSubType(fRadioactiveDecay);
  pParticleChange = &fParticleChangeForDecay;
  SelectAllVolumes();  // applies to all logical volumes
  // SelectAVolume("TargVolume"); // applies only in target
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
  flatRand = new G4RandFlat(randGen);
  gaussRand = new G4RandGauss(randGen);

  G4String BasePath = std::getenv("IonCatchersPath");
  G4String pathToFile = BasePath; pathToFile += "/Cf252FissYield.root";
  CfYieldFile = new TFile(pathToFile,"READ");
  if(!CfYieldFile) { G4cout<<"SpontaneousFission - cannot read input yield file!"<<G4endl; exit(-1); }
}

void SpontaneousFission::ProcessDescription(std::ostream& outFile) const {
  outFile<<"Spontaneous fission of Cf252 at rest inside the target/source.\n";
}

SpontaneousFission::~SpontaneousFission() {}

G4bool SpontaneousFission::IsApplicable(const G4ParticleDefinition& aParticle) {
  G4String pName = aParticle.GetParticleName();
  if(pName != "GenericIon") {
    G4cout<<"SpontaneousFission is not implemented for "<<pName<<G4endl;
    return false;
  }
  return true;
}

G4bool SpontaneousFission::IsValid(const G4ParticleDefinition& aParticle) {
  G4int pZ = aParticle.GetAtomicNumber();
  G4int pA = aParticle.GetAtomicMass();
  //  G4String targName = runInput->TargName();
  //  if(Z==28 && targName.contains("Cf")) return true;
  if(pZ==98 && pA==252) return true;
  return false;
}

void SpontaneousFission::SelectAVolume(const G4String aVolume) {
  G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  for (size_t i=0; i<theLogicalVolumes->size(); i++) {
    G4LogicalVolume* volume=(*theLogicalVolumes)[i];
    if (volume->GetName() == aVolume) {
      ValidVolumes.push_back(aVolume);
      std::sort(ValidVolumes.begin(), ValidVolumes.end()); // sort need for performing binary_search
    } else if(i == theLogicalVolumes->size()) {
      G4cout<<"SelectAVolume: "<<aVolume<<" is not a valid logical volume name"<<G4endl;
      exit(0);
    }
  }
}

void SpontaneousFission::SelectAllVolumes()  {
  G4LogicalVolumeStore* theLogicalVolumes = G4LogicalVolumeStore::GetInstance();
  ValidVolumes.clear();
  for (size_t i=0; i<theLogicalVolumes->size(); i++){
    G4LogicalVolume* volume = (*theLogicalVolumes)[i];
    ValidVolumes.push_back(volume->GetName());    
  }
  std::sort(ValidVolumes.begin(), ValidVolumes.end()); // sort needed in order to allow binary_search
  isAllVolumesMode=true;
}

// override G4VRestDiscreteProcess
G4double SpontaneousFission::GetMeanLifeTime(const G4Track& trk, G4ForceCondition*) {
  const G4ParticleDefinition* part = trk.GetDynamicParticle()->GetDefinition();
  G4int pZ = part->GetAtomicNumber();
  G4int pA = part->GetAtomicMass();
  if(pZ==98 && pA==252) return DBL_MIN; // force Cf252 to decay immediately
  if(part->GetPDGStable()) return DBL_MAX;
  G4double lifeTime = part->GetPDGLifeTime();
  if(lifeTime<0) return DBL_MAX;
  return lifeTime;
}

// override G4VRestDiscreteProcess
G4double SpontaneousFission::GetMeanFreePath(const G4Track& trk, G4double, G4ForceCondition*) {
  const G4ParticleDefinition* part = trk.GetDynamicParticle()->GetDefinition();
  G4int pZ = part->GetAtomicNumber();
  G4int pA = part->GetAtomicMass();
  if(pZ==98 && pA==252) return DBL_MIN; // force Cf252 to decay immediately
  if(part->GetPDGStable()) return DBL_MAX;
  G4double lifeTime = part->GetPDGLifeTime();
  if(lifeTime<0) return DBL_MAX;
  G4double theMass = part->GetPDGMass();
  if(theMass<DBL_MIN) {
    G4cout<<"SpontaneousFission::GetMeanFreePath - ERROR: "<<part->GetParticleName()<<" has mass null!"<<G4endl;
    exit(-1);
  }
  G4double betaGamma = trk.GetDynamicParticle()->GetTotalMomentum()/theMass;
  G4double mfp = c_light*lifeTime*betaGamma;
  if(mfp<DBL_MIN) return DBL_MIN;
  return mfp;
}

//  BuildPhysicsTable - enable print of parameters 
void SpontaneousFission::BuildPhysicsTable(const G4ParticleDefinition&) {
  if (!isInitialised) isInitialised = true;
}

G4VParticleChange* SpontaneousFission::Apply(const G4Track& trk, const G4Step&) {
  const G4ParticleDefinition* part = trk.GetDynamicParticle()->GetDefinition();

  if(IsValid(*part)) {  // check if particle is Cf252
    const G4String theVolume = trk.GetVolume()->GetLogicalVolume()->GetName();
    if (!isAllVolumesMode) { // check if particle is in target logical volume
      if (!std::binary_search(ValidVolumes.begin(), ValidVolumes.end(), theVolume))
	return Clear();
    }
    if(trk.GetTrackStatus() != fStopButAlive) { // check if particle is at rest
      G4cout<<"SpontaneousFission - ERROR: nucleus not at rest!!!"<<G4endl;
      return Clear();
    }

    // parent nucleus:
    G4double parZ = (G4double)part->GetAtomicNumber();
    G4double parA = (G4double)part->GetAtomicMass();
    G4double parM = part->GetPDGMass()/MeV;
        
    // generate input parameters:
    G4double f1Z = 0., f1A = 0;
    TH2D* hYieldAZ = (TH2D*)CfYieldFile->Get("Cf252_YxAxZ");  // for the moment just use the TH2D
    if(!hYieldAZ) { G4cout<<"SpontaneousFission::Apply - yield histogram not found!"<<G4endl; exit(-1); }
    hYieldAZ->GetRandom2(f1A, f1Z);
    G4double f2Z = parZ-f1Z;
    G4double f2A = parA-f1A-GetNeutMult(f1A);
    G4double MeanTKE = GetMeanTKE(f1A);
    G4double SigmaTKE= GetSigTKE(f1A);
    // generate the fission axis:
    G4double FAPhi = CLHEP::twopi*flatRand->fire(); FAPhi *= rad;
    G4double FAThe = std::acos(2.*flatRand->fire()-1.); FAThe *= rad;
    //    f1Z = 36.; f1A = 84.; f2Z = 36.; f2A = 84.;  // fixed fragments for testing
    //    FAPhi = 0.*rad; FAThe = 0.*rad;  // uncomment these 2 lines and the one below
    // define the two fragments:
    G4int frag1Z = lrint(f1Z); G4int frag1A = lrint(f1A);
    G4int frag2Z = lrint(f2Z); G4int frag2A = lrint(f2A);
    const G4ParticleDefinition* frag1Def = theTableOfIons->GetIon(frag1Z, frag1A, 0.);
    const G4ParticleDefinition* frag2Def = theTableOfIons->GetIon(frag2Z, frag2A, 0.);
    // generate fragment energies:
    G4double fTKE = gaussRand->fire(MeanTKE, SigmaTKE);  
    G4double frag1E = fTKE*frag2Def->GetPDGMass()/parM;
    G4double frag2E = fTKE*frag1Def->GetPDGMass()/parM;
    //    frag1E = 100.*MeV; frag2E = 100.*MeV;  // fixed fragments for testing
    // generate the fragment orientations:
    G4ThreeVector  secDir1(0.,0.,0.);
    secDir1.setRThetaPhi(1., FAThe, FAPhi);
    G4ThreeVector  secDir2(secDir1);
    secDir2 *= -1; // rotate direction of 2nd fragment with 180 degrees w.r.t. that of 1st fragment
    G4DynamicParticle* frag1 = new G4DynamicParticle(frag1Def, secDir1, frag1E);
    G4DynamicParticle* frag2 = new G4DynamicParticle(frag2Def, secDir2, frag2E);
    
    // CLEAR OUTPUT:
    fParticleChangeForDecay.Clear();
    fParticleChangeForDecay.ProposeTrackStatus(fStopAndKill);
    fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.0); // or  part->GetKineticEnergy()???
    fParticleChangeForDecay.ProposeLocalTime(trk.GetLocalTime());
    ClearNumberOfInteractionLengthLeft();
    
    // FILL THE TWO FRAGMENTS IN OUTPUT:
    G4Track* secondary1 = new G4Track(frag1, trk.GetGlobalTime(), trk.GetPosition());
    secondary1->SetGoodForTrackingFlag();
    secondary1->SetTouchableHandle(trk.GetTouchableHandle());
    fParticleChangeForDecay.AddSecondary(secondary1);
    G4Track* secondary2 = new G4Track(frag2, trk.GetGlobalTime(), trk.GetPosition());
    secondary2->SetGoodForTrackingFlag();
    secondary2->SetTouchableHandle(trk.GetTouchableHandle());
    fParticleChangeForDecay.AddSecondary(secondary2);
    // G4cout<<"SpontaneousFissionOccurence"<<std::setprecision(5)<<": Frag1(Z,A,KE)=("<<frag1Z<<","<<frag1A<<","<<frag1E/MeV<<"); Frag2(Z,A,KE)=("<<frag2Z<<","<<frag2A<<","<<frag2E/MeV<<"), FA(the,phi)=("<<FAThe/deg<<","<<FAPhi/deg<<")"<<G4endl;
  }
  
  return &fParticleChangeForDecay;
}

G4VParticleChange*  SpontaneousFission::Clear() {
  fParticleChangeForDecay.SetNumberOfSecondaries(0);
  fParticleChangeForDecay.ProposeTrackStatus(fStopAndKill) ;  // Kill the parent particle
  fParticleChangeForDecay.ProposeLocalEnergyDeposit(0.);
  ClearNumberOfInteractionLengthLeft(); // reset NumberOfInteractionLengthLeft
  return &fParticleChangeForDecay;
}

G4double SpontaneousFission::GetNeutMult(G4double fA) {
  G4double nM = 0.;
  if(fA<=72) nM = 1.45;
  else if(fA>72 && fA<=74.67) nM = 36.06 - 0.4891*fA;
  else if(fA>74.67 && fA<=75.69) nM = -152.+2.045*fA;
  else if(fA>75.69 && fA<=77.63) nM = 76.47-0.974*fA;
  else if(fA>77.63 && fA<=82.5) nM = 6.565+0.07824*fA;
  else if(fA>82.5 && fA<=94) nM = -6.99+0.0857*fA;
  else if(fA>94 && fA<=98) nM = 1.05+7.34e-08*fA;
  else if(fA>98 && fA<=102) nM = -3.891+0.0505*fA;
  else if(fA>102 && fA<=108) nM = -7.283+0.08349*fA;
  else if(fA>108 && fA<=123) nM = -10.37+0.113*fA;
  else if(fA>123 && fA<=125) nM = 34.14-0.2498*fA;
  else if(fA>125 && fA<=127) nM = 108.5-0.8442*fA;
  else if(fA>127 && fA<=130) nM = exp(48.43-0.3802*fA);
  else if(fA>130 && fA<=141) nM = -14.29+0.1129*fA;
  else if(fA>141 && fA<=162) nM = -5.976+0.0536*fA;
  else if(fA>162 && fA<=169) nM = -11.04+0.08507*fA;
  else if(fA>169 && fA<=174) nM = -36.71+0.2361*fA;
  else if(fA>174 && fA<=180) nM = 79.08-0.4331*fA;
  else if(fA>180) nM = 1.12;
  else nM = -1.;
  return nM;
}

G4double SpontaneousFission::GetMeanTKE(G4double fA) {
  G4double fE = 0.*MeV;
  if(fA<73) fE = 145.5 ;
  else if(fA>=73 && fA<88)     fE = -357.5+11.34*fA-0.06123*fA*fA;
  else if(fA>=88 && fA<100)   fE = 77.56+1.007*fA;
  else if(fA>=100 && fA<116) fE = -355.0+9.292*fA-0.03954*fA*fA;
  else if(fA>=116 && fA<127) fE = -2507.+44.56*fA-0.1837*fA*fA;
  else if(fA>=127 && fA<138) fE = -2979.+47.96*fA-0.1811*fA*fA;
  else if(fA>=138 && fA<153) fE = -341.4+8.078*fA-0.03056*fA*fA;
  else if(fA>=153 && fA<165) fE = 350.4-1.120*fA;
  else if(fA>=165 && fA<180) fE = -1496.+20.59*fA-0.06374*fA*fA;
  else if(fA>= 180) fE = 145.5 ;
  else fE = -1.;
  fE *= MeV;
  return fE;
}

G4double SpontaneousFission::GetSigTKE(G4double fA) {
  G4double fS = 0.*MeV;
  if(fA<72) fS = 10.2;
  else if(fA>=72 && fA<75)     fS = -28.24+0.533*fA;
  else if(fA>=75 && fA<80)     fS = 29.40-0.2447*fA;
  else if(fA>=80 && fA<83)     fS = 11.93-0.0257*fA;
  else if(fA>=83 && fA<88)     fS = -484.1+11.61*fA-0.06819*fA*fA;
  else if(fA>=88 && fA<110)   fS = 7.978+0.01885*fA;
  else if(fA>=110 && fA<121) fS = 176.2-3.0105*fA+0.01363*fA*fA;
  else if(fA>=121 && fA<128) fS = -660.0+10.86*fA-0.04388*fA*fA;
  else if(fA>=128 && fA<135) fS = -921.7+14.21*fA-0.05409*fA*fA;
  else if(fA>=135 && fA<146) fS = 319.5-4.281*fA+0.0148*fA*fA;
  else if(fA>=146 && fA<165) fS = 12.61-0.01775*fA;
  else if(fA>=165 && fA<170) fS = -1407+16.92*fA-0.05054*fA*fA;
  else if(fA>=170 && fA<173) fS = -400.2+4.707*fA-0.0135*fA*fA;
  else if(fA>=173 && fA<179) fS = -31.25+0.238*fA;
  else if(fA>=179 && fA<182) fS = 91.00-0.4457*fA;
  else if(fA>=182) fS = 10.2;
  else fS = -1.;
  fS *= MeV;
  return fS; 
}
