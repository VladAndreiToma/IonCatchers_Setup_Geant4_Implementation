#include "MNTDataModel.hh"
#include "G4QMDReaction.hh"

using namespace std;

MNTDataModel::MNTDataModel(RunInput* rinp, CLHEP::HepRandomEngine* rgen) : G4HadronicInteraction("MNT"),  runInput(rinp), randGen(rgen) {
  SetMinEnergy(0.*MeV);   // actual energy limits in the cross-section class
  SetMaxEnergy(100.*MeV);
  theTableOfIons = G4ParticleTable::GetParticleTable()->GetIonTable();
  G4String BasePath = std::getenv("IonCatchersPath");
  G4String pathToFile = BasePath; pathToFile += "/MNTExcitations.root";
  inFile = new TFile(pathToFile,"READ");
  if(!inFile) { G4cout<<"MNTDataModel - cannot read input file!"<<G4endl; exit(-1); }

  flatRand = new G4RandFlat(randGen);
 
  G4String TreeName = "KarpovTree"; TreeName += runInput->TargName();
  treeKin = (TTree*)inFile->Get(TreeName);
  if(!treeKin) { G4cout<<"MNTDataModel: cannot get the kinematics tree!"<<G4endl; exit(-1); }
  else {
    treeKin->SetBranchAddress("fZ",&fZ);
    treeKin->SetBranchAddress("fA",&fA);
    treeKin->SetBranchAddress("fN",&fN);
    treeKin->SetBranchAddress("fE",&fE);
    treeKin->SetBranchAddress("fS",&fS);
    treeKin->SetBranchAddress("tKin",tKin);
    treeKin->SetBranchAddress("tThe",tThe);
    treeKin->SetBranchAddress("tSig",tSig);
    for(G4int i=0; i<20000; i++) { tKin[i]=0.; tThe[i]=0.; tSig[i]=0.; }   // adjust the size to that in the input tree
  }
  
  maxKin = 4200.;   // TUNE THE maxKin SETTING FROM DATA!!!
  G4double wBinKin = 5., wBinThe = 1.;
  G4int nBinKin = G4int(maxKin/wBinKin);
  G4int nBinThe = G4int(180./wBinThe);
  histKin = new TH2D("histKin","temporary kinematics histogram",nBinKin,0.,maxKin,nBinThe,0.,180.);
}

MNTDataModel::~MNTDataModel() {
  // do not delete flatRand, histKin
  inFile->Close();
  delete inFile;
}

G4HadFinalState* MNTDataModel::ApplyYourself(const G4HadProjectile& theProjectile, G4Nucleus& theNucleus) {
  const G4ParticleDefinition* ProjecDef = theProjectile.GetDefinition();
  if(!ProjecDef) { G4cout<<"MNTDataModel::ApplyYourself - Undefined projectile! Exiting..."<<G4endl; exit(-1); }
  G4String pName = ProjecDef->GetParticleName();
  G4int ProjecZ = ProjecDef->GetAtomicNumber();
  G4int ProjecA = ProjecDef->GetAtomicMass();
  G4int TargetZ = theNucleus.GetZ_asInt();
  G4int TargetA = theNucleus.GetA_asInt();
  G4double ProjecE = theProjectile.GetKineticEnergy()/MeV;
  G4double ProjecEperN = ProjecE/ProjecA;
  // set beam energy to the nearest value available in the excitation file:
  G4int bE = NearestEnergy(ProjecE);
  G4int bEperN = 2*lround(ProjecEperN/2.);
  if(bEperN<6) bEperN = 6;
  if(TargetZ==92 && bEperN<8) bEperN = 8; // UU does not go to 6...
  if(bEperN>18) bEperN = 18;
  G4int dataEn = 0;
  if(pName=="Xe136") dataEn = bE; // Xe file uses energy
  if(pName=="U238") dataEn = bEperN; // U file uses energy per nucleon

  theParticleChange.Clear();
  theParticleChange.SetStatusChange(stopAndKill); // projectile does not survive the process

  if(IsValidProjectile(ProjecZ,ProjecA) && IsValidTarget(TargetZ)) {
    // generate first fragment (Z,A) from MNT sigma(Z,A) cross section at beam energy:
    G4String hName = "SigFrag"; hName += runInput->TargName(); hName += "_"; hName += std::to_string(dataEn);
    TH2D* hSigZA = (TH2D*)inFile->Get(hName);
    if(!hSigZA) { G4cout<<"MNTDataModel::ApplyYourself - cross section "<<hName<<" not found!"<<G4endl; exit(-1); }
    G4double tempZ = 0., tempA = 0.;
    hSigZA->GetRandom2(tempZ, tempA);
    G4int Frag1Z = lround(tempZ);
    G4int Frag1A = lround(tempA);
    const G4ParticleDefinition* Frag1Def = theTableOfIons->GetIon(Frag1Z, Frag1A, 0.);   // fragments are already evaporated

    // generate first fragment (KE,Theta) from MNT sigma(KE,Theta) cross section at above (Z,A) and beam energy:
    G4double Frag1M = Frag1Def->GetPDGMass()/MeV;
    G4double Frag1Phi = 2*CLHEP::pi*flatRand->fire(); Frag1Phi *= rad;
    G4double Frag1Kin = 0., Frag1The = 0.;
    GetKinHist(dataEn, TargetZ, Frag1Z, Frag1A);
    if(histKin->Integral()>0.) histKin->GetRandom2(Frag1Kin, Frag1The);
    else G4cout<<"MNTFailure: beamE="<<std::setprecision(5)<<ProjecEperN/MeV<<"MeV/n; Frag1(Z,A)=("<<Frag1Z<<","<<Frag1A<<") has empty kinematics histo:"<<histKin->Integral()<<G4endl;
    Frag1Kin *= MeV; Frag1The *= CLHEP::pi/180.; Frag1The *= rad;
    G4double Frag1P = sqrt(pow(Frag1Kin+Frag1M,2)-pow(Frag1M,2));   // E^2=p^2+m^2 and E=KE+m
    G4double Frag1Px = Frag1P*sin(Frag1The)*cos(Frag1Phi);
    G4double Frag1Py = Frag1P*sin(Frag1The)*sin(Frag1Phi);
    G4double Frag1Pz = Frag1P*cos(Frag1The);
    const G4ThreeVector Frag1Mom(Frag1Px, Frag1Py, Frag1Pz);

    // generate second fragment (Z,A) from MNT sigma(Z,A) cross section at beam energy:
    hSigZA->GetRandom2(tempZ, tempA);
    G4int Frag2Z = lround(tempZ);
    G4int Frag2A = lround(tempA);
    const G4ParticleDefinition* Frag2Def = theTableOfIons->GetIon(Frag2Z, Frag2A, 0.);   // fragments are already evaporated
    // generate second fragment (KE,Theta) from MNT sigma(KE,Theta) cross section at above (Z,A) and beam energy:
    G4double Frag2M = Frag2Def->GetPDGMass()/MeV;
    G4double Frag2Phi = CLHEP::pi-Frag1Phi; Frag2Phi *= rad;
    G4double Frag2Kin = 0., Frag2The = 0.;
    GetKinHist(dataEn, TargetZ, Frag2Z, Frag2A);
    if(histKin->Integral()>0.) histKin->GetRandom2(Frag2Kin, Frag2The);
    else G4cout<<"MNTFailure: beamE="<<std::setprecision(5)<<ProjecEperN/MeV<<"MeV/n; Frag2(Z,A)=("<<Frag2Z<<","<<Frag2A<<") has empty kinematics histo:"<<histKin->Integral()<<G4endl;
    Frag2Kin *= MeV; Frag2The *= CLHEP::pi/180.; Frag2The *= rad;
    G4double Frag2P = sqrt(pow(Frag2Kin+Frag2M,2)-pow(Frag2M,2));
    G4double Frag2Px = Frag2P*sin(Frag2The)*cos(Frag2Phi);
    G4double Frag2Py = Frag2P*sin(Frag2The)*sin(Frag2Phi);
    G4double Frag2Pz = Frag2P*cos(Frag2The);
    const G4ThreeVector Frag2Mom(Frag2Px, Frag2Py, Frag2Pz);

    if(Frag1Kin>0. && Frag2Kin>0.) {
      G4DynamicParticle* Frag1 = new G4DynamicParticle(Frag1Def, Frag1Mom);
      theParticleChange.AddSecondary(Frag1);
      G4DynamicParticle* Frag2 = new G4DynamicParticle(Frag2Def, Frag2Mom);
      theParticleChange.AddSecondary(Frag2);
   
      G4cout<<"MNTOccurence: beamE="<<std::setprecision(5)<<ProjecEperN/MeV<<"MeV/n->"<<dataEn<<"; Frag1(Z,A,KE,t)=("<<Frag1Z<<","<<Frag1A<<","<<Frag1Kin/MeV<<","<<Frag1The/deg<<"); Frag2(Z,A,KE,t)=("<<Frag2Z<<","<<Frag2A<<","<<Frag2Kin/MeV<<","<<Frag2The/deg<<")"<<G4endl;
    }

  } else {
    G4cout<<"MNTDefault: Target(Z,A)=("<<TargetZ<<","<<TargetA<<"), Projec(Z,A,E/A)=("<<ProjecZ<<","<<ProjecA<<","<<ProjecEperN<<")"<<G4endl;
  }

  return &theParticleChange;
}

/////////////////////////////////////////////////////////////////////
// PRIVATE HELPER METHODS FOR THE GENERATION OF VARIOUS DISTRIBUTIONS
/////////////////////////////////////////////////////////////////////

G4bool MNTDataModel::IsValidProjectile(G4int Z, G4int A){
  if(Z==92 && A==238) return true;  // U238
  if(Z==54 && A==136) return true;  // Xe136
  return false;
}

G4bool MNTDataModel::IsValidTarget(G4int Z){
  G4String targName = runInput->TargName();
  if(Z==28 && targName.contains("Ni")) return true;  // nickel
  if(Z==66 && targName.contains("Dy")) return true;  // dysprosium
  if(Z==78 && targName.contains("Pt")) return true;  // platinum
  if(Z==83 && targName.contains("Bi")) return true;  // bismuth
  if(Z==92 && targName.contains("U")) return true;  // uranium
  return false;
}

// tree branches: fZ, fA, fN, fE, fS, tKin[fN], tThe[fN], tSig[fN]
void MNTDataModel::GetKinHist(G4int beamE, G4int targZ, G4int fragZ, G4int fragA) {
  histKin->Reset();   // reset before filling
  G4int nEntries = treeKin->GetEntries();
  for(G4int iEvt=0; iEvt<nEntries; iEvt++){
    treeKin->GetEvent(iEvt);
    if(fN>=20000) { G4cout<<"MNTDataModel::GetKinHist - fN="<<fN<<" out of bounds!"<<G4endl; exit(-1); }
    if(fE != beamE) continue;   // select beam energy
    if(fZ != fragZ) continue;   // select fragment charge
    if(fA != fragA) continue;   // select fragment mass
    for(G4int iTrk=0; iTrk<fN; iTrk++) {
      if(tKin[iTrk]>=maxKin) { G4cout<<"MNTDataModel::GetKinHist - KE="<<tKin[iTrk]<<" out of bounds!"<<G4endl; exit(-1); }
      histKin->Fill(tKin[iTrk], tThe[iTrk], tSig[iTrk]);
    }
  }
}

// finds the nearest energy in the MNT input file
G4int MNTDataModel::NearestEnergy(G4double inE) {
  G4int outE = -1;
  if(inE<=835.*MeV) outE=785;
  if(inE>835.*MeV   && inE<=910.*MeV) outE=885;
  if(inE>910.*MeV   && inE<=960.*MeV) outE=935;
  if(inE>960.*MeV   && inE<=1010.*MeV) outE=985;
  if(inE>1010.*MeV && inE<=1060.*MeV) outE=1035;
  if(inE>1060.*MeV) outE=1085;
  return outE;
}
