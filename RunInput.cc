#include "RunInput.hh"

using namespace std;

RunInput::RunInput() : numberOfEvent(0), visualOn(false), runMode(-1), setup(-1), targThick(0.*um), windThick(0.*um), dWheel(0.*mm), cellPres(0.*pascal), cellTemp(0.*kelvin), beamName(""), targName(""), beamEMean(0.*MeV), beamEWidth(0.*MeV), beamPWidth(0.*mm), enerRange(0.), atimaOn(true), stepDR(0.), stepFR(0.*mm), useNeutrons(false), targNistName(""), targDens(0.*g/cm3),  beamEvent(0), beamTrack(0), beamParent(0), sourceZ(0.), sourceR(0.) {

  trkNtuple = false; // switches ON sorage of track Ntuple
  regVol.clear();

  // Load input file:
  G4String BasePath = std::getenv("IonCatchersPath");
  G4String pathToFile = BasePath; pathToFile += "/RunInput.txt";
  dataFile.open(pathToFile, ios::in | ios::binary);
  if(!dataFile) { G4cout<<"RunInput - Error: Can't open the input file."<<G4endl; exit(1); }
  // read the input parameters:
  G4int visualInt = 0, elossInt=0, useNeutInt = 0;
  dataFile>>numberOfEvent>>visualInt>>runMode>>setup>>targThick>>windThick>>dWheel>>cellPres>>cellTemp>>beamName>>targName>>beamEMean>>beamEWidth>>beamPWidth>>enerRange>>elossInt>>stepDR>>stepFR>>useNeutInt;
  // fix the units for all parameters:
  targThick*=um; windThick*=um; dWheel*=mm; cellPres*=(100*pascal); cellTemp*=kelvin; beamEMean*=MeV; beamEWidth*=MeV; beamPWidth*=mm; enerRange*=um; stepFR*=um;
  // process the boolean inputs:
  visualOn = (visualInt!=0);  atimaOn = (elossInt!=0); useNeutrons = (useNeutInt!=0);

  targNistName = "G4_";
  targDens = 0.*g/cm3;
  if(targName.contains("Pt"))         {  targNistName += "Pt"; targDens = 21.45*g/cm3; }
  else if(targName.contains("Ni"))  {  targNistName += "Ni"; targDens =   8.90*g/cm3; }
  else if(targName.contains("Dy")) { targNistName += "Dy"; targDens =   8.55*g/cm3; }
  else if(targName.contains("Bi"))  {  targNistName += "Bi"; targDens =    9.78*g/cm3; }
  else if(targName.contains("U"))  {  targNistName += "URANIUM_OXIDE"; targDens = 10.96*g/cm3; }
  else if(targName.contains("Cf"))  {  targNistName += "Cf"; targDens =  15.1*g/cm3; }
  else {G4cout<<"RunInput ERROR: target "<<targName<<" not implemented!"<<G4endl; exit(-1);}

  if(runMode<0 || runMode>2) {G4cout<<"RunInput - ERROR: undefined running run mode:"<<runMode<<G4endl; exit(0);}
  if(setup<0 || setup>2) {G4cout<<"RunInput - ERROR: undefined setup:"<<setup<<G4endl; exit(0);}
  if(setup!=1 && (beamName.contains("Cf") || targName.contains("Cf"))) {G4cout<<"RunInput - ERROR: Cf source only with INCREASE setup!"<<G4endl; exit(0);}

  // list of processes names (to be increased by hand)
  fProcList[0] = "none"; fProcList[1] = "Transportation";
  fProcList[2] = "conv"; fProcList[3] = "compt"; fProcList[4] = "phot";                           // gamma
  fProcList[5] = "eBrem"; fProcList[6] = "eIoni"; fProcList[7] = "msc"; fProcList[8] = "annihil"; // e+e-
  fProcList[9] = "hadElastic"; fProcList[10] = "neutronInelastic";                                // neutron
  fProcList[11] = "hIoni"; fProcList[12] = "nuclearStopping"; fProcList[13] = "ionIoni";          // proton,alpha (1)
  fProcList[14] = "protonInelastic"; fProcList[15] = "alphaInelastic";                            // proton,alpha (2)
  fProcList[16] = "ionIonZ>3"; fProcList[17] = "ScreenedElastic"; fProcList[18] = "ionInelastic"; // ions (1)
  fProcList[19] = "MNT"; fProcList[20] = "Decay";                                              // ions (2)
  // initial list of particle names (to be increased dynamically)
  const G4String iPartList[7] = {"gamma", "e-", "e+", "neutron", "proton", "alpha", beamName};
  for(G4int i=0; i<7; i++) fPartList.push_back(iPartList[i]);
  
  // PRINT OUT VARIOUS INPUT INFO:
  G4cout<<"\nRunInput settings:"<<G4endl;
  G4cout<<"Number of Events: "<<numberOfEvent<<G4endl;
  if(visualOn) G4cout<<"Visualisation ON"<<G4endl;
  else G4cout<<"Visualisation OFF"<<G4endl;
  G4cout<<"Running mode: "<<runMode<<G4endl;
  if(runMode==2) G4cout<<"   -> WARNING: this running mode typically produces very large output ROOT files!!!"<<G4endl;
  if(setup==0) G4cout<<"Setup: IGISOL(MNT)"<<G4endl;
  if(setup==1) G4cout<<"Setup: INCREASE(Cf source)"<<G4endl;
  if(setup==2) G4cout<<"Setup: INCREASE(MNT)"<<G4endl;
  G4cout<<"Gas pressure:"<<1.e+03*cellPres/bar<<"mbar, temperature:"<<cellTemp/kelvin<<"K"<<G4endl;
  if(setup==0) G4cout<<"IGISOL window thickness="<<windThick/um<<"microns, distance(target-window)="<<dWheel/mm<<"mm"<<G4endl;
  G4cout<<"Beam:"<<beamName<<" on Target:"<<targName<<G4endl;
  if(beamName=="Cf252" && targName=="Cf252") {
    G4cout<<"Ti degrader foil thickness: "<<targThick/um<<"um"<<G4endl;
    G4cout<<"Beam Energy: 0MeV (spontaneous fission)"<<G4endl;
  } else {
    G4cout<<"Target foil thickness: "<<targThick/um<<"um"<<G4endl;
    G4cout<<"Beam Energy Mean: "<<beamEMean/MeV<<"MeV"<<G4endl;
    G4cout<<"Beam Energy Width: "<<beamEWidth/MeV<<"MeV"<<G4endl;
    G4cout<<"Beam Position Width: "<<beamPWidth/mm<<"mm"<<G4endl;
  }
  G4cout<<"Energy Range: "<<enerRange/um<<"um"<<G4endl;
  if(atimaOn) G4cout<<"Using Atima Eloss model."<<G4endl;
  else G4cout<<"Using default Eloss model."<<G4endl;
  G4cout<<"Eloss stepping (only for mode=0,2 - see PhysicsList): dRoverRange="<<stepDR<<", finalRange="<<stepFR/um<<"um"<<G4endl;
  if(10.*stepFR>targThick) G4cout<<"   -> WARNING: big ion ionization stepping!!! finalRange/targetThickness="<<stepFR/targThick<<G4endl;
  if(useNeutrons) G4cout<<"Using neutron physics."<<G4endl;
  else G4cout<<"Not using neutron physics."<<G4endl;
  if(trkNtuple) G4cout<<"Tracking Ntuple also stored."<<G4endl;
  else G4cout<<"Tracking Ntuple not stored."<<G4endl;
}

RunInput::~RunInput() {
  for(unsigned int i=0; i<21; i++) G4cout<<"Process #"<<i<<" is "<<fProcList[i]<<G4endl; G4cout<<G4endl;
  for(unsigned int i=0; i<fPartList.size(); i++) G4cout<<"Particle #"<<i<<" is "<<fPartList[i]<<G4endl; G4cout<<G4endl;
  dataFile.close();
}

// uniquely converts process string to process integer and stores the association to a list
G4int RunInput::GetProcessID(const G4Track* track) {
  G4int procID = -1;

  G4VProcess* hitProc = const_cast<G4VProcess*>(track->GetCreatorProcess());
  G4String ProcName = "none";
  if(hitProc) ProcName = hitProc->GetProcessName();

  for(unsigned int i=0; i<21; i++)
    if(ProcName == fProcList[i]) procID = i;
  if(procID == -1) {
    G4cout<<"SteppingAction::GetProcessID - Process:"<<ProcName<<" is not on the list! Please add it."<<G4endl;
    exit(0);
  }

  return procID;
}

// uniquely converts process string to process integer and stores the association to a list
G4int RunInput::GetParticleID(const G4Track* track) {
  G4int partID = -1;

  G4ParticleDefinition* hitPart = const_cast<G4ParticleDefinition*>(track->GetDefinition());
  G4String PartName = "none";
  if(hitPart) PartName = hitPart->GetParticleName();
  unsigned int NumParts = fPartList.size();

  // find a match; if no match was found, store new process string
  for(unsigned int i=0; i<NumParts; i++)
    if(PartName == fPartList[i]) partID = i;
  if(partID == -1) {
    fPartList.push_back(PartName);
    partID = fPartList.size()-1;
  }

  return partID;
}
/*
// uniquely converts process string to process integer and stores the association to a list
G4int RunInput::GetProcessID(const G4StepPoint* StepPt) {
  G4int procID = -1;

  G4VProcess* hitProc = const_cast<G4VProcess*>(StepPt->GetProcessDefinedStep());
  G4String ProcName = "none";
  if(hitProc) ProcName = hitProc->GetProcessName();

  for(unsigned int i=0; i<21; i++)
    if(ProcName == fProcList[i]) procID = i;
  if(procID == -1) {
    G4cout<<"SteppingAction::GetProcessID - Process:"<<ProcName<<" is not on the list! Please add it."<<G4endl;
    exit(0);
  }

  return procID;
}
*/
