#include "RunAction.hh"

RunAction::RunAction(RunInput* rinp, DetectorConstruction* dc) : G4UserRunAction(),  runInput(rinp), setup(dc), fVolList(), stpNtu(-1), trkNtu(-1) {
  G4RunManager::GetRunManager()->SetPrintProgress(100000); // set printing event number per XYZ events  

  // The choice of analysis technology is done via selecting a namespace in Analysis.hh
  analysisManager = G4AnalysisManager::Instance();
  G4cout<<"Using for analysis manager: "<<analysisManager->GetType()<<G4endl;
}

RunAction::~RunAction() {
  delete G4AnalysisManager::Instance();  
}

void RunAction::BeginOfRunAction(const G4Run* /*run*/) {
  // define list of geometry volumes for output - KEEP "Test3" ON POSITION 5!!!
  // "Targ" is MNT targ or Ti foil for Cf source
  if(runInput->IncreaseCf()) {
    const G4String iVolList[11] = {"Hall", "Cell", "Targ", "Test1", "Test2", "Test3", "Test4", "Test5", "Back", "CfAu", "Cone"};
    for(G4int i=0; i<11; i++) fVolList.push_back(iVolList[i]);
  }
  if(runInput->IncreaseMNT()) {
    const G4String iVolList[8] = {"Hall", "Cell", "Targ", "Test1", "Test2", "Test3", "Test4", "Test5"};
    for(G4int i=0; i<8; i++) fVolList.push_back(iVolList[i]);
  }
  if(runInput->IgisolMNT()) {
    const G4String iVolList[8] = {"Hall", "Cell", "Targ", "Test1", "Test2", "Test3", "Window", "Pipe"};
    for(G4int i=0; i<8; i++) fVolList.push_back(iVolList[i]);
  }
  std::vector<G4String> regList = runInput->GetVolumes();
  for(unsigned int i=0; i<regList.size(); i++) fVolList.push_back(regList[i]);
  
  /* Ntuples get automatically attributed an integer identifier which value is returned from the "Create" function.
     The default start value is 0 and it is incremented by 1 for each next created ntuple. The start ntuple identifier
     value can be changed with the SetFirstNtupleId(G4int) function.
  */
  stpNtu = analysisManager->CreateNtuple("StpNTuple", "Geant4 NTuple of G4Steps"); // see EventAction.cc for entry explanation
  analysisManager->CreateNtupleIColumn(stpNtu, "evt");
  analysisManager->CreateNtupleIColumn(stpNtu, "trk");
  analysisManager->CreateNtupleIColumn(stpNtu, "pare");
  analysisManager->CreateNtupleIColumn(stpNtu, "part");
  analysisManager->CreateNtupleIColumn(stpNtu, "Z");
  analysisManager->CreateNtupleIColumn(stpNtu, "A");
  analysisManager->CreateNtupleDColumn(stpNtu, "M");
  analysisManager->CreateNtupleDColumn(stpNtu, "E");
  analysisManager->CreateNtupleDColumn(stpNtu, "T");
  analysisManager->CreateNtupleDColumn(stpNtu, "P");
  analysisManager->CreateNtupleIColumn(stpNtu, "stp");
  analysisManager->CreateNtupleDColumn(stpNtu, "x");
  analysisManager->CreateNtupleDColumn(stpNtu, "y");
  analysisManager->CreateNtupleDColumn(stpNtu, "z");
  analysisManager->CreateNtupleIColumn(stpNtu, "det");
  analysisManager->CreateNtupleDColumn(stpNtu, "edep");
  analysisManager->CreateNtupleDColumn(stpNtu, "slen");
  analysisManager->CreateNtupleDColumn(stpNtu, "Q");
  analysisManager->FinishNtuple(stpNtu);

  if(runInput->TrackNtuple()) {
    trkNtu = analysisManager->CreateNtuple("TrkNTuple", "Geant4 NTuple of G4Tracks"); // see EventAction.cc for entry explanation
    analysisManager->CreateNtupleIColumn(trkNtu, "evt");
    analysisManager->CreateNtupleIColumn(trkNtu, "trk");
    analysisManager->CreateNtupleIColumn(trkNtu, "pare");
    analysisManager->CreateNtupleIColumn(trkNtu, "part");
    analysisManager->CreateNtupleIColumn(trkNtu, "proc");
    analysisManager->CreateNtupleIColumn(trkNtu, "Z");
    analysisManager->CreateNtupleIColumn(trkNtu, "A");
    analysisManager->CreateNtupleDColumn(trkNtu, "M");
    analysisManager->CreateNtupleDColumn(trkNtu, "E");
    analysisManager->CreateNtupleDColumn(trkNtu, "T");
    analysisManager->CreateNtupleDColumn(trkNtu, "P");
    analysisManager->FinishNtuple(trkNtu);
  }

  // Open an output file
  G4String fileName = "NTuple";
  if(runInput->RunMode0()) fileName += "_Mode0";
  if(runInput->RunMode1()) fileName += "_Mode1";
  if(runInput->RunMode2()) fileName += "_Mode2";
  analysisManager->OpenFile(fileName);
}

void RunAction::EndOfRunAction(const G4Run* /*run*/) {
  // print info on xsections:
  // PrintXSections();

  // print info on EM processes:
  // PrintEMInfo();

  // save histograms & ntuple
  analysisManager->Write();
  analysisManager->CloseFile();
}

void RunAction::PrintXSections() {
  G4cout<<"RunAction::EndOfRunAction - start printing cross-section data:"<<G4endl;
  // SEE ALSO PrintEMInfo() BELOW!!!
  const G4ParticleDefinition* particle = G4ParticleTable::GetParticleTable()->FindParticle("GenericIon");
  G4VProcess* process = G4ProcessTable::GetProcessTable()->FindProcess("MNT", particle);
  const G4Element* elem  = G4NistManager::Instance()->FindOrBuildElement(runInput->TargName());
  const G4Material* mate = G4NistManager::Instance()->FindOrBuildMaterial(runInput->TargNistName());
  if(!particle || !process || !elem || !mate) {
    if(!particle) G4cout<<"RunAction Error: Particle undefined. No XSections!."<<G4endl;
    if(!process) G4cout<<"RunAction Error: Process undefined. No XSections!."<<G4endl;
    if(!elem)     G4cout<<"RunAction Error: Element undefined. No XSections!."<<G4endl;
    if(!mate)     G4cout<<"RunAction Error: Material undefined. No XSections!."<<G4endl;
  } else {
    G4double KE = 12.*238.*MeV;   // MNTCrossSections parameterizes xsection vs beam KA per nucleon
    G4double xs = G4HadronicProcessStore::Instance()->GetCrossSectionPerAtom(particle, KE,  process, elem, mate);
    xs /= millibarn;
    G4cout<<"MNT XSection at "<<KE<<"MeV for beam "<<runInput->BeamName()<<" on target "<<runInput->TargName()<<" is "<<xs<<"mb"<<G4endl;
  }
}

void RunAction::PrintEMInfo() {
  G4cout<<"RunAction::EndOfRunAction - start printing EM interaction data:"<<G4endl;
  G4EmCalculator emCal; //emCal.SetVerbose(2); uncomment the verbose to print directly from G4EmCal

  // set particle (must be ion, including proton, see cuts below):
  G4int pZ = 36, pA = 84;
  G4ParticleDefinition* particle = G4IonTable::GetIonTable()->GetIon(pZ, pA, 0.*keV);
  const G4int numE = 15;
  G4double EperA[numE] = {0.1, 0.3, 0.5, 0.7, 0.9, 1., 3., 5., 7., 9., 10., 30., 50., 70., 90.}; // energy per nucleon
  for(size_t i=0; i<numE; i++) EperA[i] *= MeV; // set unit
  G4double energy[numE] = {0.*MeV};
  for(size_t i=0; i<numE; i++) energy[i] = EperA[i]*(G4double)pA;
 
  // set materials
  const G4int numM = 2;
  G4Material* material[numM];
  material[0] = setup->GetGasDef(); G4cout<<"Material He: density="<<material[0]->GetDensity()/(mg/cm3)<<"mg/cm3"<<G4endl;
  material[1] = setup->GetTargDef(); G4cout<<"Material Ti: density="<<material[1]->GetDensity()/(mg/cm3)<<"mg/cm3"<<G4endl;

  // set range and energy cuts:
  G4ProductionCutsTable* theCoupleTable = G4ProductionCutsTable::GetProductionCutsTable();
  G4double energyCut[numM]={0.*MeV}, rangeCut[numM]={0.*MeV};
  for(size_t i=0; i<theCoupleTable->GetTableSize(); i++) {
     const G4MaterialCutsCouple* couple = theCoupleTable->GetMaterialCutsCouple(i);
     for(size_t j=0; j<numM; j++) {
       if(couple->GetMaterial() == material[j]) {
	 energyCut[j] = (*(theCoupleTable->GetEnergyCutsVector(idxG4ProtonCut)))[i]; // only for gamma,e+,e-,p
	 rangeCut[j]  = (*(theCoupleTable->GetRangeCutsVector(idxG4ProtonCut)))[i]; // only for gamma,e+,e-,p
	 G4cout<<"Material "<<material[j]->GetName()<<" - Proton energy cut:"<<G4BestUnit(energyCut[j],"Energy")
	       <<", Range cut: "<<G4BestUnit(rangeCut[j],"Length")<<G4endl;
       }
     }
  }
                
  // get processList and extract EM processes (but not MultipleScattering)
  G4ProcessVector* plist = particle->GetProcessManager()->GetProcessList();
  std::deque<G4String> emName;
  for(unsigned int j=0; j<plist->size(); j++) {
    G4String procName = (*plist)[j]->GetProcessName();
    if(((*plist)[j]->GetProcessType()==fElectromagnetic) && (procName != "msc")) emName.push_back(procName);
  }

  // PRINT OUT CROSS-SECTIONS AND MEAN FREE PATHS - add code later...

  // PRINT STOPPING POWERS:
  std::deque<G4double> dEdxRD; std::deque<G4double> dEdxFD; // restricted and full stopping powers
  for(size_t i=0; i<emName.size(); i++) { // processes
    for(size_t j=0; j<numM; j++) { // materials
      for(size_t k=0; k<numE; k++) { // energies
	G4double dedxR = emCal.ComputeDEDX(energy[k],particle,emName[i],material[j],energyCut[j]); // restricted dEdx (note final energy cut)
	G4double dedxF = emCal.ComputeDEDX(energy[k],particle,emName[i],material[j],energy[j]); // full dEdx (note final energy)
	dEdxRD.push_back(dedxR);
	dEdxFD.push_back(dedxF);
      }
    }
  }
  G4cout<<"EnerA[MeV/n]: "; for(size_t k=0; k<numE; k++) G4cout<<EperA[k]/MeV<<" "; G4cout<<G4endl;
  G4cout<<"Energy[MeV]: "; for(size_t k=0; k<numE; k++) G4cout<<energy[k]/MeV<<" "; G4cout<<G4endl;
  for(size_t i=0; i<emName.size(); i++) {
    G4cout<<"dE/dx for process "<<emName[i]<<" and particle "<<particle->GetParticleName()<<":"<<G4endl;
    for(size_t j=0; j<numM; j++) {
      G4double dens = material[j]->GetDensity();
      G4cout<<"   Restricted in material "<<material[j]->GetName()<<" [MeV/mm->MeV*cm2/mg]:"<<G4endl;
      for(size_t k=0; k<numE; k++) { G4cout<<"   "<<dEdxRD.front()*mm/MeV<<"->"<<dEdxRD.front()*mg/(MeV*cm2*dens); dEdxRD.pop_front(); } G4cout<<G4endl;
      G4cout<<"   Full in material "<<material[j]->GetName()<<"  [MeV/mm->MeV*cm2/mg]:"<<G4endl;
      for(size_t k=0; k<numE; k++) { G4cout<<"   "<<dEdxFD.front()*mm/MeV<<"->"<<dEdxFD.front()*mg/(MeV*cm2*dens); dEdxFD.pop_front();} G4cout<<G4endl;
    }
  }
    
  // PRINT RANGES: -> code below fails for material "G4_He"
  /*
  G4double EmaxTable = G4EmParameters::Instance()->MaxEnergyForCSDARange();
  for(size_t j=0; j<numE; j++) {
    G4cout<<"At "<<G4BestUnit(energy[j],"Energy")<<G4endl;
    G4double range1 = emCal.GetRangeFromRestricteDEDX(energy[j],particle,material);
    G4double range2 = range1*density;
    G4cout<<"Restricted range from: dE/dx="<<std::setw(8)<<G4BestUnit(range1,"Length")<<", dE/dx/rho="<<std::setw(8)<<G4BestUnit(range2,"Mass/Surface")<<G4endl;
    if(energy[j]<EmaxTable) {
      G4double Range1 = emCal.GetCSDARange(energy[j],particle,material);
      G4double Range2 = Range1*density;
      G4cout<<"Full  range from: dE/dx="<<std::setw(8)<<G4BestUnit(Range1,"Length")<<", dE/dx/rho="<<std::setw(8)<<G4BestUnit(Range2,"Mass/Surface")<<G4endl;
    }
  }
  */
}
