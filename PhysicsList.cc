#include "PhysicsList.hh"

#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4LossTableManager.hh"

// PARTICLES:
#include "G4Gamma.hh"
#include "G4LeptonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// PHYSICS PROCESSES:
// GAMMA
#include "G4GammaConversion.hh"
#include "G4ComptonScattering.hh"
#include "G4PhotoElectricEffect.hh"
// E+E-
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
// NEUTRONS
#include "G4HadronElasticProcess.hh"
#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"
// PROTONS
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ProtonInelasticProcess.hh"
#include "G4CascadeInterface.hh"
#include "G4ProtonInelasticCrossSection.hh"
// ALPHAS
#include "G4AlphaInelasticProcess.hh"
#include "G4BinaryLightIonReaction.hh"
#include "G4TripathiLightCrossSection.hh"
// IONS
#include "G4ionIonisation.hh"
#include "G4AtimaEnergyLossModel.hh"
#include "G4AtimaFluctuations.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4NuclearStopping.hh"
//#include "G4BraggIonModel.hh"
#include "MNTProcess.hh"
#include "MNTDataModel.hh"
#include "SpontaneousFission.hh"
//#include "G4IonInelasticProcess.hh"
//#include "G4IonsShenCrossSection.hh"
//#include "G4QMDReaction.hh"

PhysicsList::PhysicsList(RunInput* rinp, CLHEP::HepRandomEngine* rgen) : G4VModularPhysicsList(), runInput(rinp), randGen(rgen) {
  G4LossTableManager::Instance();
  SetDefaultCutValue(runInput->GetEnerRange());
}

PhysicsList::~PhysicsList() {}

void PhysicsList::ConstructParticle() {
  G4Gamma::GammaDefinition();                    // gamma
  G4LeptonConstructor::ConstructParticle();      // leptons
  G4IonConstructor::ConstructParticle();         // ions
  G4MesonConstructor::ConstructParticle();       // hadrons
  G4BaryonConstructor::ConstructParticle();
  G4ShortLivedConstructor::ConstructParticle();  // short lived
}

void PhysicsList::ConstructProcess(){
  AddTransportation(); // register the G4Transportation class with all particle classes

  G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
  auto particleIterator = GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4String particleName = particle->GetParticleName();

    if(particleName=="gamma") {
      ph->RegisterProcess(new G4ComptonScattering(),particle);
      ph->RegisterProcess(new G4PhotoElectricEffect(),particle);
      ph->RegisterProcess(new G4GammaConversion(),particle);
    } else if(particleName=="e-" || particleName=="e+") {
      ph->RegisterProcess(new G4eMultipleScattering(),particle);
      ph->RegisterProcess(new G4eBremsstrahlung(),particle);
      ph->RegisterProcess(new G4eIonisation(),particle);
      if(particleName == "e+") ph->RegisterProcess(new G4eplusAnnihilation(),particle);
    } else if(particleName=="neutron") {
      if(runInput->UseNeutrons()) {
	// Neutron elastic process
	G4HadronElasticProcess* theNeutronElastic = new G4HadronElasticProcess("neutronElastic");
	// G4DiffuseElastic* nElasticModel = new G4DiffuseElastic();
	G4NeutronHPElastic* nElasticModel = new G4NeutronHPElastic();
	nElasticModel->SetMinEnergy(0.001*MeV);
	theNeutronElastic->RegisterMe(nElasticModel);
	G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData();
	theNeutronElastic->AddDataSet(theHPElasticData);
	ph->RegisterProcess(theNeutronElastic,particle);
	// Neutron inelastic process
	G4NeutronInelasticProcess* theNeutronInelastic = new G4NeutronInelasticProcess();
	G4NeutronHPInelastic* theHPInelasticModel = new G4NeutronHPInelastic();
	theNeutronInelastic->RegisterMe(theHPInelasticModel);
	G4NeutronHPInelasticData* theHPInelasticData = new G4NeutronHPInelasticData();
	theNeutronInelastic->AddDataSet(theHPInelasticData);
	ph->RegisterProcess(theNeutronInelastic,particle);
      }
    } else if(particleName=="proton") {
      ph->RegisterProcess(new G4hMultipleScattering(),particle);
      ph->RegisterProcess(new G4hIonisation(),particle);
      G4ProtonInelasticProcess* theProtonInelastic = new G4ProtonInelasticProcess();
      theProtonInelastic->RegisterMe(new G4CascadeInterface());
      theProtonInelastic->AddDataSet(new G4ProtonInelasticCrossSection());
      ph->RegisterProcess(theProtonInelastic,particle);
    } else if(particleName=="alpha" || particleName=="He3") {
      ph->RegisterProcess(new G4hMultipleScattering(),particle);
      ph->RegisterProcess(new G4ionIonisation(),particle);
      // Inelastic
      G4AlphaInelasticProcess* theAlphaInelastic = new G4AlphaInelasticProcess();
      theAlphaInelastic->RegisterMe(new G4BinaryLightIonReaction());
      theAlphaInelastic->AddDataSet(new G4TripathiLightCrossSection());
      ph->RegisterProcess(theAlphaInelastic,particle);
    } else if (particleName == "GenericIon") {
      ph->RegisterProcess(new G4hMultipleScattering("ionmsc"), particle);
      G4ionIonisation* ionIoni = new G4ionIonisation("ionIonZ>3");
      if(runInput->UseAtima()) {
	ionIoni->SetEmModel(new G4AtimaEnergyLossModel());
	ionIoni->SetFluctModel(new G4AtimaFluctuations());
      } else {
	G4IonParametrisedLossModel* ionplm = new G4IonParametrisedLossModel();
	ionplm->DeactivateICRU73Scaling();
	ionIoni->SetEmModel(ionplm);
      }
      G4double dRoverRange=0.2, finalRange=1.0*mm;  // default G4 values
      if(runInput->RunMode1()) {
	dRoverRange=0.1; finalRange=10.*um; // large steps (to be used with gas)
      } else {
	dRoverRange=runInput->GetStepDR(); finalRange=runInput->GetStepFR(); // small steps (to be used with solid)
      }
      G4cout<<"PhysicsList - eloss stepping dRoverRange="<<dRoverRange<<", finalRange="<<finalRange/um<<"um"<<G4endl;
      ionIoni->SetStepFunction(dRoverRange, finalRange);
      ph->RegisterProcess(ionIoni,particle);   // see twiki.cern.ch/twiki/bin/view/Geant4/LoweIonParameterized
      ph->RegisterProcess(new G4NuclearStopping(), particle);                   

      if(runInput->BeamName()=="Cf252") {
	ph->RegisterProcess(new SpontaneousFission(randGen),particle);
      } else {
	MNTProcess* theMNT = new MNTProcess();
	MNTDataModel* theMNTModel = new MNTDataModel(runInput, randGen);
	// G4HadronicInteraction.cc has been modified to disable energy conservation
	std::pair<G4double, G4double> emFatal = theMNTModel->GetFatalEnergyCheckLevels();
	G4cout<<"PhysicsList - Energy conservation level:"<<emFatal.first*100.<<"%/"<<emFatal.second/GeV<<"GeV"<<G4endl;
	theMNT->RegisterMe(theMNTModel);
	ph->RegisterProcess(theMNT,particle);
      }
      
      // inelastic (reaction) with QMD model (including G4Evaporation with GEM model) and Shen XSection
      // G4IonInelasticProcess* theIonInelastic = new G4IonInelasticProcess();
      // theIonInelastic->AddDataSet(new G4IonsShenCrossSection());
      // theIonInelastic->RegisterMe(new G4QMDReaction());
      // ph->RegisterProcess(theIonInelastic,particle);
    }
    
  } // iteration over all defined particles
}

void PhysicsList::SetCuts(){
  // TRY THE FOLLOWING: fix e lower limit for cut
  // G4ProductionCutsTable::GetProductionCutsTable()->SetEnergyRange(100*eV, 1*GeV);

  DumpCutValuesTable();
}
