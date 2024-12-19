#include "MNTCrossSection.hh"

using namespace std;

#include "G4CrossSectionFactory.hh"
G4_DECLARE_XS_FACTORY(MNTCrossSection);

MNTCrossSection::MNTCrossSection() : G4VCrossSectionDataSet(Default_Name()) {}

MNTCrossSection::~MNTCrossSection() {}

void MNTCrossSection::CrossSectionDescription(std::ostream& outFile) const {
  outFile<<"MNTCrossSection provides a parameterization of the Karpov cross-section for nucleus-nucleus process in the energy range 5-20MeV.\n";
}

G4bool MNTCrossSection::IsIsoApplicable(const G4DynamicParticle* particle, G4int Z, G4int /*A*/, const G4Element* /*elm*/, const G4Material* /*mat*/) {
  //   const G4double partEner = particle->GetKineticEnergy();
  const G4ParticleDefinition* partDef = particle->GetParticleDefinition();
  G4String pName = partDef->GetParticleName();
  G4bool isNi = (Z==28);
  G4bool isPt = (Z==78);
  G4bool isBi = (Z==83);
  G4bool isUr = (Z==92);
  if(pName=="Xe136" && (isPt || isBi)) return true;
  if(pName=="U238" && (isNi || isBi || isUr)) return true;
  return false;
}

G4double MNTCrossSection::GetIsoCrossSection(const G4DynamicParticle* aPart, G4int Z, G4int /*A*/,
						const G4Isotope* /*iso*/, const G4Element* /*elm*/, const G4Material* /*mat*/) {
  G4double sigma = 0.*millibarn;
  G4double bE = aPart->GetKineticEnergy()/MeV;
  const G4ParticleDefinition* partDef = aPart->GetParticleDefinition();
  G4double bEperN = bE/partDef->GetAtomicMass();
  G4String pName = partDef->GetParticleName();
  
  // cross sections for U238 beam:
  if(pName=="U238") {
    G4double par0 = 0., par1 = 0., par2 = 0., par3 = 0.;
    if(Z==28) {               // for target Ni
      if(bEperN>5.25 && bEperN<21.) {   // fit safe region
	par0 = -8840.0*millibarn; par1 = 2315.*millibarn; par2 = -133.8*millibarn; par3 = 2.759*millibarn;
      }
    } else if(Z==66) {    // for target Dy
      if(bEperN>6.28 && bEperN<21.)  {   // fit safe region
	par0 = -13158.0*millibarn; par1 = 3258.0*millibarn; par2 = -193.0*millibarn; par3 = 4.045*millibarn;
      }
    } else if(Z==83) {   // for target Bi
      if(bEperN>6.28 && bEperN<21.) {   // fit safe region
	par0 = -15040.0*millibarn; par1 = 3545.0*millibarn;  par2 = -210.9*millibarn; par3 = 4.441*millibarn;
      }
    } else if(Z==92) {   // for target U
      if(bEperN>5.57 && bEperN<21.) {   // fit safe region
	par0 = -8787.0*millibarn; par1 = 2062.6*millibarn;  par2 = -96.27*millibarn; par3 = 1.658*millibarn;
      }
    } else { G4cout<<"MNTCrossSection: beam/target "<<pName<<"/"<<Z<<" not implemented!"<<G4endl; exit(-1); }
    sigma = par0 + par1*bEperN + par2*bEperN*bEperN + par3*bEperN*bEperN*bEperN;
  }
  // cross sections for Xe136 beam:
  if(pName=="Xe136") {
    if(Z==78) { // for target Pt
      const G4int nVal = 7;
      G4double En[nVal] = {700, 785, 885, 935, 985, 1035, 1085};
      G4double Sig[nVal] = {0, 1206, 1971, 2288, 2578, 2838, 3067};
      for(G4int iv=0; iv<nVal; iv++) { En[iv] *= MeV; Sig[iv] *= millibarn; }
      G4double inter=0., slope=0.;
      for(G4int iv=0; iv<(nVal-1); iv++) {
	if(bE>En[iv] && bE<=En[iv+1]) {
	  slope = (Sig[iv+1]-Sig[iv])/(En[iv+1]-En[iv]);
	  inter = Sig[iv]-(slope*En[iv]);
	}
      }
      sigma = inter+slope*bE;
    }
    if(Z==83) { // for target Bi
      const G4int nVal = 6;
      G4double En[nVal] = {736, 740, 790, 840, 890, 940};
      G4double Sig[nVal] = {0, 994, 1468, 1872, 2229, 2281};
      for(G4int iv=0; iv<nVal; iv++) { En[iv] *= MeV; Sig[iv] *= millibarn; }
      G4double inter=0., slope=0.;
      for(G4int iv=0; iv<(nVal-1); iv++) {
	if(bE>En[iv] && bE<=En[iv+1]) {
	  slope = (Sig[iv+1]-Sig[iv])/(En[iv+1]-En[iv]);
	  inter = Sig[iv]-(slope*En[iv]);
	}
      }
      sigma = inter+slope*bE;
    }
  }
  
  if(sigma<0) {
    G4cout<<"MNTCrossSection: negative xsection="<<sigma/millibarn<<" mb at "<<bEperN<<" MeV/n"<<G4endl;
    exit(-1);
  }
  // sigma *= 1000.; // GENERATE LOTS OF INTERACTIONS FOR TESTING PURPOSES!!!

  // G4cout<<"MNTCrossSection:"<<sigma/barn<<"b at E="<<bEperN<<"MeV/n for U238 on (Z,A)=("<<Z<<","<<A<<")"<<G4endl;
  return sigma;
}
