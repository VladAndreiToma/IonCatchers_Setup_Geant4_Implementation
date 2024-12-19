#ifndef MNTCrossSection_h
#define MNTCrossSection_h 1

#include <iostream>
#include <fstream>

#include "G4VCrossSectionDataSet.hh"
#include "G4ParticleDefinition.hh"
#include "G4DynamicParticle.hh"
#include "G4SystemOfUnits.hh"
#include "G4Material.hh"

// Author: Paul Constantin
// RunInput pointer cannot be passed to MNTCrossSection because of inheritance errors
class MNTCrossSection : public G4VCrossSectionDataSet {
public:
  MNTCrossSection();
  virtual ~MNTCrossSection();
  
  static const char* Default_Name() {return "MNTXS";}
  virtual void CrossSectionDescription(std::ostream&) const;

  /*
    In G4VCrossSectionDataSet:
    - both IsElementApplicable() and IsIsoApplicable() return false and both GetElementCrossSection() and GetIsoCrossSection() throw exceptions.
    Hence, implement only the one which you want the framework to use:
    - if you have isotope-wise cross-sections, implement just IsIsoApplicable and GetIsoCrossSection;
    - if you have element-wise cross-sections, implement just IsElementApplicable and GetElementCrossSection;
    - the arguments G4Element and G4Material are needed to access low-energy neutron cross sections, but are not required for others.
  */
//   virtual G4bool IsElementApplicable(const G4DynamicParticle*, G4int Z, const G4Material* mat=0);
  virtual G4bool IsIsoApplicable(const G4DynamicParticle*, G4int Z, G4int A, const G4Element* elm=0, const G4Material* mat=0);

//   virtual G4double GetElementCrossSection(const G4DynamicParticle*, G4int Z, const G4Material* mat=0);
  virtual G4double GetIsoCrossSection(const G4DynamicParticle*, G4int Z, G4int A, const G4Isotope* iso=0, const G4Element* elm=0, const G4Material* mat=0);
};

#endif
