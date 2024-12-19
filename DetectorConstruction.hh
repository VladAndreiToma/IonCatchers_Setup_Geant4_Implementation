#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4LogicalVolume.hh"
#include "G4SubtractionSolid.hh"
#include "G4TessellatedSolid.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4AutoDelete.hh"
#include "G4VisExtent.hh"

#include "CADMesh.hh"
#include "RunInput.hh"
#include "ElectricField.hh"

// Author: Paul Constantin
// Detector construction class to define:
// Construct(): geometry, materials, and volumes.
// ConstructSDandField(): sensitive deterctors and the field.
class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  DetectorConstruction(RunInput*);
  virtual ~DetectorConstruction();
  
  virtual G4VPhysicalVolume* Construct();    
  virtual void ConstructSDandField();
  G4Material* GetGasDef() { return m_gas; }
  G4Material* GetTargDef() { return m_targ; }
  
private:
  RunInput* runInput;
  G4String BasePath;
  G4Cache<ElectricField*> fEmField; // the cell electric field

  G4NistManager* nist;
  G4Material* m_vac;
  G4Material* m_gas;
  G4Material* m_targ;
  
  G4Point3D GetStructCenter(CADMesh*, G4double*);
};

#endif
