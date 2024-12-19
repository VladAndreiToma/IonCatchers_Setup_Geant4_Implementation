#include "ElectricField.hh"

#include "G4UniformElectricField.hh"
#include "G4UniformMagField.hh"
#include "G4MagneticField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4EquationOfMotion.hh"
#include "G4EqMagElectricField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4MagIntegratorDriver.hh"
#include "G4ChordFinder.hh"

#include "G4ExplicitEuler.hh"
#include "G4ImplicitEuler.hh"
#include "G4SimpleRunge.hh"
#include "G4SimpleHeum.hh"
#include "G4ClassicalRK4.hh"
#include "G4HelixExplicitEuler.hh"
#include "G4HelixImplicitEuler.hh"
#include "G4HelixSimpleRunge.hh"
#include "G4CashKarpRKF45.hh"
#include "G4RKG3_Stepper.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

// stepper type=4(ClassicalRK4, the default stepper), minimal step of 10 microns
ElectricField::ElectricField() : fFieldManager(0), fChordFinder(0), fEquation(0), fEMfield(0), fElFieldValue(), fStepper(0), fIntgrDriver(0), fStepperType(4), fMinStep(0.010*mm) {
  fEMfield = new G4UniformElectricField(G4ThreeVector(0.,0.,0.));
  fEquation = new G4EqMagElectricField(fEMfield);
  fFieldManager = GetGlobalFieldManager();
  UpdateField();
}

// stepper type=4(ClassicalRK4, the default stepper), minimal step of 10 microns
ElectricField::ElectricField(G4ThreeVector fieldVector) : fFieldManager(0), fChordFinder(0), fEquation(0), fEMfield(0), fElFieldValue(), fStepper(0), fIntgrDriver(0), fStepperType(4), fMinStep(0.010*mm) {
  fEMfield = new G4UniformElectricField(fieldVector);
  fEquation = new G4EqMagElectricField(fEMfield);
  fFieldManager = GetGlobalFieldManager();
  UpdateField();
}

ElectricField::~ElectricField() {
  delete fChordFinder;
  delete fStepper;
  delete fEquation;
  delete fEMfield;
}

void ElectricField::UpdateField(){
  // Register this field to 'global' Field Manager and
  // Create Stepper and Chord Finder with predefined type, minstep (resp.)
  SetStepper();
  G4cout<<"The minimal step is equal to "<<fMinStep/mm<<" mm"<<G4endl;
  fFieldManager->SetDetectorField(fEMfield);
  if(fChordFinder) delete fChordFinder;
  fIntgrDriver = new G4MagInt_Driver(fMinStep, fStepper, fStepper->GetNumberOfVariables());
  fChordFinder = new G4ChordFinder(fIntgrDriver);
  fFieldManager->SetChordFinder(fChordFinder);
}

void ElectricField::SetStepper(){
  // Set stepper according to the stepper type
  G4int nvar = 8;
  if(fStepper) delete fStepper;

  switch(fStepperType) {
    case 0:
      fStepper = new G4ExplicitEuler( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4ExplicitEuler is called"<<G4endl;
      break;
    case 1:
      fStepper = new G4ImplicitEuler( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4ImplicitEuler is called"<<G4endl;
      break;
    case 2:
      fStepper = new G4SimpleRunge( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4SimpleRunge is called"<<G4endl;
      break;
    case 3:
      fStepper = new G4SimpleHeum( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4SimpleHeum is called"<<G4endl;
      break;
    case 4:
      fStepper = new G4ClassicalRK4( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4ClassicalRK4 (default) is called"<<G4endl;
      break;
    case 5:
      fStepper = new G4CashKarpRKF45( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4CashKarpRKF45 is called"<<G4endl;
      break;
    case 6:
      fStepper = 0; // new G4RKG3_Stepper( fEquation, nvar );
      G4cout<<"ElectricField::SetStepper - G4RKG3_Stepper is not currently working for Electric Field"<<G4endl;
      break;
    case 7:
      fStepper = 0; // new G4HelixExplicitEuler( fEquation );
      G4cout<<"ElectricField::SetStepper - G4HelixExplicitEuler is not valid for Electric Field"<<G4endl;
      break;
    case 8:
      fStepper = 0; // new G4HelixImplicitEuler( fEquation );
      G4cout<<"ElectricField::SetStepper - G4HelixImplicitEuler is not valid for Electric Field"<<G4endl;
      break;
    case 9:
      fStepper = 0; // new G4HelixSimpleRunge( fEquation );
      G4cout<<"ElectricField::SetStepper - G4HelixSimpleRunge is not valid for Electric Field"<<G4endl;
      break;
    default: fStepper = 0;
  }
}

void ElectricField::SetFieldValue(G4double fieldValue){
  // Set the value of the Global Field to fieldValue along Z
  G4ThreeVector fieldVector( 0.0, 0.0, fieldValue );
  SetFieldValue(fieldVector);
}

void ElectricField::SetFieldValue(G4ThreeVector fieldVector){
  if(fEMfield) delete fEMfield;
  // Set the value of the Global Field value to fieldVector
  // Find the Field Manager for the global field
  G4FieldManager* fieldMgr= GetGlobalFieldManager();
  if(fieldVector != G4ThreeVector(0.,0.,0.)) {
    fEMfield = new G4UniformElectricField(fieldVector);
  } else {
    fEMfield = 0; // if the new field's value is zero, then it is best to insure that it is not used for propagation.
  }
  fieldMgr->SetDetectorField(fEMfield);
  fEquation->SetFieldObj(fEMfield);  // must now point to the new field
}

G4FieldManager* ElectricField::GetGlobalFieldManager() {
  return G4TransportationManager::GetTransportationManager()->GetFieldManager();
}
