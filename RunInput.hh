#ifndef RunInput_h
#define RunInput_h 1

#include <string>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4VProcess.hh"
#include "G4Track.hh"

using namespace std;

// Author: Paul Constantin
// Loads input parameters and distributes data around the module 
class RunInput {
public:
  RunInput();
  ~RunInput();

  G4int GetProcessID(const G4Track*);
  G4int GetParticleID(const G4Track*);
  //  G4int GetProcessID(const G4StepPoint*);
  
  void RegisterVolume(G4String v) { regVol.push_back(v); }
  std::vector<G4String> GetVolumes() { return regVol; }
  
  G4int GetNumEvents()           { return numberOfEvent; }
  void SetNumEvents(G4int v)  { numberOfEvent = v; }
  G4bool IsVisualOn()                { return visualOn; }
  G4bool RunMode0()                { return (runMode==0); }
  G4bool RunMode1()                { return (runMode==1); }
  G4bool RunMode2()                { return (runMode==2); }
  G4bool IgisolMNT()                  { return (setup==0); }
  G4bool IncreaseCf()                 { return (setup==1); }
  G4bool IncreaseMNT()             { return (setup==2); }
  G4double GetTargThick()        { return targThick; }
  G4double GetWindowThick()    { return windThick; }
  G4double GetDistTH()             { return dWheel; }
  G4double GetCellPres()          { return cellPres; }
  G4double GetCellTemp()        { return cellTemp; }
  G4String BeamName()           { return beamName;}
  G4String TargName()             { return targName;}
  G4String TargNistName()       { return targNistName;}
  G4double TargDens()             { return targDens;}
  G4double GetBeamEMean()   { return beamEMean; }
  G4double GetBeamEWidth()  { return beamEWidth; }
  G4double GetBeamPWidth()  { return beamPWidth; }
  G4double GetEnerRange()     { return enerRange; }
  G4bool UseAtima()                 { return atimaOn;}
  G4double GetStepDR()           { return stepDR; }
  G4double GetStepFR()           { return stepFR; }
  G4bool UseNeutrons()            { return useNeutrons; }
  void SetBeamEvent(G4int v)  { beamEvent = v; }
  void SetBeamTrack(G4int v)   { beamTrack = v; }
  void SetBeamParent(G4int v) { beamParent = v; }
  G4int GetBeamEvent()           { return beamEvent; }
  G4int GetBeamTrack()           { return beamTrack; }
  G4int GetBeamParent()          { return beamParent; }
  void SetSourceZ(G4double v) { sourceZ = v; }
  G4double GetSourceZ()          { return sourceZ; }
  void SetSourceR(G4double v) { sourceR = v; }
  G4double GetSourceR()          { return sourceR; }
  G4bool TrackNtuple()              { return trkNtuple; }

private:
  G4int numberOfEvent;
  G4bool visualOn;
  G4int runMode;
  G4int setup;
  G4double targThick;
  G4double windThick;
  G4double dWheel;
  G4double cellPres;
  G4double cellTemp;
  G4String beamName;
  G4String targName;
  G4double beamEMean;
  G4double beamEWidth;
  G4double beamPWidth;
  G4double enerRange;
  G4bool atimaOn;
  G4double stepDR;
  G4double stepFR;
  G4bool useNeutrons;
  G4String targNistName;
  G4double targDens;
  G4int beamEvent;
  G4int beamTrack;
  G4int beamParent;
  G4double sourceZ;
  G4double sourceR;
  G4bool trkNtuple;
  
  std::vector<G4String> regVol;
  std::ifstream dataFile;
  G4String fProcList[21];
  std::vector<G4String> fPartList;
};

#endif
