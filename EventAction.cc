#include "EventAction.hh"

EventAction::EventAction(RunInput* rinp, RunAction* run) : G4UserEventAction(), runInput(rinp), fRunAction(run), eventID(0),
							   sTrk(), sPrt(), sPre(), sZ(), sA(), sM(), sE(), sT(), sP(), sStp(), sx(), sy(), sz(), sDet(), sEDep(), sSLen(), sQ(),
							   tTrk(), tPrt(), tPre(), tPrc(), tZ(), tA(), tM(), tE(), tT(), tP() {
  analysisManager = G4AnalysisManager::Instance();
}

EventAction::~EventAction() {
  delete G4AnalysisManager::Instance();  
}

void EventAction::BeginOfEventAction(const G4Event* event) {
  eventID = event->GetEventID();
  // clear step vectors:
  sTrk.clear(); sPrt.clear(); sPre.clear();
  sZ.clear(); sA.clear(); sM.clear();
  sE.clear(); sT.clear(); sP.clear();
  sStp.clear(); sDet.clear();
  sx.clear(); sy.clear(); sz.clear();
  sEDep.clear(); sSLen.clear();sQ.clear();
  if(runInput->TrackNtuple()) {
    // clear track vectors:
    tTrk.clear(); tPrt.clear(); tPre.clear(); tPrc.clear();
    tZ.clear(); tA.clear(); tM.clear();
    tE.clear(); tT.clear(); tP.clear();
  }
}

void EventAction::EndOfEventAction(const G4Event* event) {
  unsigned int stpEntries = sStp.size();  // relies on same number of entries in each step vector
  unsigned int trkEntries = tTrk.size();   // relies on same number of entries in each track vector

  // fill output NTuple from vectors if:
  // either (1) event from released ion
  // or (2) event from beam and has a released ion
  G4bool fillNTuple = false;
  G4int trackID =  -1, parentID =  -1;
  if(runInput->RunMode1()) {   // event from released ion
    fillNTuple = true;
    eventID = runInput->GetBeamEvent();  // get event from input file
    trackID = runInput->GetBeamTrack();  // get track from input file
    parentID = runInput->GetBeamParent();  // get track from input file
    if(runInput->TrackNtuple())
      for(unsigned int iT=0; iT<trkEntries; iT++) FillTrkNtuple(iT, eventID, trackID, parentID);
  } else {   // event from beam
    eventID = event->GetEventID();
    if(runInput->TrackNtuple())
      for(unsigned int iT=0; iT<trkEntries; iT++) FillTrkNtuple(iT, eventID, tTrk[iT], tPre[iT]);
    for(unsigned int iS=0; iS<stpEntries; iS++) {
      G4bool reactionFrag = (sPrt[iS]>6);   // fill step Ntuple if current event has a released reaction fragment from:
      G4bool releasedGas = (sDet[iS]==5); // target (INCREASE+MNT), window (IGISOL+MNT), or Ti foil (INCREASE+Cf252)
      if(reactionFrag && releasedGas) {
	fillNTuple = true;
	break;
      }
    }  
  }

  if(runInput->RunMode2()) {
    for(unsigned int iS=0; iS<stpEntries; iS++) FillStpNtuple(iS, eventID, sTrk[iS], sPre[iS]);
    return;
  }

  if(fillNTuple) {
    if(runInput->RunMode1()) { //  event from released ion
      for(unsigned int iS=0; iS<stpEntries; iS++) FillStpNtuple(iS, eventID, trackID, parentID);
    } else { //  event from beam
      G4int currTrk = -1;
      G4bool storeStp = false;
      for(unsigned int iS=0; iS<stpEntries; iS++) {
	if(currTrk != sTrk[iS]) {
	  storeStp = true;
	  currTrk = sTrk[iS];
	}
	if(sDet[iS]==5) storeStp = false;   // fill only until release
	if(storeStp) FillStpNtuple(iS, eventID, sTrk[iS], sPre[iS]);
      }
    }
  }

  //  if(stpEntries>0) PrintStp();
  //  if(trkEntries>0) PrintTrk();
}

void EventAction::FillStpNtuple(unsigned int i, G4int evt, G4int trk, G4int par) {
  G4int stpNtu = fRunAction->GetStpNtuple();
  if(stpNtu<0) { G4cout<<"EventAction::EndOfEventAction - ERROR: invalid step ntuple!"<<G4endl; exit(-1); }
  analysisManager->FillNtupleIColumn(stpNtu, 0, evt);            // event ID
  analysisManager->FillNtupleIColumn(stpNtu, 1, trk);             // track ID
  analysisManager->FillNtupleIColumn(stpNtu, 2, par);            // parent ID
  analysisManager->FillNtupleIColumn(stpNtu, 3,  sPrt[i]);       // particle number
  analysisManager->FillNtupleIColumn(stpNtu, 4,  sZ[i]);         // particle PDG charge
  analysisManager->FillNtupleIColumn(stpNtu, 5,  sA[i]);         // particle atomic mass
  analysisManager->FillNtupleDColumn(stpNtu, 6,  sM[i]);       // particle mass
  analysisManager->FillNtupleDColumn(stpNtu, 7,  sE[i]);        // step energy (at prePoint)
  analysisManager->FillNtupleDColumn(stpNtu, 8,  sT[i]);        // step theta (at prePoint)
  analysisManager->FillNtupleDColumn(stpNtu, 9,  sP[i]);        // step theta (at prePoint)
  analysisManager->FillNtupleIColumn(stpNtu, 10, sStp[i]);     // step ID
  analysisManager->FillNtupleDColumn(stpNtu, 11, sx[i]);       // prePoint x-coordinate
  analysisManager->FillNtupleDColumn(stpNtu, 12, sy[i]);       // prePoint y-coordinate
  analysisManager->FillNtupleDColumn(stpNtu, 13, sz[i]);       // prePoint z-coordinate
  analysisManager->FillNtupleIColumn(stpNtu, 14, sDet[i]);     // detector number
  analysisManager->FillNtupleDColumn(stpNtu, 15, sEDep[i]); // step deposited energy
  analysisManager->FillNtupleDColumn(stpNtu, 16, sSLen[i]);  // step length
  analysisManager->FillNtupleDColumn(stpNtu, 17, sQ[i]);       // prepoint charge
  analysisManager->AddNtupleRow(stpNtu);
}

void EventAction::FillTrkNtuple(unsigned int i, G4int evt, G4int trk, G4int par) {
  G4int trkNtu = fRunAction->GetTrkNtuple();
  if(trkNtu<0) { G4cout<<"EventAction::EndOfEventAction - ERROR: invalid track ntuple!"<<G4endl; exit(-1); }
  analysisManager->FillNtupleIColumn(trkNtu, 0, evt);      // event ID
  analysisManager->FillNtupleIColumn(trkNtu, 1, trk);       // track ID
  analysisManager->FillNtupleIColumn(trkNtu, 2, par);      // parent ID
  analysisManager->FillNtupleIColumn(trkNtu, 3,  tPrt[i]); // particle number
  analysisManager->FillNtupleIColumn(trkNtu, 4,  tPrc[i]); // proces number
  analysisManager->FillNtupleIColumn(trkNtu, 5,  tZ[i]);    // particle PDG charge
  analysisManager->FillNtupleIColumn(trkNtu, 6,  tA[i]);    // particle atomic mass
  analysisManager->FillNtupleDColumn(trkNtu, 7,  tM[i]);  // particle mass
  analysisManager->FillNtupleDColumn(trkNtu, 8,  tE[i]);   // track energy
  analysisManager->FillNtupleDColumn(trkNtu, 9,  tT[i]);   // track theta
  analysisManager->FillNtupleDColumn(trkNtu, 10,  tP[i]); // track theta
  analysisManager->AddNtupleRow(trkNtu);
}

void EventAction::PrintStp() {
  G4int fPrec = 5;
  for(unsigned int i=0; i<sStp.size(); i++)
    G4cout<<"STEP: Det="<<sDet[i]<<", Part="<<sPrt[i]<<", Track/Step="<<"/"<<sTrk[i]<<"/"<<sStp[i]<<std::setprecision(fPrec)
	  <<", HE="<<G4BestUnit(sE[i],"Energy")<<", DE="<<G4BestUnit(sEDep[i],"Energy")<<", SL="<<G4BestUnit(sSLen[i],"Length")
	  <<", Z="<<sZ[i]<<", Q="<<sQ[i]<<", A="<<sA[i]<<", x/y/z="<<G4BestUnit(sx[i],"Length")
	  <<"/"<<G4BestUnit(sy[i],"Length")<<"/"<<G4BestUnit(sz[i],"Length")<<G4endl;
}

void EventAction::PrintTrk() {
  G4int fPrec = 5;
  for(unsigned int i=0; i<tTrk.size(); i++)
    G4cout<<"TRACK: Trk="<<tTrk[i]<<", Part="<<tPrt[i]<<", Parent="<<tPre[i]<<", Z="<<tZ[i]<<", A="<<tA[i]<<std::setprecision(fPrec)
	  <<", E="<<G4BestUnit(tE[i],"Energy")<<", T="<<G4BestUnit(tT[i],"Angle")<<G4endl;
}
