// ==================================================================
// HPGeSensitiveDetector.cc
// Professional Sensitive Detector Implementation
// ==================================================================

#include "HPGeSensitiveDetector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHistory.hh"
#include "G4VProcess.hh"  // ADDED: Include full definition of G4VProcess
#include "G4ios.hh"

// ==================================================================
// Constructor
// ==================================================================
HPGeSensitiveDetector::HPGeSensitiveDetector(const G4String& name, G4int detectorID)
: G4VSensitiveDetector(name),
  fDetectorID(detectorID),
  fTotalEdep(0.),
  fEventTime(0.),
  fFirstHit(true),
  fPeakEnergy(0.),
  fMultiplicity(0) {
}

// ==================================================================
// Destructor
// ==================================================================
HPGeSensitiveDetector::~HPGeSensitiveDetector() {}

// ==================================================================
// Initialize - called at beginning of each event
// ==================================================================
void HPGeSensitiveDetector::Initialize(G4HCofThisEvent*) {
    Reset();
}

// ==================================================================
// Process hits - called for each step in sensitive volume
// ==================================================================
G4bool HPGeSensitiveDetector::ProcessHits(G4Step* step, G4TouchableHistory*) {
    
    // Get energy deposited in this step
    G4double edep = step->GetTotalEnergyDeposit();
    
    if (edep <= 0.) return false;
    
    // Accumulate total energy
    fTotalEdep += edep;
    
    // Record time of first hit
    if (fFirstHit) {
        fEventTime = step->GetPreStepPoint()->GetGlobalTime();
        fFirstHit = false;
    }
    
    // Track peak energy deposit
    if (edep > fPeakEnergy) {
        fPeakEnergy = edep;
    }
    
    // Create detailed hit information
    HPGeHitInfo hit;
    hit.energy = edep;
    hit.time = step->GetPreStepPoint()->GetGlobalTime();
    hit.position = step->GetPreStepPoint()->GetPosition();
    
    // Track information
    G4Track* track = step->GetTrack();
    hit.trackID = track->GetTrackID();
    hit.parentID = track->GetParentID();
    hit.particleName = track->GetDefinition()->GetParticleName();
    
    // Process information
    const G4VProcess* process = step->GetPostStepPoint()->GetProcessDefinedStep();
    if (process) {
        hit.processName = process->GetProcessName();
    } else {
        hit.processName = "Primary";
    }
    
    fHits.push_back(hit);
    
    return true;
}

// ==================================================================
// End of event
// ==================================================================
void HPGeSensitiveDetector::EndOfEvent(G4HCofThisEvent*) {
    // Calculate multiplicity (number of distinct interactions)
    // For simplicity, count number of hits
    fMultiplicity = fHits.size();
    
    // Data is now ready to be accessed by EventAction
}

// ==================================================================
// Reset for next event
// ==================================================================
void HPGeSensitiveDetector::Reset() {
    fTotalEdep = 0.;
    fEventTime = 0.;
    fFirstHit = true;
    fPeakEnergy = 0.;
    fMultiplicity = 0;
    fHits.clear();
}

// ==================================================================
// Get energy deposits
// ==================================================================
std::vector<G4double> HPGeSensitiveDetector::GetEnergyDeposits() const {
    std::vector<G4double> energies;
    for (const auto& hit : fHits) {
        energies.push_back(hit.energy);
    }
    return energies;
}

// ==================================================================
// Get peak energy
// ==================================================================
G4double HPGeSensitiveDetector::GetPeakEnergy() const {
    return fPeakEnergy;
}

// ==================================================================
// Get multiplicity
// ==================================================================
G4int HPGeSensitiveDetector::GetMultiplicity() const {
    return fMultiplicity;
}
