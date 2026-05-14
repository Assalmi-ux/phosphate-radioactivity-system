// ==================================================================
// HPGeSensitiveDetector.hh
// Professional Sensitive Detector for HPGe
// 
// Features:
// - Energy deposition tracking
// - Time information
// - Event-by-event statistics
// - Thread-safe data collection
// ==================================================================

#ifndef HPGeSensitiveDetector_h
#define HPGeSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4TouchableHistory.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <map>

// ==================================================================
// Hit information structure
// ==================================================================
struct HPGeHitInfo {
    G4double energy;          // Energy deposited
    G4double time;            // Global time of hit
    G4ThreeVector position;   // Position of hit
    G4int trackID;            // Track ID
    G4int parentID;           // Parent track ID
    G4String particleName;    // Particle type
    G4String processName;     // Process that created the hit
};

// ==================================================================
// HPGeSensitiveDetector class
// ==================================================================
class HPGeSensitiveDetector : public G4VSensitiveDetector {
public:
    HPGeSensitiveDetector(const G4String& name, G4int detectorID);
    virtual ~HPGeSensitiveDetector();
    
    // Mandatory methods from G4VSensitiveDetector
    virtual void Initialize(G4HCofThisEvent* hce);
    virtual G4bool ProcessHits(G4Step* step, G4TouchableHistory* history);
    virtual void EndOfEvent(G4HCofThisEvent* hce);
    
    // Access accumulated data for this event
    G4double GetTotalEnergyDeposit() const { return fTotalEdep; }
    G4double GetEventTime() const { return fEventTime; }
    G4int GetDetectorID() const { return fDetectorID; }
    G4int GetNumberOfHits() const { return fHits.size(); }
    
    // Access detailed hit information
    const std::vector<HPGeHitInfo>& GetHits() const { return fHits; }
    
    // Get energy spectrum within event (for Compton tracking)
    std::vector<G4double> GetEnergyDeposits() const;
    
    // Reset for next event
    void Reset();
    
    // Statistics
    G4double GetPeakEnergy() const;  // Highest energy deposit in event
    G4int GetMultiplicity() const;    // Number of separate interactions
    
private:
    G4int fDetectorID;
    G4double fTotalEdep;             // Total energy deposited in event
    G4double fEventTime;             // Time of first hit
    G4bool fFirstHit;
    
    std::vector<HPGeHitInfo> fHits;  // Detailed hit information
    
    // Statistics for event
    G4double fPeakEnergy;
    G4int fMultiplicity;
};

#endif
