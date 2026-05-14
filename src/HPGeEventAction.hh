// ==================================================================
// HPGeEventAction.hh
// Professional Event Action with Efficiency Calculation
// ==================================================================

#ifndef HPGeEventAction_h
#define HPGeEventAction_h 1

#include "G4UserEventAction.hh"
#include "G4Event.hh"
#include "HPGeSensitiveDetector.hh"
#include "HPGeTCSManager.hh"
#include <vector>
#include <fstream>
#include <map>
#include <string>

// ==================================================================
// Event data structure for output
// ==================================================================
struct EventData {
    G4int eventID;
    G4bool isCoincidence;
    G4int nDetectors;                      
    std::vector<G4int> detectorIDs;
    std::vector<G4double> energies;        
    std::vector<G4double> trueEnergies;    
    std::vector<G4double> times;
    G4bool isFullEnergy;                   
    G4double comptonFraction;              
};

// ==================================================================
// Efficiency tracking structure
// ==================================================================
struct EfficiencyData {
    G4double energy;                       
    G4long totalEmitted;                   
    G4long totalDetected;                  
    G4long fullEnergyDetected;             
    G4long comptonDetected;                
    G4long escapePeakDetected;             
    
    G4double totalEfficiency;              
    G4double peakEfficiency;               
    G4double peakToTotal;                  
    
    EfficiencyData() : 
        energy(0), totalEmitted(0), totalDetected(0),
        fullEnergyDetected(0), comptonDetected(0), escapePeakDetected(0),
        totalEfficiency(0), peakEfficiency(0), peakToTotal(0) {}
};

// ==================================================================
// HPGeEventAction class
// ==================================================================
class HPGeEventAction : public G4UserEventAction {
public:
    // Constructor taking the vector of detectors (Matches your .cc)
    HPGeEventAction(const std::vector<HPGeSensitiveDetector*>& detectors);
    virtual ~HPGeEventAction();
    
    // Mandatory methods
    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);
    
    // --- Configuration Methods (Called by main.cc) ---
    void SetOutputFile(const std::string& filename);
    
    // Getters and Setters for Physics/Analysis parameters
    inline void SetCoincidenceWindow(G4double window) { fCoincidenceWindow = window; }
    inline void SetEnergyResolution(G4double resolution) { fEnergyResolution = resolution; }
    inline void SetEnergyThreshold(G4double threshold) { fEnergyThreshold = threshold; }
    inline void SetFullEnergyWindow(G4double window) { fFullEnergyWindow = window; }
    
    // TCS Support (These were causing your errors!)
    inline void EnableTCS(G4bool enable) { fTCSEnabled = enable; }
    inline G4bool IsTCSEnabled() const { return fTCSEnabled; } // Fixed Error
    inline void SetTCSManager(HPGeTCSManager* manager) { fTCSManager = manager; }
    inline HPGeTCSManager* GetTCSManager() const { return fTCSManager; }

    // *** NEW: Set actual emission counts from PrimaryGenerator ***
    void SetEmissionCounts(const std::map<G4double, G4long>& counts);

    // Efficiency Setup
    void SetExpectedEnergies(const std::vector<G4double>& energies,
                            const std::vector<G4double>& intensities);
    
    const std::map<G4double, EfficiencyData>& GetEfficiencyData() const {
        return fEfficiencyMap;
    }
    
    // Output Methods
    void PrintStatistics() const;
    void ExportEfficiencyData(const std::string& filename) const;
    
private:
    // --- Helper Internal Methods (Matches your .cc) ---
    G4double ApplyEnergyResolution(G4double energy);
    G4bool IsFullEnergyPeak(G4double measuredEnergy, G4double trueEnergy);
    G4bool CheckCoincidence(const std::vector<G4double>& times);
    G4double FindClosestExpectedEnergy(G4double energy);
    void UpdateEfficiency(G4double energy, G4bool detected, G4bool fullEnergy);
    void WriteEvent(const EventData& eventData);
    
    // --- Member Variables ---
    std::vector<HPGeSensitiveDetector*> fDetectors;
    std::ofstream fOutputFile;
    
    // Configuration parameters
    G4double fCoincidenceWindow;     
    G4double fEnergyResolution;      
    G4double fEnergyThreshold;       
    G4double fFullEnergyWindow;      
    
    // Expected energies for efficiency calculation
    std::vector<G4double> fExpectedEnergies;
    std::vector<G4double> fExpectedIntensities;
    
    // Efficiency tracking
    std::map<G4double, EfficiencyData> fEfficiencyMap;
    
    // Statistics
    G4long fTotalEvents;
    G4long fDetectedEvents;
    G4long fCoincidenceEvents;
    G4long fFullEnergyEvents;
    G4long fComptonEvents;
    
    // TCS Support
    G4bool fTCSEnabled;
    HPGeTCSManager* fTCSManager;
};

#endif
