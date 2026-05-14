// ==================================================================
// HPGeTCSManager.hh
// True Coincidence Summing Manager - Enhanced Version
// 
// NEW Features:
// - Dynamic solid angle calculation based on geometry
// - Support for efficiency-based COI calculation
// - Geometry-aware TCS corrections
// ==================================================================

#ifndef HPGeTCSManager_h
#define HPGeTCSManager_h 1

#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include <vector>
#include <map>
#include <string>

// ==================================================================
// Geometry configuration for solid angle calculation
// ==================================================================
struct GeometryConfig {
    G4String sourceType;              // "Point", "Disk", "Cylinder", "Marinelli_1000mL", "Marinelli_450mL"
    G4double detectorRadius;          // Detector crystal radius
    G4double detectorHeight;          // Detector crystal height
    G4double sourceDistance;          // Distance from detector face (for point sources)
    
    // Disk/Cylinder source parameters
    G4double sourceRadius;            // Radius of disk or cylinder source
    G4double sourceHeight;            // Height/thickness of disk or cylinder source
    
    // Marinelli-specific parameters
    G4double marinelliInnerRadius;    // Inner radius of Marinelli beaker
    G4double marinelliOuterRadius;    // Outer radius of Marinelli beaker
    G4double marinelliHeight;         // Height of active volume
    G4double marinelliInnerHeight;    // Height of inner well
    
    GeometryConfig() :
        sourceType("Point"),
        detectorRadius(30.0*mm),
        detectorHeight(60.0*mm),
        sourceDistance(50.0*mm),
        sourceRadius(25.0*mm),
        sourceHeight(10.0*mm),
        marinelliInnerRadius(42.5*mm),
        marinelliOuterRadius(65.0*mm),
        marinelliHeight(77.0*mm),
        marinelliInnerHeight(68.0*mm) {}
};

// ==================================================================
// Cascade structure for decay schemes
// ==================================================================
struct GammaCascade {
    std::vector<G4double> energies;        // Energies in cascade
    std::vector<G4double> times;           // Time delays between emissions
    G4double probability;                  // Cascade probability
    G4String parentIsotope;                // Parent isotope
    
    GammaCascade() : probability(0.0) {}
};

// ==================================================================
// TCS Event structure
// ==================================================================
struct TCSEvent {
    G4int eventID;
    std::vector<G4double> primaryEnergies;    // Original gamma energies
    std::vector<G4double> detectedEnergies;   // Detected energies (may be summed)
    std::vector<G4double> times;              // Detection times
    std::vector<G4bool> isSummed;             // Flag for summed peaks
    G4double sumEnergy;                       // Total sum energy if coincident
    G4bool isCoincidence;                     // True if coincidence detected
};

// ==================================================================
// TCS Statistics structure (Enhanced)
// ==================================================================
struct TCSStatistics {
    G4double energy;                       // Gamma energy
    G4long totalEmitted;                   // Total gammas emitted
    G4long singleDetected;                 // Single gamma detected
    G4long coincidenceDetected;            // Detected in coincidence
    G4long sumPeakDetected;                // Detected as sum peak
    G4double COI;                          // Coincidence correction factor
    
    // New: For efficiency-based COI calculation
    G4long detectedWithTCS;                // Counts with TCS enabled (Mode A)
    G4long detectedNoTCS;                  // Counts with TCS disabled (Mode B)
    G4double efficiencyTCS;                // Efficiency with TCS (ε_TCS)
    G4double efficiencyNoTCS;              // Efficiency without TCS (ε_NoTCS)
    G4double COI_efficiency;               // COI = ε_NoTCS / ε_TCS
    
    TCSStatistics() : 
        energy(0), totalEmitted(0), singleDetected(0),
        coincidenceDetected(0), sumPeakDetected(0), COI(1.0),
        detectedWithTCS(0), detectedNoTCS(0),
        efficiencyTCS(0), efficiencyNoTCS(0), COI_efficiency(1.0) {}
};

// ==================================================================
// HPGeTCSManager class
// ==================================================================
class HPGeTCSManager {
public:
    HPGeTCSManager();
    ~HPGeTCSManager();
    
    // Configuration
    void SetCoincidenceWindow(G4double window) { fCoincidenceWindow = window; }
    void SetDetectorSolidAngle(G4double angle) { fSolidAngle = angle; }
    void SetGeometryConfig(const GeometryConfig& config);
    
    // Calculate solid angle from geometry
    G4double CalculateSolidAngle() const;
    G4double CalculateSolidAnglePoint() const;
    G4double CalculateSolidAngleDisk() const;
    G4double CalculateSolidAngleCylinder() const;
    G4double CalculateSolidAngleMarinelli() const;
    G4double GetSolidAngle() const { return fSolidAngle; }
    
    // Cascade definition
    void AddCascade(const GammaCascade& cascade);
    void LoadIsotopeScheme(const G4String& isotope);
    void ClearCascades() { fCascades.clear(); }
    
    // Process events
    TCSEvent ProcessEvent(const std::vector<G4double>& energies,
                         const std::vector<G4double>& times);
    
    // Calculate corrections
    G4double CalculateCOI(G4double energy);
    G4double GetCoincidenceProbability(G4double energy1, G4double energy2);
    
    // Statistics
    void UpdateStatistics(const TCSEvent& event);
    void UpdateEfficiencyStatistics(G4double energy, G4long detected, G4bool isTCSMode);
    const std::map<G4double, TCSStatistics>& GetStatistics() const { 
        return fStatistics; 
    }
    void PrintStatistics() const;
    void ExportTCSData(const G4String& filename) const;
    void ExportEfficiencyData(const G4String& filename) const;
    
    // Check if energies are from cascade
    G4bool IsFromCascade(G4double energy1, G4double energy2) const;
    std::vector<G4double> GetCascadeEnergies(const G4String& isotope) const;
    
private:
    // Find cascades containing specific energy
    std::vector<GammaCascade> FindCascadesWithEnergy(G4double energy) const;
    
    // Check if two events are coincident
    G4bool AreCoincident(G4double time1, G4double time2) const;
    
    // Calculate geometric efficiency for cascade
    G4double CalculateGeometricEfficiency(const GammaCascade& cascade) const;
    
    // Member variables
    G4double fCoincidenceWindow;           // Time window for coincidence
    G4double fSolidAngle;                  // Detector solid angle (calculated or set)
    GeometryConfig fGeometryConfig;        // Geometry configuration
    
    std::vector<GammaCascade> fCascades;   // Defined cascades
    std::map<G4double, TCSStatistics> fStatistics;  // TCS statistics
    
    // Pre-calculated coincidence probabilities
    std::map<std::pair<G4double, G4double>, G4double> fCoincidenceProbabilities;
    
    // Standard isotope schemes
    void LoadCo60Scheme();
    void LoadY88Scheme();
    void LoadEu152Scheme();
    void LoadBa133Scheme();
    void LoadNa22Scheme();
};

#endif
