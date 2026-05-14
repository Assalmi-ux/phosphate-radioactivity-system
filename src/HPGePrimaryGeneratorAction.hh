// ==================================================================
// HPGePrimaryGeneratorAction.hh - WITH MARINELLI & CASCADE SUPPORT
// Professional Primary Generator - ISOCS Style
// 
// Features:
// - Isotropic 4π emission
// - Multiple gamma lines with proper intensities
// - Standard isotope library
// - Custom energy support
// - Realistic source geometries (point, disk, volume, MARINELLI)
// - Activity-based normalization
// - TRUE COINCIDENCE SUMMING: Cascade emission support
// ==================================================================

#ifndef HPGePrimaryGeneratorAction_h
#define HPGePrimaryGeneratorAction_h 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include <vector>
#include <string>
#include <map>

// ==================================================================
// Gamma line structure
// ==================================================================
struct GammaLine {
    G4double energy;       // Energy in keV
    G4double intensity;    // Relative intensity (per 100 decays)
    G4String name;         // Line identification
    
    GammaLine(G4double e = 0, G4double i = 0, G4String n = "") :
        energy(e), intensity(i), name(n) {}
};

// ==================================================================
// HPGePrimaryGeneratorAction class
// ==================================================================
class HPGePrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
public:
    HPGePrimaryGeneratorAction();
    virtual ~HPGePrimaryGeneratorAction();
    
    // Mandatory method from G4VUserPrimaryGeneratorAction
    virtual void GeneratePrimaries(G4Event* anEvent);
    
    // Source configuration
    void SetSourcePosition(G4ThreeVector pos) { fSourcePosition = pos; }
    void SetSourceType(const std::string& type);  // "point", "disk", "volume", "marinelli", "cartridge", "filter"
    void SetSourceRadius(G4double radius) { fSourceRadius = radius; }
    void SetSourceHeight(G4double height) { fSourceHeight = height; }
    
    // MARINELLI-specific configuration
    void EnableMarinelli(bool enable) { fMarinelliEnabled = enable; }
    void SetMarinelliGeometry(G4double outerD, G4double innerD, 
                             G4double totalH, G4double wellH);
    
    // Gamma line configuration
    void AddGammaLine(G4double energy, G4double intensity, G4String name = "");
    void ClearGammaLines() { fGammaLines.clear(); }
    void SetIsotope(const std::string& isotope);
    
    // Getters
    G4ThreeVector GetSourcePosition() const { return fSourcePosition; }
    const std::vector<GammaLine>& GetGammaLines() const { return fGammaLines; }
    G4int GetNumberOfGammaLines() const { return fGammaLines.size(); }
    G4String GetCurrentIsotope() const { return fCurrentIsotope; }
    
    // Statistics
    G4long GetTotalPrimariesGenerated() const { return fTotalPrimaries; }
    std::map<G4double, G4long> GetEnergyStatistics() const { return fEnergyStats; }
    
    // Print information
    void PrintConfiguration() const;
    
    // TCS/Cascade mode control
    void SetCascadeMode(G4bool enable) { fCascadeMode = enable; }
    G4bool GetCascadeMode() const { return fCascadeMode; }
    
private:
    // Sampling methods
    G4ThreeVector SampleSourcePosition();
    G4ThreeVector SampleIsotropicDirection();
    G4double SelectGammaEnergy();
    
    // MARINELLI sampling methods
    G4ThreeVector SampleMarinelliPosition();
    
    // Initialize standard isotopes
    void InitializeIsotopeLibrary();
    
    // CASCADE EMISSION - NEW for TCS support
    void EmitSingleGamma(G4Event* anEvent, 
                        G4ThreeVector position,
                        G4double energy);
    
    void EmitCascade(G4Event* anEvent,
                    G4ThreeVector position,
                    std::vector<G4double> energies,
                    std::vector<G4double> times = {});
    
    // Member variables
    G4ParticleGun* fParticleGun;
    
    // Source geometry
    G4ThreeVector fSourcePosition;
    std::string fSourceType;
    G4double fSourceRadius;
    G4double fSourceHeight;
    
    // MARINELLI geometry parameters
    bool fMarinelliEnabled;
    G4double fMarinelliOuterD;    // Outer diameter
    G4double fMarinelliInnerD;    // Inner diameter (well)
    G4double fMarinelliTotalH;    // Total height
    G4double fMarinelliWellH;     // Well height
    
    // Gamma lines
    std::vector<GammaLine> fGammaLines;
    
    // Current isotope name - NEW for TCS support
    G4String fCurrentIsotope;
    
    // TCS/Cascade mode flag - TRUE = emit cascades, FALSE = single gammas
    G4bool fCascadeMode;
    
    // Statistics
    G4long fTotalPrimaries;
    std::map<G4double, G4long> fEnergyStats;  // Energy -> count
    
    // Cumulative distribution for energy sampling
    std::vector<G4double> fCumulativeIntensities;
    void UpdateCumulativeDistribution();
};

#endif
