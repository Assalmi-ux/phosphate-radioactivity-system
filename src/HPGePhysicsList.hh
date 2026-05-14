// ==================================================================
// HPGePhysicsList.hh
// Professional Physics List for High-Precision Gamma Spectroscopy
// 
// Features:
// - Livermore low-energy models for accurate transport
// - Fluorescence, Auger, PIXE
// - Rayleigh scattering
// - Atomic deexcitation
// - Optimal cuts for germanium detector
// ==================================================================

#ifndef HPGePhysicsList_h
#define HPGePhysicsList_h 1

#include "G4VModularPhysicsList.hh"

// ==================================================================
// HPGePhysicsList class
// ==================================================================
class HPGePhysicsList: public G4VModularPhysicsList {
public:
    HPGePhysicsList(G4bool fastMode = false);
    virtual ~HPGePhysicsList();

    // Mandatory method from G4VModularPhysicsList
    virtual void SetCuts();
    
    // Optional: set custom cuts
    void SetCutForGamma(G4double cut);
    void SetCutForElectron(G4double cut);
    void SetCutForPositron(G4double cut);
    
    // Fast mode control
    void SetFastMode(G4bool fast) { fFastMode = fast; }
    G4bool IsFastMode() const { return fFastMode; }
    
private:
    G4double fCutForGamma;
    G4double fCutForElectron;
    G4double fCutForPositron;
    G4bool fFastMode;
};

#endif
