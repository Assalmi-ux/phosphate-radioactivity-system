// ==================================================================
// HPGePhysicsList.cc
// Professional Physics List for Precision Gamma Spectroscopy
// ==================================================================

#include "HPGePhysicsList.hh"


// Geant4 physics lists
#include "G4EmLivermorePhysics.hh"
#include "G4EmPenelopePhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4EmParameters.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4SystemOfUnits.hh"

// ==================================================================
// Constructor - with optional fast mode
// ==================================================================
HPGePhysicsList::HPGePhysicsList(G4bool fastMode) 
: G4VModularPhysicsList(),
  fCutForGamma(0.1*mm),
  fCutForElectron(0.1*mm),
  fCutForPositron(0.1*mm),
  fFastMode(fastMode) {
    
    SetVerboseLevel(1);

    // ================================================================
    // EM Physics - Choose based on fast mode
    // ================================================================
    if (fFastMode) {
        // FAST MODE: Use Penelope (faster than Livermore, still good)
        RegisterPhysics(new G4EmPenelopePhysics());
        
        // Increase cuts for speed
        fCutForGamma = 0.1*mm;
        fCutForElectron = 1.0*mm;
        fCutForPositron = 1.0*mm;
        
        G4cout << "\n*** FAST MODE ENABLED ***" << G4endl;
        G4cout << "Using Penelope physics with relaxed cuts" << G4endl;
    } else {
        // PRECISION MODE: Use Livermore (best for low energy)
        RegisterPhysics(new G4EmLivermorePhysics());
    }
    
    // ================================================================
    // Extra EM processes
    // ================================================================
    RegisterPhysics(new G4EmExtraPhysics());

    // ================================================================
    // Decay Physics
    // ================================================================
    RegisterPhysics(new G4DecayPhysics());

    // ================================================================
    // Radioactive Decay
    // ================================================================
    RegisterPhysics(new G4RadioactiveDecayPhysics());

    // ================================================================
    // Configure EM parameters
    // ================================================================
    G4EmParameters* param = G4EmParameters::Instance();
    param->SetDefaults();
    
    // Energy range
    param->SetMinEnergy(100*eV);
    param->SetMaxEnergy(10*GeV);
    
    // Number of bins - fewer for fast mode
    if (fFastMode) {
        param->SetNumberOfBinsPerDecade(10);  // Coarser binning
    } else {
        param->SetNumberOfBinsPerDecade(20);  // Fine binning
    }
    
    param->SetVerbose(0);
    
    // ================================================================
    // Physics processes - disable some in fast mode
    // ================================================================
    if (fFastMode) {
        param->SetFluo(false);             // Disable fluorescence
        param->SetAuger(false);            // Disable Auger
        param->SetAugerCascade(false);     // Disable Auger cascade
        param->SetPixe(false);             // Disable PIXE
        param->SetDeexcitationIgnoreCut(false);
    } else {
        param->SetFluo(true);
        param->SetAuger(true);
        param->SetAugerCascade(true);
        param->SetPixe(true);
        param->SetDeexcitationIgnoreCut(true);
    }
    
    param->SetMscRangeFactor(0.04);
    
    G4cout << "\n=== Physics List Configuration ===" << G4endl;
    G4cout << "Mode: " << (fFastMode ? "FAST" : "PRECISION") << G4endl;
    G4cout << "EM Model: " << (fFastMode ? "Penelope" : "Livermore") << G4endl;
    G4cout << "Energy range: " << param->MinKinEnergy()/eV << " eV - " 
           << param->MaxKinEnergy()/GeV << " GeV" << G4endl;
    G4cout << "Fluorescence: " << (param->Fluo() ? "ON" : "OFF") << G4endl;
    G4cout << "Auger: " << (param->Auger() ? "ON" : "OFF") << G4endl;
    G4cout << "PIXE: " << (param->Pixe() ? "ON" : "OFF") << G4endl;
    G4cout << "Production cuts: " << fCutForGamma/mm << " mm" << G4endl;
    G4cout << "===================================" << G4endl;
}

// ==================================================================
// Destructor
// ==================================================================
HPGePhysicsList::~HPGePhysicsList() {}

// ==================================================================
// Set production cuts
// ==================================================================
void HPGePhysicsList::SetCuts() {
    // Set default cuts
    SetCutsWithDefault();

    // ================================================================
    // Set specific cuts for high precision
    // ================================================================
    // Production cuts define the threshold for producing secondary particles
    // Smaller cuts = more accurate but slower simulation
    // 
    // For HPGe detectors:
    // - gamma: 0.1 mm (good precision without excessive secondaries)
    // - e-/e+: 0.1 mm (track all significant electrons)
    // ================================================================
    
    SetCutValue(fCutForGamma, "gamma");
    SetCutValue(fCutForElectron, "e-");
    SetCutValue(fCutForPositron, "e+");
    SetCutValue(0.1*mm, "proton");
    
    G4cout << "\n=== Production Cuts ===" << G4endl;
    G4cout << "Gamma:      " << fCutForGamma/mm << " mm" << G4endl;
    G4cout << "Electron:   " << fCutForElectron/mm << " mm" << G4endl;
    G4cout << "Positron:   " << fCutForPositron/mm << " mm" << G4endl;
    G4cout << "Proton:     " << 0.1 << " mm" << G4endl;
    G4cout << "=======================" << G4endl;
    
    if (verboseLevel > 0) {
        DumpCutValuesTable();
    }
}

// ==================================================================
// Set custom cuts
// ==================================================================
void HPGePhysicsList::SetCutForGamma(G4double cut) {
    fCutForGamma = cut;
    SetCutValue(cut, "gamma");
}

void HPGePhysicsList::SetCutForElectron(G4double cut) {
    fCutForElectron = cut;
    SetCutValue(cut, "e-");
}

void HPGePhysicsList::SetCutForPositron(G4double cut) {
    fCutForPositron = cut;
    SetCutValue(cut, "e+");
}
