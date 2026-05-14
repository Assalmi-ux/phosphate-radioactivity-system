// ==================================================================
// HPGePrimaryGeneratorAction.cc - WITH MARINELLI & CASCADE SUPPORT
// Professional Primary Generator Implementation with TCS Support
// FIXED VERSION: Custom mode now supports cascade emission!
// ==================================================================

#include "HPGePrimaryGeneratorAction.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include <numeric>
#include <algorithm>

// ==================================================================
// Constructor
// ==================================================================
HPGePrimaryGeneratorAction::HPGePrimaryGeneratorAction()
: G4VUserPrimaryGeneratorAction(),
  fSourcePosition(0., 0., 5.*cm),
  fSourceType("point"),
  fSourceRadius(1.*cm),
  fSourceHeight(2.*cm),
  fMarinelliEnabled(false),
  fMarinelliOuterD(0),
  fMarinelliInnerD(0),
  fMarinelliTotalH(0),
  fMarinelliWellH(0),
  fCurrentIsotope("Cs137"),
  fCascadeMode(true),  // Default: cascade mode ON
  fTotalPrimaries(0) {
    
    fParticleGun = new G4ParticleGun(1);
    
    // Set default particle type to gamma
    G4ParticleDefinition* particle = G4Gamma::GammaDefinition();
    fParticleGun->SetParticleDefinition(particle);
    
    // Initialize isotope library
    InitializeIsotopeLibrary();
    
    // Default: Cs-137 source
    SetIsotope("Cs137");
}

// ==================================================================
// Destructor
// ==================================================================
HPGePrimaryGeneratorAction::~HPGePrimaryGeneratorAction() {
    delete fParticleGun;
}

// ==================================================================
// Helper: Emit a single gamma
// ==================================================================
void HPGePrimaryGeneratorAction::EmitSingleGamma(G4Event* anEvent, 
                                                  G4ThreeVector position,
                                                  G4double energy) {
    fParticleGun->SetParticlePosition(position);
    fParticleGun->SetParticleEnergy(energy);
    fParticleGun->SetParticleMomentumDirection(SampleIsotropicDirection());
    fParticleGun->SetParticleTime(0.*ns);
    fParticleGun->GeneratePrimaryVertex(anEvent);
}

// ==================================================================
// Helper: Emit cascade gammas
// ==================================================================
void HPGePrimaryGeneratorAction::EmitCascade(G4Event* anEvent,
                                             G4ThreeVector position,
                                             std::vector<G4double> energies,
                                             std::vector<G4double> times) {
    
    // If no times specified, assume prompt (all at t=0)
    if (times.empty()) {
        times.resize(energies.size(), 0.*ns);
    }
    
    for (size_t i = 0; i < energies.size(); i++) {
        fParticleGun->SetParticlePosition(position);
        fParticleGun->SetParticleEnergy(energies[i]);
        fParticleGun->SetParticleMomentumDirection(SampleIsotropicDirection());
        fParticleGun->SetParticleTime(times[i]);
        fParticleGun->GeneratePrimaryVertex(anEvent);
    }
}

// ==================================================================
// Generate primaries - WITH CASCADE SUPPORT FOR TCS
// ==================================================================
void HPGePrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {
    
    // Sample source position (same for all gammas in cascade)
    G4ThreeVector position = SampleSourcePosition();
    
    // ================================================================
    // CASCADE ISOTOPES - Behavior depends on fCascadeMode
    // fCascadeMode = true:  Emit multiple gammas per event (TCS enabled)
    // fCascadeMode = false: Emit single gamma per event (no TCS)
    // ================================================================
    
    if (fCurrentIsotope == "Co60") {
        // Co60: Two gammas in prompt cascade (99.85% probability)
        // 1173.228 keV → 1332.492 keV
        
        if (fCascadeMode) {
            // CASCADE MODE: Emit BOTH gammas simultaneously
            if (G4UniformRand() < 0.9985) {
                std::vector<G4double> cascadeEnergies = {
                    1173.228*keV, 
                    1332.492*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                
                // Update statistics for BOTH energies
                fEnergyStats[1173.228*keV]++;
                fEnergyStats[1332.492*keV]++;
                fTotalPrimaries += 2;
            }
            else {
                // Rare case: only one gamma
                G4double energy = SelectGammaEnergy();
                EmitSingleGamma(anEvent, position, energy);
                fEnergyStats[energy]++;
                fTotalPrimaries++;
            }
        }
        else {
            // NON-CASCADE MODE: Emit SINGLE gamma per event
            // This simulates measuring each gamma line independently
            // Select based on emission probability
            G4double energy = SelectGammaEnergy();
            EmitSingleGamma(anEvent, position, energy);
            fEnergyStats[energy]++;
            fTotalPrimaries++;
        }
    }
    
    else if (fCurrentIsotope == "Y88") {
        // Y-88: Two gammas in cascade (93.7% probability)
        // 898.042 keV → 1836.063 keV
        
        if (fCascadeMode) {
            // CASCADE MODE
            if (G4UniformRand() < 0.937) {
                std::vector<G4double> cascadeEnergies = {
                    898.042*keV,
                    1836.063*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                
                fEnergyStats[898.042*keV]++;
                fEnergyStats[1836.063*keV]++;
                fTotalPrimaries += 2;
            }
            else {
                G4double energy = SelectGammaEnergy();
                EmitSingleGamma(anEvent, position, energy);
                fEnergyStats[energy]++;
                fTotalPrimaries++;
            }
        }
        else {
            // NON-CASCADE MODE: Single gamma
            G4double energy = SelectGammaEnergy();
            EmitSingleGamma(anEvent, position, energy);
            fEnergyStats[energy]++;
            fTotalPrimaries++;
        }
    }
    
    else if (fCurrentIsotope == "Na22") {
        // Na-22: 1274.537 keV + two 511 keV annihilation gammas
        // 90.44% cascade probability
        
        if (fCascadeMode) {
            // CASCADE MODE
            if (G4UniformRand() < 0.9044) {
                // Emit all three gammas (positron annihilation + nuclear gamma)
                std::vector<G4double> cascadeEnergies = {
                    511.0*keV,    // Annihilation 1
                    511.0*keV,    // Annihilation 2
                    1274.537*keV  // Na-22 gamma
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                
                fEnergyStats[511.0*keV] += 2;  // Two 511 keV
                fEnergyStats[1274.537*keV]++;
                fTotalPrimaries += 3;
            }
            else {
                G4double energy = SelectGammaEnergy();
                EmitSingleGamma(anEvent, position, energy);
                fEnergyStats[energy]++;
                fTotalPrimaries++;
            }
        }
        else {
            // NON-CASCADE MODE: Single gamma
            G4double energy = SelectGammaEnergy();
            EmitSingleGamma(anEvent, position, energy);
            fEnergyStats[energy]++;
            fTotalPrimaries++;
        }
    }
    
    else if (fCurrentIsotope == "Eu152") {
        // Eu-152: Complex decay scheme with multiple cascades
        // Simplified: Main cascades only
        
        if (fCascadeMode) {
            // CASCADE MODE
            G4double rand = G4UniformRand();
            
            if (rand < 0.2654) {
                // 121-344 keV cascade (26.54%)
                std::vector<G4double> cascadeEnergies = {
                    121.782*keV,
                    344.279*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                fEnergyStats[121.782*keV]++;
                fEnergyStats[344.279*keV]++;
                fTotalPrimaries += 2;
            }
            else if (rand < 0.2654 + 0.0751) {
                // 121-244 keV cascade (7.51%)
                std::vector<G4double> cascadeEnergies = {
                    121.782*keV,
                    244.697*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                fEnergyStats[121.782*keV]++;
                fEnergyStats[244.697*keV]++;
                fTotalPrimaries += 2;
            }
            else if (rand < 0.2654 + 0.0751 + 0.0227) {
                // 344-411 keV cascade (2.27%)
                std::vector<G4double> cascadeEnergies = {
                    344.279*keV,
                    411.117*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                fEnergyStats[344.279*keV]++;
                fEnergyStats[411.117*keV]++;
                fTotalPrimaries += 2;
            }
            else {
                // Single gamma emission (remaining probability)
                G4double energy = SelectGammaEnergy();
                EmitSingleGamma(anEvent, position, energy);
                fEnergyStats[energy]++;
                fTotalPrimaries++;
            }
        }
        else {
            // NON-CASCADE MODE: Single gamma
            G4double energy = SelectGammaEnergy();
            EmitSingleGamma(anEvent, position, energy);
            fEnergyStats[energy]++;
            fTotalPrimaries++;
        }
    }
    
    else if (fCurrentIsotope == "Ba133") {
        // Ba-133: Two main cascades
        
        if (fCascadeMode) {
            // CASCADE MODE
            G4double rand = G4UniformRand();
            
            if (rand < 0.329) {
                // 81-356 keV cascade (32.9%)
                std::vector<G4double> cascadeEnergies = {
                    80.998*keV,
                    356.013*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                fEnergyStats[80.998*keV]++;
                fEnergyStats[356.013*keV]++;
                fTotalPrimaries += 2;
            }
            else if (rand < 0.329 + 0.183) {
                // 81-302 keV cascade (18.3%)
                std::vector<G4double> cascadeEnergies = {
                    80.998*keV,
                    302.851*keV
                };
                EmitCascade(anEvent, position, cascadeEnergies);
                fEnergyStats[80.998*keV]++;
                fEnergyStats[302.851*keV]++;
                fTotalPrimaries += 2;
            }
            else {
                // Single gamma
                G4double energy = SelectGammaEnergy();
                EmitSingleGamma(anEvent, position, energy);
                fEnergyStats[energy]++;
                fTotalPrimaries++;
            }
        }
        else {
            // NON-CASCADE MODE: Single gamma
            G4double energy = SelectGammaEnergy();
            EmitSingleGamma(anEvent, position, energy);
            fEnergyStats[energy]++;
            fTotalPrimaries++;
        }
    }
    
    // ================================================================
    // CUSTOM MODE - Enhanced with CASCADE DETECTION
    // FIXED: Now respects fCascadeMode flag!
    // ================================================================
    
    else if (fCurrentIsotope == "Custom") {
        // Check if custom energies match known cascade patterns
        // This allows Custom mode to emit cascades for TCS analysis
        
        if (fCascadeMode && fGammaLines.size() >= 2) {
            // CASCADE MODE ENABLED - Check for cascade patterns
            bool isCo60Pair = false;
            bool isY88Pair = false;
            
            // Check if first two energies match Co-60 cascade
            if (fGammaLines.size() >= 2) {
                G4double e1 = fGammaLines[0].energy;
                G4double e2 = fGammaLines[1].energy;
                
                if ((std::abs(e1 - 1173.228) < 1.0 && std::abs(e2 - 1332.492) < 1.0) ||
                    (std::abs(e1 - 1332.492) < 1.0 && std::abs(e2 - 1173.228) < 1.0)) {
                    isCo60Pair = true;
                }
            }
            
            // Check if energies 2-3 (indices) match Y-88 cascade
            if (fGammaLines.size() >= 4) {
                G4double e3 = fGammaLines[2].energy;
                G4double e4 = fGammaLines[3].energy;
                
                if ((std::abs(e3 - 898.042) < 1.0 && std::abs(e4 - 1836.063) < 1.0) ||
                    (std::abs(e3 - 1836.063) < 1.0 && std::abs(e4 - 898.042) < 1.0)) {
                    isY88Pair = true;
                }
            }
            
            // Emit cascades if recognized patterns found
            if (isCo60Pair && isY88Pair) {
                // Both Co-60 and Y-88 present - randomly choose one per event
                G4double rand = G4UniformRand();
                
                if (rand < 0.5) {
                    // Emit Co-60 cascade (1173 + 1332 keV)
                    if (G4UniformRand() < 0.9985) {
                        std::vector<G4double> cascadeEnergies = {
                            1173.228*keV, 
                            1332.492*keV
                        };
                        EmitCascade(anEvent, position, cascadeEnergies);
                        fEnergyStats[1173.228*keV]++;
                        fEnergyStats[1332.492*keV]++;
                        fTotalPrimaries += 2;
                    } else {
                        // Rare single gamma
                        G4double energy = SelectGammaEnergy();
                        EmitSingleGamma(anEvent, position, energy);
                        fEnergyStats[energy]++;
                        fTotalPrimaries++;
                    }
                } else {
                    // Emit Y-88 cascade (898 + 1836 keV)
                    if (G4UniformRand() < 0.937) {
                        std::vector<G4double> cascadeEnergies = {
                            898.042*keV,
                            1836.063*keV
                        };
                        EmitCascade(anEvent, position, cascadeEnergies);
                        fEnergyStats[898.042*keV]++;
                        fEnergyStats[1836.063*keV]++;
                        fTotalPrimaries += 2;
                    } else {
                        // Rare single gamma
                        G4double energy = SelectGammaEnergy();
                        EmitSingleGamma(anEvent, position, energy);
                        fEnergyStats[energy]++;
                        fTotalPrimaries++;
                    }
                }
            }
            else if (isCo60Pair) {
                // Only Co-60 cascade detected
                if (G4UniformRand() < 0.9985) {
                    std::vector<G4double> cascadeEnergies = {
                        1173.228*keV, 
                        1332.492*keV
                    };
                    EmitCascade(anEvent, position, cascadeEnergies);
                    fEnergyStats[1173.228*keV]++;
                    fEnergyStats[1332.492*keV]++;
                    fTotalPrimaries += 2;
                } else {
                    G4double energy = SelectGammaEnergy();
                    EmitSingleGamma(anEvent, position, energy);
                    fEnergyStats[energy]++;
                    fTotalPrimaries++;
                }
            }
            else if (isY88Pair) {
                // Only Y-88 cascade detected
                if (G4UniformRand() < 0.937) {
                    std::vector<G4double> cascadeEnergies = {
                        898.042*keV,
                        1836.063*keV
                    };
                    EmitCascade(anEvent, position, cascadeEnergies);
                    fEnergyStats[898.042*keV]++;
                    fEnergyStats[1836.063*keV]++;
                    fTotalPrimaries += 2;
                } else {
                    G4double energy = SelectGammaEnergy();
                    EmitSingleGamma(anEvent, position, energy);
                    fEnergyStats[energy]++;
                    fTotalPrimaries++;
                }
            }
            else {
                // No recognized cascade pattern - emit single gamma
                static bool warningPrinted = false;
                if (!warningPrinted) {
                    G4cout << "\n*** WARNING: Custom energies don't match known cascade patterns ***" << G4endl;
                    G4cout << "    Emitting single gammas. TCS will not work properly." << G4endl;
                    G4cout << "    For TCS, use 'isotope = Co60' or 'isotope = Y88' instead." << G4endl;
                    G4cout << "    Or ensure custom energies match: 1173.228, 1332.492 (Co60)" << G4endl;
                    G4cout << "                                    898.042, 1836.063 (Y88)\n" << G4endl;
                    warningPrinted = true;
                }
                G4double energy = SelectGammaEnergy();
                EmitSingleGamma(anEvent, position, energy);
                fEnergyStats[energy]++;
                fTotalPrimaries++;
            }
        }
        else {
            // NON-CASCADE MODE or less than 2 energies - emit single gamma
            // This is the FIX: Custom mode now respects fCascadeMode = false!
            G4double energy = SelectGammaEnergy();
            EmitSingleGamma(anEvent, position, energy);
            fEnergyStats[energy]++;
            fTotalPrimaries++;
        }
    }
    
    // ================================================================
    // OTHER SINGLE-GAMMA ISOTOPES - Emit one gamma per event
    // ================================================================
    
    else {
        // Cs-137, Mn-54, Am-241, and unknown isotopes
        // These have no significant cascade gammas
        G4double energy = SelectGammaEnergy();
        EmitSingleGamma(anEvent, position, energy);
        fEnergyStats[energy]++;
        fTotalPrimaries++;
    }
}

// ==================================================================
// Sample source position
// ==================================================================
G4ThreeVector HPGePrimaryGeneratorAction::SampleSourcePosition() {
    
    if (fSourceType == "point") {
        // Point source
        return fSourcePosition;
    }
    else if (fSourceType == "disk") {
        // Disk source - uniform distribution in disk with thickness
        // Sample radius: r = R * sqrt(U) for uniform area distribution
        G4double r = fSourceRadius * std::sqrt(G4UniformRand());
        G4double phi = twopi * G4UniformRand();
        
        // Sample height within disk thickness (centered at source position)
        G4double h = fSourceHeight * (G4UniformRand() - 0.5);  // -H/2 to +H/2
        
        G4double x = fSourcePosition.x() + r * std::cos(phi);
        G4double y = fSourcePosition.y() + r * std::sin(phi);
        G4double z = fSourcePosition.z() + h;
        
        return G4ThreeVector(x, y, z);
    }
    else if (fSourceType == "volume") {
        // Cylindrical volume source
        // Sample radius: r = R * sqrt(U) for uniform volume distribution in cylinder
        G4double r = fSourceRadius * std::sqrt(G4UniformRand());
        G4double phi = twopi * G4UniformRand();
        G4double h = fSourceHeight * (G4UniformRand() - 0.5);  // -H/2 to +H/2
        
        G4double x = fSourcePosition.x() + r * std::cos(phi);
        G4double y = fSourcePosition.y() + r * std::sin(phi);
        G4double z = fSourcePosition.z() + h;
        
        return G4ThreeVector(x, y, z);
    }
    else if (fSourceType == "marinelli") {
        // MARINELLI beaker - complex geometry (annular + bottom)
        return SampleMarinelliPosition();
    }
    else if (fSourceType == "cartridge") {
        // Cartridge source (Orano LEA Type D) - uniform in active cylindrical volume
        // Active matrix is a cylinder centered at fSourcePosition
        G4double r = fSourceRadius * std::sqrt(G4UniformRand());
        G4double phi = twopi * G4UniformRand();
        G4double h = fSourceHeight * (G4UniformRand() - 0.5);  // -H/2 to +H/2
        
        G4double x = fSourcePosition.x() + r * std::cos(phi);
        G4double y = fSourcePosition.y() + r * std::sin(phi);
        G4double z = fSourcePosition.z() + h;
        
        return G4ThreeVector(x, y, z);
    }
    else if (fSourceType == "filter") {
        // Filter source (Orano LEA Type M) - uniform in thin disk active surface
        // Active area is a thin disk centered at fSourcePosition
        G4double r = fSourceRadius * std::sqrt(G4UniformRand());
        G4double phi = twopi * G4UniformRand();
        G4double h = fSourceHeight * (G4UniformRand() - 0.5);  // thin disk
        
        G4double x = fSourcePosition.x() + r * std::cos(phi);
        G4double y = fSourcePosition.y() + r * std::sin(phi);
        G4double z = fSourcePosition.z() + h;
        
        return G4ThreeVector(x, y, z);
    }
    
    return fSourcePosition;
}

// ==================================================================
// Sample Marinelli position - UNIFORM in BOTH volumes
// ==================================================================
G4ThreeVector HPGePrimaryGeneratorAction::SampleMarinelliPosition() {
    
    if (!fMarinelliEnabled) {
        G4cout << "WARNING: Marinelli not enabled! Using point source." << G4endl;
        return fSourcePosition;
    }
    
    G4double wall = 2.0*mm;  // Wall thickness
    G4double airGap = 2.0*mm;
    
    // Calculate geometry
    G4double outerR = fMarinelliOuterD / 2.0;
    G4double innerR = fMarinelliInnerD / 2.0;
    
    
    // Annular region dimensions (around detector well, Z+)
    G4double annularRInner = innerR + wall;
    G4double annularROuter = outerR - wall;
    G4double annularHeight = fMarinelliWellH - wall;
    
    // Bottom region dimensions (solid cylinder, Z-)
    G4double bottomRadius = outerR - wall;
    G4double bottomHeight = fMarinelliTotalH - fMarinelliWellH - wall;
    
    // Calculate volumes
    G4double annularVolume = pi * (annularROuter*annularROuter - annularRInner*annularRInner) * annularHeight;
    G4double bottomVolume = pi * bottomRadius*bottomRadius * bottomHeight;
    G4double totalVolume = annularVolume + bottomVolume;
    
    // Probability of sampling from each volume (proportional to volume)
    G4double pAnnular = annularVolume / totalVolume;
    
    // Z position reference
    G4double beakerTop = airGap + fMarinelliWellH;
    
    // Randomly choose which volume to sample from
    G4double randomChoice = G4UniformRand();
    
    if (randomChoice < pAnnular) {
        // ===== ANNULAR REGION (puits) - Around detector =====
        // Sample radius in annular region (uniform in area)
        G4double u = G4UniformRand();
        G4double r = std::sqrt(u * (annularROuter*annularROuter - annularRInner*annularRInner) + annularRInner*annularRInner);
        G4double phi = twopi * G4UniformRand();
        
        // Sample height in annular region (uniform in height)
        G4double h = annularHeight * G4UniformRand();
        
        // Position in annular region
        // Starts at airGap and goes up to (airGap + annularHeight)
        G4double x = r * std::cos(phi);
        G4double y = r * std::sin(phi);
        G4double z = airGap + h;
        
        return G4ThreeVector(x, y, z);
    }
    else {
        // ===== BOTTOM REGION (cylindre plein) - Below detector =====
        // Sample radius in full cylinder (uniform in area)
        G4double r = bottomRadius * std::sqrt(G4UniformRand());
        G4double phi = twopi * G4UniformRand();
        
        // Sample height in bottom region (uniform in height)
        G4double h = bottomHeight * G4UniformRand();
        
        // Position in bottom region
        // Starts at (airGap - bottomHeight) and goes up to airGap
        G4double x = r * std::cos(phi);
        G4double y = r * std::sin(phi);
        G4double z = airGap - bottomHeight + h;
        
        return G4ThreeVector(x, y, z);
    }
}

// ==================================================================
// Set Marinelli geometry
// ==================================================================
void HPGePrimaryGeneratorAction::SetMarinelliGeometry(G4double outerD, G4double innerD, 
                                                      G4double totalH, G4double wellH) {
    fMarinelliOuterD = outerD;
    fMarinelliInnerD = innerD;
    fMarinelliTotalH = totalH;
    fMarinelliWellH = wellH;
    fMarinelliEnabled = true;
    
    G4cout << "\n=========================================" << G4endl;
    G4cout << "MARINELLI Beaker Geometry Configured:" << G4endl;
    G4cout << "  Outer diameter: " << outerD/cm << " cm" << G4endl;
    G4cout << "  Inner diameter: " << innerD/cm << " cm" << G4endl;
    G4cout << "  Total height: " << totalH/cm << " cm" << G4endl;
    G4cout << "  Well height: " << wellH/cm << " cm" << G4endl;
    
    // Calculate and display volumes
    G4double wall = 2.0*mm;
    G4double outerR = outerD / 2.0;
    G4double innerR = innerD / 2.0;
    G4double annularRInner = innerR + wall;
    G4double annularROuter = outerR - wall;
    G4double annularHeight = wellH - wall;
    G4double bottomRadius = outerR - wall;
    G4double bottomHeight = totalH - wellH - wall;
    
    G4double annularVol = pi * (annularROuter*annularROuter - annularRInner*annularRInner) * annularHeight;
    G4double bottomVol = pi * bottomRadius*bottomRadius * bottomHeight;
    
    G4cout << "  Annular volume: " << annularVol/cm3 << " cm³ (" 
           << 100.0*annularVol/(annularVol+bottomVol) << "%)" << G4endl;
    G4cout << "  Bottom volume: " << bottomVol/cm3 << " cm³ (" 
           << 100.0*bottomVol/(annularVol+bottomVol) << "%)" << G4endl;
    G4cout << "  Total volume: " << (annularVol+bottomVol)/cm3 << " cm³" << G4endl;
    G4cout << "=========================================" << G4endl;
}

// ==================================================================
// Sample isotropic direction (4π)
// ==================================================================
G4ThreeVector HPGePrimaryGeneratorAction::SampleIsotropicDirection() {
    // Uniform sampling on unit sphere
    // cos(theta) uniform in [-1, 1]
    // phi uniform in [0, 2π]
    
    G4double cosTheta = 2.0 * G4UniformRand() - 1.0;  // [-1, 1]
    G4double sinTheta = std::sqrt(1.0 - cosTheta * cosTheta);
    G4double phi = twopi * G4UniformRand();
    
    G4double x = sinTheta * std::cos(phi);
    G4double y = sinTheta * std::sin(phi);
    G4double z = cosTheta;
    
    return G4ThreeVector(x, y, z);
}

// ==================================================================
// Select gamma energy from spectrum
// ==================================================================
G4double HPGePrimaryGeneratorAction::SelectGammaEnergy() {
    if (fGammaLines.empty()) {
        G4cout << "WARNING: No gamma lines defined! Using 661.7 keV" << G4endl;
        return 661.7 * keV;
    }
    
    if (fGammaLines.size() == 1) {
        return fGammaLines[0].energy * keV;
    }
    
    // Sample from cumulative distribution
    G4double rand = G4UniformRand();
    
    for (size_t i = 0; i < fCumulativeIntensities.size(); i++) {
        if (rand < fCumulativeIntensities[i]) {
            return fGammaLines[i].energy * keV;
        }
    }
    
    // Fallback
    return fGammaLines.back().energy * keV;
}

// ==================================================================
// Add gamma line
// ==================================================================
void HPGePrimaryGeneratorAction::AddGammaLine(G4double energy, G4double intensity, G4String name) {
    GammaLine line(energy, intensity, name);
    fGammaLines.push_back(line);
    UpdateCumulativeDistribution();
}

// ==================================================================
// Update cumulative distribution for sampling
// ==================================================================
void HPGePrimaryGeneratorAction::UpdateCumulativeDistribution() {
    if (fGammaLines.empty()) return;
    
    // Calculate total intensity
    G4double totalIntensity = 0.0;
    for (const auto& line : fGammaLines) {
        totalIntensity += line.intensity;
    }
    
    // Build cumulative distribution
    fCumulativeIntensities.clear();
    G4double cumulative = 0.0;
    
    for (const auto& line : fGammaLines) {
        cumulative += line.intensity / totalIntensity;
        fCumulativeIntensities.push_back(cumulative);
    }
    
    // Ensure last value is exactly 1.0
    fCumulativeIntensities.back() = 1.0;
}

// ==================================================================
// Set source type
// ==================================================================
void HPGePrimaryGeneratorAction::SetSourceType(const std::string& type) {
    if (type == "point" || type == "disk" || type == "volume" || type == "marinelli" ||
        type == "cartridge" || type == "filter") {
        fSourceType = type;
    } else {
        G4cout << "WARNING: Unknown source type '" << type << "'. Using 'point'." << G4endl;
        fSourceType = "point";
    }
}

// ==================================================================
// Set isotope - Initialize standard isotopes
// ==================================================================
void HPGePrimaryGeneratorAction::SetIsotope(const std::string& isotope) {
    fCurrentIsotope = isotope;  // Store current isotope for cascade logic
    ClearGammaLines();
    
    if (isotope == "Custom") {
        // Custom energies will be added manually
        G4cout << "Custom energy mode - gamma lines must be added manually" << G4endl;
        G4cout << "NOTE: If custom energies match Co-60 or Y-88 patterns," << G4endl;
        G4cout << "      cascade emission will be automatically enabled." << G4endl;
        return;
    }
    
    if (isotope == "Cs137") {
        AddGammaLine(661.657, 85.1, "Cs-137 (Ba-137m)");
    }
    else if (isotope == "Co60") {
        AddGammaLine(1173.228, 99.85, "Co60 gamma-1");
        AddGammaLine(1332.492, 99.9826, "Co60 gamma-2");
    }
    else if (isotope == "Na22") {
       // AddGammaLine(511.0, 180.7, "e+ annihilation (2×)");
        AddGammaLine(1274.537, 99.94, "Na-22");
    }
    else if (isotope == "Am241") {
        AddGammaLine(59.5409, 35.9, "Am-241");
        //AddGammaLine(26.3446, 2.4, "Np L X-ray");
    }
    else if (isotope == "Ba133") {
        AddGammaLine(30.625, 33.1, "Ba-133 (Cs K-alpha1)");
        AddGammaLine(30.973, 17.8, "Ba-133 (Cs K-alpha2)");
        AddGammaLine(34.987, 5.7, "Ba-133 (Cs K-beta1)");
        AddGammaLine(53.1622, 2.14, "Ba-133 (1)");
        AddGammaLine(79.6142, 2.62, "Ba-133 (2)");
        AddGammaLine(80.9979, 32.9, "Ba-133 (3)");
        AddGammaLine(160.612, 0.645, "Ba-133 (4)");
        AddGammaLine(223.234, 0.45, "Ba-133 (5)");
        AddGammaLine(276.398, 7.16, "Ba-133 (6)");
        AddGammaLine(302.853, 18.33, "Ba-133 (7)");
        AddGammaLine(356.017, 62.05, "Ba-133 (8)");
        AddGammaLine(383.849, 8.94, "Ba-133 (9)");
    }
    else if (isotope == "Eu152") {
        AddGammaLine(121.782, 28.58, "Eu-152 (1)");
        AddGammaLine(244.697, 7.61, "Eu-152 (2)");
        AddGammaLine(344.278, 26.57, "Eu-152 (3)");
        AddGammaLine(411.116, 2.237, "Eu-152 (4)");
        AddGammaLine(443.965, 2.821, "Eu-152 (5)");
        AddGammaLine(778.904, 12.96, "Eu-152 (6)");
        AddGammaLine(867.380, 4.245, "Eu-152 (7)");
        AddGammaLine(964.079, 14.61, "Eu-152 (8)");
        AddGammaLine(1085.869, 10.21, "Eu-152 (9)");
        AddGammaLine(1112.074, 13.67, "Eu-152 (10)");
        AddGammaLine(1408.013, 21.01, "Eu-152 (11)");
    }
    else if (isotope == "Mn54") {
        AddGammaLine(834.848, 99.976, "Mn-54");
    }
    else if (isotope == "Cd109") {
        AddGammaLine(88.034, 3.644, "Cd-109");
    }
    else if (isotope == "Sn113") {
        AddGammaLine(391.698, 64.97, "Sn-113");
    }
    else if (isotope == "Y88") {
        AddGammaLine(898.042, 93.7, "Y-88 (1)");
        AddGammaLine(1836.063, 99.2, "Y-88 (2)");
    }
    else if (isotope == "Co57") {
        AddGammaLine(122.06065, 85.60, "Co-57 (1)");
        AddGammaLine(136.4736, 10.68, "Co-57 (2)");
    }
    else if (isotope == "Ce139") {
        AddGammaLine(165.8575, 79.90, "Ce-139");
    }        
    else if (isotope == "Sr85") {
        AddGammaLine(514.0048, 96, "Sr-85");
    }    
    else if (isotope == "Zn65") {
        AddGammaLine(1115.52, 50.04, "Zn-65");
    } 
            
    else {
        G4cout << "WARNING: Unknown isotope '" << isotope << "'. Using Cs-137." << G4endl;
        AddGammaLine(661.657, 85.1, "Cs-137 (default)");
    }
}

// ==================================================================
// Initialize isotope library
// ==================================================================
void HPGePrimaryGeneratorAction::InitializeIsotopeLibrary() {
    // This method can be extended with more isotopes
    // Currently, isotopes are loaded on-demand via SetIsotope()
}

// ==================================================================
// Print configuration
// ==================================================================
void HPGePrimaryGeneratorAction::PrintConfiguration() const {
    G4cout << "\n=== Primary Generator Configuration ===" << G4endl;
    G4cout << "Current isotope: " << fCurrentIsotope << G4endl;
    
    if (fSourceType != "marinelli") {
        G4cout << "Source position: (" << fSourcePosition.x()/cm << ", "
               << fSourcePosition.y()/cm << ", "
               << fSourcePosition.z()/cm << ") cm" << G4endl;
    }
    
    G4cout << "Source type: " << fSourceType << G4endl;
    
    if (fSourceType == "disk") {
        G4cout << "Source radius: " << fSourceRadius/cm << " cm" << G4endl;
    } else if (fSourceType == "volume") {
        G4cout << "Source radius: " << fSourceRadius/cm << " cm" << G4endl;
        G4cout << "Source height: " << fSourceHeight/cm << " cm" << G4endl;
    } else if (fSourceType == "marinelli" && fMarinelliEnabled) {
        G4cout << "Marinelli geometry configured (see above)" << G4endl;
    }
    
    G4cout << "\nGamma lines:" << G4endl;
    for (const auto& line : fGammaLines) {
        G4cout << "  " << line.energy << " keV (" << line.intensity << "%)";
        if (!line.name.empty()) {
            G4cout << " - " << line.name;
        }
        G4cout << G4endl;
    }
    
    // Note about cascade emission
    if (fCurrentIsotope == "Co60" || fCurrentIsotope == "Y88" || 
        fCurrentIsotope == "Na22" || fCurrentIsotope == "Eu152" || 
        fCurrentIsotope == "Ba133") {
        G4cout << "\n*** CASCADE EMISSION ENABLED for " << fCurrentIsotope << " ***" << G4endl;
        G4cout << "Multiple gammas will be emitted per decay event for TCS analysis" << G4endl;
    }
    else if (fCurrentIsotope == "Custom") {
        G4cout << "\n*** CUSTOM MODE with CASCADE AUTO-DETECTION ***" << G4endl;
        G4cout << "If energies match Co-60 or Y-88 patterns, cascades will be emitted" << G4endl;
    }
    
    G4cout << "=======================================" << G4endl;
}
