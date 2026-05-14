// ==================================================================
// HPGeTCSManager.cc
// True Coincidence Summing Manager Implementation
// ENHANCED: Dynamic Solid Angle + Efficiency-Based COI
// ==================================================================

#include "HPGeTCSManager.hh"
#include "G4SystemOfUnits.hh"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cmath>

// ==================================================================
// Constructor
// ==================================================================
HPGeTCSManager::HPGeTCSManager() 
    : fCoincidenceWindow(100.0*ns),
      fSolidAngle(0.0) {  // Initialize to 0 - will be calculated
    // Default geometry (Point source at 50mm)
    fGeometryConfig.sourceType = "Point";
    fGeometryConfig.detectorRadius = 30.0*mm;
    fGeometryConfig.detectorHeight = 60.0*mm;
    fGeometryConfig.sourceDistance = 50.0*mm;
    
    // Calculate initial solid angle
    fSolidAngle = CalculateSolidAngle();
}

// ==================================================================
// Destructor
// ==================================================================
HPGeTCSManager::~HPGeTCSManager() {
}

// ==================================================================
// Set geometry configuration and recalculate solid angle
// ==================================================================
void HPGeTCSManager::SetGeometryConfig(const GeometryConfig& config) {
    fGeometryConfig = config;
    
    // Recalculate solid angle based on new geometry
    fSolidAngle = CalculateSolidAngle();
    
    // Recalculate coincidence probabilities for all cascades
    fCoincidenceProbabilities.clear();
    std::vector<GammaCascade> cascades = fCascades;  // Copy
    fCascades.clear();
    for (const auto& cascade : cascades) {
        AddCascade(cascade);  // This will recalculate with new solid angle
    }
    
    G4cout << "\n=== Geometry Configuration Updated ===" << G4endl;
    G4cout << "Source Type: " << fGeometryConfig.sourceType << G4endl;
    G4cout << "Calculated Solid Angle: " << fSolidAngle*100 << "% (Ω/4π)" << G4endl;
    G4cout << "=====================================" << G4endl;
}

// ==================================================================
// Calculate solid angle based on current geometry
// ==================================================================
G4double HPGeTCSManager::CalculateSolidAngle() const {
    if (fGeometryConfig.sourceType == "Point") {
        return CalculateSolidAnglePoint();
    } else if (fGeometryConfig.sourceType == "Disk") {
        return CalculateSolidAngleDisk();
    } else if (fGeometryConfig.sourceType == "Cylinder") {
        return CalculateSolidAngleCylinder();
    } else if (fGeometryConfig.sourceType.contains("Marinelli")) {
        return CalculateSolidAngleMarinelli();
    } else {
        // Default fallback
        G4cout << "WARNING: Unknown source type '" << fGeometryConfig.sourceType 
               << "'. Using point source approximation." << G4endl;
        return CalculateSolidAnglePoint();
    }
}

// ==================================================================
// Calculate solid angle for point source
// Formula: Ω/4π = 0.5 * (1 - d/√(d² + R²))
// where d = distance from detector face, R = detector radius
// ==================================================================
G4double HPGeTCSManager::CalculateSolidAnglePoint() const {
    G4double d = fGeometryConfig.sourceDistance;
    G4double R = fGeometryConfig.detectorRadius;
    
    if (d <= 0) {
        G4cout << "WARNING: Invalid source distance. Using default solid angle 0.1" << G4endl;
        return 0.1;
    }
    
    // Solid angle for disk detector and point source on axis
    G4double solidAngle = 0.5 * (1.0 - d / std::sqrt(d*d + R*R));
    
    return solidAngle;
}

// ==================================================================
// Calculate solid angle for disk source (extended source)
// Uses numerical integration over the disk volume
// The disk has radius R_src, thickness h_src, centered at distance d
// ==================================================================
G4double HPGeTCSManager::CalculateSolidAngleDisk() const {
    G4double d = fGeometryConfig.sourceDistance;      // Distance to disk center
    G4double R_det = fGeometryConfig.detectorRadius;  // Detector radius
    G4double R_src = fGeometryConfig.sourceRadius;    // Disk radius
    G4double h_src = fGeometryConfig.sourceHeight;    // Disk thickness
    
    if (d <= 0 || R_src <= 0) {
        G4cout << "WARNING: Invalid disk source parameters. Using point approximation." << G4endl;
        return CalculateSolidAnglePoint();
    }
    
    // Numerical integration over the disk volume
    // Using cylindrical coordinates (r, phi, z)
    // Integrate solid angle weighted by volume element
    
    const int nR = 20;   // Radial divisions
    const int nZ = 10;   // Z divisions (thickness)
    
    G4double totalOmega = 0.0;
    G4double totalVolume = 0.0;
    
    G4double dr = R_src / nR;
    G4double dz = h_src / nZ;
    
    for (int ir = 0; ir < nR; ir++) {
        G4double r = (ir + 0.5) * dr;  // Radial position
        G4double dV_r = 2.0 * CLHEP::pi * r * dr;  // Ring volume element (without dz)
        
        for (int iz = 0; iz < nZ; iz++) {
            G4double z_offset = (iz + 0.5) * dz - h_src / 2.0;  // Z offset from center
            G4double z_point = d + z_offset;  // Actual Z distance to detector
            
            if (z_point <= 0) continue;  // Skip points behind detector
            
            // Volume element
            G4double dV = dV_r * dz;
            
            // For off-axis point at radius r and distance z_point,
            // use approximate solid angle (average over azimuthal angle)
            // This is an approximation valid for r << z_point
            G4double d_eff = std::sqrt(z_point * z_point + r * r);
            G4double omega_point = 0.5 * (1.0 - d_eff / std::sqrt(d_eff * d_eff + R_det * R_det));
            
            // Correction factor for off-axis position (cosine effect)
            G4double cos_factor = z_point / d_eff;
            omega_point *= cos_factor;
            
            totalOmega += omega_point * dV;
            totalVolume += dV;
        }
    }
    
    G4double avgOmega = (totalVolume > 0) ? totalOmega / totalVolume : 0.1;
    
    G4cout << "Disk source solid angle calculation:" << G4endl;
    G4cout << "  Source radius: " << R_src/mm << " mm" << G4endl;
    G4cout << "  Source thickness: " << h_src/mm << " mm" << G4endl;
    G4cout << "  Distance to center: " << d/mm << " mm" << G4endl;
    G4cout << "  Average solid angle: " << avgOmega * 100 << "% (Ω/4π)" << G4endl;
    
    return avgOmega;
}

// ==================================================================
// Calculate solid angle for cylindrical volume source
// Uses numerical integration over the cylinder volume
// ==================================================================
G4double HPGeTCSManager::CalculateSolidAngleCylinder() const {
    G4double d = fGeometryConfig.sourceDistance;      // Distance to cylinder center
    G4double R_det = fGeometryConfig.detectorRadius;  // Detector radius
    G4double R_src = fGeometryConfig.sourceRadius;    // Cylinder radius
    G4double h_src = fGeometryConfig.sourceHeight;    // Cylinder height
    
    if (d <= 0 || R_src <= 0 || h_src <= 0) {
        G4cout << "WARNING: Invalid cylinder source parameters. Using point approximation." << G4endl;
        return CalculateSolidAnglePoint();
    }
    
    // Numerical integration over the cylinder volume
    // Using cylindrical coordinates (r, phi, z)
    
    const int nR = 20;   // Radial divisions
    const int nZ = 30;   // Z divisions (more for cylinder height)
    
    G4double totalOmega = 0.0;
    G4double totalVolume = 0.0;
    
    G4double dr = R_src / nR;
    G4double dz = h_src / nZ;
    
    for (int ir = 0; ir < nR; ir++) {
        G4double r = (ir + 0.5) * dr;  // Radial position
        G4double dV_r = 2.0 * CLHEP::pi * r * dr;  // Ring volume element (without dz)
        
        for (int iz = 0; iz < nZ; iz++) {
            G4double z_offset = (iz + 0.5) * dz - h_src / 2.0;  // Z offset from center
            G4double z_point = d + z_offset;  // Actual Z distance to detector
            
            if (z_point <= 0) continue;  // Skip points behind detector
            
            // Volume element
            G4double dV = dV_r * dz;
            
            // For off-axis point at radius r and distance z_point,
            // calculate effective distance and solid angle
            G4double d_eff = std::sqrt(z_point * z_point + r * r);
            G4double omega_point = 0.5 * (1.0 - d_eff / std::sqrt(d_eff * d_eff + R_det * R_det));
            
            // Correction factor for off-axis position (cosine effect)
            G4double cos_factor = z_point / d_eff;
            omega_point *= cos_factor;
            
            totalOmega += omega_point * dV;
            totalVolume += dV;
        }
    }
    
    G4double avgOmega = (totalVolume > 0) ? totalOmega / totalVolume : 0.1;
    
    G4cout << "Cylinder source solid angle calculation:" << G4endl;
    G4cout << "  Source radius: " << R_src/mm << " mm" << G4endl;
    G4cout << "  Source height: " << h_src/mm << " mm" << G4endl;
    G4cout << "  Distance to center: " << d/mm << " mm" << G4endl;
    G4cout << "  Average solid angle: " << avgOmega * 100 << "% (Ω/4π)" << G4endl;
    
    return avgOmega;
}

// ==================================================================
// Calculate solid angle for Marinelli beaker
// This is a complex geometry - using numerical integration approach
// The Marinelli surrounds the detector, giving much higher solid angle
// ==================================================================
G4double HPGeTCSManager::CalculateSolidAngleMarinelli() const {
    G4double R_det = fGeometryConfig.detectorRadius;
    G4double H_det = fGeometryConfig.detectorHeight;
    G4double R_inner = fGeometryConfig.marinelliInnerRadius;
    G4double R_outer = fGeometryConfig.marinelliOuterRadius;
    G4double H_marinelli = fGeometryConfig.marinelliHeight;
    G4double H_inner = fGeometryConfig.marinelliInnerHeight;
    
    // Marinelli geometry creates high solid angle due to 4π surround geometry
    // Approximate calculation based on:
    // 1. Annular region around detector (high solid angle ~0.8-0.9)
    // 2. Top region above detector (moderate solid angle ~0.2-0.3)
    
    // Volume fractions
    G4double V_annular = CLHEP::pi * (R_outer*R_outer - R_inner*R_inner) * H_inner;
    G4double V_top = CLHEP::pi * R_outer*R_outer * (H_marinelli - H_inner);
    G4double V_total = V_annular + V_top;
    
    // Weighted solid angle contributions
    G4double omega_annular = 0.85;  // High solid angle for annular region
    G4double omega_top = 0.25;      // Lower solid angle for top region
    
    // Volume-weighted average
    G4double omega_marinelli = (omega_annular * V_annular + omega_top * V_top) / V_total;
    
    // Adjust for Marinelli type
    if (fGeometryConfig.sourceType == "Marinelli_1000mL") {
        // 1000mL has larger volume, slightly higher solid angle
        omega_marinelli *= 1.0;
    } else if (fGeometryConfig.sourceType == "Marinelli_450mL") {
        // 450mL is smaller, slightly lower solid angle
        omega_marinelli *= 0.95;
    }
    
    return omega_marinelli;
}

// ==================================================================
// Add cascade to the manager
// ==================================================================
void HPGeTCSManager::AddCascade(const GammaCascade& cascade) {
    fCascades.push_back(cascade);
    
    // Pre-calculate coincidence probabilities for this cascade
    for (size_t i = 0; i < cascade.energies.size(); ++i) {
        for (size_t j = i + 1; j < cascade.energies.size(); ++j) {
            G4double E1 = cascade.energies[i];
            G4double E2 = cascade.energies[j];
            
            // Probability of detecting both gammas
            // P = Prob(emission) * GeometricEff(1) * GeometricEff(2)
            G4double prob = cascade.probability * fSolidAngle * fSolidAngle;
            
            // Store both orderings
            fCoincidenceProbabilities[std::make_pair(E1, E2)] = prob;
            fCoincidenceProbabilities[std::make_pair(E2, E1)] = prob;
        }
    }
}

// ==================================================================
// Process event for TCS
// ==================================================================
TCSEvent HPGeTCSManager::ProcessEvent(const std::vector<G4double>& energies,
                                     const std::vector<G4double>& times) {
    TCSEvent event;
    event.primaryEnergies = energies;
    event.times = times;
    event.isCoincidence = false;
    event.sumEnergy = 0.0;
    
    // Initialize detected energies and summed flags
    event.detectedEnergies = energies;
    event.isSummed.resize(energies.size(), false);
    
    if (energies.empty()) return event;

    // Check for coincidences
    std::vector<bool> processed(energies.size(), false);
    
    for (size_t i = 0; i < energies.size(); ++i) {
        if (processed[i]) continue;
        
        for (size_t j = i + 1; j < energies.size(); ++j) {
            if (processed[j]) continue;
            
            // Check time coincidence
            if (AreCoincident(times[i], times[j])) {
                // Check if these energies are from a known cascade
                if (IsFromCascade(energies[i], energies[j])) {
                    event.isCoincidence = true;
                    
                    // Create sum peak
                    G4double sumE = energies[i] + energies[j];
                    event.sumEnergy = sumE;
                    
                    // Mark original peaks as summed
                    event.isSummed[i] = true;
                    event.isSummed[j] = true;
                    
                    // Replace with sum peak in detected energies
                    event.detectedEnergies[i] = sumE;
                    event.detectedEnergies[j] = 0.0;  // Zero out
                    
                    processed[i] = true;
                    processed[j] = true;
                }
            }
        }
    }
    
    // Remove zero energies
    event.detectedEnergies.erase(
        std::remove(event.detectedEnergies.begin(), 
                   event.detectedEnergies.end(), 0.0),
        event.detectedEnergies.end());
    
    return event;
}

// ==================================================================
// Calculate COI (Coincidence Correction Factor) - Theoretical
// ==================================================================
G4double HPGeTCSManager::CalculateCOI(G4double energy) {
    // COI = 1 / (1 - P_coincidence)
    // where P_coincidence is the probability of losing counts to summing
    
    G4double totalProbability = 0.0;
    
    // Find all cascades containing this energy
    std::vector<GammaCascade> cascades = FindCascadesWithEnergy(energy);
    
    for (const auto& cascade : cascades) {
        // For each other gamma in the cascade
        for (const auto& otherE : cascade.energies) {
            if (std::abs(otherE - energy) > 0.1*keV) {  // Not the same gamma
                totalProbability += GetCoincidenceProbability(energy, otherE);
            }
        }
    }
    
    // Ensure COI >= 1.0
    G4double COI = 1.0;
    if (totalProbability < 0.99) {  // Avoid division by very small numbers
        COI = 1.0 / (1.0 - totalProbability);
    }
    
    return COI;
}

// ==================================================================
// Get coincidence probability for two energies
// ==================================================================
G4double HPGeTCSManager::GetCoincidenceProbability(G4double energy1, 
                                                  G4double energy2) {
    auto key = std::make_pair(energy1, energy2);
    auto it = fCoincidenceProbabilities.find(key);
    
    if (it != fCoincidenceProbabilities.end()) {
        return it->second;
    }
    
    return 0.0;
}

// ==================================================================
// Update statistics with TCS event
// ==================================================================
void HPGeTCSManager::UpdateStatistics(const TCSEvent& event) {
    
    // Loop over PRIMARY (Emitted) energies to update "Total Emitted"
    for (size_t i = 0; i < event.primaryEnergies.size(); ++i) {
        G4double energy = event.primaryEnergies[i];
        
        // Identify the closest canonical cascade energy
        bool isKnown = false;
        for (const auto& cascade : fCascades) {
            for (const auto& cascadeE : cascade.energies) {
                if (std::abs(energy - cascadeE) < 10.0*keV) {
                    energy = cascadeE;
                    isKnown = true;
                    break;
                }
            }
            if (isKnown) break;
        }
        
        if (!isKnown) continue;

        if (fStatistics.find(energy) == fStatistics.end()) {
            fStatistics[energy] = TCSStatistics();
            fStatistics[energy].energy = energy;
        }
        
        // INCREMENT EMITTED
        fStatistics[energy].totalEmitted++;
    }
    
    // Now check DETECTED energies
    if (!event.detectedEnergies.empty()) {
         for (size_t i = 0; i < event.detectedEnergies.size(); ++i) {
             G4double detE = event.detectedEnergies[i];
             
             // Identify if this detected energy matches a known statistic
             bool found = false;
             G4double statsE = 0;
             
             // Check against active statistics
             for (auto& pair : fStatistics) {
                 if (std::abs(detE - pair.first) < 10.0*keV) {
                     statsE = pair.first;
                     found = true;
                     break;
                 }
             }
             
             if (found) {
                 fStatistics[statsE].singleDetected++;
             }
         }
    }
    
    // Handle Sum Peaks & Coincidences
    if (event.isCoincidence) {
        if (event.sumEnergy > 0) {
             bool isKnownSum = false;
             G4double knownSumE = 0;
             
             for (const auto& cascade : fCascades) {
                 if (cascade.energies.size() >= 2) {
                      G4double s = cascade.energies[0] + cascade.energies[1];
                      if (std::abs(event.sumEnergy - s) < 20.*keV) {
                          isKnownSum = true;
                          knownSumE = s;
                          break;
                      }
                 }
             }
             
             if (isKnownSum) {
                 if (fStatistics.find(knownSumE) == fStatistics.end()) {
                    fStatistics[knownSumE] = TCSStatistics();
                    fStatistics[knownSumE].energy = knownSumE;
                 }
                 fStatistics[knownSumE].sumPeakDetected++;
             }
        }
    }
    
    // Update COI factors
    for (auto& pair : fStatistics) {
        pair.second.COI = CalculateCOI(pair.first);
    }
}

// ==================================================================
// Update efficiency statistics (NEW - for Mode A/B comparison)
// ==================================================================
void HPGeTCSManager::UpdateEfficiencyStatistics(G4double energy, 
                                               G4long detected, 
                                               G4bool isTCSMode) {
    // Find or create statistics entry
    if (fStatistics.find(energy) == fStatistics.end()) {
        fStatistics[energy] = TCSStatistics();
        fStatistics[energy].energy = energy;
    }
    
    if (isTCSMode) {
        // Mode A: TCS enabled (reality)
        fStatistics[energy].detectedWithTCS = detected;
    } else {
        // Mode B: TCS disabled (fake singles)
        fStatistics[energy].detectedNoTCS = detected;
    }
    
    // Calculate efficiencies if we have emitted counts
    if (fStatistics[energy].totalEmitted > 0) {
        fStatistics[energy].efficiencyTCS = 
            static_cast<G4double>(fStatistics[energy].detectedWithTCS) / 
            static_cast<G4double>(fStatistics[energy].totalEmitted);
        
        fStatistics[energy].efficiencyNoTCS = 
            static_cast<G4double>(fStatistics[energy].detectedNoTCS) / 
            static_cast<G4double>(fStatistics[energy].totalEmitted);
        
        // Calculate efficiency-based COI
        if (fStatistics[energy].efficiencyTCS > 0) {
            fStatistics[energy].COI_efficiency = 
                fStatistics[energy].efficiencyNoTCS / fStatistics[energy].efficiencyTCS;
        } else {
            fStatistics[energy].COI_efficiency = 1.0;
        }
    }
}

// ==================================================================
// Check if energies are from cascade
// ==================================================================
G4bool HPGeTCSManager::IsFromCascade(G4double energy1, G4double energy2) const {
    for (const auto& cascade : fCascades) {
        bool hasE1 = false, hasE2 = false;
        
        for (const auto& E : cascade.energies) {
            if (std::abs(E - energy1) < 20.0*keV) hasE1 = true;
            if (std::abs(E - energy2) < 20.0*keV) hasE2 = true;
        }
        
        if (hasE1 && hasE2) return true;
    }
    
    return false;
}

// ==================================================================
// Get cascade energies for isotope
// ==================================================================
std::vector<G4double> HPGeTCSManager::GetCascadeEnergies(const G4String& isotope) const {
    std::vector<G4double> energies;
    
    for (const auto& cascade : fCascades) {
        if (cascade.parentIsotope == isotope) {
            for (const auto& E : cascade.energies) {
                if (std::find(energies.begin(), energies.end(), E) == energies.end()) {
                    energies.push_back(E);
                }
            }
        }
    }
    
    return energies;
}

// ==================================================================
// Find cascades containing specific energy
// ==================================================================
std::vector<GammaCascade> HPGeTCSManager::FindCascadesWithEnergy(G4double energy) const {
    std::vector<GammaCascade> result;
    
    for (const auto& cascade : fCascades) {
        for (const auto& E : cascade.energies) {
            if (std::abs(E - energy) < 1.0*keV) {
                result.push_back(cascade);
                break;
            }
        }
    }
    
    return result;
}

// ==================================================================
// Check if two events are coincident
// ==================================================================
G4bool HPGeTCSManager::AreCoincident(G4double time1, G4double time2) const {
    return std::abs(time1 - time2) < fCoincidenceWindow;
}

// ==================================================================
// Calculate geometric efficiency for cascade
// ==================================================================
G4double HPGeTCSManager::CalculateGeometricEfficiency(const GammaCascade& cascade) const {
    return std::pow(fSolidAngle, cascade.energies.size());
}

// ==================================================================
// Print statistics
// ==================================================================
void HPGeTCSManager::PrintStatistics() const {
    std::cout << "\n=============================================" << std::endl;
    std::cout << "True Coincidence Summing Statistics:" << std::endl;
    std::cout << "Solid Angle: " << fSolidAngle*100 << "% (Ω/4π)" << std::endl;
    std::cout << "=============================================" << std::endl;
    std::cout << std::setw(10) << "Energy" 
              << std::setw(12) << "Emitted"
              << std::setw(12) << "Single"
              << std::setw(12) << "Coincid."
              << std::setw(12) << "Sum Peak"
              << std::setw(10) << "COI" << std::endl;
    std::cout << "---------------------------------------------" << std::endl;
    
    for (const auto& pair : fStatistics) {
        const TCSStatistics& stats = pair.second;
        
        if (stats.totalEmitted == 0 && stats.singleDetected == 0 && 
            stats.coincidenceDetected == 0 && stats.sumPeakDetected == 0) continue;
        
        std::cout << std::setw(10) << stats.energy/keV
                  << std::setw(12) << stats.totalEmitted
                  << std::setw(12) << stats.singleDetected
                  << std::setw(12) << stats.coincidenceDetected
                  << std::setw(12) << stats.sumPeakDetected
                  << std::setw(10) << std::fixed << std::setprecision(3) 
                  << stats.COI << std::endl;
    }
    std::cout << "=============================================" << std::endl;
}

// ==================================================================
// Export TCS data to file
// ==================================================================
void HPGeTCSManager::ExportTCSData(const G4String& filename) const {
    std::ofstream file(filename);
    
    file << "# True Coincidence Summing Analysis Report" << std::endl;
    file << "# Coincidence Window: " << fCoincidenceWindow/ns << " ns" << std::endl;
    file << "# Solid Angle: " << fSolidAngle*100 << "% (Ω/4π)" << std::endl;
    file << "# Source Type: " << fGeometryConfig.sourceType << std::endl;
    file << "#" << std::endl;
    file << "# Energy(keV) Emitted Single Coincidence SumPeak COI" << std::endl;
    
    for (const auto& pair : fStatistics) {
        const TCSStatistics& stats = pair.second;
        file << stats.energy/keV << " "
             << stats.totalEmitted << " "
             << stats.singleDetected << " "
             << stats.coincidenceDetected << " "
             << stats.sumPeakDetected << " "
             << std::fixed << std::setprecision(4) << stats.COI << std::endl;
    }
    
    file.close();
}

// ==================================================================
// Export efficiency-based COI data (NEW)
// ==================================================================
void HPGeTCSManager::ExportEfficiencyData(const G4String& filename) const {
    std::ofstream file(filename);
    
    file << "# Efficiency-Based COI Analysis" << std::endl;
    file << "# Solid Angle: " << fSolidAngle*100 << "% (Ω/4π)" << std::endl;
    file << "# Source Type: " << fGeometryConfig.sourceType << std::endl;
    file << "#" << std::endl;
    file << "# Energy(keV) Emitted DetTCS DetNoTCS EffTCS EffNoTCS COI_theory COI_eff" << std::endl;
    
    for (const auto& pair : fStatistics) {
        const TCSStatistics& stats = pair.second;
        if (stats.totalEmitted > 0) {
            file << std::fixed << std::setprecision(2) << stats.energy/keV << " "
                 << stats.totalEmitted << " "
                 << stats.detectedWithTCS << " "
                 << stats.detectedNoTCS << " "
                 << std::setprecision(6) << stats.efficiencyTCS << " "
                 << stats.efficiencyNoTCS << " "
                 << std::setprecision(4) << stats.COI << " "
                 << stats.COI_efficiency << std::endl;
        }
    }
    
    file.close();
}

// ==================================================================
// Load Co60 decay scheme
// ==================================================================
void HPGeTCSManager::LoadCo60Scheme() {
    GammaCascade cascade;
    cascade.parentIsotope = "Co60";
    cascade.energies.push_back(1173.228*keV);
    cascade.energies.push_back(1332.492*keV);
    cascade.times.push_back(0.0);
    cascade.times.push_back(0.0);
    cascade.probability = 0.9985;   
    AddCascade(cascade);
}

// ==================================================================
// Load Y-88 decay scheme
// ==================================================================
void HPGeTCSManager::LoadY88Scheme() {
    GammaCascade cascade1;
    cascade1.parentIsotope = "Y88";
    cascade1.energies.push_back(898.042*keV);
    cascade1.energies.push_back(1836.063*keV);
    cascade1.times.push_back(0.0);
    cascade1.times.push_back(0.0);
    cascade1.probability = 0.937;
    AddCascade(cascade1);
}

// ==================================================================
// Load Eu-152 decay scheme
// ==================================================================
void HPGeTCSManager::LoadEu152Scheme() {
    GammaCascade cascade1;
    cascade1.parentIsotope = "Eu152";
    cascade1.energies.push_back(344.279*keV);
    cascade1.energies.push_back(411.117*keV);
    cascade1.times.push_back(0.0);
    cascade1.times.push_back(0.0);
    cascade1.probability = 0.0227;
    AddCascade(cascade1);
    
    GammaCascade cascade2;
    cascade2.parentIsotope = "Eu152";
    cascade2.energies.push_back(121.782*keV);
    cascade2.energies.push_back(344.279*keV);
    cascade2.times.push_back(0.0);
    cascade2.times.push_back(0.0);
    cascade2.probability = 0.2654;
    AddCascade(cascade2);
    
    GammaCascade cascade3;
    cascade3.parentIsotope = "Eu152";
    cascade3.energies.push_back(121.782*keV);
    cascade3.energies.push_back(244.697*keV);
    cascade3.times.push_back(0.0);
    cascade3.times.push_back(0.0);
    cascade3.probability = 0.0751;
    AddCascade(cascade3);
}

// ==================================================================
// Load Ba-133 decay scheme
// ==================================================================
void HPGeTCSManager::LoadBa133Scheme() {
    GammaCascade cascade1;
    cascade1.parentIsotope = "Ba133";
    cascade1.energies.push_back(80.998*keV);
    cascade1.energies.push_back(356.013*keV);
    cascade1.times.push_back(0.0);
    cascade1.times.push_back(0.0);
    cascade1.probability = 0.329;
    AddCascade(cascade1);
    
    GammaCascade cascade2;
    cascade2.parentIsotope = "Ba133";
    cascade2.energies.push_back(80.998*keV);
    cascade2.energies.push_back(302.851*keV);
    cascade2.times.push_back(0.0);
    cascade2.times.push_back(0.0);
    cascade2.probability = 0.183;
    AddCascade(cascade2);
}

// ==================================================================
// Load Na-22 decay scheme
// ==================================================================
void HPGeTCSManager::LoadNa22Scheme() {
    GammaCascade cascade;
    cascade.parentIsotope = "Na22";
    cascade.energies.push_back(511.0*keV);
    cascade.energies.push_back(1274.537*keV);
    cascade.times.push_back(0.0);
    cascade.times.push_back(0.0);
    cascade.probability = 0.9044;
    AddCascade(cascade);
}

// ==================================================================
// Load isotope scheme by name
// ==================================================================
void HPGeTCSManager::LoadIsotopeScheme(const G4String& isotope) {
    fCascades.erase(
        std::remove_if(fCascades.begin(), fCascades.end(),
            [&isotope](const GammaCascade& c) { 
                return c.parentIsotope == isotope; 
            }),
        fCascades.end());
    
    if (isotope == "Co60") {
        LoadCo60Scheme();
    } else if (isotope == "Y88") {
        LoadY88Scheme();
    } else if (isotope == "Eu152") {
        LoadEu152Scheme();
    } else if (isotope == "Ba133") {
        LoadBa133Scheme();
    } else if (isotope == "Na22") {
        LoadNa22Scheme();
    }
}
