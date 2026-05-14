// ==================================================================
// HPGeEventAction.cc
// Professional Event Action with Accurate Efficiency Calculation
// ISOCS/LabSOCS Style Implementation
// FIXED: Removed extra closing brace causing compilation error
// ==================================================================

#include "HPGeEventAction.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"
#include <cmath>
#include <algorithm>
#include <iomanip>

// ==================================================================
// Constructor
// ==================================================================
HPGeEventAction::HPGeEventAction(const std::vector<HPGeSensitiveDetector*>& detectors)
: G4UserEventAction(),
  fDetectors(detectors),
  fCoincidenceWindow(100.*ns),
  fEnergyResolution(0.002),      // 0.2% FWHM at 1332 keV
  fEnergyThreshold(10.*keV),     // 10 keV detection threshold
  fFullEnergyWindow(0.03),       // ±3% window for full energy peak
  fTotalEvents(0),
  fDetectedEvents(0),
  fCoincidenceEvents(0),
  fFullEnergyEvents(0),
  fComptonEvents(0),
  fTCSEnabled(false),
  fTCSManager(nullptr) {
}

// ==================================================================
// Destructor
// ==================================================================
HPGeEventAction::~HPGeEventAction() {
    if (fOutputFile.is_open()) {
        PrintStatistics();
        fOutputFile.close();
    }
}

// ==================================================================
// Begin of event action
// ==================================================================
void HPGeEventAction::BeginOfEventAction(const G4Event*) {
    // Reset all detectors
    for (auto detector : fDetectors) {
        detector->Reset();
    }
}

// ==================================================================
// End of event action - CRITICAL FOR EFFICIENCY CALCULATION
// ==================================================================
void HPGeEventAction::EndOfEventAction(const G4Event* event) {
    
    fTotalEvents++;
    
    // Collect energy deposits from all detectors
    std::vector<G4double> measuredEnergies;
    std::vector<G4double> trueEnergies;
    std::vector<G4double> times;
    std::vector<G4int> detIDs;
    
    G4bool anyDetection = false;
    
    for (auto detector : fDetectors) {
        G4double totalEdep = detector->GetTotalEnergyDeposit();
        
        if (totalEdep > fEnergyThreshold) {
            anyDetection = true;
            
            // Check if this energy matches a known cascade sum
            bool isSumPeak = false;
            
            if (fTCSEnabled && fTCSManager && !fExpectedEnergies.empty()) {
                // Build list of possible cascade sums
                std::vector<std::pair<G4double, G4double>> cascadePairs;
                
                // For 2-energy cascades, build all possible pairs
                if (fExpectedEnergies.size() >= 2) {
                    for (size_t i = 0; i < fExpectedEnergies.size(); i++) {
                        for (size_t j = i+1; j < fExpectedEnergies.size(); j++) {
                            cascadePairs.push_back({fExpectedEnergies[i], fExpectedEnergies[j]});
                        }
                    }
                }
                
                // Check if total energy matches any cascade sum
                // Use LARGE tolerance to account for escape peaks and partial deposits
                G4double tolerance = 100.0*keV;  // ±400 keV - realistic for HPGe
                
                for (const auto& pair : cascadePairs) {
                    G4double sumEnergy = pair.first + pair.second;
                    
                    if (std::abs(totalEdep - sumEnergy) < tolerance) {
                        // This IS a sum peak! Split into components for TCS analysis
                        isSumPeak = true;
                        
                        measuredEnergies.push_back(ApplyEnergyResolution(pair.first));
                        measuredEnergies.push_back(ApplyEnergyResolution(pair.second));
                        trueEnergies.push_back(pair.first);
                        trueEnergies.push_back(pair.second);
                        times.push_back(detector->GetEventTime());
                        times.push_back(detector->GetEventTime());
                        detIDs.push_back(detector->GetDetectorID());
                        detIDs.push_back(detector->GetDetectorID());
                        break;
                    }
                }
            }
            
            // If NOT a sum peak, record the energy as-is
            if (!isSumPeak) {
                measuredEnergies.push_back(ApplyEnergyResolution(totalEdep));
                trueEnergies.push_back(totalEdep);
                times.push_back(detector->GetEventTime());
                detIDs.push_back(detector->GetDetectorID());
            }
        }
    }
    
    // Check if event was detected
    if (anyDetection) {
        fDetectedEvents++;
    }
    
    // ==================================================================
    // TCS PROCESSING - FIXED LOGIC
    // ==================================================================
    
    std::vector<G4double> finalEnergies = measuredEnergies;
    std::vector<G4double> finalTimes = times;
    TCSEvent tcsEvent;
    tcsEvent.isCoincidence = false; 
    
    if (fTCSEnabled && fTCSManager) {
        // IMPORTANT: We pass fExpectedEnergies as the "primary" (emitted) energies.
        // We pass 'trueEnergies' (detected) only if detection occurred.
        
        // Construct the event object explicitly for the manager
        tcsEvent.primaryEnergies = fExpectedEnergies; // The "Truth" of what was emitted
        tcsEvent.detectedEnergies = trueEnergies;     // What we actually saw
        tcsEvent.times = times;
        
        // Only perform coincidence logic if we actually detected something
        if (anyDetection && !trueEnergies.empty()) {
             // Recalculate summing based on detection
             TCSEvent processedEvent = fTCSManager->ProcessEvent(trueEnergies, times);
             tcsEvent.isCoincidence = processedEvent.isCoincidence;
             tcsEvent.sumEnergy = processedEvent.sumEnergy;
             tcsEvent.isSummed = processedEvent.isSummed;
             
             // Update our final display energies if a sum occurred
             if (tcsEvent.isCoincidence) {
                 finalEnergies.clear();
                 finalTimes.clear();
                 for (size_t i = 0; i < processedEvent.detectedEnergies.size(); ++i) {
                    if (processedEvent.detectedEnergies[i] > fEnergyThreshold) {
                        finalEnergies.push_back(ApplyEnergyResolution(processedEvent.detectedEnergies[i]));
                        finalTimes.push_back(times.empty() ? 0.0 : times[0]);
                    }
                }
             }
        }
        
        // UPDATE STATISTICS: This runs for EVERY event now.
        fTCSManager->UpdateStatistics(tcsEvent);
    }
    
    // Check for coincidences (both TCS and multi-detector)
    G4bool isCoincidence = false;
    
    // TCS coincidence: both gammas detected in same detector (sum peak)
    if (fTCSEnabled && fTCSManager && tcsEvent.isCoincidence) {
        isCoincidence = true;
    }
    // Multi-detector coincidence: gammas in different detectors
    else if (CheckCoincidence(finalTimes) && finalEnergies.size() >= 2) {
        isCoincidence = true;
    }
    
    // Increment counter once per event if coincidence detected
    if (isCoincidence && anyDetection) {
        fCoincidenceEvents++;
    }
    
    // Determine if this is a full energy event
    G4bool isFullEnergy = false;
    G4double closestExpectedEnergy = 0.;
    
    if (!finalEnergies.empty()) {
        // Sum all energies detected (for multi-detector)
        G4double totalDetectedEnergy = 0.;
        for (auto e : trueEnergies) {
            totalDetectedEnergy += e;
        }
        
        // Find closest expected energy
        closestExpectedEnergy = FindClosestExpectedEnergy(totalDetectedEnergy);
        
        // Check if within full energy window
        if (closestExpectedEnergy > 0.) {
            isFullEnergy = IsFullEnergyPeak(totalDetectedEnergy, closestExpectedEnergy);
            
            if (isFullEnergy) {
                fFullEnergyEvents++;
            } else {
                fComptonEvents++;
            }
            
            // Update efficiency statistics
            UpdateEfficiency(closestExpectedEnergy, true, isFullEnergy);
        }
    } else {
        // Event not detected - still need to update efficiency for all expected energies
        for (const auto& energy : fExpectedEnergies) {
            UpdateEfficiency(energy, false, false);
        }
    }
    
    // Write to output file
    if (fOutputFile.is_open() && anyDetection) {
        EventData eventData;
        eventData.eventID = event->GetEventID();
        eventData.isCoincidence = isCoincidence || (fTCSEnabled && tcsEvent.isCoincidence);
        eventData.nDetectors = finalEnergies.size();
        eventData.detectorIDs = detIDs;
        eventData.energies = finalEnergies;
        eventData.trueEnergies = trueEnergies;
        eventData.times = finalTimes;
        eventData.isFullEnergy = isFullEnergy;
        
        if (closestExpectedEnergy > 0.) {
            G4double totalEnergy = 0.;
            for (auto e : trueEnergies) totalEnergy += e;
            eventData.comptonFraction = 1.0 - (totalEnergy / closestExpectedEnergy);
        } else {
            eventData.comptonFraction = 0.;
        }
        
        WriteEvent(eventData);
    }
}

// ==================================================================
// Apply energy resolution (Gaussian broadening)
// ==================================================================
G4double HPGeEventAction::ApplyEnergyResolution(G4double energy) {
    const G4double referenceE = 1332.*keV;
    G4double fwhm = fEnergyResolution * energy * sqrt(referenceE / energy);
    G4double sigma = fwhm / 2.355;
    G4double measuredE = G4RandGauss::shoot(energy, sigma);
    return (measuredE > 0.) ? measuredE : 0.;
}

// ==================================================================
// Check if event is full energy peak
// ==================================================================
G4bool HPGeEventAction::IsFullEnergyPeak(G4double measuredEnergy, G4double trueEnergy) {
    G4double windowLow = trueEnergy * (1.0 - fFullEnergyWindow);
    G4double windowHigh = trueEnergy * (1.0 + fFullEnergyWindow);
    return (measuredEnergy >= windowLow && measuredEnergy <= windowHigh);
}

// ==================================================================
// Check for coincidence
// ==================================================================
G4bool HPGeEventAction::CheckCoincidence(const std::vector<G4double>& times) {
    if (times.size() < 2) return false;
    for (size_t i = 0; i < times.size() - 1; i++) {
        for (size_t j = i + 1; j < times.size(); j++) {
            if (std::abs(times[i] - times[j]) < fCoincidenceWindow) {
                return true;
            }
        }
    }
    return false;
}

// ==================================================================
// Find closest expected energy
// ==================================================================
G4double HPGeEventAction::FindClosestExpectedEnergy(G4double energy) {
    if (fExpectedEnergies.empty()) return 0.;
    
    G4double closest = fExpectedEnergies[0];
    G4double minDiff = std::abs(energy - closest);
    
    for (const auto& expectedE : fExpectedEnergies) {
        G4double diff = std::abs(energy - expectedE);
        if (diff < minDiff) {
            minDiff = diff;
            closest = expectedE;
        }
    }
    if (minDiff / closest < 0.20) {
        return closest;
    }
    return 0.;
}

// ==================================================================
// Update efficiency statistics
// ==================================================================
void HPGeEventAction::UpdateEfficiency(G4double energy, G4bool detected, G4bool fullEnergy) {
    if (fEfficiencyMap.find(energy) == fEfficiencyMap.end()) {
        EfficiencyData data;
        data.energy = energy;
        fEfficiencyMap[energy] = data;
    }
    
    EfficiencyData& data = fEfficiencyMap[energy];
    //data.totalEmitted++; // Always increment per event
    data.totalEmitted = fTotalEvents; //TotalEmmitted=beamOn
    
    if (detected) {
        data.totalDetected++;
        if (fullEnergy) {
            data.fullEnergyDetected++;
        } else {
            data.comptonDetected++;
        }
    }
    
    if (data.totalEmitted > 0) {
        data.totalEfficiency = (G4double)data.totalDetected / (G4double)data.totalEmitted;
        data.peakEfficiency = (G4double)data.fullEnergyDetected / (G4double)data.totalEmitted;
    }
    if (data.totalDetected > 0) {
        data.peakToTotal = (G4double)data.fullEnergyDetected / (G4double)data.totalDetected;
    }
}

// ==================================================================
// Write event to output file
// ==================================================================
void HPGeEventAction::WriteEvent(const EventData& eventData) {
    if (!fOutputFile.is_open()) return;
    
    fOutputFile << eventData.eventID << " "
                << (eventData.isCoincidence ? 1 : 0) << " "
                << eventData.nDetectors << " ";
    
    for (size_t i = 0; i < eventData.energies.size(); i++) {
        fOutputFile << eventData.detectorIDs[i] << " "
                    << eventData.energies[i]/keV << " "
                    << eventData.times[i]/ns << " ";
    }
    
    fOutputFile << (eventData.isFullEnergy ? 1 : 0) << " "
                << eventData.comptonFraction << std::endl;
}

// ==================================================================
// Set output file
// ==================================================================
void HPGeEventAction::SetOutputFile(const std::string& filename) {
    if (fOutputFile.is_open()) {
        fOutputFile.close();
    }
    fOutputFile.open(filename);
    if (fOutputFile.is_open()) {
        fOutputFile << "# HPGe Simulation Output - Professional Version" << std::endl;
        fOutputFile << "# Format: EventID Coincidence NDetectors [DetID Energy(keV) Time(ns)] ... FullEnergy ComptonFraction" << std::endl;
    }
}

// ==================================================================
// Set expected energies
// ==================================================================
void HPGeEventAction::SetExpectedEnergies(const std::vector<G4double>& energies,
                                          const std::vector<G4double>& intensities) {
    fExpectedEnergies = energies;
    fExpectedIntensities = intensities;
    
    for (const auto& energy : energies) {
        EfficiencyData data;
        data.energy = energy;
        fEfficiencyMap[energy] = data;
    }
    
    G4cout << "\n=== Expected Energies for Efficiency Calculation ===" << G4endl;
    for (size_t i = 0; i < energies.size(); i++) {
        G4cout << "  " << energies[i]/keV << " keV";
        if (i < intensities.size()) {
            G4cout << " (intensity: " << intensities[i] << "%)";
        }
        G4cout << G4endl;
    }
}

// ==================================================================
// Set actual emission counts from PrimaryGenerator - CRITICAL FIX!
// ==================================================================
void HPGeEventAction::SetEmissionCounts(const std::map<G4double, G4long>& counts) {
    G4cout << "\n=== Updating Emission Counts from PrimaryGenerator ===" << G4endl;
    
    // Update totalEmitted based on actual emissions from PrimaryGenerator
    for (const auto& pair : counts) {
        G4double energy = pair.first;
        G4long emittedCount = pair.second;
        
        // Find matching energy in efficiency map (with tolerance)
        G4bool found = false;
        for (auto& effPair : fEfficiencyMap) {
            G4double storedEnergy = effPair.first;
            
            // Match energies with tolerance (within 1 keV)
            if (std::abs(energy - storedEnergy) < 1.0*keV) {
                effPair.second.totalEmitted = emittedCount;
                
                // Recalculate efficiencies with correct emission count
                auto& data = effPair.second;
                if (data.totalEmitted > 0) {
                    data.totalEfficiency = (G4double)data.totalDetected / (G4double)data.totalEmitted;
                    data.peakEfficiency = (G4double)data.fullEnergyDetected / (G4double)data.totalEmitted;
                }
                
                G4cout << "  " << storedEnergy/keV << " keV: " 
                       << emittedCount << " emitted" << G4endl;
                found = true;
                break;
            }
        }
        
        if (!found) {
            G4cout << "  WARNING: No matching energy in efficiency map for " 
                   << energy/keV << " keV" << G4endl;
        }
    }
    
    G4cout << "=== Emission Counts Updated ===" << G4endl;
}

// ==================================================================
// Print statistics
// ==================================================================
void HPGeEventAction::PrintStatistics() const {
    G4cout << "\n" << std::string(70, '=') << G4endl;
    G4cout << "=== SIMULATION STATISTICS ===" << G4endl;
    G4cout << std::string(70, '=') << G4endl;
    
    G4cout << "\nEvent Statistics:" << G4endl;
    G4cout << "  Total events:           " << fTotalEvents << G4endl;
    G4cout << "  Detected events:        " << fDetectedEvents 
           << " (" << (fTotalEvents > 0 ? 100.0*fDetectedEvents/fTotalEvents : 0.0) << "%)" << G4endl;
    
    if (!fEfficiencyMap.empty()) {
        G4cout << "\n" << std::string(70, '-') << G4endl;
        G4cout << "EFFICIENCY DATA (ISOCS-Style Calculation)" << G4endl;
        G4cout << std::string(70, '-') << G4endl;
        G4cout << std::setw(10) << "Energy"
               << std::setw(12) << "Emitted"
               << std::setw(12) << "Detected"
               << std::setw(12) << "Full Peak"
               << std::setw(12) << "Total Eff"
               << std::setw(12) << "Peak Eff" << G4endl;
        G4cout << std::setw(10) << "(keV)"
               << std::setw(12) << "(count)"
               << std::setw(12) << "(count)"
               << std::setw(12) << "(count)"
               << std::setw(12) << "(%)"
               << std::setw(12) << "(%)" << G4endl;
        G4cout << std::string(70, '-') << G4endl;
        
        for (const auto& pair : fEfficiencyMap) {
            const EfficiencyData& data = pair.second;
            G4cout << std::setw(10) << std::fixed << std::setprecision(2) << data.energy/keV
                   << std::setw(12) << data.totalEmitted
                   << std::setw(12) << data.totalDetected
                   << std::setw(12) << data.fullEnergyDetected
                   << std::setw(12) << std::setprecision(4) << data.totalEfficiency*100
                   << std::setw(12) << std::setprecision(4) << data.peakEfficiency*100 << G4endl;
        }
        G4cout << std::string(70, '-') << G4endl;
    }
    
    G4cout << std::string(70, '=') << G4endl;
    
    // Print TCS statistics if enabled
    if (fTCSEnabled && fTCSManager) {
        fTCSManager->PrintStatistics();
    }
}

void HPGeEventAction::ExportEfficiencyData(const std::string& filename) const {
    std::ofstream file(filename);
    if (!file.is_open()) return;
    
    file << "# HPGe Efficiency Data" << std::endl;
    for (const auto& pair : fEfficiencyMap) {
        const EfficiencyData& data = pair.second;
        file << std::fixed << std::setprecision(3) << data.energy/keV << " "
             << data.totalEmitted << " "
             << data.totalDetected << " "
             << data.fullEnergyDetected << " "
             << std::setprecision(6) << data.totalEfficiency*100 << " "
             << data.peakEfficiency*100 << std::endl;
    }
    file.close();
}
