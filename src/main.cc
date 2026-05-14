// ==================================================================
// main.cc - HPGe Monte Carlo with Enhanced TCS Support
// VERSION 2.0 - Dynamic Solid Angle & Efficiency-Based COI
// ==================================================================

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4SDManager.hh"

#include "HPGeDetectorConstruction.hh"
#include "HPGePhysicsList.hh"
#include "HPGePrimaryGeneratorAction.hh"
#include "HPGeEventAction.hh"
#include "HPGeSensitiveDetector.hh"
#include "HPGeTCSManager.hh"

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

// ==================================================================
// Configuration loader with Marinelli support and TCS enhancement
// ==================================================================
void LoadConfiguration(const G4String& configFile,
                      HPGeDetectorConstruction* detector,
                      HPGePrimaryGeneratorAction* primary,
                      HPGeEventAction* eventAction) {
    
    std::ifstream file(configFile);
    if (!file.is_open()) {
        G4cout << "⚠ Config file not found: " << configFile << G4endl;
        G4cout << "Using default configuration..." << G4endl;
        return;
    }
    
    G4cout << "\n=== Loading Configuration ===" << G4endl;
    G4cout << "File: " << configFile << G4endl;
    
    G4double srcX = 0.0, srcY = 0.0, srcZ = 5.0;
    G4bool isCustomEnergy = false;
    std::vector<G4double> customEnergies;
    std::vector<G4double> customIntensities;
    std::string isotopeName;
    std::string sourceType = "point";
    
    // Marinelli configuration
    G4bool useMarinelli = false;
    std::string marinelliType = "1000ml";
    std::string marinelliFillMaterial = "water";
    
    // Disk source configuration
    G4bool useDiskSource = false;
    G4double diskRadius = 25.0;     // mm
    G4double diskThickness = 5.0;   // mm
    G4double diskPositionZ = -50.0; // mm (negative = below detector)
    std::string diskMaterial = "water";
    
    // Cylinder source configuration
    G4bool useCylinderSource = false;
    G4double cylinderInnerRadius = 30.0;   // mm
    G4double cylinderOuterRadius = 32.0;   // mm
    G4double cylinderHeight = 50.0;        // mm
    G4double cylinderWallThickness = 2.0;  // mm
    G4double cylinderBottomThickness = 2.0;// mm
    G4double cylinderPositionZ = -60.0;    // mm
    std::string cylinderWallMaterial = "polypropylene";
    std::string cylinderFillMaterial = "water";
    
    // Cartridge source configuration (Orano LEA Type D)
    G4bool useCartridgeSource = false;
    G4double cartridgeHousingDiameter = 56.0;   // mm
    G4double cartridgeHousingHeight = 26.5;     // mm
    G4double cartridgeActiveDiameter = 52.0;    // mm
    G4double cartridgeActiveThickness = 20.5;   // mm
    G4double cartridgeWallThickness = 2.0;      // mm
    G4double cartridgePositionZ = -40.0;        // mm
    std::string cartridgeHousingMaterial = "polycarbonate";
    std::string cartridgeActiveMaterial = "activated_carbon";
    
    // Filter source configuration (Orano LEA Type M)
    G4bool useFilterSource = false;
    G4double filterOuterDiameter = 53.0;    // mm (M53 default)
    G4double filterActiveDiameter = 47.0;   // mm
    G4double filterThickness = 0.5;         // mm
    G4double filterSealThickness = 0.05;    // mm
    G4double filterPositionZ = -30.0;       // mm
    std::string filterMaterial = "water";
    std::string filterSealMaterial = "polyester";
    std::string filterType = "M53";
    
    // TCS configuration
    G4bool enableTCS = false;
    G4double tcsCoincidenceWindow = 100.0; // ns
    
    // Fast mode configuration
    G4bool fastMode = false;
    
    // Point source configuration
    G4double pointSourceDistance = 50.0; // mm
    
    // Detector geometry parameters (defaults = Canberra Mirion S/N:21206)
    G4bool geometryFromConfig = false;
    G4bool extendedGeometryFromConfig = false;
    G4double cfg_crystalDiameter = 60.5;     // mm (Canberra default)
    G4double cfg_crystalLength = 45.4;       // mm
    G4double cfg_holeDiameter = 9.5;         // mm
    G4double cfg_holeDepth = 19.0;           // mm
    G4double cfg_geDeadLayer = 0.5;          // mm (outer electrode)
    G4double cfg_geDeadLayerFront = 0.0004;  // mm (0.4 um front dead layer)
    G4double cfg_liDeadLayer = 0.0003;       // mm (0.3 um inner electrode)
    G4double cfg_liDeadLayerFront = 0.0003;  // mm
    G4double cfg_alWindowThickness = 0.6;    // mm (carbon epoxy window)
    G4double cfg_alCupThickness = 1.5;       // mm (endcap side wall)
    G4double cfg_vacuumGap = 6.0;            // mm (Ge to endcap)
    
    // Extended Canberra-specific parameters
    G4double cfg_endcapDiameter = 76.2;       // mm
    G4double cfg_endcapLength = 8.0;          // mm
    std::string cfg_endcapMaterial = "Al";
    std::string cfg_windowMaterial = "carbon_epoxy";
    G4double cfg_geToEndcapDistance = 6.0;     // mm
    G4double cfg_holderOuterDiameter = 69.0;   // mm
    G4double cfg_holderTotalLength = 95.2;     // mm
    std::string cfg_holderMaterial = "Al";
    G4double cfg_hdpeThickness = 6.5;          // mm
    std::string cfg_geometryMode = "canberra"; // canberra or fluka
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string key, equals, value;
        if (!(iss >> key >> equals >> value)) continue;
        
        // Source position
        if (key == "source_x") {
            srcX = std::stod(value);
        }
        else if (key == "source_y") {
            srcY = std::stod(value);
        }
        else if (key == "source_z") {
            srcZ = std::stod(value);
        }
        else if (key == "point_source_distance") {
            pointSourceDistance = std::stod(value);
        }
        
        // Isotope configuration
        else if (key == "isotope") {
            if (value == "Custom") {
                isCustomEnergy = true;
                G4cout << "  Custom energy mode activated" << G4endl;
            } else if (value == "Multi") {
                G4cout << "  Multi-isotope mode detected" << G4endl;
            } else {
                isotopeName = value;
                primary->SetIsotope(value);
                G4cout << "  Isotope: " << value << G4endl;
            }
        }
        
        // Custom energies
        else if (key == "custom_energies") {
            std::stringstream ss(value);
            std::string token;
            while (std::getline(ss, token, ',')) {
                try {
                    G4double energy = std::stod(token);
                    customEnergies.push_back(energy);
                } catch (...) {
                    G4cout << "  ⚠ Invalid energy value: " << token << G4endl;
                }
            }
        }
        
        // Custom intensities
        else if (key == "custom_intensities") {
            std::stringstream ss(value);
            std::string token;
            while (std::getline(ss, token, ',')) {
                try {
                    G4double intensity = std::stod(token);
                    customIntensities.push_back(intensity);
                } catch (...) {
                    G4cout << "  ⚠ Invalid intensity value: " << token << G4endl;
                }
            }
        }
        
        // Source geometry
        else if (key == "source_type") {
            sourceType = value;
            primary->SetSourceType(value);
            G4cout << "  Source type: " << value << G4endl;
            if (value == "marinelli") {
                useMarinelli = true;
            }
            else if (value == "disk") {
                useDiskSource = true;
            }
            else if (value == "volume") {
                useCylinderSource = true;
            }
            else if (value == "cartridge") {
                useCartridgeSource = true;
            }
            else if (value == "filter") {
                useFilterSource = true;
            }
        }
        
        // Disk source configuration
        else if (key == "disk_radius") {
            diskRadius = std::stod(value);
            G4cout << "  Disk radius: " << diskRadius << " mm" << G4endl;
        }
        else if (key == "disk_thickness") {
            diskThickness = std::stod(value);
            G4cout << "  Disk thickness: " << diskThickness << " mm" << G4endl;
        }
        else if (key == "disk_position_z") {
            diskPositionZ = std::stod(value);
            G4cout << "  Disk position Z: " << diskPositionZ << " mm" << G4endl;
        }
        else if (key == "disk_material") {
            diskMaterial = value;
            G4cout << "  Disk material: " << diskMaterial << G4endl;
        }
        
        // Cylinder source configuration
        else if (key == "cylinder_inner_radius") {
            cylinderInnerRadius = std::stod(value);
            G4cout << "  Cylinder inner radius: " << cylinderInnerRadius << " mm" << G4endl;
        }
        else if (key == "cylinder_outer_radius") {
            cylinderOuterRadius = std::stod(value);
            G4cout << "  Cylinder outer radius: " << cylinderOuterRadius << " mm" << G4endl;
        }
        else if (key == "cylinder_height") {
            cylinderHeight = std::stod(value);
            G4cout << "  Cylinder height: " << cylinderHeight << " mm" << G4endl;
        }
        else if (key == "cylinder_wall_thickness") {
            cylinderWallThickness = std::stod(value);
            G4cout << "  Cylinder wall thickness: " << cylinderWallThickness << " mm" << G4endl;
        }
        else if (key == "cylinder_bottom_thickness") {
            cylinderBottomThickness = std::stod(value);
            G4cout << "  Cylinder bottom thickness: " << cylinderBottomThickness << " mm" << G4endl;
        }
        else if (key == "cylinder_position_z") {
            cylinderPositionZ = std::stod(value);
            G4cout << "  Cylinder position Z: " << cylinderPositionZ << " mm" << G4endl;
        }
        else if (key == "cylinder_wall_material") {
            cylinderWallMaterial = value;
            G4cout << "  Cylinder wall material: " << cylinderWallMaterial << G4endl;
        }
        else if (key == "cylinder_fill_material") {
            cylinderFillMaterial = value;
            G4cout << "  Cylinder fill material: " << cylinderFillMaterial << G4endl;
        }
        
        // Cartridge source configuration (Orano LEA Type D)
        else if (key == "cartridge_housing_diameter") {
            cartridgeHousingDiameter = std::stod(value);
            G4cout << "  Cartridge housing diameter: " << cartridgeHousingDiameter << " mm" << G4endl;
        }
        else if (key == "cartridge_housing_height") {
            cartridgeHousingHeight = std::stod(value);
            G4cout << "  Cartridge housing height: " << cartridgeHousingHeight << " mm" << G4endl;
        }
        else if (key == "cartridge_active_diameter") {
            cartridgeActiveDiameter = std::stod(value);
            G4cout << "  Cartridge active diameter: " << cartridgeActiveDiameter << " mm" << G4endl;
        }
        else if (key == "cartridge_active_thickness") {
            cartridgeActiveThickness = std::stod(value);
            G4cout << "  Cartridge active thickness: " << cartridgeActiveThickness << " mm" << G4endl;
        }
        else if (key == "cartridge_wall_thickness") {
            cartridgeWallThickness = std::stod(value);
            G4cout << "  Cartridge wall thickness: " << cartridgeWallThickness << " mm" << G4endl;
        }
        else if (key == "cartridge_position_z") {
            cartridgePositionZ = std::stod(value);
            G4cout << "  Cartridge position Z: " << cartridgePositionZ << " mm" << G4endl;
        }
        else if (key == "cartridge_housing_material") {
            cartridgeHousingMaterial = value;
            G4cout << "  Cartridge housing material: " << value << G4endl;
        }
        else if (key == "cartridge_active_material") {
            cartridgeActiveMaterial = value;
            G4cout << "  Cartridge active material: " << value << G4endl;
        }
        
        // Filter source configuration (Orano LEA Type M)
        else if (key == "filter_outer_diameter") {
            filterOuterDiameter = std::stod(value);
            G4cout << "  Filter outer diameter: " << filterOuterDiameter << " mm" << G4endl;
        }
        else if (key == "filter_active_diameter") {
            filterActiveDiameter = std::stod(value);
            G4cout << "  Filter active diameter: " << filterActiveDiameter << " mm" << G4endl;
        }
        else if (key == "filter_thickness") {
            filterThickness = std::stod(value);
            G4cout << "  Filter thickness: " << filterThickness << " mm" << G4endl;
        }
        else if (key == "filter_seal_thickness") {
            filterSealThickness = std::stod(value);
            G4cout << "  Filter seal thickness: " << filterSealThickness << " mm" << G4endl;
        }
        else if (key == "filter_position_z") {
            filterPositionZ = std::stod(value);
            G4cout << "  Filter position Z: " << filterPositionZ << " mm" << G4endl;
        }
        else if (key == "filter_material") {
            filterMaterial = value;
            G4cout << "  Filter material: " << value << G4endl;
        }
        else if (key == "filter_seal_material") {
            filterSealMaterial = value;
            G4cout << "  Filter seal material: " << value << G4endl;
        }
        else if (key == "filter_type") {
            filterType = value;
            G4cout << "  Filter type: " << value << G4endl;
        }
        
        // Marinelli configuration
        else if (key == "marinelli_type") {
            marinelliType = value;
            G4cout << "  Marinelli type: " << value << G4endl;
        }
        else if (key == "marinelli_fill_material") {
            marinelliFillMaterial = value;
            G4cout << "  Marinelli fill material: " << value << G4endl;
        }
        
        // Event action parameters (only if eventAction is valid)
        else if (eventAction && key == "coincidence_window") {
            G4double window = std::stod(value);
            eventAction->SetCoincidenceWindow(window*ns);
            G4cout << "  Coincidence window: " << window << " ns" << G4endl;
        }
        else if (eventAction && key == "energy_resolution") {
            G4double res = std::stod(value);
            eventAction->SetEnergyResolution(res/100.);
            G4cout << "  Energy resolution: " << res << " % FWHM" << G4endl;
        }
        else if (eventAction && key == "energy_threshold") {
            G4double threshold = std::stod(value);
            eventAction->SetEnergyThreshold(threshold*keV);
            G4cout << "  Energy threshold: " << threshold << " keV" << G4endl;
        }
        
        // TCS parameters (kept for backward compatibility, but cascade mode is now automatic)
        else if (key == "enable_tcs") {
            enableTCS = (value == "true" || value == "1");
            G4cout << "  TCS flag in config: " << (enableTCS ? "true" : "false") << G4endl;
            // NOTE: Cascade mode is now set automatically based on isotope type
            // enable_tcs only affects the TCSManager, not the emission mode
        }
        else if (key == "tcs_coincidence_window") {
            tcsCoincidenceWindow = std::stod(value);
            G4cout << "  TCS coincidence window: " << tcsCoincidenceWindow << " ns" << G4endl;
        }
        else if (key == "fast_mode") {
            fastMode = (value == "true" || value == "1");
            G4cout << "  Fast mode: " << (fastMode ? "ENABLED" : "disabled") << G4endl;
        }
        
        // ===== NEW: Detector geometry parameters =====
        else if (key == "crystal_diameter") {
            cfg_crystalDiameter = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Crystal diameter: " << cfg_crystalDiameter << " mm" << G4endl;
        }
        else if (key == "crystal_length") {
            cfg_crystalLength = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Crystal length: " << cfg_crystalLength << " mm" << G4endl;
        }
        else if (key == "hole_diameter") {
            cfg_holeDiameter = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Hole diameter: " << cfg_holeDiameter << " mm" << G4endl;
        }
        else if (key == "hole_depth") {
            cfg_holeDepth = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Hole depth: " << cfg_holeDepth << " mm" << G4endl;
        }
        else if (key == "ge_dead_layer") {
            cfg_geDeadLayer = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Ge dead layer (lateral): " << cfg_geDeadLayer << " mm" << G4endl;
        }
        else if (key == "ge_dead_layer_front") {
            cfg_geDeadLayerFront = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Ge dead layer (front): " << cfg_geDeadLayerFront << " mm" << G4endl;
        }
        else if (key == "li_dead_layer") {
            cfg_liDeadLayer = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Li dead layer (lateral): " << cfg_liDeadLayer << " mm" << G4endl;
        }
        else if (key == "li_dead_layer_front") {
            cfg_liDeadLayerFront = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Li dead layer (front): " << cfg_liDeadLayerFront << " mm" << G4endl;
        }
        else if (key == "al_window_thickness") {
            cfg_alWindowThickness = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Al window thickness: " << cfg_alWindowThickness << " mm" << G4endl;
        }
        else if (key == "al_cup_thickness") {
            cfg_alCupThickness = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Al cup thickness: " << cfg_alCupThickness << " mm" << G4endl;
        }
        else if (key == "vacuum_gap") {
            cfg_vacuumGap = std::stod(value);
            geometryFromConfig = true;
            G4cout << "  Vacuum gap: " << cfg_vacuumGap << " mm" << G4endl;
        }
        
        // Extended Canberra geometry parameters
        else if (key == "endcap_diameter") {
            cfg_endcapDiameter = std::stod(value);
            extendedGeometryFromConfig = true;
            G4cout << "  Endcap diameter: " << cfg_endcapDiameter << " mm" << G4endl;
        }
        else if (key == "endcap_length") {
            cfg_endcapLength = std::stod(value);
            extendedGeometryFromConfig = true;
            G4cout << "  Endcap length: " << cfg_endcapLength << " mm" << G4endl;
        }
        else if (key == "endcap_material") {
            cfg_endcapMaterial = value;
            extendedGeometryFromConfig = true;
            G4cout << "  Endcap material: " << cfg_endcapMaterial << G4endl;
        }
        else if (key == "window_material") {
            cfg_windowMaterial = value;
            extendedGeometryFromConfig = true;
            G4cout << "  Window material: " << cfg_windowMaterial << G4endl;
        }
        else if (key == "ge_to_endcap_distance") {
            cfg_geToEndcapDistance = std::stod(value);
            extendedGeometryFromConfig = true;
            G4cout << "  Ge to endcap distance: " << cfg_geToEndcapDistance << " mm" << G4endl;
        }
        else if (key == "holder_outer_diameter") {
            cfg_holderOuterDiameter = std::stod(value);
            extendedGeometryFromConfig = true;
            G4cout << "  Holder outer diameter: " << cfg_holderOuterDiameter << " mm" << G4endl;
        }
        else if (key == "holder_total_length") {
            cfg_holderTotalLength = std::stod(value);
            extendedGeometryFromConfig = true;
            G4cout << "  Holder total length: " << cfg_holderTotalLength << " mm" << G4endl;
        }
        else if (key == "holder_material") {
            cfg_holderMaterial = value;
            extendedGeometryFromConfig = true;
            G4cout << "  Holder material: " << cfg_holderMaterial << G4endl;
        }
        else if (key == "hdpe_thickness") {
            cfg_hdpeThickness = std::stod(value);
            extendedGeometryFromConfig = true;
            G4cout << "  HDPE thickness: " << cfg_hdpeThickness << " mm" << G4endl;
        }
        else if (key == "geometry_mode") {
            cfg_geometryMode = value;
            G4cout << "  Geometry mode: " << cfg_geometryMode << G4endl;
        }
    }
    
    // Apply geometry mode (canberra or fluka)
    if (cfg_geometryMode == "fluka" || cfg_geometryMode == "FLUKA") {
        detector->SetUseCanberraGeometry(false);
        G4cout << "  Geometry mode: FLUKA (simplified)" << G4endl;
    } else {
        detector->SetUseCanberraGeometry(true);
        G4cout << "  Geometry mode: Canberra (real geometry)" << G4endl;
    }
    
    // Apply detector geometry
    if (extendedGeometryFromConfig || geometryFromConfig) {
        detector->SetDetectorGeometryExtended(
            cfg_crystalDiameter * mm,
            cfg_crystalLength * mm,
            cfg_holeDiameter * mm,
            cfg_holeDepth * mm,
            cfg_geDeadLayer * mm,
            cfg_geDeadLayerFront * mm,
            cfg_liDeadLayer * mm,
            cfg_liDeadLayerFront * mm,
            cfg_alWindowThickness * mm,
            cfg_alCupThickness * mm,
            cfg_vacuumGap * mm,
            cfg_endcapDiameter * mm,
            cfg_endcapLength * mm,
            cfg_endcapMaterial,
            cfg_windowMaterial,
            cfg_geToEndcapDistance * mm,
            cfg_holderOuterDiameter * mm,
            cfg_holderTotalLength * mm,
            cfg_holderMaterial,
            cfg_hdpeThickness * mm
        );
    }
    
    // Configure Marinelli if enabled
    if (useMarinelli) {
        MarinelliConfig marinelliConfig;
        marinelliConfig.SetType(marinelliType);
        marinelliConfig.fillMaterial = marinelliFillMaterial;
        marinelliConfig.isActive = true;
        
        detector->SetMarinelliConfig(marinelliConfig);
        detector->EnableMarinelli(true);
        
        primary->SetMarinelliGeometry(
            marinelliConfig.outerDiameter,
            marinelliConfig.innerDiameter,
            marinelliConfig.totalHeight,
            marinelliConfig.innerHeight
        );
        primary->EnableMarinelli(true);
        
        G4cout << "\n  ✓ Marinelli beaker enabled:" << G4endl;
        G4cout << "    Type: " << marinelliType << G4endl;
        G4cout << "    Fill material: " << marinelliFillMaterial << G4endl;
    }
    
    // Configure Disk source if enabled
    if (useDiskSource) {
        DiskSourceConfig diskConfig;
        diskConfig.isActive = true;
        diskConfig.radius = diskRadius * mm;
        diskConfig.thickness = diskThickness * mm;
        diskConfig.positionZ = diskPositionZ * mm;
        diskConfig.material = diskMaterial;
        
        detector->SetDiskSourceConfig(diskConfig);
        detector->EnableDiskSource(true);
        
        // Configure primary generator for disk source
        primary->SetSourceRadius(diskRadius * mm);
        primary->SetSourceHeight(diskThickness * mm);  // CRITICAL: Set thickness for Z sampling
        primary->SetSourcePosition(G4ThreeVector(0, 0, diskPositionZ * mm));
        
        G4cout << "\n  ✓ Disk source enabled:" << G4endl;
        G4cout << "    Radius: " << diskRadius << " mm" << G4endl;
        G4cout << "    Thickness: " << diskThickness << " mm" << G4endl;
        G4cout << "    Position Z: " << diskPositionZ << " mm" << G4endl;
        G4cout << "    Material: " << diskMaterial << G4endl;
        G4cout << "    === PRIMARY GENERATOR CONFIG ===" << G4endl;
        G4cout << "    Sampling radius: " << diskRadius << " mm" << G4endl;
        G4cout << "    Sampling height: " << diskThickness << " mm" << G4endl;
        G4cout << "    Sampling center Z: " << diskPositionZ << " mm" << G4endl;
        G4cout << "    Z range: [" << (diskPositionZ - diskThickness/2.0) << ", " 
               << (diskPositionZ + diskThickness/2.0) << "] mm" << G4endl;
    }
    
    // Configure Cylinder source if enabled
    if (useCylinderSource) {
        CylinderSourceConfig cylConfig;
        cylConfig.isActive = true;
        cylConfig.innerRadius = cylinderInnerRadius * mm;
        cylConfig.outerRadius = cylinderOuterRadius * mm;
        cylConfig.height = cylinderHeight * mm;
        cylConfig.wallThickness = cylinderWallThickness * mm;
        cylConfig.bottomThickness = cylinderBottomThickness * mm;
        cylConfig.positionZ = cylinderPositionZ * mm;
        cylConfig.wallMaterial = cylinderWallMaterial;
        cylConfig.fillMaterial = cylinderFillMaterial;
        
        detector->SetCylinderSourceConfig(cylConfig);
        detector->EnableCylinderSource(true);
        
        // Configure primary generator for cylinder source
        // IMPORTANT: Sample radius = outer radius - wall thickness
        // This matches the geometry in HPGeDetectorConstruction where innerR = outerR - wall
        G4double sampleRadius = cylinderOuterRadius - cylinderWallThickness;  // Correct calculation
        G4double sampleHeight = cylinderHeight - cylinderBottomThickness;
        G4double sampleCenterZ = cylinderPositionZ + cylinderBottomThickness + sampleHeight / 2.0;
        
        primary->SetSourceRadius(sampleRadius * mm);
        primary->SetSourceHeight(sampleHeight * mm);
        primary->SetSourcePosition(G4ThreeVector(0, 0, sampleCenterZ * mm));
        
        G4cout << "\n  ✓ Cylinder source enabled:" << G4endl;
        G4cout << "    Inner radius: " << cylinderInnerRadius << " mm" << G4endl;
        G4cout << "    Outer radius: " << cylinderOuterRadius << " mm" << G4endl;
        G4cout << "    Height: " << cylinderHeight << " mm" << G4endl;
        G4cout << "    Wall thickness: " << cylinderWallThickness << " mm" << G4endl;
        G4cout << "    Bottom thickness: " << cylinderBottomThickness << " mm" << G4endl;
        G4cout << "    Position Z: " << cylinderPositionZ << " mm" << G4endl;
        G4cout << "    Wall material: " << cylinderWallMaterial << G4endl;
        G4cout << "    Fill material: " << cylinderFillMaterial << G4endl;
        G4cout << "    === PRIMARY GENERATOR CONFIG ===" << G4endl;
        G4cout << "    Sampling radius: " << sampleRadius << " mm" << G4endl;
        G4cout << "    Sampling height: " << sampleHeight << " mm" << G4endl;
        G4cout << "    Sampling center Z: " << sampleCenterZ << " mm" << G4endl;
        G4cout << "    Z range: [" << (sampleCenterZ - sampleHeight/2.0) << ", " 
               << (sampleCenterZ + sampleHeight/2.0) << "] mm" << G4endl;
    }
    
    // Configure Cartridge source if enabled (Orano LEA Type D)
    if (useCartridgeSource) {
        CartridgeSourceConfig cartConfig;
        cartConfig.isActive = true;
        cartConfig.housingDiameter = cartridgeHousingDiameter * mm;
        cartConfig.housingHeight = cartridgeHousingHeight * mm;
        cartConfig.activeDiameter = cartridgeActiveDiameter * mm;
        cartConfig.activeThickness = cartridgeActiveThickness * mm;
        cartConfig.wallThickness = cartridgeWallThickness * mm;
        cartConfig.positionZ = cartridgePositionZ * mm;
        cartConfig.housingMaterial = cartridgeHousingMaterial;
        cartConfig.activeMaterial = cartridgeActiveMaterial;
        
        detector->SetCartridgeSourceConfig(cartConfig);
        detector->EnableCartridgeSource(true);
        
        // Configure primary generator for cartridge source
        G4double activeR = cartridgeActiveDiameter / 2.0;
        G4double activeCenterZ = cartridgePositionZ + cartridgeHousingHeight / 2.0;
        
        primary->SetSourceRadius(activeR * mm);
        primary->SetSourceHeight(cartridgeActiveThickness * mm);
        primary->SetSourcePosition(G4ThreeVector(0, 0, activeCenterZ * mm));
        
        G4cout << "\n  Cartridge source (Type D) enabled:" << G4endl;
        G4cout << "    Housing: " << cartridgeHousingDiameter << " x " << cartridgeHousingHeight << " mm" << G4endl;
        G4cout << "    Active matrix: " << cartridgeActiveDiameter << " x " << cartridgeActiveThickness << " mm" << G4endl;
        G4cout << "    Position Z: " << cartridgePositionZ << " mm" << G4endl;
    }
    
    // Configure Filter source if enabled (Orano LEA Type M)
    if (useFilterSource) {
        FilterSourceConfig filtConfig;
        filtConfig.isActive = true;
        filtConfig.outerDiameter = filterOuterDiameter * mm;
        filtConfig.activeDiameter = filterActiveDiameter * mm;
        filtConfig.filterThickness = filterThickness * mm;
        filtConfig.sealThickness = filterSealThickness * mm;
        filtConfig.positionZ = filterPositionZ * mm;
        filtConfig.filterMaterial = filterMaterial;
        filtConfig.sealMaterial = filterSealMaterial;
        filtConfig.filterType = filterType;
        
        detector->SetFilterSourceConfig(filtConfig);
        detector->EnableFilterSource(true);
        
        // Configure primary generator for filter source
        G4double activeR = filterActiveDiameter / 2.0;
        
        primary->SetSourceRadius(activeR * mm);
        primary->SetSourceHeight(filterThickness * mm);
        primary->SetSourcePosition(G4ThreeVector(0, 0, filterPositionZ * mm));
        
        G4cout << "\n  Filter source (Type " << filterType << ") enabled:" << G4endl;
        G4cout << "    Outer diameter: " << filterOuterDiameter << " mm" << G4endl;
        G4cout << "    Active diameter: " << filterActiveDiameter << " mm" << G4endl;
        G4cout << "    Position Z: " << filterPositionZ << " mm" << G4endl;
    }
    
    // Apply custom energies if specified
    if (isCustomEnergy && !customEnergies.empty()) {
        primary->ClearGammaLines();
        
        if (customIntensities.empty() || customIntensities.size() != customEnergies.size()) {
            G4cout << "  Using equal intensities for all energies" << G4endl;
            customIntensities.resize(customEnergies.size(), 100.0);
        }
        
        G4cout << "  Custom gamma lines:" << G4endl;
        for (size_t i = 0; i < customEnergies.size(); i++) {
            primary->AddGammaLine(customEnergies[i], customIntensities[i]);
            G4cout << "    " << customEnergies[i] << " keV (intensity: " 
                   << customIntensities[i] << "%)" << G4endl;
        }
    }
    
    // ============================================================================
    // AUTOMATIC CASCADE MODE DETECTION
    // Mode A (predefined isotope): Cascade mode based on isotope type
    // Mode B (Custom): Always single gamma emission (no cascade)
    // ============================================================================
    
    if (isCustomEnergy) {
        // Custom mode = Mode B = Independent single gamma emission
        // Cascade mode is ALWAYS disabled for Custom
        primary->SetCascadeMode(false);
        G4cout << "\n  *** CUSTOM MODE: Independent single-gamma emission ***" << G4endl;
        G4cout << "  Each energy line emitted independently (no cascades)" << G4endl;
    }
    else if (!isotopeName.empty()) {
        // Predefined isotope = Mode A = Realistic emission
        // Cascade mode depends on whether isotope has cascades
        
        // List of cascade isotopes
        bool hasCascade = (isotopeName == "Co60" || 
                          isotopeName == "Y88" || 
                          isotopeName == "Na22" || 
                          isotopeName == "Eu152" || 
                          isotopeName == "Ba133");
        
        primary->SetCascadeMode(hasCascade);
        
        if (hasCascade) {
            G4cout << "\n  *** CASCADE ISOTOPE DETECTED: " << isotopeName << " ***" << G4endl;
            G4cout << "  Cascade emission ENABLED (multiple gammas per decay)" << G4endl;
        } else {
            G4cout << "\n  *** SINGLE-GAMMA ISOTOPE: " << isotopeName << " ***" << G4endl;
            G4cout << "  Single gamma emission (no cascades)" << G4endl;
        }
    }
    
    // Set source position (only for point source - NOT for Marinelli, Disk, or Cylinder)
    // These special geometries have their own position configuration above
    if (!useMarinelli && !useDiskSource && !useCylinderSource && !useCartridgeSource && !useFilterSource) {
        primary->SetSourcePosition(G4ThreeVector(srcX*cm, srcY*cm, srcZ*cm));
        detector->SetSourcePosition(srcX*cm, srcY*cm, srcZ*cm);
        G4cout << "  Source position: (" << srcX << ", " << srcY << ", " << srcZ << ") cm" << G4endl;
    }
    
    // Configure expected energies for efficiency calculation
    const auto& gammaLines = primary->GetGammaLines();
    std::vector<G4double> expectedEnergies;
    std::vector<G4double> expectedIntensities;
    
    for (const auto& gammaLine : gammaLines) {
        expectedEnergies.push_back(gammaLine.energy * keV);
        expectedIntensities.push_back(gammaLine.intensity);
    }
    
    if (eventAction && !expectedEnergies.empty()) {
        eventAction->SetExpectedEnergies(expectedEnergies, expectedIntensities);
    }
    
    // ============================================================================
    // ENHANCED TCS CONFIGURATION with Dynamic Solid Angle
    // ============================================================================
    if (enableTCS && eventAction) {
        HPGeTCSManager* tcsManager = new HPGeTCSManager();
        
        G4cout << "\n=== Configuring Enhanced TCS System ===" << G4endl;
        
        // ---------------------------------------------------------------
        // NEW: Configure Geometry for Dynamic Solid Angle Calculation
        // ---------------------------------------------------------------
        GeometryConfig geomConfig;
        
        // Set detector dimensions from detector construction
        const std::vector<HPGeConfig>& detConfigs = detector->GetDetectorConfigs();
        if (!detConfigs.empty()) {
            const HPGeConfig& detConfig = detConfigs[0];
            geomConfig.detectorRadius = detConfig.crystalDiameter / 2.0;
            geomConfig.detectorHeight = detConfig.crystalLength;
            
            G4cout << "Detector Configuration:" << G4endl;
            G4cout << "  Radius: " << geomConfig.detectorRadius/mm << " mm" << G4endl;
            G4cout << "  Height: " << geomConfig.detectorHeight/mm << " mm" << G4endl;
        } else {
            // Default values if not available
            geomConfig.detectorRadius = 30.0*mm;
            geomConfig.detectorHeight = 60.0*mm;
            G4cout << "  Using default detector dimensions" << G4endl;
        }
        
        // Configure based on source type
        if (sourceType == "point") {
            geomConfig.sourceType = "Point";
            
            // Calculate distance from source position to detector face
            G4ThreeVector sourcePos = detector->GetSourcePosition();
            G4double detectorZ = 0.0;  // Assuming detector is at z=0
            geomConfig.sourceDistance = std::abs(sourcePos.z() - detectorZ);
            
            // Use configured point source distance if available
            if (pointSourceDistance > 0) {
                geomConfig.sourceDistance = pointSourceDistance * mm;
            }
            
            G4cout << "\nPoint Source Configuration:" << G4endl;
            G4cout << "  Distance from detector: " << geomConfig.sourceDistance/mm << " mm" << G4endl;
            
        } else if (sourceType == "marinelli") {
            
            // Get Marinelli configuration from detector construction
            MarinelliConfig marinelliConfig = detector->GetMarinelliConfig();
            
            if (marinelliType == "1000ml") {
                geomConfig.sourceType = "Marinelli_1000mL";
                
                // 1000 mL Marinelli dimensions
                geomConfig.marinelliInnerRadius = marinelliConfig.innerDiameter / 2.0;
                geomConfig.marinelliOuterRadius = marinelliConfig.outerDiameter / 2.0;
                geomConfig.marinelliHeight = marinelliConfig.totalHeight;
                geomConfig.marinelliInnerHeight = marinelliConfig.innerHeight;
                
            } else if (marinelliType == "450ml") {
                geomConfig.sourceType = "Marinelli_450mL";
                
                // 450 mL Marinelli dimensions
                geomConfig.marinelliInnerRadius = marinelliConfig.innerDiameter / 2.0;
                geomConfig.marinelliOuterRadius = marinelliConfig.outerDiameter / 2.0;
                geomConfig.marinelliHeight = marinelliConfig.totalHeight;
                geomConfig.marinelliInnerHeight = marinelliConfig.innerHeight;
            }
            
            G4cout << "\nMarinelli Beaker Configuration:" << G4endl;
            G4cout << "  Type: " << marinelliType << G4endl;
            G4cout << "  Inner Radius: " << geomConfig.marinelliInnerRadius/mm << " mm" << G4endl;
            G4cout << "  Outer Radius: " << geomConfig.marinelliOuterRadius/mm << " mm" << G4endl;
            G4cout << "  Height: " << geomConfig.marinelliHeight/mm << " mm" << G4endl;
            G4cout << "  Inner Height: " << geomConfig.marinelliInnerHeight/mm << " mm" << G4endl;
            
        } else if (sourceType == "disk") {
            // Disk source configuration for TCS
            geomConfig.sourceType = "Disk";
            
            // Get disk configuration
            DiskSourceConfig diskConfig = detector->GetDiskSourceConfig();
            
            // Calculate effective distance (from disk center to detector)
            geomConfig.sourceDistance = std::abs(diskConfig.positionZ);
            geomConfig.sourceRadius = diskConfig.radius;
            geomConfig.sourceHeight = diskConfig.thickness;
            
            G4cout << "\nDisk Source Configuration (TCS):" << G4endl;
            G4cout << "  Radius: " << diskConfig.radius/mm << " mm" << G4endl;
            G4cout << "  Thickness: " << diskConfig.thickness/mm << " mm" << G4endl;
            G4cout << "  Distance: " << geomConfig.sourceDistance/mm << " mm" << G4endl;
            
        } else if (sourceType == "volume") {
            // Cylinder volume source configuration for TCS
            geomConfig.sourceType = "Cylinder";
            
            // Get cylinder configuration
            CylinderSourceConfig cylConfig = detector->GetCylinderSourceConfig();
            
            // Calculate sample dimensions
            G4double sampleRadius = cylConfig.outerRadius - cylConfig.wallThickness;
            G4double sampleHeight = cylConfig.height - cylConfig.bottomThickness;
            G4double sampleCenterZ = cylConfig.positionZ + cylConfig.bottomThickness + sampleHeight / 2.0;
            
            geomConfig.sourceDistance = std::abs(sampleCenterZ);
            geomConfig.sourceRadius = sampleRadius;
            geomConfig.sourceHeight = sampleHeight;
            
            G4cout << "\nCylinder Source Configuration (TCS):" << G4endl;
            G4cout << "  Sample Radius: " << sampleRadius/mm << " mm" << G4endl;
            G4cout << "  Sample Height: " << sampleHeight/mm << " mm" << G4endl;
            G4cout << "  Sample Center Z: " << sampleCenterZ/mm << " mm" << G4endl;
            G4cout << "  Distance: " << geomConfig.sourceDistance/mm << " mm" << G4endl;
        }
        
        // Apply geometry configuration to TCS Manager
        // This will automatically calculate the solid angle
        tcsManager->SetGeometryConfig(geomConfig);
        
        G4cout << "\n✓ Solid Angle Calculated: " << tcsManager->GetSolidAngle()*100 
               << "% (Ω/4π)" << G4endl;
        
        // ---------------------------------------------------------------
        // Load isotope decay scheme
        // ---------------------------------------------------------------
        if (!isotopeName.empty()) {
            tcsManager->LoadIsotopeScheme(isotopeName);
            G4cout << "✓ TCS decay scheme loaded for: " << isotopeName << G4endl;
        }
        
        // Set coincidence window
        tcsManager->SetCoincidenceWindow(tcsCoincidenceWindow * ns);
        
        // Pass TCS manager to event action
        eventAction->SetTCSManager(tcsManager);
        eventAction->EnableTCS(true);
        
        G4cout << "✓ True Coincidence Summing (TCS) enabled" << G4endl;
        G4cout << "  Coincidence window: " << tcsCoincidenceWindow << " ns" << G4endl;
        G4cout << "======================================" << G4endl;
    }
    
    file.close();
    G4cout << "=== Configuration Loaded ===" << G4endl;
}

// ==================================================================
// Quick config pre-scan for fast_mode (before physics list creation)
// ==================================================================
G4bool PreScanFastMode(const G4String& configFile) {
    std::ifstream file(configFile);
    if (!file.is_open()) return false;
    
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream iss(line);
        std::string key, equals, value;
        if (!(iss >> key >> equals >> value)) continue;
        
        if (key == "fast_mode" && (value == "true" || value == "1")) {
            file.close();
            return true;
        }
    }
    file.close();
    return false;
}

// ==================================================================
// Main function (ENHANCED with TCS v2.0)
// ==================================================================
int main(int argc, char** argv) {
    
    G4cout << "\n" << std::string(70, '=') << G4endl;
    G4cout << "HPGe Monte Carlo Simulation - Enhanced TCS v2.0" << G4endl;
    G4cout << "Features:" << G4endl;
    G4cout << "  • Dynamic Solid Angle Calculation" << G4endl;
    G4cout << "  • Efficiency-Based COI Method" << G4endl;
    G4cout << "  • Marinelli Beaker Support" << G4endl;
    G4cout << "  • ISOCS/LabSOCS Style Analysis" << G4endl;
    G4cout << "  • Fast Mode Option" << G4endl;
    G4cout << std::string(70, '=') << G4endl;
    
    // Get configuration file name early
    G4String configFile = "config.txt";
    if (argc > 2) {
        configFile = argv[2];
    }
    
    // Pre-scan config for fast_mode setting
    G4bool fastMode = PreScanFastMode(configFile);
    if (fastMode) {
        G4cout << "\n*** FAST MODE DETECTED IN CONFIG ***" << G4endl;
    }
    
    // UI Executive for interactive mode
    G4UIExecutive* ui = nullptr;
    if (argc == 1) {
        ui = new G4UIExecutive(argc, argv);
    }
    
    // Run Manager
    G4RunManager* runManager = new G4RunManager;
    
    // Mandatory Initialization Classes
    HPGeDetectorConstruction* detector = new HPGeDetectorConstruction();
    runManager->SetUserInitialization(detector);
    
    // Physics list with fast mode option
    HPGePhysicsList* physics = new HPGePhysicsList(fastMode);
    runManager->SetUserInitialization(physics);
    
    HPGePrimaryGeneratorAction* primary = new HPGePrimaryGeneratorAction();
    runManager->SetUserAction(primary);
    
    // ============================================================================
    // CRITICAL: Load Configuration BEFORE Initialize()
    // This allows Marinelli configuration to affect geometry construction
    // ============================================================================
    
    // Create temporary event action for pre-initialization configuration
    std::vector<HPGeSensitiveDetector*> tempSDs;
    HPGeEventAction* tempEventAction = new HPGeEventAction(tempSDs);
    
    // Load configuration FIRST - this will enable Marinelli if configured
    G4cout << "\n=== Pre-Initialize Configuration ===" << G4endl;
    LoadConfiguration(configFile, detector, primary, tempEventAction);
    
    delete tempEventAction;  // Clean up temporary event action
    
    // NOW Initialize geometry and physics
    G4cout << "\n=== Initializing Geant4 ===" << G4endl;
    runManager->Initialize();
    G4cout << "=== Initialization Complete ===" << G4endl;
    
    // Retrieve sensitive detectors (only available after Initialize)
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    std::vector<HPGeSensitiveDetector*> sensitiveDetectors;
    
    for (size_t i = 0; i < detector->GetDetectorConfigs().size(); i++) {
        G4String sdName = "HPGe_SD_" + std::to_string(i);
        HPGeSensitiveDetector* sd = dynamic_cast<HPGeSensitiveDetector*>(
            sdManager->FindSensitiveDetector(sdName)
        );
        if (sd) {
            sensitiveDetectors.push_back(sd);
        } else {
            G4cout << "✗ ERROR: Could not find sensitive detector: " << sdName << G4endl;
        }
    }
    
    if (sensitiveDetectors.empty()) {
        G4cout << "✗ FATAL: No sensitive detectors found!" << G4endl;
        return 1;
    }
    
    // Create the real Event Action with actual sensitive detectors
    HPGeEventAction* eventAction = new HPGeEventAction(sensitiveDetectors);
    eventAction->SetOutputFile("output.dat");
    runManager->SetUserAction(eventAction);
    
    // Reload configuration to set event action parameters and TCS
    G4cout << "\n=== Post-Initialize Configuration ===" << G4endl;
    LoadConfiguration(configFile, detector, primary, eventAction);
    
    // ============================================================================
    // Print simulation summary
    // ============================================================================
    G4cout << "\n" << std::string(70, '=') << G4endl;
    G4cout << "SIMULATION SETUP SUMMARY" << G4endl;
    G4cout << std::string(70, '=') << G4endl;
    G4cout << "Detectors:       " << sensitiveDetectors.size() << G4endl;
    G4cout << "Output file:     output.dat" << G4endl;
    G4cout << "Efficiency file: efficiency.dat (will be created)" << G4endl;
    
    // Detector information
    for (size_t i = 0; i < detector->GetDetectorConfigs().size(); i++) {
        G4double volume = detector->GetActiveVolume(i);
        G4double mass = detector->GetActiveMass(i);
        G4double distance = detector->GetCrystalToSourceDistance(i);
        
        G4cout << "\nDetector " << i << ":" << G4endl;
        G4cout << "  Active volume:  " << volume/cm3 << " cm³" << G4endl;
        G4cout << "  Active mass:    " << mass/g << " g" << G4endl;
        G4cout << "  Distance:       " << distance/cm << " cm" << G4endl;
    }
    
    // Marinelli status
    if (detector->IsMarinelliActive()) {
        G4cout << "\n✓ Marinelli beaker ACTIVE" << G4endl;
        auto marinelliConfig = detector->GetMarinelliConfig();
        G4cout << "  Type: " << marinelliConfig.beakerType << G4endl;
        G4cout << "  Fill material: " << marinelliConfig.fillMaterial << G4endl;
        G4cout << "  Z position: " << marinelliConfig.GetZPosition()/mm << " mm" << G4endl;
    }
    
    primary->PrintConfiguration();
    
    G4cout << std::string(70, '=') << G4endl;
    
    // Visualization
    G4VisManager* visManager = new G4VisExecutive;
    visManager->Initialize();
    
    // UI Manager
    G4UImanager* UImanager = G4UImanager::GetUIpointer();
    
    if (!ui) {
        // Batch mode
        G4String command = "/control/execute ";
        G4String fileName = argv[1];
        UImanager->ApplyCommand(command + fileName);
    }
    else {
        // Interactive mode
        UImanager->ApplyCommand("/control/execute vis.mac");
        ui->SessionStart();
        delete ui;
    }
    
    // *** CRITICAL FIX: Update emission counts from primary generator ***
    // This ensures totalEmitted reflects actual emissions, not just total events
    eventAction->SetEmissionCounts(primary->GetEnergyStatistics());
    
    eventAction->ExportEfficiencyData("efficiency.dat");

    
    // Cleanup
    delete visManager;
    delete runManager;
    
    return 0;
}
