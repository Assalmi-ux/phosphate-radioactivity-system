// ==================================================================
// HPGeDetectorConstruction.cc (Enhanced with Canberra Mirion Real Geometry)
// ==================================================================

#include "HPGeDetectorConstruction.hh"
#include "HPGeSensitiveDetector.hh"

#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"

#include <cmath>

// ==================================================================
// Constructor
// ==================================================================
HPGeDetectorConstruction::HPGeDetectorConstruction()
: G4VUserDetectorConstruction(),
  fSourcePosition(0., 0., 5.*cm),
  fWorldSize(1.*m),
  fGe(nullptr), fAl(nullptr), fC(nullptr), fVacuum(nullptr), fAir(nullptr), 
  fLi(nullptr), fBoron(nullptr), fLead(nullptr), fWater(nullptr), 
  fSteel(nullptr), fPerspex(nullptr), fPolypropylene(nullptr), fSoil(nullptr),
  fGlass(nullptr), fPVC(nullptr), fHDPE(nullptr), fCarbonEpoxy(nullptr),
  fPolycarbonate(nullptr), fPolyester(nullptr), fActivatedCarbon(nullptr),
  fUseCanberraGeometry(true),
  fMarinelliLogical(nullptr), fWorldPhysical(nullptr),
  fDiskSourceLogical(nullptr), fCylinderSourceLogical(nullptr), fCylinderWallLogical(nullptr),
  fCartridgeSourceLogical(nullptr), fCartridgeHousingLogical(nullptr),
  fFilterSourceLogical(nullptr), fFilterSealLogical(nullptr) {
    
    DefineMaterials();
    
    // Default: use Canberra Mirion S/N:21206 real geometry
    HPGeConfig canberraConfig;
    SetCanberraGeometryParameters(canberraConfig);
    fDetectorConfigs.push_back(canberraConfig);
    
    // Initialize Marinelli config
    fMarinelliConfig.isActive = false;
    
    // Initialize Disk source config
    fDiskSourceConfig.isActive = false;
    
    // Initialize Cylinder source config
    fCylinderSourceConfig.isActive = false;
    
    // Initialize Cartridge source config (Type D)
    fCartridgeSourceConfig.isActive = false;
    
    // Initialize Filter source config (Type M53)
    fFilterSourceConfig.isActive = false;
}

// ==================================================================
// Destructor
// ==================================================================
HPGeDetectorConstruction::~HPGeDetectorConstruction() {}

// ==================================================================
// Define materials
// ==================================================================
void HPGeDetectorConstruction::DefineMaterials() {
    G4NistManager* nist = G4NistManager::Instance();
    
    // Germanium
    fGe = nist->FindOrBuildMaterial("G4_Ge");
    
    // Aluminum
    fAl = nist->FindOrBuildMaterial("G4_Al");
    
    fC = nist->FindOrBuildMaterial("G4_C");
    // Air
    fAir = nist->FindOrBuildMaterial("G4_AIR");
    
    // Vacuum
    fVacuum = new G4Material("Vacuum", 1., 1.008*g/mole, 1.e-25*g/cm3,
                            kStateGas, 0.1*kelvin, 1.e-19*pascal);
    
    // Lithium
    fLi = nist->FindOrBuildMaterial("G4_Li");
    
    // Boron
    fBoron = nist->FindOrBuildMaterial("G4_B");
    
    // Lead
    fLead = nist->FindOrBuildMaterial("G4_Pb");
    
    // Water
    fWater = nist->FindOrBuildMaterial("G4_WATER");
    
    // Stainless Steel
    fSteel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    
    // Polymethylmethacrylate (Perspex)
    G4Element* H = nist->FindOrBuildElement("H");
    G4Element* C = nist->FindOrBuildElement("C");
    G4Element* O = nist->FindOrBuildElement("O");
    
    fPerspex = new G4Material("Perspex", 1.19*g/cm3, 3);
    fPerspex->AddElement(H, 8);
    fPerspex->AddElement(C, 5);
    fPerspex->AddElement(O, 2);
    
    // Polypropylene (C3H6) for Marinelli beaker
    fPolypropylene = new G4Material("Polypropylene", 0.9*g/cm3, 2);
    fPolypropylene->AddElement(C, 3);
    fPolypropylene->AddElement(H, 6);
    
    // Soil (simplified composition)
    G4Element* Si = nist->FindOrBuildElement("Si");
    G4Element* Al = nist->FindOrBuildElement("Al");
    G4Element* Fe = nist->FindOrBuildElement("Fe");
    
    fSoil = new G4Material("Soil", 1.5*g/cm3, 5);
    fSoil->AddElement(Si, 45);
    fSoil->AddElement(O, 40);
    fSoil->AddElement(Al, 8);
    fSoil->AddElement(Fe, 5);
    fSoil->AddElement(C, 2);
    
    // Glass (Borosilicate - simplified for containers)
    G4Element* B = nist->FindOrBuildElement("B");
    G4Element* Na = nist->FindOrBuildElement("Na");
    
    fGlass = new G4Material("Glass", 2.23*g/cm3, 5);
    fGlass->AddElement(Si, 80);
    fGlass->AddElement(O, 12);
    fGlass->AddElement(B, 4);
    fGlass->AddElement(Na, 3);
    fGlass->AddElement(Al, 1);
    
    // PVC (for some containers)
    G4Element* Cl = nist->FindOrBuildElement("Cl");
    
    fPVC = new G4Material("PVC", 1.4*g/cm3, 3);
    fPVC->AddElement(C, 2);
    fPVC->AddElement(H, 3);
    fPVC->AddElement(Cl, 1);
    
    // HDPE (High-density polyethylene)
    fHDPE = new G4Material("HDPE", 0.95*g/cm3, 2);
    fHDPE->AddElement(C, 2);
    fHDPE->AddElement(H, 4);
    
    // Carbon Epoxy (window material for Canberra detector)
    G4Element* N = nist->FindOrBuildElement("N");
    fCarbonEpoxy = new G4Material("CarbonEpoxy", 1.55*g/cm3, 4);
    fCarbonEpoxy->AddElement(C, 76);
    fCarbonEpoxy->AddElement(H, 8);
    fCarbonEpoxy->AddElement(O, 12);
    fCarbonEpoxy->AddElement(N, 4);
    
    // Polycarbonate (C16H14O3) - housing for Orano LEA sources
    fPolycarbonate = new G4Material("Polycarbonate", 1.20*g/cm3, 3);
    fPolycarbonate->AddElement(C, 16);
    fPolycarbonate->AddElement(H, 14);
    fPolycarbonate->AddElement(O, 3);
    
    // Polyester (PET, C10H8O4) - sealing for paper filter sources
    fPolyester = new G4Material("Polyester", 1.38*g/cm3, 3);
    fPolyester->AddElement(C, 10);
    fPolyester->AddElement(H, 8);
    fPolyester->AddElement(O, 4);
    
    // Activated Carbon (amorphous carbon, porous)
    fActivatedCarbon = new G4Material("ActivatedCarbon", 0.45*g/cm3, 1);
    fActivatedCarbon->AddElement(C, 1);
    
    G4cout << "\n=== Materials Defined ===" << G4endl;
    G4cout << "Germanium density: " << fGe->GetDensity()/(g/cm3) << " g/cm³" << G4endl;
    G4cout << "Polypropylene density: " << fPolypropylene->GetDensity()/(g/cm3) << " g/cm³" << G4endl;
    G4cout << "Glass density: " << fGlass->GetDensity()/(g/cm3) << " g/cm³" << G4endl;
    G4cout << "PVC density: " << fPVC->GetDensity()/(g/cm3) << " g/cm³" << G4endl;
    G4cout << "HDPE density: " << fHDPE->GetDensity()/(g/cm3) << " g/cm³" << G4endl;
}

// ==================================================================
// Set FLUKA geometry parameters
// ==================================================================
void HPGeDetectorConstruction::SetFLUKAGeometryParameters(HPGeConfig& config) {
    
    G4double cryst_dia = 65 * mm;
    G4double cryst_ht = 45.3 * mm;
    G4double hole_dia = 12 * mm;
    
    G4double dead_lat = 1 * mm;
    G4double dead_Front = 0.9 * mm;
    G4double dead_inL = 1 * mm;
    G4double dead_inF = 2 * mm;
    
    G4double end_thk = 1.5 * mm;
    G4double end_thk1 = 0.6 * mm;
    G4double vac_gap = 7.5 * mm;
    
    G4double boron_thk = 0.03 * mm;
    G4double hole_depth = 15 * mm;
    
    config.detectorName = "HPGe_FLUKA";
    config.detectorID = 0;
    
    config.crystalDiameter = cryst_dia;
    config.crystalLength = cryst_ht;
    config.holeDiameter = hole_dia;
    config.holeDepth = hole_depth;
    
    config.geDeadLayer = dead_lat / 2.0;
    config.geDeadLayerFront = dead_Front;
    config.liDeadLayer = dead_inL / 2.0;
    config.liDeadLayerFront = dead_inF;
    
    config.alWindowThickness = end_thk1;
    config.alCupThickness = end_thk;
    config.vacuumGap = vac_gap;
    
    config.boronThickness = boron_thk;
    
    config.posX = 0.0;
    config.posY = 0.0;
    config.posZ = 0.0;
    
    G4cout << "\n=== FLUKA Geometry Parameters Loaded ===" << G4endl;
    G4cout << "Crystal: Ø" << cryst_dia/mm << " mm × " << cryst_ht/mm << " mm" << G4endl;
}

// ==================================================================
// Set Canberra Mirion S/N:21206 geometry parameters
// ==================================================================
void HPGeDetectorConstruction::SetCanberraGeometryParameters(HPGeConfig& config) {
    
    config.detectorName = "HPGe_Canberra_21206";
    config.detectorID = 0;
    
    // Crystal dimensions from characterization sheet
    config.crystalDiameter = 60.5 * mm;
    config.crystalLength = 45.4 * mm;
    config.holeDiameter = 9.5 * mm;
    config.holeDepth = 19.0 * mm;
    config.crystalCornerRadius = 1.0 * mm;
    config.crystalBottomOuterStep = 28.0 * mm;
    config.crystalBottomInnerStep = 20.0 * mm;
    
    // Dead layers from characterization sheet
    config.geDeadLayer = 0.5 * mm;             // Outer Electrode Thickness: 0.5 mm
    config.geDeadLayerFront = 0.0004 * mm;     // Front Dead Layer: 0.4 um
    config.liDeadLayer = 0.0003 * mm;          // Inner Electrode: 0.3 um eq. Ge
    config.liDeadLayerFront = 0.0003 * mm;
    
    // Endcap/Housing from characterization sheet
    config.alWindowThickness = 0.6 * mm;       // Window from assembly drawing
    config.alCupThickness = 1.5 * mm;
    config.vacuumGap = 6.0 * mm;              // Ge front to center endcap outside
    config.endcapDiameter = 76.2 * mm;         // Assembly outer diameter
    config.endcapLength = 8.0 * mm;
    config.endcapMaterial = "Al";
    config.windowMaterial = "carbon_epoxy";     // Carbon Epoxy window
    config.geToEndcapDistance = 6.0 * mm;
    
    // Holder from characterization sheet
    config.holderOuterDiameter1 = 69.0 * mm;
    config.holderOuterDiameter2 = 64.6 * mm;
    config.holderOuterDiameter3 = 63.0 * mm;
    config.holderInnerDepth1 = 66.7 * mm;
    config.holderInnerDepth2 = 69.8 * mm;
    config.holderTotalLength = 95.2 * mm;
    config.holderTopHeight = 8.6 * mm;
    config.holderStepHeight = 19.4 * mm;
    config.holderMaterial = "Al";
    
    // HDPE spacer
    config.hdpeThickness = 6.5 * mm;
    
    config.boronThickness = 0.03 * mm;
    
    config.posX = 0.0;
    config.posY = 0.0;
    config.posZ = 0.0;
    
    G4cout << "\n=== Canberra Mirion S/N:21206 Geometry Parameters Loaded ===" << G4endl;
    G4cout << "Crystal: Ø" << config.crystalDiameter/mm << " mm × " 
           << config.crystalLength/mm << " mm" << G4endl;
    G4cout << "Hole: Ø" << config.holeDiameter/mm << " mm, depth " 
           << config.holeDepth/mm << " mm" << G4endl;
    G4cout << "Front Dead Layer: " << config.geDeadLayerFront/mm * 1000.0 << " um" << G4endl;
    G4cout << "Outer Electrode: " << config.geDeadLayer/mm << " mm" << G4endl;
    G4cout << "Inner Electrode: " << config.liDeadLayer/mm * 1000.0 << " um eq. Ge" << G4endl;
    G4cout << "Endcap: Ø" << config.endcapDiameter/mm << " mm, " 
           << config.endcapMaterial << G4endl;
    G4cout << "Window: " << config.windowMaterial << ", " 
           << config.alWindowThickness/mm << " mm" << G4endl;
    G4cout << "Ge to endcap: " << config.geToEndcapDistance/mm << " mm" << G4endl;
    G4cout << "Holder: " << config.holderMaterial << ", Ø" 
           << config.holderOuterDiameter1/mm << " mm" << G4endl;
}

// ==================================================================
// Construct world and detectors
// ==================================================================
G4VPhysicalVolume* HPGeDetectorConstruction::Construct() {
    
    G4double worldSize = 2.0 * m;
    G4Box* solidWorld = new G4Box("World", worldSize/2, worldSize/2, worldSize/2);
    G4LogicalVolume* logicWorld = new G4LogicalVolume(solidWorld, fAir, "World");
    fWorldPhysical = new G4PVPlacement(0, G4ThreeVector(), logicWorld,
                                       "World", 0, false, 0);
    
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    G4cout << "\n=== Constructing HPGe Detector ===" << G4endl;
    
    ConstructLeadCastle(logicWorld);
    
    for(size_t i = 0; i < fDetectorConfigs.size(); i++) {
        const HPGeConfig& config = fDetectorConfigs[i];
        
        G4cout << "\nDetector " << i << ": " << config.detectorName << G4endl;
        
        G4LogicalVolume* detLogic;
        if (fUseCanberraGeometry) {
            detLogic = ConstructHPGeDetectorCanberra(config);
        } else {
            detLogic = ConstructHPGeDetectorFLUKA(config);
        }
        
        G4RotationMatrix* rotation = new G4RotationMatrix();
        rotation->rotateX(180.*deg);
        
        G4double crystalH = config.crystalLength / 2.0;
        G4double windowH = config.alWindowThickness / 2.0;
        G4double windowZ = crystalH + config.vacuumGap + windowH;
        
        G4double zShift = windowZ;
        
        G4ThreeVector detPos(config.posX, config.posY, config.posZ + zShift);
        new G4PVPlacement(rotation, detPos, detLogic, config.detectorName, 
                         logicWorld, false, i);
        
        G4cout << "  Detector positioned with window at Z=0" << G4endl;
        
        G4double activeVol = CalculateActiveVolume(config);
        G4cout << "  Active volume: " << activeVol/cm3 << " cm³" << G4endl;
    }
    
    // Construct Marinelli beaker if enabled
    if (fMarinelliConfig.isActive) {
        ConstructMarinelliBeaker(logicWorld);
    }
    
    // Construct Disk source if enabled
    if (fDiskSourceConfig.isActive) {
        ConstructDiskSource(logicWorld);
    }
    
    // Construct Cylinder source if enabled
    if (fCylinderSourceConfig.isActive) {
        ConstructCylinderSource(logicWorld);
    }
    
    // Construct Cartridge source (Type D) if enabled
    if (fCartridgeSourceConfig.isActive) {
        ConstructCartridgeSource(logicWorld);
    }
    
    // Construct Filter source (Type M53) if enabled
    if (fFilterSourceConfig.isActive) {
        ConstructFilterSource(logicWorld);
    }
    
    G4cout << "\n=== Detector Construction Complete ===" << G4endl;
    
    return fWorldPhysical;
}

// ==================================================================
// Construct Marinelli beaker
// ==================================================================
void HPGeDetectorConstruction::ConstructMarinelliBeaker(G4LogicalVolume* worldLogical) {
    
    G4cout << "\n=== Constructing Marinelli Beaker (FINAL CORRECTED - ROTATED 180°) ===" << G4endl;
    G4cout << "Type: " << fMarinelliConfig.beakerType << G4endl;
    G4cout << "Fill material: " << fMarinelliConfig.fillMaterial << G4endl;
    
    G4double D = fMarinelliConfig.outerDiameter;
    G4double d = fMarinelliConfig.innerDiameter;
    G4double H = fMarinelliConfig.totalHeight;
    G4double h = fMarinelliConfig.innerHeight;
    G4double wall = fMarinelliConfig.wallThickness;
    
    G4cout << "  Outer diameter: " << D/mm << " mm" << G4endl;
    G4cout << "  Inner diameter: " << d/mm << " mm" << G4endl;
    G4cout << "  Total height: " << H/mm << " mm" << G4endl;
    G4cout << "  Inner height (well depth): " << h/mm << " mm" << G4endl;
    
    // ========================================================================
    // CALCUL DES POSITIONS (INVERSÉ - ROTATION 180°)
    // ========================================================================
    // Référence: Fenêtre détecteur à Z = 0, cristal dans Z+
    // Gap entre échantillon et fenêtre: 2 mm
    
    G4double airGap = -3.0*mm;
    
    // APRÈS ROTATION 180°:
    // - Le puits (annulaire) sera en Z positif (entoure le détecteur)
    // - Le fond plein (top) sera en Z négatif (en bas)
    
    // Position du haut du beaker (maintenant en haut après rotation)
    G4double beakerTop = airGap + h;
    
    // Position du centre du beaker pour le placement
    G4double beakerCenterZ = beakerTop - H/2.0;
    
    G4cout << "  Beaker top Z: " << beakerTop/mm << " mm" << G4endl;
    G4cout << "  Beaker center Z: " << beakerCenterZ/mm << " mm" << G4endl;
    G4cout << "  Beaker bottom Z: " << (beakerTop - H)/mm << " mm" << G4endl;
    G4cout << "  Sample starts at Z: " << (beakerTop - wall)/mm << " mm" << G4endl;
    G4cout << "  Sample bottom (at window level) at Z: " << airGap/mm << " mm" << G4endl;
    
    // ========================================================================
    // POLYPROPYLENE SHELL (murs du beaker) - AVEC ROTATION 180°
    // ========================================================================
    
    // Corps extérieur complet
    G4Tubs* solidOuterBody = new G4Tubs("MarinelliOuter",
                                        0,
                                        D/2.0,
                                        H/2.0,
                                        0, twopi);
    
    // Cavité interne pour le puits (maintenant vers le haut après rotation)
    G4Tubs* solidInnerNeck = new G4Tubs("MarinelliInnerNeck",
                                        0,
                                        d/2.0,
                                        h/2.0,
                                        0, twopi);
    
    // Cavité extérieure (espace échantillon)
    G4double outerCavityR = D/2.0 - wall;
    G4double outerCavityH = H - wall;
    
    G4Tubs* solidOuterCavity = new G4Tubs("MarinelliOuterCavity",
                                          d/2.0 + wall,
                                          outerCavityR,
                                          outerCavityH/2.0,
                                          0, twopi);
    
    // Soustraire le puits interne (maintenant du côté positif)
    G4SubtractionSolid* solidBeakerShell1 = new G4SubtractionSolid("BeakerShell1",
                                                                    solidOuterBody,
                                                                    solidInnerNeck,
                                                                    0,
                                                                    G4ThreeVector(0, 0, +(H/2.0 - h/2.0))); // INVERSÉ!
    
    // Soustraire la cavité extérieure
    G4SubtractionSolid* solidBeakerShell = new G4SubtractionSolid("BeakerShell",
                                                                   solidBeakerShell1,
                                                                   solidOuterCavity,
                                                                   0,
                                                                   G4ThreeVector(0, 0, -wall/2.0)); // INVERSÉ!
    
    G4LogicalVolume* logicBeakerShell = new G4LogicalVolume(solidBeakerShell,
                                                            fPolypropylene,
                                                            "MarinelliShell");
    
    // Cyan semi-transparent
    G4VisAttributes* beakerVis = new G4VisAttributes(G4Colour(0.5, 0.8, 0.9, 0.3));
    beakerVis->SetForceSolid(true);
    logicBeakerShell->SetVisAttributes(beakerVis);
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, beakerCenterZ),
                     logicBeakerShell,
                     "MarinelliBeaker_phys",
                     worldLogical,
                     false,
                     0);
    
    // ========================================================================
    // RÉGION ÉCHANTILLON (SOURCE) - DEUX PARTIES - INVERSÉES
    // ========================================================================
    
    // Matériau de remplissage
    G4Material* fillMat = fWater;
    if (fMarinelliConfig.fillMaterial == "soil") {
        fillMat = fSoil;
    } else if (fMarinelliConfig.fillMaterial == "air") {
        fillMat = fAir;
    }
    
    // Couleur selon matériau
    G4Colour sampleColor;
    if (fMarinelliConfig.fillMaterial == "water") {
        sampleColor = G4Colour(0.2, 0.5, 1.0, 0.5);  // Bleu
    } else if (fMarinelliConfig.fillMaterial == "soil") {
        sampleColor = G4Colour(0.6, 0.4, 0.2, 0.6);  // Marron
    } else {
        sampleColor = G4Colour(0.8, 0.8, 0.8, 0.3);  // Gris
    }
    
    // ------------------------------------------------------------------------
    // PARTIE 1: RÉGION ANNULAIRE (autour du puits - EN HAUT, Z POSITIF)
    // ------------------------------------------------------------------------
    
    G4double annularRInner = d/2.0 + wall;
    G4double annularROuter = D/2.0 - wall;
    G4double annularHeight = h - wall;  // Hauteur du puits moins paroi
    
    G4Tubs* solidAnnularSample = new G4Tubs("MarinelliSample_Annular",
                                            annularRInner,
                                            annularROuter,
                                            annularHeight/2.0,
                                            0, twopi);
    
    G4LogicalVolume* logicAnnularSample = new G4LogicalVolume(solidAnnularSample,
                                                              fillMat,
                                                              "MarinelliSample_Annular");
    
    G4VisAttributes* annularVis = new G4VisAttributes(sampleColor);
    annularVis->SetForceSolid(true);
    logicAnnularSample->SetVisAttributes(annularVis);
    
    // Position: centré dans la région annulaire (maintenant en HAUT, Z+)
    G4double annularZ = beakerTop - wall - annularHeight/2.0;
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, annularZ),
                     logicAnnularSample,
                     "MarinelliSample_Annular_phys",
                     worldLogical,
                     false,
                     0);
    
    // ------------------------------------------------------------------------
    // PARTIE 2: RÉGION TOP (cylindre plein - EN BAS, Z NÉGATIF)
    // ------------------------------------------------------------------------
    
    G4double topRadius = D/2.0 - wall;
    G4double topHeight = H - h - wall;  // Hauteur restante
    
    G4Tubs* solidTopSample = new G4Tubs("MarinelliSample_Top",
                                        0,              // PAS DE TROU - cylindre plein!
                                        topRadius,
                                        topHeight/2.0,
                                        0, twopi);
    
    G4LogicalVolume* logicTopSample = new G4LogicalVolume(solidTopSample,
                                                          fillMat,
                                                          "MarinelliSample_Top");
    
    G4VisAttributes* topVis = new G4VisAttributes(sampleColor);
    topVis->SetForceSolid(true);
    logicTopSample->SetVisAttributes(topVis);
    
    // Position: en dessous de la région annulaire (maintenant en BAS, Z-)
    G4double topZ = beakerTop - wall - annularHeight - topHeight/2.0;
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, topZ),
                     logicTopSample,
                     "MarinelliSample_Top_phys",
                     worldLogical,
                     false,
                     1);
    
    // Stocker le volume annulaire comme volume principal
    fMarinelliLogical = logicAnnularSample;
    
    // ========================================================================
    // RÉSUMÉ DE LA CONSTRUCTION
    // ========================================================================
    
    G4cout << "\n  ✓ Marinelli beaker construit avec succès (ROTATION 180°)" << G4endl;
    G4cout << "\n  GÉOMÉTRIE (APRÈS ROTATION):" << G4endl;
    G4cout << "    Fenêtre détecteur: Z = 0 mm" << G4endl;
    G4cout << "    Cristal détecteur: Z > 0 (vers le haut)" << G4endl;
    G4cout << "    Gap air: " << airGap/mm << " mm" << G4endl;
    
    G4cout << "\n  PARTIE 1 - Région annulaire (puits entoure détecteur - Z+):" << G4endl;
    G4cout << "    R_intérieur: " << annularRInner/mm << " mm" << G4endl;
    G4cout << "    R_extérieur: " << annularROuter/mm << " mm" << G4endl;
    G4cout << "    Hauteur: " << annularHeight/mm << " mm" << G4endl;
    G4cout << "    Position Z: " << annularZ/mm << " mm (centre)" << G4endl;
    G4cout << "    De Z = " << (annularZ - annularHeight/2.0)/mm 
           << " mm à Z = " << (annularZ + annularHeight/2.0)/mm << " mm" << G4endl;
    
    G4cout << "\n  PARTIE 2 - Région bottom (cylindre plein - Z-):" << G4endl;
    G4cout << "    Rayon: " << topRadius/mm << " mm (PAS DE TROU)" << G4endl;
    G4cout << "    Hauteur: " << topHeight/mm << " mm" << G4endl;
    G4cout << "    Position Z: " << topZ/mm << " mm (centre)" << G4endl;
    G4cout << "    De Z = " << (topZ - topHeight/2.0)/mm 
           << " mm à Z = " << (topZ + topHeight/2.0)/mm << " mm" << G4endl;
    G4cout << "    Bas de l'échantillon: Z = " << (topZ - topHeight/2.0)/mm << " mm" << G4endl;
    
    G4double annularVol = pi*(annularROuter*annularROuter - annularRInner*annularRInner)*annularHeight;
    G4double topVol = pi*topRadius*topRadius*topHeight;
    G4cout << "\n  VOLUMES:" << G4endl;
    G4cout << "    Annulaire: " << annularVol/cm3 << " cm³" << G4endl;
    G4cout << "    Bottom: " << topVol/cm3 << " cm³" << G4endl;
    G4cout << "    Total: " << (annularVol + topVol)/cm3 << " cm³" << G4endl;
}

// ==================================================================
// Set Marinelli configuration
// ==================================================================
void HPGeDetectorConstruction::SetMarinelliConfig(const MarinelliConfig& config) {
    fMarinelliConfig = config;
}

// ==================================================================
// Enable Marinelli
// ==================================================================
void HPGeDetectorConstruction::EnableMarinelli(G4bool enable) {
    fMarinelliConfig.isActive = enable;
}

// ==================================================================
// Set Disk source configuration
// ==================================================================
void HPGeDetectorConstruction::SetDiskSourceConfig(const DiskSourceConfig& config) {
    fDiskSourceConfig = config;
}

// ==================================================================
// Enable Disk source
// ==================================================================
void HPGeDetectorConstruction::EnableDiskSource(G4bool enable) {
    fDiskSourceConfig.isActive = enable;
}

// ==================================================================
// Set Cylinder source configuration
// ==================================================================
void HPGeDetectorConstruction::SetCylinderSourceConfig(const CylinderSourceConfig& config) {
    fCylinderSourceConfig = config;
}

// ==================================================================
// Enable Cylinder source
// ==================================================================
void HPGeDetectorConstruction::EnableCylinderSource(G4bool enable) {
    fCylinderSourceConfig.isActive = enable;
}

// ==================================================================
// Set Cartridge source configuration
// ==================================================================
void HPGeDetectorConstruction::SetCartridgeSourceConfig(const CartridgeSourceConfig& config) {
    fCartridgeSourceConfig = config;
}

// ==================================================================
// Enable Cartridge source
// ==================================================================
void HPGeDetectorConstruction::EnableCartridgeSource(G4bool enable) {
    fCartridgeSourceConfig.isActive = enable;
}

// ==================================================================
// Set Filter source configuration
// ==================================================================
void HPGeDetectorConstruction::SetFilterSourceConfig(const FilterSourceConfig& config) {
    fFilterSourceConfig = config;
}

// ==================================================================
// Enable Filter source
// ==================================================================
void HPGeDetectorConstruction::EnableFilterSource(G4bool enable) {
    fFilterSourceConfig.isActive = enable;
}

// ==================================================================
// Helper: Get material by name
// ==================================================================
G4Material* HPGeDetectorConstruction::GetMaterialByName(const G4String& name) {
    if (name == "water" || name == "Water") return fWater;
    if (name == "soil" || name == "Soil") return fSoil;
    if (name == "air" || name == "Air") return fAir;
    if (name == "polypropylene" || name == "Polypropylene" || name == "PP") return fPolypropylene;
    if (name == "glass" || name == "Glass") return fGlass;
    if (name == "pvc" || name == "PVC") return fPVC;
    if (name == "hdpe" || name == "HDPE") return fHDPE;
    if (name == "perspex" || name == "Perspex" || name == "PMMA") return fPerspex;
    if (name == "steel" || name == "Steel") return fSteel;
    if (name == "aluminum" || name == "Aluminum" || name == "Al") return fAl;
    if (name == "polycarbonate" || name == "Polycarbonate" || name == "PC") return fPolycarbonate;
    if (name == "polyester" || name == "Polyester" || name == "PET") return fPolyester;
    if (name == "activated_carbon" || name == "ActivatedCarbon") return fActivatedCarbon;
    if (name == "carbon_epoxy" || name == "CarbonEpoxy") return fCarbonEpoxy;
    
    G4cout << "WARNING: Unknown material '" << name << "'. Using water." << G4endl;
    return fWater;
}

// ==================================================================
// Construct Disk source
// ==================================================================
void HPGeDetectorConstruction::ConstructDiskSource(G4LogicalVolume* worldLogical) {
    
    G4cout << "\n=== Constructing Disk Source ===" << G4endl;
    G4cout << "Radius: " << fDiskSourceConfig.radius/mm << " mm" << G4endl;
    G4cout << "Thickness: " << fDiskSourceConfig.thickness/mm << " mm" << G4endl;
    G4cout << "Position Z: " << fDiskSourceConfig.positionZ/mm << " mm" << G4endl;
    G4cout << "Material: " << fDiskSourceConfig.material << G4endl;
    
    // Get fill material
    G4Material* fillMat = GetMaterialByName(fDiskSourceConfig.material);
    
    // Create disk solid
    G4Tubs* solidDiskSource = new G4Tubs("DiskSource",
                                          0,
                                          fDiskSourceConfig.radius,
                                          fDiskSourceConfig.thickness / 2.0,
                                          0, twopi);
    
    fDiskSourceLogical = new G4LogicalVolume(solidDiskSource,
                                              fillMat,
                                              "DiskSource_Logical");
    
    // Visualization - Green semi-transparent
    G4Colour diskColor;
    if (fDiskSourceConfig.material == "water") {
        diskColor = G4Colour(0.2, 0.5, 1.0, 0.6);  // Blue
    } else if (fDiskSourceConfig.material == "soil") {
        diskColor = G4Colour(0.6, 0.4, 0.2, 0.6);  // Brown
    } else {
        diskColor = G4Colour(0.0, 1.0, 0.5, 0.6);  // Green
    }
    
    G4VisAttributes* diskVis = new G4VisAttributes(diskColor);
    diskVis->SetForceSolid(true);
    fDiskSourceLogical->SetVisAttributes(diskVis);
    
    // Position: Z position is center of disk
    G4ThreeVector diskPos(0, 0, fDiskSourceConfig.positionZ);
    
    new G4PVPlacement(0,
                     diskPos,
                     fDiskSourceLogical,
                     "DiskSource_phys",
                     worldLogical,
                     false,
                     0);
    
    // Calculate and display volume
    G4double diskVolume = pi * fDiskSourceConfig.radius * fDiskSourceConfig.radius * fDiskSourceConfig.thickness;
    
    G4cout << "  ✓ Disk source constructed" << G4endl;
    G4cout << "  Volume: " << diskVolume/cm3 << " cm³" << G4endl;
    G4cout << "  Center position: Z = " << fDiskSourceConfig.positionZ/mm << " mm" << G4endl;
    G4cout << "  Extends from Z = " << (fDiskSourceConfig.positionZ - fDiskSourceConfig.thickness/2.0)/mm 
           << " mm to Z = " << (fDiskSourceConfig.positionZ + fDiskSourceConfig.thickness/2.0)/mm << " mm" << G4endl;
}

// ==================================================================
// Construct Cylinder source (with container walls)
// ==================================================================
void HPGeDetectorConstruction::ConstructCylinderSource(G4LogicalVolume* worldLogical) {
    
    G4cout << "\n=== Constructing Cylinder Source ===" << G4endl;
    G4cout << "Inner radius: " << fCylinderSourceConfig.innerRadius/mm << " mm" << G4endl;
    G4cout << "Outer radius: " << fCylinderSourceConfig.outerRadius/mm << " mm" << G4endl;
    G4cout << "Height: " << fCylinderSourceConfig.height/mm << " mm" << G4endl;
    G4cout << "Wall thickness: " << fCylinderSourceConfig.wallThickness/mm << " mm" << G4endl;
    G4cout << "Bottom thickness: " << fCylinderSourceConfig.bottomThickness/mm << " mm" << G4endl;
    G4cout << "Position Z (bottom): " << fCylinderSourceConfig.positionZ/mm << " mm" << G4endl;
    G4cout << "Wall material: " << fCylinderSourceConfig.wallMaterial << G4endl;
    G4cout << "Fill material: " << fCylinderSourceConfig.fillMaterial << G4endl;
    
    // Get materials
    G4Material* wallMat = GetMaterialByName(fCylinderSourceConfig.wallMaterial);
    G4Material* fillMat = GetMaterialByName(fCylinderSourceConfig.fillMaterial);
    
    G4double wall = fCylinderSourceConfig.wallThickness;
    G4double bottom = fCylinderSourceConfig.bottomThickness;
    G4double outerR = fCylinderSourceConfig.outerRadius;
    G4double innerR = outerR - wall;  // Inner radius is outer - wall
    G4double totalH = fCylinderSourceConfig.height;
    G4double sampleH = totalH - bottom;  // Sample height
    
    // Center Z position (from bottom position + half height)
    G4double centerZ = fCylinderSourceConfig.positionZ + totalH / 2.0;
    
    // ========================================================================
    // CONTAINER WALLS (using CSG subtraction)
    // ========================================================================
    
    // Outer cylinder (full container)
    G4Tubs* solidContainerOuter = new G4Tubs("CylinderContainerOuter",
                                              0,
                                              outerR,
                                              totalH / 2.0,
                                              0, twopi);
    
    // Inner cavity (to subtract)
    G4Tubs* solidContainerCavity = new G4Tubs("CylinderContainerCavity",
                                               0,
                                               innerR,
                                               sampleH / 2.0,
                                               0, twopi);
    
    // Subtract cavity from outer (offset upward to leave bottom)
    G4SubtractionSolid* solidContainerWalls = new G4SubtractionSolid("CylinderContainerWalls",
                                                                      solidContainerOuter,
                                                                      solidContainerCavity,
                                                                      0,
                                                                      G4ThreeVector(0, 0, bottom / 2.0));
    
    fCylinderWallLogical = new G4LogicalVolume(solidContainerWalls,
                                                wallMat,
                                                "CylinderWalls_Logical");
    
    // Visualization for container walls - gray/transparent
    G4Colour wallColor;
    if (fCylinderSourceConfig.wallMaterial == "glass") {
        wallColor = G4Colour(0.7, 0.8, 0.9, 0.3);  // Light blue transparent
    } else if (fCylinderSourceConfig.wallMaterial == "polypropylene" || 
               fCylinderSourceConfig.wallMaterial == "hdpe") {
        wallColor = G4Colour(0.8, 0.8, 0.8, 0.4);  // White/gray transparent
    } else {
        wallColor = G4Colour(0.6, 0.6, 0.6, 0.5);  // Gray
    }
    
    G4VisAttributes* wallVis = new G4VisAttributes(wallColor);
    wallVis->SetForceSolid(true);
    fCylinderWallLogical->SetVisAttributes(wallVis);
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, centerZ),
                     fCylinderWallLogical,
                     "CylinderWalls_phys",
                     worldLogical,
                     false,
                     0);
    
    // ========================================================================
    // SAMPLE VOLUME (inside container)
    // ========================================================================
    
    G4Tubs* solidSample = new G4Tubs("CylinderSample",
                                      0,
                                      innerR,
                                      sampleH / 2.0,
                                      0, twopi);
    
    fCylinderSourceLogical = new G4LogicalVolume(solidSample,
                                                  fillMat,
                                                  "CylinderSample_Logical");
    
    // Visualization for sample - colored based on material
    G4Colour sampleColor;
    if (fCylinderSourceConfig.fillMaterial == "water") {
        sampleColor = G4Colour(0.2, 0.5, 1.0, 0.5);  // Blue
    } else if (fCylinderSourceConfig.fillMaterial == "soil") {
        sampleColor = G4Colour(0.6, 0.4, 0.2, 0.6);  // Brown
    } else {
        sampleColor = G4Colour(0.5, 0.8, 0.5, 0.5);  // Light green
    }
    
    G4VisAttributes* sampleVis = new G4VisAttributes(sampleColor);
    sampleVis->SetForceSolid(true);
    fCylinderSourceLogical->SetVisAttributes(sampleVis);
    
    // Sample center position (offset from container center by bottom/2)
    G4double sampleCenterZ = centerZ + bottom / 2.0;
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, sampleCenterZ),
                     fCylinderSourceLogical,
                     "CylinderSample_phys",
                     worldLogical,
                     false,
                     0);
    
    // ========================================================================
    // SUMMARY
    // ========================================================================
    
    G4double wallVolume = pi * (outerR*outerR - innerR*innerR) * totalH + pi * innerR*innerR * bottom;
    G4double sampleVolume = pi * innerR * innerR * sampleH;
    
    G4cout << "\n  ✓ Cylinder source constructed" << G4endl;
    G4cout << "  Container volume: " << wallVolume/cm3 << " cm³" << G4endl;
    G4cout << "  Sample volume: " << sampleVolume/cm3 << " cm³" << G4endl;
    G4cout << "  Container extends from Z = " << fCylinderSourceConfig.positionZ/mm 
           << " mm to Z = " << (fCylinderSourceConfig.positionZ + totalH)/mm << " mm" << G4endl;
    G4cout << "  Sample extends from Z = " << (fCylinderSourceConfig.positionZ + bottom)/mm 
           << " mm to Z = " << (fCylinderSourceConfig.positionZ + totalH)/mm << " mm" << G4endl;
}

// ==================================================================
// Construct Cartridge Source (Orano LEA Type D - Activated Carbon)
// ==================================================================
void HPGeDetectorConstruction::ConstructCartridgeSource(G4LogicalVolume* worldLogical) {
    
    G4cout << "\n=== Constructing Cartridge Source (Orano LEA Type D) ===" << G4endl;
    G4cout << "Housing diameter: " << fCartridgeSourceConfig.housingDiameter/mm << " mm" << G4endl;
    G4cout << "Housing height: " << fCartridgeSourceConfig.housingHeight/mm << " mm" << G4endl;
    G4cout << "Active diameter: " << fCartridgeSourceConfig.activeDiameter/mm << " mm" << G4endl;
    G4cout << "Active thickness: " << fCartridgeSourceConfig.activeThickness/mm << " mm" << G4endl;
    G4cout << "Position Z: " << fCartridgeSourceConfig.positionZ/mm << " mm" << G4endl;
    
    G4double housingR = fCartridgeSourceConfig.housingDiameter / 2.0;
    G4double housingH = fCartridgeSourceConfig.housingHeight;
    G4double activeR = fCartridgeSourceConfig.activeDiameter / 2.0;
    G4double activeH = fCartridgeSourceConfig.activeThickness;
    G4double wallThk = fCartridgeSourceConfig.wallThickness;
    G4double centerZ = fCartridgeSourceConfig.positionZ + housingH / 2.0;
    
    G4Material* housingMat = GetMaterialByName(fCartridgeSourceConfig.housingMaterial);
    G4Material* activeMat = GetMaterialByName(fCartridgeSourceConfig.activeMaterial);
    
    // Housing: outer cylinder with cavity
    G4Tubs* solidHousingOuter = new G4Tubs("CartridgeHousingOuter",
                                            0, housingR,
                                            housingH / 2.0,
                                            0, twopi);
    
    G4Tubs* solidHousingCavity = new G4Tubs("CartridgeHousingCavity",
                                             0, activeR,
                                             activeH / 2.0,
                                             0, twopi);
    
    G4SubtractionSolid* solidHousing = new G4SubtractionSolid("CartridgeHousing",
                                                               solidHousingOuter,
                                                               solidHousingCavity,
                                                               0,
                                                               G4ThreeVector(0, 0, (housingH - activeH) / 2.0 - wallThk));
    
    fCartridgeHousingLogical = new G4LogicalVolume(solidHousing,
                                                    housingMat,
                                                    "CartridgeHousing_Logical");
    
    G4VisAttributes* housingVis = new G4VisAttributes(G4Colour(0.85, 0.85, 0.80, 0.4));
    housingVis->SetForceSolid(true);
    fCartridgeHousingLogical->SetVisAttributes(housingVis);
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, centerZ),
                     fCartridgeHousingLogical,
                     "CartridgeHousing_phys",
                     worldLogical,
                     false,
                     0);
    
    // Active matrix: activated carbon disk
    G4Tubs* solidActive = new G4Tubs("CartridgeActive",
                                      0, activeR,
                                      activeH / 2.0,
                                      0, twopi);
    
    fCartridgeSourceLogical = new G4LogicalVolume(solidActive,
                                                   activeMat,
                                                   "CartridgeActive_Logical");
    
    G4VisAttributes* activeVis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 0.7));
    activeVis->SetForceSolid(true);
    fCartridgeSourceLogical->SetVisAttributes(activeVis);
    
    G4double activeCenterZ = centerZ + (housingH - activeH) / 2.0 - wallThk;
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, activeCenterZ),
                     fCartridgeSourceLogical,
                     "CartridgeActive_phys",
                     worldLogical,
                     false,
                     0);
    
    G4double activeVolume = pi * activeR * activeR * activeH;
    G4double housingVolume = pi * housingR * housingR * housingH;
    
    G4cout << "  Active matrix volume: " << activeVolume/cm3 << " cm3" << G4endl;
    G4cout << "  Housing volume: " << housingVolume/cm3 << " cm3" << G4endl;
    G4cout << "  Center position: Z = " << centerZ/mm << " mm" << G4endl;
}

// ==================================================================
// Construct Filter Source (Orano LEA Type M - Paper Filter)
// ==================================================================
void HPGeDetectorConstruction::ConstructFilterSource(G4LogicalVolume* worldLogical) {
    
    G4cout << "\n=== Constructing Filter Source (Orano LEA Type M) ===" << G4endl;
    G4cout << "Outer diameter: " << fFilterSourceConfig.outerDiameter/mm << " mm" << G4endl;
    G4cout << "Active diameter: " << fFilterSourceConfig.activeDiameter/mm << " mm" << G4endl;
    G4cout << "Filter thickness: " << fFilterSourceConfig.filterThickness/mm << " mm" << G4endl;
    G4cout << "Seal thickness: " << fFilterSourceConfig.sealThickness/mm << " mm" << G4endl;
    G4cout << "Position Z: " << fFilterSourceConfig.positionZ/mm << " mm" << G4endl;
    
    G4double outerR = fFilterSourceConfig.outerDiameter / 2.0;
    G4double activeR = fFilterSourceConfig.activeDiameter / 2.0;
    G4double filterH = fFilterSourceConfig.filterThickness;
    G4double sealH = fFilterSourceConfig.sealThickness;
    G4double totalH = filterH + 2.0 * sealH;  // Sealed between two polyester sheets
    G4double centerZ = fFilterSourceConfig.positionZ;
    
    G4Material* filterMat = GetMaterialByName(fFilterSourceConfig.filterMaterial);
    G4Material* sealMat = GetMaterialByName(fFilterSourceConfig.sealMaterial);
    
    // Polyester seal (outer envelope covering full diameter)
    G4Tubs* solidSeal = new G4Tubs("FilterSeal",
                                    0, outerR,
                                    totalH / 2.0,
                                    0, twopi);
    
    // Cavity for filter paper
    G4Tubs* solidSealCavity = new G4Tubs("FilterSealCavity",
                                          0, activeR,
                                          filterH / 2.0,
                                          0, twopi);
    
    G4SubtractionSolid* solidSealFinal = new G4SubtractionSolid("FilterSealFinal",
                                                                  solidSeal,
                                                                  solidSealCavity,
                                                                  0,
                                                                  G4ThreeVector(0, 0, 0));
    
    fFilterSealLogical = new G4LogicalVolume(solidSealFinal,
                                              sealMat,
                                              "FilterSeal_Logical");
    
    G4VisAttributes* sealVis = new G4VisAttributes(G4Colour(0.9, 0.9, 0.85, 0.3));
    sealVis->SetForceSolid(true);
    fFilterSealLogical->SetVisAttributes(sealVis);
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, centerZ),
                     fFilterSealLogical,
                     "FilterSeal_phys",
                     worldLogical,
                     false,
                     0);
    
    // Active filter paper (radionuclides deposited on paper)
    G4Tubs* solidFilter = new G4Tubs("FilterActive",
                                      0, activeR,
                                      filterH / 2.0,
                                      0, twopi);
    
    fFilterSourceLogical = new G4LogicalVolume(solidFilter,
                                                filterMat,
                                                "FilterActive_Logical");
    
    G4VisAttributes* filterVis = new G4VisAttributes(G4Colour(1.0, 0.9, 0.0, 0.7));
    filterVis->SetForceSolid(true);
    fFilterSourceLogical->SetVisAttributes(filterVis);
    
    new G4PVPlacement(0,
                     G4ThreeVector(0, 0, centerZ),
                     fFilterSourceLogical,
                     "FilterActive_phys",
                     worldLogical,
                     false,
                     0);
    
    G4double activeArea = pi * activeR * activeR;
    G4double totalArea = pi * outerR * outerR;
    
    G4cout << "  Active area: " << activeArea/cm2 << " cm2" << G4endl;
    G4cout << "  Total area: " << totalArea/cm2 << " cm2" << G4endl;
    G4cout << "  Total thickness (with seal): " << totalH/mm << " mm" << G4endl;
    G4cout << "  Center position: Z = " << centerZ/mm << " mm" << G4endl;
}

// ==================================================================
// Construct lead castle
// ==================================================================
void HPGeDetectorConstruction::ConstructLeadCastle(G4LogicalVolume* motherVolume) {
    
    G4double lead_x = 14.0 * cm;
    G4double lead_inner_x = 10.5 * cm;
    
    G4Box* solidLeadOuter = new G4Box("LeadOuter", 
                                      lead_x, lead_x, (53.3 + 29.0)*cm/2.0);
    
    G4Box* solidLeadCavity = new G4Box("LeadCavity",
                                       lead_inner_x, lead_inner_x, (37.5 + 13.0)*cm/2.0);
    
    G4SubtractionSolid* solidLeadCastle = new G4SubtractionSolid("LeadCastle",
                                                                 solidLeadOuter,
                                                                 solidLeadCavity);
    
    G4LogicalVolume* logicLeadCastle = new G4LogicalVolume(solidLeadCastle, fLead, "LeadCastle");
    
    G4VisAttributes* leadVis = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3, 0.5));
    leadVis->SetForceSolid(true);
    logicLeadCastle->SetVisAttributes(leadVis);
    
    G4double zPos = (-53.3 + 29.0) * cm / 2.0;
    new G4PVPlacement(0, G4ThreeVector(0, 0, zPos), logicLeadCastle,
                     "LeadCastle_phys", motherVolume, false, 0);
    
    G4cout << "  Lead castle constructed" << G4endl;
}

// ==================================================================
// Construct HPGe detector
// ==================================================================
G4LogicalVolume* HPGeDetectorConstruction::ConstructHPGeDetectorFLUKA(const HPGeConfig& config) {
    
    G4double crystalR = config.crystalDiameter / 2.0;
    G4double crystalH = config.crystalLength / 2.0;
    G4double holeR = config.holeDiameter / 2.0;
    G4double holeD = config.holeDepth;
    
    G4double geDeadLat = config.geDeadLayer;
    G4double geDeadFront = config.geDeadLayerFront;
    G4double liDeadLat = config.liDeadLayer;
    G4double liDeadFront = config.liDeadLayerFront;
    G4double boronThk = config.boronThickness;
    
    G4double alHolderThickness = 1.0 * mm;
    
    G4double containerR = crystalR + boronThk + alHolderThickness + config.vacuumGap + config.alCupThickness + 5.*mm;
    G4double containerH = crystalH + config.vacuumGap + config.alWindowThickness + 10.*mm;
    
    G4Tubs* solidContainer = new G4Tubs("Container", 0., containerR, containerH, 0., twopi);
    G4LogicalVolume* logicContainer = new G4LogicalVolume(solidContainer, fAir, "Container");
    logicContainer->SetVisAttributes(G4VisAttributes::GetInvisible());
 
 
     // 1. Define Z-boundaries based on the hole starting from the back (-crystalH)
G4double holeStartZ = -crystalH;
G4double holeEndZ   = -crystalH + holeD;   
    // Active Ge volumes
    G4double activeOuterR = crystalR - geDeadLat;
    G4double activeInnerR = holeR + liDeadLat;
    
    G4double activeFrontZ = -crystalH + geDeadFront;
    //G4double activeBackZ = crystalH - geDeadFront;
    
    G4double activeBackZ = crystalH - geDeadFront;
    
    // Upper part (the solid tip) starts where the hole ends
G4double activeUpperFrontZ = holeEndZ;
G4double activeUpperBackZ  = crystalH - geDeadFront;
    
G4double activeUpperHalfHeight = (activeUpperBackZ - activeUpperFrontZ) / 2.0;
G4double activeUpperCenterZ    = (activeUpperFrontZ + activeUpperBackZ) / 2.0;
    
    G4Tubs* solidActiveUpper = new G4Tubs("ActiveGe_Upper",
                                          0,
                                          activeOuterR,
                                          activeUpperHalfHeight,
                                          0., twopi);
    

    
    G4double activeLowerFrontZ = activeFrontZ;
    //G4double activeLowerBackZ = activeBackZ - holeD;
    G4double activeLowerBackZ = holeEndZ;
    
    G4double activeLowerHalfHeight = (activeLowerBackZ - activeLowerFrontZ) / 2.0;
    G4double activeLowerCenterZ = (activeLowerFrontZ + activeLowerBackZ) / 2.0;
    
    G4Tubs* solidActiveLower = new G4Tubs("ActiveGe_Lower",
                                          activeInnerR,
                                          activeOuterR,
                                          activeLowerHalfHeight,
                                          0., twopi);
    
    G4UnionSolid* solidActiveGe = new G4UnionSolid("ActiveGe_Union",
                                                   solidActiveLower,
                                                   solidActiveUpper,
                                                   0,
                                                   G4ThreeVector(0, 0, activeUpperCenterZ - activeLowerCenterZ));
    
    G4LogicalVolume* logicActiveGe = new G4LogicalVolume(solidActiveGe, fGe, "ActiveGe");
    
    G4VisAttributes* activeGeVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.8));
    activeGeVis->SetForceSolid(true);
    logicActiveGe->SetVisAttributes(activeGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, activeLowerCenterZ), logicActiveGe,
                     "ActiveGe_phys", logicContainer, false, 0);
    
    fActiveGeLogicals.push_back(logicActiveGe);
    
    // Boron outer
    G4double boronOuterR = crystalR + boronThk;
    
    G4Tubs* solidBoronOuter = new G4Tubs("BoronOuter",
                                         crystalR,
                                         boronOuterR,
                                         crystalH,
                                         0., twopi);
    
    G4LogicalVolume* logicBoronOuter = new G4LogicalVolume(solidBoronOuter, fBoron, "BoronOuter");
    
    G4VisAttributes* boronVis = new G4VisAttributes(G4Colour(0.6, 0.4, 0.2, 0.4));
    boronVis->SetForceSolid(true);
    logicBoronOuter->SetVisAttributes(boronVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicBoronOuter,
                     "BoronOuter_phys", logicContainer, false, 0);
    
    // Aluminum holder
    G4double holderInnerR = boronOuterR;
    G4double holderOuterR = holderInnerR + alHolderThickness;
    
    G4Tubs* solidAlHolder = new G4Tubs("AlHolder",
                                       holderInnerR,
                                       holderOuterR,
                                       crystalH,
                                       0., twopi);
    
    G4LogicalVolume* logicAlHolder = new G4LogicalVolume(solidAlHolder, fAl, "AlHolder");
    
    G4VisAttributes* holderVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.6));
    holderVis->SetForceSolid(true);
    logicAlHolder->SetVisAttributes(holderVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicAlHolder,
                     "AlHolder_phys", logicContainer, false, 0);
    
    // Dead layers
    G4Tubs* solidDeadGeLat = new G4Tubs("DeadGeLat",
                                        activeOuterR,
                                        crystalR,
                                        crystalH - geDeadFront,
                                        0., twopi);
    
    G4LogicalVolume* logicDeadGeLat = new G4LogicalVolume(solidDeadGeLat, fGe, "DeadGeLat");
    
    G4VisAttributes* deadGeVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.3));
    deadGeVis->SetForceSolid(true);
    logicDeadGeLat->SetVisAttributes(deadGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, geDeadFront/2.0), logicDeadGeLat,
                     "DeadGeLat_phys", logicContainer, false, 0);
    
    G4double deadFrontZ = crystalH - geDeadFront/2.0;
    
    G4Tubs* solidDeadGeFront = new G4Tubs("DeadGeFront",
                                          0,
                                          activeOuterR,
                                          geDeadFront/2.0,
                                          0., twopi);
    
    G4LogicalVolume* logicDeadGeFront = new G4LogicalVolume(solidDeadGeFront, fGe, "DeadGeFront");
    logicDeadGeFront->SetVisAttributes(deadGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, deadFrontZ), logicDeadGeFront,
                     "DeadGeFront_phys", logicContainer, false, 0);
    
    G4double deadBottomZ = -crystalH + geDeadFront/2.0;
    
    G4Tubs* solidDeadGeBottom = new G4Tubs("DeadGeBottom",
                                           0,
                                           crystalR,
                                           geDeadFront/2.0,
                                           0., twopi);
    
    G4LogicalVolume* logicDeadGeBottom = new G4LogicalVolume(solidDeadGeBottom, fGe, "DeadGeBottom");
    logicDeadGeBottom->SetVisAttributes(deadGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, deadBottomZ), logicDeadGeBottom,
                     "DeadGeBottom_phys", logicContainer, false, 0);
    
    // Li dead layers
    //G4double liZ = -crystalH + holeD/2.0;
    G4double liZ = holeStartZ + holeD/2.0;
    
    G4Tubs* solidLiDeadLat = new G4Tubs("LiDeadLat",
                                        holeR,
                                        holeR + liDeadLat,
                                        holeD/2.0,
                                        0., twopi);
    
    G4LogicalVolume* logicLiDeadLat = new G4LogicalVolume(solidLiDeadLat, fLi, "LiDeadLat");
    
    G4VisAttributes* liVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.0, 0.4));
    liVis->SetForceSolid(true);
    logicLiDeadLat->SetVisAttributes(liVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, liZ), logicLiDeadLat,
                     "LiDeadLat_phys", logicContainer, false, 0);
    
    //G4double liFrontZ = -crystalH + liDeadFront/2.0;
    G4double liFrontZ = holeEndZ - liDeadFront/2.0;
    
    G4Tubs* solidLiDeadFront = new G4Tubs("LiDeadFront",
                                          0,
                                          holeR,
                                          liDeadFront/2.0,
                                          0., twopi);
    
    G4LogicalVolume* logicLiDeadFront = new G4LogicalVolume(solidLiDeadFront, fLi, "LiDeadFront");
    logicLiDeadFront->SetVisAttributes(liVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, liFrontZ), logicLiDeadFront,
                     "LiDeadFront_phys", logicContainer, false, 0);
    
    // Vacuum gap
    G4double vacR = holderOuterR + config.vacuumGap;
    //G4double vacZ = crystalH + config.vacuumGap/2.0;
    
    G4Tubs* solidVacuum = new G4Tubs("Vacuum",
                                     holderOuterR,
                                     vacR,
                                     config.vacuumGap/2.0 + crystalH,
                                     0., twopi);
    
    G4LogicalVolume* logicVacuum = new G4LogicalVolume(solidVacuum, fVacuum, "Vacuum");
    logicVacuum->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicVacuum,
                     "Vacuum_phys", logicContainer, false, 0);
    
    // Al cup
    G4double cupInnerR = vacR;
    G4double cupOuterR = cupInnerR + config.alCupThickness;
    
    G4Tubs* solidAlCup = new G4Tubs("AlCup",
                                    cupInnerR,
                                    cupOuterR,
                                    crystalH + config.vacuumGap/2.0,
                                    0., twopi);
    
    G4LogicalVolume* logicAlCup = new G4LogicalVolume(solidAlCup, fAl, "AlCup");
    
    G4VisAttributes* alVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.8));
    alVis->SetForceSolid(true);
    logicAlCup->SetVisAttributes(alVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, config.vacuumGap/2.0), logicAlCup,
                     "AlCup_phys", logicContainer, false, 0);
    
    // Al window
    G4double windowZ = crystalH + config.vacuumGap + config.alWindowThickness/2.0;
    
    G4Tubs* solidAlWindow = new G4Tubs("AlWindow",
                                       0,
                                       cupInnerR,
                                       config.alWindowThickness/2.0,
                                       0., twopi);
    
    G4LogicalVolume* logicAlWindow = new G4LogicalVolume(solidAlWindow, fC, "AlWindow");
    logicAlWindow->SetVisAttributes(alVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, windowZ), logicAlWindow,
                     "AlWindow_phys", logicContainer, false, 0);
    
    return logicContainer;
}

// ==================================================================
// Construct HPGe detector with Canberra Mirion real geometry
// Based on S/N: 21206 characterization sheet
// ==================================================================
G4LogicalVolume* HPGeDetectorConstruction::ConstructHPGeDetectorCanberra(const HPGeConfig& config) {
    
    G4cout << "\n=== Constructing Canberra Mirion Real Geometry ===" << G4endl;
    
    // ================================================================
    // CRYSTAL DIMENSIONS
    // ================================================================
    G4double crystalR = config.crystalDiameter / 2.0;
    G4double crystalH = config.crystalLength / 2.0;  // half-height
    G4double holeR = config.holeDiameter / 2.0;
    G4double holeD = config.holeDepth;
    
    // Dead layers
    G4double geDeadLat = config.geDeadLayer;           // 0.5 mm outer electrode
    G4double geDeadFront = config.geDeadLayerFront;    // 0.4 um front dead layer
    G4double liDeadLat = config.liDeadLayer;            // 0.3 um inner electrode
    G4double liDeadFront = config.liDeadLayerFront;
    G4double boronThk = config.boronThickness;
    
    // ================================================================
    // ENDCAP DIMENSIONS (from assembly drawing)
    // ================================================================
    G4double endcapR = config.endcapDiameter / 2.0;          // 76.2/2 = 38.1 mm
    G4double endcapWallThickness = config.alCupThickness;     // 1.5 mm side wall
    G4double windowThickness = config.alWindowThickness;      // 0.6 mm window
    G4double geToEndcap = config.geToEndcapDistance;           // 6.0 mm
    
    // ================================================================
    // HOLDER DIMENSIONS (from characterization sheet)
    // ================================================================
    G4double holderR1 = config.holderOuterDiameter1 / 2.0;   // 69.0/2 = 34.5 mm
    G4double holderR3 = config.holderOuterDiameter3 / 2.0;   // 63.0/2 = 31.5 mm
    G4double holderTopH = config.holderTopHeight;              // 8.6 mm
    G4double holderTotalL = config.holderTotalLength;          // 95.2 mm
    
    // HDPE spacer
    G4double hdpeThk = config.hdpeThickness;                   // 6.5 mm
    
    // ================================================================
    // CONTAINER: encompasses all detector components
    // ================================================================
    G4double containerR = endcapR + 5.*mm;
    G4double containerH = crystalH + geToEndcap + windowThickness + hdpeThk + 20.*mm;
    
    G4Tubs* solidContainer = new G4Tubs("Container", 0., containerR, containerH, 0., twopi);
    G4LogicalVolume* logicContainer = new G4LogicalVolume(solidContainer, fAir, "Container");
    logicContainer->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    // ================================================================
    // 1. ACTIVE GERMANIUM CRYSTAL
    //    Hole starts from back face (-crystalH), extends holeD upward
    // ================================================================
    G4double holeStartZ = -crystalH;
    G4double holeEndZ = -crystalH + holeD;
    
    G4double activeOuterR = crystalR - geDeadLat;
    G4double activeInnerR = holeR + liDeadLat;
    
    G4double activeFrontZ = -crystalH + geDeadFront;
    
    // Upper part: solid tip above the hole
    G4double activeUpperFrontZ = holeEndZ;
    G4double activeUpperBackZ = crystalH - geDeadFront;
    G4double activeUpperHalfHeight = (activeUpperBackZ - activeUpperFrontZ) / 2.0;
    G4double activeUpperCenterZ = (activeUpperFrontZ + activeUpperBackZ) / 2.0;
    
    G4Tubs* solidActiveUpper = new G4Tubs("ActiveGe_Upper",
                                          0, activeOuterR,
                                          activeUpperHalfHeight,
                                          0., twopi);
    
    // Lower part: annular region around the hole
    G4double activeLowerFrontZ = activeFrontZ;
    G4double activeLowerBackZ = holeEndZ;
    G4double activeLowerHalfHeight = (activeLowerBackZ - activeLowerFrontZ) / 2.0;
    G4double activeLowerCenterZ = (activeLowerFrontZ + activeLowerBackZ) / 2.0;
    
    G4Tubs* solidActiveLower = new G4Tubs("ActiveGe_Lower",
                                          activeInnerR, activeOuterR,
                                          activeLowerHalfHeight,
                                          0., twopi);
    
    // Union of upper and lower active regions
    G4UnionSolid* solidActiveGe = new G4UnionSolid("ActiveGe_Union",
                                                   solidActiveLower,
                                                   solidActiveUpper,
                                                   0,
                                                   G4ThreeVector(0, 0, activeUpperCenterZ - activeLowerCenterZ));
    
    G4LogicalVolume* logicActiveGe = new G4LogicalVolume(solidActiveGe, fGe, "ActiveGe");
    
    G4VisAttributes* activeGeVis = new G4VisAttributes(G4Colour(0.0, 0.5, 1.0, 0.8));
    activeGeVis->SetForceSolid(true);
    logicActiveGe->SetVisAttributes(activeGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, activeLowerCenterZ), logicActiveGe,
                     "ActiveGe_phys", logicContainer, false, 0);
    
    fActiveGeLogicals.push_back(logicActiveGe);
    
    // ================================================================
    // 2. OUTER DEAD LAYER (Ge - 0.5 mm outer electrode)
    // ================================================================
    // Lateral dead layer
    G4Tubs* solidDeadGeLat = new G4Tubs("DeadGeLat",
                                        activeOuterR, crystalR,
                                        crystalH - geDeadFront,
                                        0., twopi);
    
    G4LogicalVolume* logicDeadGeLat = new G4LogicalVolume(solidDeadGeLat, fGe, "DeadGeLat");
    G4VisAttributes* deadGeVis = new G4VisAttributes(G4Colour(0.5, 0.5, 0.5, 0.3));
    deadGeVis->SetForceSolid(true);
    logicDeadGeLat->SetVisAttributes(deadGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, geDeadFront/2.0), logicDeadGeLat,
                     "DeadGeLat_phys", logicContainer, false, 0);
    
    // Front dead layer (0.4 um)
    G4double deadFrontZ = crystalH - geDeadFront/2.0;
    G4Tubs* solidDeadGeFront = new G4Tubs("DeadGeFront",
                                          0, activeOuterR,
                                          geDeadFront/2.0,
                                          0., twopi);
    
    G4LogicalVolume* logicDeadGeFront = new G4LogicalVolume(solidDeadGeFront, fGe, "DeadGeFront");
    logicDeadGeFront->SetVisAttributes(deadGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, deadFrontZ), logicDeadGeFront,
                     "DeadGeFront_phys", logicContainer, false, 0);
    
    // Bottom dead layer
    G4double deadBottomZ = -crystalH + geDeadFront/2.0;
    G4Tubs* solidDeadGeBottom = new G4Tubs("DeadGeBottom",
                                           0, crystalR,
                                           geDeadFront/2.0,
                                           0., twopi);
    
    G4LogicalVolume* logicDeadGeBottom = new G4LogicalVolume(solidDeadGeBottom, fGe, "DeadGeBottom");
    logicDeadGeBottom->SetVisAttributes(deadGeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, deadBottomZ), logicDeadGeBottom,
                     "DeadGeBottom_phys", logicContainer, false, 0);
    
    // ================================================================
    // 3. INNER DEAD LAYER (Li - 0.3 um inner electrode)
    // ================================================================
    G4double liZ = holeStartZ + holeD/2.0;
    
    G4Tubs* solidLiDeadLat = new G4Tubs("LiDeadLat",
                                        holeR, holeR + liDeadLat,
                                        holeD/2.0,
                                        0., twopi);
    
    G4LogicalVolume* logicLiDeadLat = new G4LogicalVolume(solidLiDeadLat, fLi, "LiDeadLat");
    G4VisAttributes* liVis = new G4VisAttributes(G4Colour(0.8, 0.8, 0.0, 0.4));
    liVis->SetForceSolid(true);
    logicLiDeadLat->SetVisAttributes(liVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, liZ), logicLiDeadLat,
                     "LiDeadLat_phys", logicContainer, false, 0);
    
    // Inner electrode front cap
    G4double liFrontZ = holeEndZ - liDeadFront/2.0;
    G4Tubs* solidLiDeadFront = new G4Tubs("LiDeadFront",
                                          0, holeR,
                                          liDeadFront/2.0,
                                          0., twopi);
    
    G4LogicalVolume* logicLiDeadFront = new G4LogicalVolume(solidLiDeadFront, fLi, "LiDeadFront");
    logicLiDeadFront->SetVisAttributes(liVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, liFrontZ), logicLiDeadFront,
                     "LiDeadFront_phys", logicContainer, false, 0);
    
    // ================================================================
    // 4. BORON OUTER COATING
    // ================================================================
    G4double boronOuterR = crystalR + boronThk;
    
    G4Tubs* solidBoronOuter = new G4Tubs("BoronOuter",
                                         crystalR, boronOuterR,
                                         crystalH,
                                         0., twopi);
    
    G4LogicalVolume* logicBoronOuter = new G4LogicalVolume(solidBoronOuter, fBoron, "BoronOuter");
    G4VisAttributes* boronVis = new G4VisAttributes(G4Colour(0.6, 0.4, 0.2, 0.4));
    boronVis->SetForceSolid(true);
    logicBoronOuter->SetVisAttributes(boronVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicBoronOuter,
                     "BoronOuter_phys", logicContainer, false, 0);
    
    // ================================================================
    // 5. ALUMINUM HOLDER (from characterization sheet)
    //    Stepped cylinder: top ring Ø69.0, body Ø63.0
    // ================================================================
    G4double holderWallThickness = 1.0 * mm;
    G4double holderInnerR = boronOuterR;
    G4double holderOuterR = holderInnerR + holderWallThickness;
    
    G4Tubs* solidAlHolder = new G4Tubs("AlHolder",
                                       holderInnerR, holderOuterR,
                                       crystalH,
                                       0., twopi);
    
    G4LogicalVolume* logicAlHolder = new G4LogicalVolume(solidAlHolder, fAl, "AlHolder");
    G4VisAttributes* holderVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.6));
    holderVis->SetForceSolid(true);
    logicAlHolder->SetVisAttributes(holderVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicAlHolder,
                     "AlHolder_phys", logicContainer, false, 0);
    
    // ================================================================
    // 6. VACUUM GAP
    // ================================================================
    G4double vacInnerR = holderOuterR;
    G4double vacOuterR = endcapR - endcapWallThickness;
    
    G4Tubs* solidVacuum = new G4Tubs("Vacuum",
                                     vacInnerR, vacOuterR,
                                     crystalH + geToEndcap/2.0,
                                     0., twopi);
    
    G4LogicalVolume* logicVacuum = new G4LogicalVolume(solidVacuum, fVacuum, "Vacuum");
    logicVacuum->SetVisAttributes(G4VisAttributes::GetInvisible());
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, 0), logicVacuum,
                     "Vacuum_phys", logicContainer, false, 0);
    
    // ================================================================
    // 7. ALUMINUM ENDCAP (side walls)
    //    Ø76.2 mm outer diameter from assembly drawing
    // ================================================================
    G4double cupInnerR = endcapR - endcapWallThickness;
    G4double cupOuterR = endcapR;
    
    G4Tubs* solidAlEndcap = new G4Tubs("AlEndcap",
                                       cupInnerR, cupOuterR,
                                       crystalH + geToEndcap/2.0,
                                       0., twopi);
    
    G4LogicalVolume* logicAlEndcap = new G4LogicalVolume(solidAlEndcap, fAl, "AlEndcap");
    G4VisAttributes* alVis = new G4VisAttributes(G4Colour(0.7, 0.7, 0.7, 0.8));
    alVis->SetForceSolid(true);
    logicAlEndcap->SetVisAttributes(alVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, geToEndcap/2.0), logicAlEndcap,
                     "AlEndcap_phys", logicContainer, false, 0);
    
    // ================================================================
    // 8. CARBON EPOXY WINDOW (or Al window depending on config)
    //    Positioned at the front of the endcap
    // ================================================================
    G4double windowZ = crystalH + geToEndcap + windowThickness/2.0;
    
    G4Tubs* solidWindow = new G4Tubs("DetectorWindow",
                                     0, cupInnerR,
                                     windowThickness/2.0,
                                     0., twopi);
    
    // Select window material based on config
    G4Material* windowMat = fCarbonEpoxy;  // Default: Carbon Epoxy
    if (config.windowMaterial == "Al" || config.windowMaterial == "aluminum") {
        windowMat = fAl;
    } else if (config.windowMaterial == "carbon" || config.windowMaterial == "C") {
        windowMat = fC;
    }
    
    G4LogicalVolume* logicWindow = new G4LogicalVolume(solidWindow, windowMat, "DetectorWindow");
    
    G4VisAttributes* windowVis = new G4VisAttributes(G4Colour(0.2, 0.2, 0.2, 0.7));
    windowVis->SetForceSolid(true);
    logicWindow->SetVisAttributes(windowVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, windowZ), logicWindow,
                     "DetectorWindow_phys", logicContainer, false, 0);
    
    // ================================================================
    // 9. HDPE SPACER (between holder and cold finger area)
    // ================================================================
    G4double hdpeInnerR = 0.;
    G4double hdpeOuterR = holderInnerR;
    G4double hdpeZ = -crystalH - hdpeThk/2.0;
    
    G4Tubs* solidHDPE = new G4Tubs("HDPE_Spacer",
                                    hdpeInnerR, hdpeOuterR,
                                    hdpeThk/2.0,
                                    0., twopi);
    
    G4LogicalVolume* logicHDPE = new G4LogicalVolume(solidHDPE, fHDPE, "HDPE_Spacer");
    
    G4VisAttributes* hdpeVis = new G4VisAttributes(G4Colour(0.9, 0.9, 0.8, 0.5));
    hdpeVis->SetForceSolid(true);
    logicHDPE->SetVisAttributes(hdpeVis);
    
    new G4PVPlacement(0, G4ThreeVector(0, 0, hdpeZ), logicHDPE,
                     "HDPE_Spacer_phys", logicContainer, false, 0);
    
    // ================================================================
    // SUMMARY
    // ================================================================
    G4cout << "  Crystal active volume computed" << G4endl;
    G4cout << "  Outer electrode (Ge dead layer): " << geDeadLat/mm << " mm" << G4endl;
    G4cout << "  Front dead layer: " << geDeadFront/mm * 1000.0 << " um" << G4endl;
    G4cout << "  Inner electrode: " << liDeadLat/mm * 1000.0 << " um" << G4endl;
    G4cout << "  Window: " << config.windowMaterial << " (" 
           << windowThickness/mm << " mm)" << G4endl;
    G4cout << "  Endcap: " << config.endcapMaterial << " (Ø" 
           << config.endcapDiameter/mm << " mm)" << G4endl;
    G4cout << "  Holder: " << config.holderMaterial << " (Ø" 
           << config.holderOuterDiameter1/mm << " mm)" << G4endl;
    G4cout << "  HDPE spacer: " << hdpeThk/mm << " mm" << G4endl;
    
    return logicContainer;
}

// ==================================================================
// Construct SD and field
// ==================================================================
void HPGeDetectorConstruction::ConstructSDandField() {
    G4SDManager* sdManager = G4SDManager::GetSDMpointer();
    
    G4cout << "\n=== Creating Sensitive Detectors ===" << G4endl;
    
    for (size_t i = 0; i < fDetectorConfigs.size(); i++) {
        G4String sdName = "HPGe_SD_" + std::to_string(i);
        HPGeSensitiveDetector* sd = new HPGeSensitiveDetector(sdName, 
                                                              fDetectorConfigs[i].detectorID);
        sdManager->AddNewDetector(sd);
        SetSensitiveDetector("ActiveGe", sd);
        
        G4cout << "  Created: " << sdName << " (ID: " 
               << fDetectorConfigs[i].detectorID << ")" << G4endl;
    }
}

// ==================================================================
// Calculate active volume
// ==================================================================
G4double HPGeDetectorConstruction::CalculateActiveVolume(const HPGeConfig& config) const {
    G4double crystalR = config.crystalDiameter / 2.0;
    G4double crystalH = config.crystalLength / 2.0;
    G4double holeR = config.holeDiameter / 2.0;
    G4double holeD = config.holeDepth;
    
    G4double geDeadLat = config.geDeadLayer;
    G4double geDeadFront = config.geDeadLayerFront;
    G4double liDeadLat = config.liDeadLayer;
    G4double liDeadFront = config.liDeadLayerFront;
    
    G4double activeOuterR = crystalR - geDeadLat;
    G4double activeInnerR = holeR + liDeadLat;
    
    G4double activeFrontZ = -crystalH + geDeadFront;
    G4double activeBackZ = crystalH - geDeadFront;
    
    G4double upperHeight = holeD;
    G4double upperVolume = pi * activeOuterR * activeOuterR * upperHeight;
    
    G4double lowerHeight = activeBackZ - activeFrontZ - holeD;
    G4double lowerVolume = pi * (activeOuterR*activeOuterR - activeInnerR*activeInnerR) * lowerHeight;
    
    G4double totalVolume = upperVolume + lowerVolume;
    
    return totalVolume;
}

// ==================================================================
// Get active volume
// ==================================================================
G4double HPGeDetectorConstruction::GetActiveVolume(G4int detectorID) const {
    if (detectorID >= 0 && detectorID < (G4int)fDetectorConfigs.size()) {
        return CalculateActiveVolume(fDetectorConfigs[detectorID]);
    }
    return 0.0;
}

// ==================================================================
// Get active mass
// ==================================================================
G4double HPGeDetectorConstruction::GetActiveMass(G4int detectorID) const {
    G4double volume = GetActiveVolume(detectorID);
    G4double density = fGe->GetDensity();
    return volume * density;
}

// ==================================================================
// Get crystal to source distance
// ==================================================================
G4double HPGeDetectorConstruction::GetCrystalToSourceDistance(G4int detectorID) const {
    if (detectorID >= 0 && detectorID < (G4int)fDetectorConfigs.size()) {
        const HPGeConfig& config = fDetectorConfigs[detectorID];
        G4ThreeVector detPos(config.posX, config.posY, config.posZ);
        
        G4double distance = (fSourcePosition - detPos).mag();
        return distance;
    }
    return 0.0;
}

// ==================================================================
// Add detector configuration
// ==================================================================
void HPGeDetectorConstruction::AddDetector(const HPGeConfig& config) {
    fDetectorConfigs.push_back(config);
}

// ==================================================================
// Set source position
// ==================================================================
void HPGeDetectorConstruction::SetSourcePosition(G4double x, G4double y, G4double z) {
    fSourcePosition = G4ThreeVector(x, y, z);
}

// ==================================================================
// Construct source boundary visualization
// ==================================================================
void HPGeDetectorConstruction::ConstructSourceBoundary(G4LogicalVolume* worldLogical) {
    
    G4cout << "\n=== Constructing Source Boundary Visualization ===" << G4endl;
    G4cout << "Type: " << fSourceBoundaryConfig.type << G4endl;
    
    G4VSolid* solidBoundary = nullptr;
    G4ThreeVector position = fSourceBoundaryConfig.position;
    
    if (fSourceBoundaryConfig.type == "point") {
        // Small sphere to mark point source
        G4double markerRadius = 2.*mm;
        solidBoundary = new G4Sphere("SourceBoundary_Point", 
                                     0, markerRadius, 
                                     0, twopi, 
                                     0, pi);
        G4cout << "  Point source marker at: (" 
               << position.x()/cm << ", "
               << position.y()/cm << ", "
               << position.z()/cm << ") cm" << G4endl;
    }
    else if (fSourceBoundaryConfig.type == "disk") {
        // Thin disk
        G4double radius = fSourceBoundaryConfig.radius;
        G4double thickness = 0.5*mm;  // Very thin for visualization
        solidBoundary = new G4Tubs("SourceBoundary_Disk",
                                   0, radius,
                                   thickness/2.0,
                                   0, twopi);
        G4cout << "  Disk source boundary:" << G4endl;
        G4cout << "    Radius: " << radius/cm << " cm" << G4endl;
        G4cout << "    Position: (" << position.x()/cm << ", "
               << position.y()/cm << ", "
               << position.z()/cm << ") cm" << G4endl;
    }
    else if (fSourceBoundaryConfig.type == "volume") {
        // Cylindrical volume wireframe
        G4double radius = fSourceBoundaryConfig.radius;
        G4double height = fSourceBoundaryConfig.height;
        solidBoundary = new G4Tubs("SourceBoundary_Volume",
                                   radius - 0.5*mm,  // Thin shell
                                   radius,
                                   height/2.0,
                                   0, twopi);
        G4cout << "  Volume source boundary:" << G4endl;
        G4cout << "    Radius: " << radius/cm << " cm" << G4endl;
        G4cout << "    Height: " << height/cm << " cm" << G4endl;
        G4cout << "    Position: (" << position.x()/cm << ", "
               << position.y()/cm << ", "
               << position.z()/cm << ") cm" << G4endl;
    }
    else if (fSourceBoundaryConfig.type == "marinelli") {
        // Marinelli boundary is already visible as the beaker itself
        G4cout << "  Marinelli boundaries visible via beaker geometry" << G4endl;
        return;
    }
    
    if (solidBoundary) {
        fSourceBoundaryLogical = new G4LogicalVolume(solidBoundary,
                                                     fVacuum,
                                                     "SourceBoundary_Logical");
        
        // Set visualization attributes - bright yellow wireframe
        G4VisAttributes* boundaryVis = new G4VisAttributes(G4Colour(1.0, 1.0, 0.0, 0.8));
        boundaryVis->SetForceWireframe(true);
        boundaryVis->SetLineWidth(3);
        boundaryVis->SetForceSolid(false);
        fSourceBoundaryLogical->SetVisAttributes(boundaryVis);
        
        new G4PVPlacement(0,
                         position,
                         fSourceBoundaryLogical,
                         "SourceBoundary_phys",
                         worldLogical,
                         false,
                         0);
        
        G4cout << "  ✓ Source boundary visualization constructed" << G4endl;
        G4cout << "  Color: Bright Yellow Wireframe" << G4endl;
    }
}

// ==================================================================
// Set source boundary configuration
// ==================================================================
void HPGeDetectorConstruction::SetSourceBoundaryVisualization(
    const std::string& type,
    const G4ThreeVector& position,
    G4double radius,
    G4double height) {
    
    fSourceBoundaryConfig.visible = true;
    fSourceBoundaryConfig.type = type;
    fSourceBoundaryConfig.position = position;
    fSourceBoundaryConfig.radius = radius;
    fSourceBoundaryConfig.height = height;
    
    G4cout << "\n=== Source Boundary Visualization Configured ===" << G4endl;
    G4cout << "Type: " << type << G4endl;
}

// ==================================================================
// Set detector geometry from config.txt (NEW)
// Must be called BEFORE Initialize()
// ==================================================================
void HPGeDetectorConstruction::SetDetectorGeometry(
    G4double crystalDiameter, G4double crystalLength,
    G4double holeDiameter, G4double holeDepth,
    G4double geDeadLayer, G4double geDeadLayerFront,
    G4double liDeadLayer, G4double liDeadLayerFront,
    G4double alWindowThickness, G4double alCupThickness,
    G4double vacuumGap) {
    
    if (fDetectorConfigs.empty()) {
        G4cout << "⚠ Warning: No detector configuration to modify" << G4endl;
        return;
    }
    
    // Modify the first (FLUKA) configuration
    HPGeConfig& config = fDetectorConfigs[0];
    
    config.crystalDiameter = crystalDiameter;
    config.crystalLength = crystalLength;
    config.holeDiameter = holeDiameter;
    config.holeDepth = holeDepth;
    config.geDeadLayer = geDeadLayer;
    config.geDeadLayerFront = geDeadLayerFront;
    config.liDeadLayer = liDeadLayer;
    config.liDeadLayerFront = liDeadLayerFront;
    config.alWindowThickness = alWindowThickness;
    config.alCupThickness = alCupThickness;
    config.vacuumGap = vacuumGap;
    
    G4cout << "\n=== Detector Geometry Updated from Config ===" << G4endl;
    G4cout << "  Crystal: Ø" << crystalDiameter/mm << " mm × " << crystalLength/mm << " mm" << G4endl;
    G4cout << "  Hole: Ø" << holeDiameter/mm << " mm × " << holeDepth/mm << " mm" << G4endl;
    G4cout << "  Dead layers (Ge): front=" << geDeadLayerFront/mm << " mm, lateral=" << geDeadLayer/mm << " mm" << G4endl;
    G4cout << "  Dead layers (Li): front=" << liDeadLayerFront/mm << " mm, lateral=" << liDeadLayer/mm << " mm" << G4endl;
    G4cout << "  Al window: " << alWindowThickness/mm << " mm, Al cup: " << alCupThickness/mm << " mm" << G4endl;
    G4cout << "  Vacuum gap: " << vacuumGap/mm << " mm" << G4endl;
}

// ==================================================================
// Set extended detector geometry (Canberra-specific parameters)
// ==================================================================
void HPGeDetectorConstruction::SetDetectorGeometryExtended(
    G4double crystalDiameter, G4double crystalLength,
    G4double holeDiameter, G4double holeDepth,
    G4double geDeadLayer, G4double geDeadLayerFront,
    G4double liDeadLayer, G4double liDeadLayerFront,
    G4double alWindowThickness, G4double alCupThickness,
    G4double vacuumGap,
    G4double endcapDiameter, G4double endcapLength,
    const G4String& endcapMaterial, const G4String& windowMaterial,
    G4double geToEndcapDistance,
    G4double holderOuterDiameter1, G4double holderTotalLength,
    const G4String& holderMaterial,
    G4double hdpeThickness) {
    
    // Apply basic geometry first
    SetDetectorGeometry(crystalDiameter, crystalLength, holeDiameter, holeDepth,
                       geDeadLayer, geDeadLayerFront, liDeadLayer, liDeadLayerFront,
                       alWindowThickness, alCupThickness, vacuumGap);
    
    if (fDetectorConfigs.empty()) return;
    
    HPGeConfig& config = fDetectorConfigs[0];
    
    config.endcapDiameter = endcapDiameter;
    config.endcapLength = endcapLength;
    config.endcapMaterial = endcapMaterial;
    config.windowMaterial = windowMaterial;
    config.geToEndcapDistance = geToEndcapDistance;
    config.holderOuterDiameter1 = holderOuterDiameter1;
    config.holderTotalLength = holderTotalLength;
    config.holderMaterial = holderMaterial;
    config.hdpeThickness = hdpeThickness;
    
    // Note: fUseCanberraGeometry is now controlled externally via
    // SetUseCanberraGeometry() based on geometry_mode config parameter.
    // Default remains true (Canberra mode).
    
    G4cout << "\n=== Extended Geometry Updated ===" << G4endl;
    G4cout << "  Mode: " << (fUseCanberraGeometry ? "Canberra" : "FLUKA") << G4endl;
    G4cout << "  Endcap: Ø" << endcapDiameter/mm << " mm, " << endcapMaterial << G4endl;
    G4cout << "  Window: " << windowMaterial << G4endl;
    G4cout << "  Ge to endcap: " << geToEndcapDistance/mm << " mm" << G4endl;
    G4cout << "  Holder: Ø" << holderOuterDiameter1/mm << " mm, " << holderMaterial << G4endl;
    G4cout << "  HDPE: " << hdpeThickness/mm << " mm" << G4endl;
}
