// ==================================================================
// HPGeDetectorConstruction.hh (Enhanced with Canberra Mirion Real Geometry)
// ==================================================================

#ifndef HPGeDetectorConstruction_h
#define HPGeDetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4Sphere.hh"
#include <vector>
#include <string>

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Material;

// ==================================================================
// Configuration structure for HPGe detector
// ==================================================================
struct HPGeConfig {
    G4String detectorName;
    G4int detectorID;
    
    // Crystal dimensions
    G4double crystalDiameter;
    G4double crystalLength;
    G4double holeDiameter;
    G4double holeDepth;
    G4double crystalCornerRadius;      // Corner rounding radius (R1)
    G4double crystalBottomOuterStep;   // Bottom outer step diameter
    G4double crystalBottomInnerStep;   // Bottom inner step diameter
    
    // Dead layer thicknesses
    G4double geDeadLayer;              // Outer electrode (lateral)
    G4double geDeadLayerFront;         // Front dead layer
    G4double liDeadLayer;              // Inner electrode (lateral)
    G4double liDeadLayerFront;         // Inner electrode (front)
    
    // Housing parameters
    G4double alWindowThickness;        // Endcap window thickness
    G4double alCupThickness;           // Endcap side wall thickness
    G4double vacuumGap;                // Gap between crystal and endcap
    
    // Endcap geometry (from characterization sheet)
    G4double endcapDiameter;           // Endcap outer diameter
    G4double endcapLength;             // Endcap total length
    G4String endcapMaterial;           // Endcap material (Al, etc.)
    G4String windowMaterial;           // Window material (carbon_epoxy, Al, etc.)
    G4double geToEndcapDistance;       // Ge front to center endcap outside
    
    // Holder geometry (from characterization sheet)
    G4double holderOuterDiameter1;     // Top ring outer diameter
    G4double holderOuterDiameter2;     // Middle step outer diameter
    G4double holderOuterDiameter3;     // Body outer diameter
    G4double holderInnerDepth1;        // First inner depth
    G4double holderInnerDepth2;        // Second inner depth
    G4double holderTotalLength;        // Total holder length
    G4double holderTopHeight;          // Top ring height
    G4double holderStepHeight;         // Step section height
    G4String holderMaterial;           // Holder material (Al, etc.)
    
    // HDPE spacer
    G4double hdpeThickness;            // HDPE spacer thickness
    
    // Boron coating
    G4double boronThickness;
    
    // Position
    G4double posX, posY, posZ;
    
    // Default constructor with Canberra Mirion S/N:21206 values
    HPGeConfig() :
        detectorName("HPGe"),
        detectorID(0),
        // Crystal: from characterization sheet
        crystalDiameter(60.5*mm),
        crystalLength(45.4*mm),
        holeDiameter(9.5*mm),
        holeDepth(19.0*mm),
        crystalCornerRadius(1.0*mm),
        crystalBottomOuterStep(28.0*mm),
        crystalBottomInnerStep(20.0*mm),
        // Dead layers: from characterization sheet S/N:21206
        geDeadLayer(0.5*mm),           // Outer Electrode Thickness: 0.5 mm
        geDeadLayerFront(0.0004*mm),   // Front Dead Layer: 0.4 um
        liDeadLayer(0.0003*mm),        // Inner Electrode: 0.3 um eq. Ge
        liDeadLayerFront(0.0003*mm),   // Inner Electrode: 0.3 um eq. Ge
        // Endcap
        alWindowThickness(0.6*mm),     // Window thickness from assembly
        alCupThickness(1.5*mm),        // Endcap side wall
        vacuumGap(6.0*mm),            // Ge front to endcap outside: 6.0 mm
        endcapDiameter(76.2*mm),       // Assembly outer diameter
        endcapLength(8.0*mm),          // Endcap wall height from assembly
        endcapMaterial("Al"),
        windowMaterial("carbon_epoxy"),
        geToEndcapDistance(6.0*mm),     // Ge front to center endcap outside
        // Holder: from characterization sheet
        holderOuterDiameter1(69.0*mm),
        holderOuterDiameter2(64.6*mm),
        holderOuterDiameter3(63.0*mm),
        holderInnerDepth1(66.7*mm),
        holderInnerDepth2(69.8*mm),
        holderTotalLength(95.2*mm),
        holderTopHeight(8.6*mm),
        holderStepHeight(19.4*mm),
        holderMaterial("Al"),
        // HDPE
        hdpeThickness(6.5*mm),
        // Boron coating
        boronThickness(0.03*mm),
        posX(0.), posY(0.), posZ(0.) {}
};

// ==================================================================
// Marinelli beaker configuration structure
// ==================================================================
struct MarinelliConfig {
    G4String beakerType;
    G4double outerDiameter;
    G4double innerDiameter;
    G4double totalHeight;
    G4double innerHeight;
    G4double wallThickness;
    G4String fillMaterial;
    G4bool isActive;
    
    // Default constructor for 1000ml
    MarinelliConfig() :
        beakerType("1000ml"),
        outerDiameter(130.*mm),
        innerDiameter(85.*mm),
        totalHeight(152.*mm),
        innerHeight(77.*mm),
        wallThickness(2.*mm),
        fillMaterial("water"),
        isActive(false) {}
    
    // Set dimensions based on type
    void SetType(const std::string& type) {
        if (type == "450ml") {
            beakerType = "450ml";
            outerDiameter = 114.*mm;
            innerDiameter = 77.*mm;
            totalHeight = 101.*mm;
            innerHeight = 68.*mm;
        } else {
            beakerType = "1000ml";
            outerDiameter = 130.*mm;
            innerDiameter = 85.*mm;
            totalHeight = 152.*mm;
            innerHeight = 77.*mm;
        }
    }
    
    // Calculate Z position
    G4double GetZPosition() const {
        return -2.*mm - (innerHeight / 2.0);
    }
};

// ==================================================================
// Disk source configuration structure
// ==================================================================
struct DiskSourceConfig {
    G4bool isActive;
    G4double radius;              // Disk radius
    G4double thickness;           // Disk thickness (height)
    G4double positionZ;           // Z position (negative = below detector)
    G4String material;            // Source material (water, soil, etc.)
    
    DiskSourceConfig() :
        isActive(false),
        radius(25.*mm),
        thickness(5.*mm),
        positionZ(-50.*mm),
        material("water") {}
};

// ==================================================================
// Cylindrical volume source configuration structure
// ==================================================================
struct CylinderSourceConfig {
    G4bool isActive;
    G4double innerRadius;         // Inner radius (sample volume)
    G4double outerRadius;         // Outer radius (including walls)
    G4double height;              // Total height
    G4double wallThickness;       // Wall thickness
    G4double bottomThickness;     // Bottom wall thickness
    G4double positionZ;           // Z position of bottom (negative = below detector)
    G4String wallMaterial;        // Container material (polypropylene, glass, etc.)
    G4String fillMaterial;        // Sample material (water, soil, etc.)
    
    CylinderSourceConfig() :
        isActive(false),
        innerRadius(30.*mm),
        outerRadius(32.*mm),
        height(50.*mm),
        wallThickness(2.*mm),
        bottomThickness(2.*mm),
        positionZ(-60.*mm),
        wallMaterial("polypropylene"),
        fillMaterial("water") {}
};

// ==================================================================
// Cartridge source configuration (Orano LEA Type D)
// ==================================================================
struct CartridgeSourceConfig {
    G4bool isActive;
    G4double housingDiameter;     // Outer housing diameter
    G4double housingHeight;       // Total housing height
    G4double activeDiameter;      // Active matrix diameter
    G4double activeThickness;     // Active matrix thickness
    G4double wallThickness;       // Housing wall thickness
    G4double positionZ;           // Z position (negative = below detector)
    G4String housingMaterial;     // Housing material (polycarbonate)
    G4String activeMaterial;      // Active material (activated_carbon)
    
    CartridgeSourceConfig() :
        isActive(false),
        housingDiameter(56.0*mm),     // Type D: Ø56 mm
        housingHeight(26.5*mm),       // Type D: 26.5 mm
        activeDiameter(52.0*mm),      // Active matrix: Ø52 mm
        activeThickness(20.5*mm),     // Active matrix: 20.5 mm
        wallThickness(2.0*mm),        // Polycarbonate wall
        positionZ(-40.*mm),
        housingMaterial("polycarbonate"),
        activeMaterial("activated_carbon") {}
};

// ==================================================================
// Filter source configuration (Orano LEA Type M)
// ==================================================================
struct FilterSourceConfig {
    G4bool isActive;
    G4double outerDiameter;       // Overall filter diameter
    G4double activeDiameter;      // Active surface diameter
    G4double filterThickness;     // Paper filter thickness
    G4double sealThickness;       // Polyester seal sheet thickness
    G4double positionZ;           // Z position (negative = below detector)
    G4String filterMaterial;      // Filter material (water approximation for deposited isotopes)
    G4String sealMaterial;        // Seal material (polyester)
    G4String filterType;          // Filter type code (M43, M50, M53, etc.)
    
    FilterSourceConfig() :
        isActive(false),
        outerDiameter(53.0*mm),       // Type M53: Ø53 mm
        activeDiameter(47.0*mm),      // Active surface: Ø47 mm
        filterThickness(0.5*mm),      // Paper filter ~0.5 mm
        sealThickness(0.05*mm),       // Polyester seal sheet ~50 um
        positionZ(-30.*mm),
        filterMaterial("water"),
        sealMaterial("polyester"),
        filterType("M53") {}
    
    void SetType(const std::string& type) {
        filterType = type;
        if (type == "M43" || type == "M43-51") {
            activeDiameter = 43.0*mm;
            outerDiameter = 51.0*mm;
        } else if (type == "M45" || type == "M45-51") {
            activeDiameter = 45.0*mm;
            outerDiameter = 51.0*mm;
        } else if (type == "M47" || type == "M47-53") {
            activeDiameter = 47.0*mm;
            outerDiameter = 53.0*mm;
        } else if (type == "M50") {
            activeDiameter = 50.0*mm;
            outerDiameter = 63.0*mm;
        } else if (type == "M53") {
            activeDiameter = 47.0*mm;
            outerDiameter = 53.0*mm;
        } else if (type == "M60") {
            activeDiameter = 53.0*mm;
            outerDiameter = 60.0*mm;
        } else if (type == "M120") {
            activeDiameter = 120.0*mm;
            outerDiameter = 130.0*mm;
        }
    }
};

// ==================================================================
// Source boundary visualization configuration
// ==================================================================
struct SourceBoundaryConfig {
    G4bool visible;
    std::string type;
    G4ThreeVector position;
    G4double radius;
    G4double height;
    
    SourceBoundaryConfig() :
        visible(false), type("none"),
        position(0, 0, 0), radius(0), height(0) {}
};

// ==================================================================
// HPGe Detector Construction Class
// ==================================================================
class HPGeDetectorConstruction : public G4VUserDetectorConstruction {
public:
    HPGeDetectorConstruction();
    virtual ~HPGeDetectorConstruction();
    
    // Main construction methods
    virtual G4VPhysicalVolume* Construct();
    virtual void ConstructSDandField();
    
    // Configuration methods
    void AddDetector(const HPGeConfig& config);
    void SetSourcePosition(G4double x, G4double y, G4double z);
    void SetFLUKAGeometryParameters(HPGeConfig& config);
    void SetCanberraGeometryParameters(HPGeConfig& config);
    
    // Marinelli beaker configuration
    void SetMarinelliConfig(const MarinelliConfig& config);
    void EnableMarinelli(G4bool enable);
    MarinelliConfig GetMarinelliConfig() const { return fMarinelliConfig; }
    G4bool IsMarinelliActive() const { return fMarinelliConfig.isActive; }
    
    // Disk source configuration
    void SetDiskSourceConfig(const DiskSourceConfig& config);
    void EnableDiskSource(G4bool enable);
    DiskSourceConfig GetDiskSourceConfig() const { return fDiskSourceConfig; }
    G4bool IsDiskSourceActive() const { return fDiskSourceConfig.isActive; }
    
    // Cylinder source configuration
    void SetCylinderSourceConfig(const CylinderSourceConfig& config);
    void EnableCylinderSource(G4bool enable);
    CylinderSourceConfig GetCylinderSourceConfig() const { return fCylinderSourceConfig; }
    G4bool IsCylinderSourceActive() const { return fCylinderSourceConfig.isActive; }
    
    // Cartridge source configuration (Orano LEA Type D)
    void SetCartridgeSourceConfig(const CartridgeSourceConfig& config);
    void EnableCartridgeSource(G4bool enable);
    CartridgeSourceConfig GetCartridgeSourceConfig() const { return fCartridgeSourceConfig; }
    G4bool IsCartridgeSourceActive() const { return fCartridgeSourceConfig.isActive; }
    
    // Filter source configuration (Orano LEA Type M)
    void SetFilterSourceConfig(const FilterSourceConfig& config);
    void EnableFilterSource(G4bool enable);
    FilterSourceConfig GetFilterSourceConfig() const { return fFilterSourceConfig; }
    G4bool IsFilterSourceActive() const { return fFilterSourceConfig.isActive; }
    
    // Detector geometry configuration from config.txt
    void SetDetectorGeometry(G4double crystalDiameter, G4double crystalLength,
                            G4double holeDiameter, G4double holeDepth,
                            G4double geDeadLayer, G4double geDeadLayerFront,
                            G4double liDeadLayer, G4double liDeadLayerFront,
                            G4double alWindowThickness, G4double alCupThickness,
                            G4double vacuumGap);
    
    // Extended geometry with Canberra-specific parameters
    void SetDetectorGeometryExtended(
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
        G4double hdpeThickness);
    
    // Source boundary visualization
    void ConstructSourceBoundary(G4LogicalVolume* worldLogical);
    void SetSourceBoundaryVisualization(const std::string& type,
                                       const G4ThreeVector& position,
                                       G4double radius = 0,
                                       G4double height = 0);
    
    // Getters
    G4double GetActiveVolume(G4int detectorID = 0) const;
    G4double GetActiveMass(G4int detectorID = 0) const;
    G4ThreeVector GetSourcePosition() const { return fSourcePosition; }
    const std::vector<HPGeConfig>& GetDetectorConfigs() const { return fDetectorConfigs; }
    G4double GetCrystalToSourceDistance(G4int detectorID = 0) const;
    G4LogicalVolume* GetMarinelliLogicalVolume() const { return fMarinelliLogical; }
    
private:
    // Material definition
    void DefineMaterials();
    
    // Detector construction methods
    G4LogicalVolume* ConstructHPGeDetector(const HPGeConfig& config);
    G4LogicalVolume* ConstructHPGeDetectorFLUKA(const HPGeConfig& config);
    G4LogicalVolume* ConstructHPGeDetectorCanberra(const HPGeConfig& config);
    void ConstructLeadCastle(G4LogicalVolume* motherVolume);
    
    // Marinelli beaker construction
    void ConstructMarinelliBeaker(G4LogicalVolume* worldLogical);
    
    // Disk source construction
    void ConstructDiskSource(G4LogicalVolume* worldLogical);
    
    // Cylinder source construction
    void ConstructCylinderSource(G4LogicalVolume* worldLogical);
    
    // Cartridge source construction (Orano LEA Type D)
    void ConstructCartridgeSource(G4LogicalVolume* worldLogical);
    
    // Filter source construction (Orano LEA Type M)
    void ConstructFilterSource(G4LogicalVolume* worldLogical);
    
    // Helper: Get material by name
    G4Material* GetMaterialByName(const G4String& name);
    
    // Calculation methods
    G4double CalculateActiveVolume(const HPGeConfig& config) const;
    
    // Member variables
    G4ThreeVector fSourcePosition;
    G4double fWorldSize;
    
    // Materials
    G4Material* fGe;
    G4Material* fAl;
    G4Material* fC;
    G4Material* fVacuum;
    G4Material* fAir;
    G4Material* fLi;
    G4Material* fBoron;
    G4Material* fLead;
    G4Material* fWater;
    G4Material* fSteel;
    G4Material* fPerspex;
    G4Material* fPolypropylene;
    G4Material* fSoil;
    G4Material* fGlass;
    G4Material* fPVC;
    G4Material* fHDPE;
    G4Material* fCarbonEpoxy;      // Window material for Canberra detector
    G4Material* fPolycarbonate;     // Housing for Orano LEA sources
    G4Material* fPolyester;         // Seal for paper filter sources
    G4Material* fActivatedCarbon;   // Active matrix for cartridge sources
    
    // Geometry mode flag
    G4bool fUseCanberraGeometry;
    
    // Detector configurations
    std::vector<HPGeConfig> fDetectorConfigs;
    std::vector<G4LogicalVolume*> fActiveGeLogicals;
    
    // Marinelli beaker
    MarinelliConfig fMarinelliConfig;
    G4LogicalVolume* fMarinelliLogical;
    G4VPhysicalVolume* fWorldPhysical;
    
    // Disk source
    DiskSourceConfig fDiskSourceConfig;
    G4LogicalVolume* fDiskSourceLogical;
    
    // Cylinder source
    CylinderSourceConfig fCylinderSourceConfig;
    G4LogicalVolume* fCylinderSourceLogical;
    G4LogicalVolume* fCylinderWallLogical;
    
    // Cartridge source (Orano LEA Type D)
    CartridgeSourceConfig fCartridgeSourceConfig;
    G4LogicalVolume* fCartridgeSourceLogical;
    G4LogicalVolume* fCartridgeHousingLogical;
    
    // Filter source (Orano LEA Type M)
    FilterSourceConfig fFilterSourceConfig;
    G4LogicalVolume* fFilterSourceLogical;
    G4LogicalVolume* fFilterSealLogical;
    
    // Source boundary visualization
    SourceBoundaryConfig fSourceBoundaryConfig;
    G4LogicalVolume* fSourceBoundaryLogical;
};

#endif
