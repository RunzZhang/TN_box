//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// --------------------------------------------------------------
//   GEANT 4 - Underground Dark Matter Detector Advanced Example
//
//      For information related to this code contact: Alex Howard
//      e-mail: alexander.howard@cern.ch
// --------------------------------------------------------------
// Comments
//
//                  Underground Advanced
//               by A. Howard and H. Araujo 
//                    (27th November 2001)
//
// DetectorConstruction program
// --------------------------------------------------------------

#include "DMXDetectorConstruction.hh"
#include "DMXDetectorMessenger.hh"

#include "DMXScintSD.hh"
#include "DMXPmtSD.hh"


#include "G4Material.hh"
#include "G4MaterialTable.hh"
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4UnitsTable.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpBoundaryProcess.hh"

#include "G4FieldManager.hh"
#include "G4UniformElectricField.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4EqMagElectricField.hh"
#include "G4ClassicalRK4.hh"
#include "G4ChordFinder.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UserLimits.hh"

#include "G4RunManager.hh"
#include "G4SystemOfUnits.hh"
//Logan
#include "G4NCrystal/G4NCrystal.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::DMXDetectorConstruction()  
{
  // create commands for interactive definition of time cuts:
  detectorMessenger = new DMXDetectorMessenger(this);

  theUserLimitsForRoom     = 0; 
  theUserLimitsForDetector = 0; 
  // default time cut = infinite
  //  - note also number of steps cut in stepping action = MaxNoSteps
  theMaxTimeCuts      = DBL_MAX;
  theMaxStepSize      = DBL_MAX;
  theDetectorStepSize = DBL_MAX;
  theRoomTimeCut      = 1000. * nanosecond;
  theMinEkine         = 250.0*eV; // minimum kinetic energy required in volume
  theRoomMinEkine     = 250.0*eV; // minimum kinetic energy required in volume
  
  //Zero the G4Cache objects to contain logical volumes
  LXeSD.Put(0);
  pmtSD.Put(0);
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
DMXDetectorConstruction::~DMXDetectorConstruction() 
{
  delete theUserLimitsForRoom;
  delete theUserLimitsForDetector;
  delete detectorMessenger;
}



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
void DMXDetectorConstruction::DefineMaterials() 
{

#include "DMXDetectorMaterial.icc"

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....
G4VPhysicalVolume* DMXDetectorConstruction::Construct() {

  DefineMaterials();
  // Get nist material manager
  //G4NistManager* nist = G4NistManager::Instance();
  
  //Logan
  G4Material * mat_aluminium = G4NCrystal::createMaterial("Al_sg225.ncmat");
  G4Material* sapphireNCrystal_mat = G4NCrystal::createMaterial("Al2O3_sg167_Corundum.ncmat;bragg=0"); 
  G4Material* HDPENCrystal_mat = G4NCrystal::createMaterial("Polyethylene_CH2.ncmat;density=0.96gcm3");
  G4Material* vacuumNCrystal_mat = G4NCrystal::createMaterial("void.ncmat");

  // Envelope parameters
  //
  G4double S_l = 7.5*cm, S_w=5*cm, V_l = 3*cm,N_l = 10*cm, P_w = 25*cm, P_l = 6*cm, P_p = (S_l+P_l)+N_l, A_l = 30.48*cm;
  //python for quick test: S_l = 75, S_w =5, V_l=3, N_l =10, P_w =25, P_l =6, P_p= S_l+P_l+N_l, A_l = 30.48, SBC_l = 200, SD_T = 2, SD_D =100
  //sapphire window length and width, S_l and S_w; V_l empty source length; N_l, vacuum notch length; P_w lenth of HPDE subtracted by the vaccum source space on x or y direction;
   //P_l is distance between sapphire and source,P_p is the middle HDPE + sapphire +10cm(vacuum notch length=10cm) offset
  //A_l test argon length
  G4double SBC_l = 2*m, SD_T = 2*cm, SD_D = 1*m;
  // SBC_length is the height of SBC chamber and SD_T is sensitive detector thickness, which is test argon thickness in z axis, and SD_D is the sensitive detector distance which is the 2nd SD distance to HDPE box

  //G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 3.5*m;
  G4double world_sizeZ  = 3.5*m;
  //G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       world_sizeXY, world_sizeXY, world_sizeZ);     //its size
      
    logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        vacuum_mat,           //its material
                        "logicWorld");            //its name
                                   
    physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "physWorld",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking

                     

  // Filter
  G4Box* solidPoly = new G4Box("solidPoly", 0.5*(2*P_w+3*V_l+2*P_p), 0.5*(2*P_w+3*V_l+2*P_p), 0.5*(P_w+V_l+P_p)); //0.5*(P_w+V_l+P_l+S_l)
  G4Box* vacuumNotch = new G4Box("vacuumNotch", 0.5*S_w, 0.5*S_w, 0.5*(P_p-P_l-S_l)); //0.5*(P_w+V_l+P_l+S_l) //0.5*(P_p-P_l-S_l)
//  G4Tubs* vacuumNotch = new G4Tubs("vacuumNotch",0,           // inner radius
//                                   0.5*S_w,      // outer radius
//                                   0.5*(P_p-P_l-S_l),  // half-length along Z
//                                   0,           // start angle
//                                   2*M_PI);     // full cylinder
  G4ThreeVector trans1 = G4ThreeVector(0.,0.,0.5*(V_l+2*P_p-N_l)); // make z+ side the vaccum notch touch the Poly edge
  G4ThreeVector trans2 = G4ThreeVector(0,0.5*(V_l+2*P_p-N_l),0.); // make z+ side the vaccum notch touch the Poly edge
  G4RotationMatrix rot1 = G4RotationMatrix(0.,0.,0.);
  G4RotationMatrix rot2;
  rot2.rotateX(90*deg);
  G4Transform3D transform1 = G4Transform3D(rot1,trans1);
  G4Transform3D transform2 = G4Transform3D(rot2,trans2);
  G4VSolid* solidS1 = new G4SubtractionSolid("solidS1" ,solidPoly, vacuumNotch, transform1);
  G4VSolid* solidS2 = new G4SubtractionSolid("solidS2" ,solidS1, vacuumNotch, transform2);
  logicS = new G4LogicalVolume(solidS2, HDPENCrystal_mat, "logicS"); //HDPENCrystal_mat
  physS = new G4PVPlacement(0, G4ThreeVector(0.,0.,-0.5*(P_p-P_w)), logicS, "physS", logicWorld, false, 0);
 
  // Sapphire Window
  G4Box* solidSap = new G4Box("solidSap", 0.5*S_w, 0.5*S_w, 0.5*S_l);
  logicSap = new G4LogicalVolume(solidSap, sapphireNCrystal_mat, "logicSap");  //sapphireNCrystal_mat
   physSap1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(S_l+V_l+2*P_l)+0.5*(P_p-P_w)), logicSap, "physSap1", logicS, false, 0);



   physSap2 = new G4PVPlacement(&rot2, G4ThreeVector(0,0.5*(S_l+V_l+2*P_l),+0.5*(P_p-P_w)), logicSap, "physSap2", logicS, false, 0);

//   G4Tubs* solidSap = new G4Tubs("solidSap",
//                                   0,           // inner radius
//                                   0.5*S_w,      // outer radius
//                                   0.5*S_l,  // half-length along Z
//                                   0,           // start angle
//                                   2*M_PI);     // full cylinder
//  logicSap = new G4LogicalVolume(solidSap, sapphireNCrystal_mat, "logicSap");  //sapphireNCrystal_mat
//   physSap = new G4PVPlacement(0, G4ThreeVector(0.,0.,(0.5*(V_l+S_l)+P_l)-0.5*(P_p-P_w)), logicSap, "physSap", logicS, false, 0);


  // Sapphire Window for test cross section
  //part 1
  //G4Box* solidSap = new G4Box("solidSap", 0.5*S_w, 0.5*S_w, 1*mm);
  //logicSap = new G4LogicalVolume(solidSap, sapphireNCrystal_mat, "logicSap");  //sapphireNCrystal_mat
  //physSap = new G4PVPlacement(0, G4ThreeVector(1*m,0.,(0.5*(V_l+S_l)+P_l)-0.5*(P_p-P_w)), logicSap, "physSap", logicWorld, false, 0);
  //part 2
  //G4Box* solidSap2 = new G4Box("solidSap2", 0.5*S_w, 0.5*S_w, 0.5*S_l);
  //logicSap2 = new G4LogicalVolume(solidSap2, sapphireNCrystal_mat, "logicSap2");  //sapphireNCrystal_mat
  //physSap2 = new G4PVPlacement(0, G4ThreeVector(1*m,0.,(0.5*(V_l+S_l)+P_l)-0.5*(P_p-P_w)), logicSap2, "physSap2", logicWorld, false, 0);  // transform to make test cross section


  // Test Argon you also need to comment out "physAr" in the last section in this code
//  G4Box* solidAr = new G4Box("solidAr", 0.5*A_l, 0.5*A_l, 0.5*A_l);
//  logicAr = new G4LogicalVolume(solidAr, LAr_mat, "logicAr");  //sapphire_mat
//  physAr = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(P_w+V_l+P_p)+0.5*(P_p-P_w)+1*m+0.5*A_l), logicAr, "physAr", logicWorld, false, 0);

  G4Box* solidVac2 = new G4Box("solidVac2", 0.5*SBC_l, 0.5*SBC_l, 0.5*SD_T);
  logicVac2 = new G4LogicalVolume(solidVac2, vacuumNCrystal_mat, "logicVac2");  //sapphire_mat
  physVac2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(V_l+P_p)+0.5*(P_p-P_w)+0.5*SD_T), logicVac2, "physVac2", logicWorld, false, 0);

  G4Box* solidVac3 = new G4Box("solidVac3", 0.5*SBC_l, 0.5*SBC_l, 0.5*SD_T);
  logicVac3 = new G4LogicalVolume(solidVac3, vacuumNCrystal_mat, "logicVac3");  //sapphire_mat
//  physVac3 = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(V_l*P_p)+0.5*(P_p-P_w)+SD_D+0.5*SD_T), logicVac3, "physVac3", logicWorld, false, 0);

  physVac3 = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.5*(V_l*P_p)+0.5*(P_p-P_w)+10*cm+0.5*SD_T), logicVac3, "physVac3", logicWorld, false, 0);

  // SD before
 // G4Box* solidSD1 = new G4Box("solidSD1", S_l, S_l, 1*mm); 
  //logicSD1 = new G4LogicalVolume(solidSD1, vacuum_mat, "logicSD1");                    
  //physSD1 = new G4PVPlacement(0, G4ThreeVector(0.,0.,5*cm), logicSD1, "physSD1", logicWorld, false, 0);   
  
  // Empty Inside
  G4Box* solidSD2 = new G4Box("solidSD2", 0.5*V_l, 0.5*V_l, 0.5*V_l);

  logicSD2 = new G4LogicalVolume(solidSD2, vacuumNCrystal_mat, "logicSD2");
  physSD2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,+0.5*(P_p-P_w)), logicSD2, "physSD2", logicS, false, 0);  //0.5*(P_l+S_l-P_w)
  //physSD3 = new G4PVPlacement(0, G4ThreeVector(0.,0.,V_l+P_p -0.5*(2*P_p-P_w+V_l)), logicSD2, "physSD3", logicS, false, 0);  //0.5*(P_l+S_l-P_w)
  //cylinder empty
//  G4Tubs* solidSD2 = new G4Tubs("solidSD2",
//                                   0,           // inner radius
//                                   0.5*V_l,      // outer radius
//                                   0.5*V_l,  // half-length along Z
//                                   0,           // start angle
//                                   2*M_PI);     // full cylinder
//  logicSD2 = new G4LogicalVolume(solidSD2, vacuumNCrystal_mat, "logicSD2");
//  physSD2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,-0.5*(P_p-P_w)), logicSD2, "physSD2", logicS, false, 0);  //0.5*(P_l+S_l-P_w)

 /* // Sapphire Window
  G4Box* solidWindow = new G4Box("solidWindow", 0.25*(S_l-V_l), 0.5*V_l, 0.5*V_l); 
  logicWindow = new G4LogicalVolume(solidWindow, vacuum_mat, "logicWindow");                    
  physWindow = new G4PVPlacement(0, G4ThreeVector(-0.25*(S_l+V_l),0.,0.), logicWindow, "physWindow", logicWorld, false, 0); */

  //
  //always return the physical World
  //
  return physWorld;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXDetectorConstruction::ConstructSDandField()
{
  // ......................................................................
  // sensitive detectors ..................................................
  // ......................................................................

  if (LXeSD.Get() == 0) 
    {    
      G4String name="/DMXDet/LXeSD";
      DMXScintSD* aSD = new DMXScintSD(name);
      LXeSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(LXeSD.Get()); 
  if(logicS){
//      SetSensitiveDetector(logicS,LXeSD.Get());
//      SetSensitiveDetector(logicSD2,LXeSD.Get());
      SetSensitiveDetector(logicWorld,LXeSD.Get());
//      SetSensitiveDetector(logicSap,LXeSD.Get());
//      SetSensitiveDetector(logicAr,LXeSD.Get());
      }
  /*if (LXe_log)    
    SetSensitiveDetector(LXe_log,LXeSD.Get());

  if (pmtSD.Get() == 0)
    {
      G4String name="/DMXDet/pmtSD";
      DMXPmtSD* aSD = new DMXPmtSD(name);
      pmtSD.Put(aSD);
    }
  G4SDManager::GetSDMpointer()->AddNewDetector(pmtSD.Get()); 
  if (phcath_log)
    SetSensitiveDetector(phcath_log,pmtSD.Get());*/

  return;
}
 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetRoomEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for ROOM
  theRoomMinEkine = val;
  if (theUserLimitsForRoom != 0) 
    {
      theUserLimitsForRoom->SetUserMinEkine(val); 
      G4cout << " Changing Room energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMinEkine
void DMXDetectorConstruction::SetEnergyCut(G4double val)
{
  // set minimum charged particle energy cut - NB: for Xenon Detector
  theMinEkine = val;
  if (theUserLimitsForDetector != 0) 
    {
      theUserLimitsForDetector->SetUserMinEkine(val);
      G4cout << "Changing Detector energy cut to: " << G4BestUnit(val,"Energy")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetRoomTimeCut(G4double val)
{
  // set room time cut:
  theRoomTimeCut = val;
  if (theUserLimitsForRoom != 0) 
    {
      theUserLimitsForRoom->SetUserMaxTime(val);
      G4cout << " Changing Room Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

// specific method to G4UserLimits:= SetUserMaxTime
void DMXDetectorConstruction::SetTimeCut(G4double val)
{
  // set detector time cut:
  theMaxTimeCuts = val;
  if (theUserLimitsForDetector != 0) 
    {
      theUserLimitsForDetector->SetUserMaxTime(val);
      G4cout << " Changing Detector Time cut to: " << G4BestUnit(val,"Time")
	     << G4endl;
    }
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....



