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
// ScintSD (scintillator sensitive detector definition) program
// --------------------------------------------------------------

#include "DMXScintSD.hh"

#include "DMXScintHit.hh"
#include "DMXDetectorConstruction.hh"

#include "G4VPhysicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4SDManager.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Ions.hh"
#include "G4ios.hh"
//Logan 
#include "G4VProcess.hh"
#include "G4AnalysisManager.hh"
//#include "G4root.hh"
//

#include "G4OpticalPhoton.hh"

//Logan
#include <iostream>
#include <fstream>
#include <cstring>
//


int n=1;
int n1=1;
int HasHit=0;
int HasHitAr=0;
G4double eki = 0;
G4double eke = 0;
G4ThreeVector Position(0.,0.,0.);
G4ThreeVector MomentumDire(0.,0.,0.);
G4ThreeVector MomentumDiri(0.,0.,0.);

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXScintSD::DMXScintSD(G4String name) 
  :G4VSensitiveDetector(name)
{
  G4String HCname="scintillatorCollection";
  collectionName.insert(HCname);

  Info1.open("General.csv");
  Info.open("Informacion.csv");
  Info << "Event,"<<"Hit,"<<"KE_i,"<<"KE_e";
  Info1 << "Event,"<<"KE_e,"<<"x,"<<"y,"<<"z";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXScintSD::~DMXScintSD(){ Info.close(); }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::Initialize(G4HCofThisEvent*)
{
  scintillatorCollection = new DMXScintHitsCollection
    (SensitiveDetectorName,collectionName[0]);

  HitID = -1;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

G4bool DMXScintSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{

  //need to know if this is an optical photon and exclude it:
  if(aStep->GetTrack()->GetDefinition()
     == G4OpticalPhoton::OpticalPhotonDefinition()) return false;
  /*if(aStep->GetTrack()->GetDefinition() == G4Electron::ElectronDefinition()){
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return false;}
  if(aStep->GetTrack()->GetDefinition() == G4Proton::ProtonDefinition()){
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return false;}
  if(aStep->GetTrack()->GetDefinition() != G4Neutron::NeutronDefinition()){
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    return false;}*/
  /*if (aStep->GetTrack()->GetTrackID() != 1){
    aStep->GetTrack()->SetTrackStatus(fKillTrackAndSecondaries);
    }*/

  
  G4double edep = aStep->GetTotalEnergyDeposit();
  G4double ek = aStep->GetPostStepPoint()->GetKineticEnergy();
  G4ParticleDefinition* particleType = aStep->GetTrack()->GetDefinition();
  G4String particleName = particleType->GetParticleName();

  G4double posx = aStep->GetPreStepPoint()->GetPosition().x();
  G4double posy = aStep->GetPreStepPoint()->GetPosition().y();
  G4double posz = aStep->GetPreStepPoint()->GetPosition().z();
  G4double momx = aStep->GetPreStepPoint()->GetMomentumDirection().x();
  G4double momy = aStep->GetPreStepPoint()->GetMomentumDirection().y();
  G4double momz = aStep->GetPreStepPoint()->GetMomentumDirection().z();

  G4String Volume = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
  G4bool FirstStep = aStep-> IsFirstStepInVolume();
  G4bool LastStep = aStep-> IsLastStepInVolume();
  G4double StepLength = aStep->GetStepLength();
  G4bool hit = false;

  G4String IntProcessName = aStep -> GetPostStepPoint() -> GetProcessDefinedStep() -> GetProcessName();

  G4double stepl = 0.;
  if (particleType->GetPDGCharge() != 0.)
    stepl = aStep->GetStepLength();
  
  // if ((edep==0.)&&(stepl==0.)) return false;      


  // fill in hit
  DMXScintHit* newHit = new DMXScintHit();
  newHit->SetEdep(edep);
  newHit->SetPos(aStep->GetPostStepPoint()->GetPosition());
  newHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
  newHit->SetParticle(particleName);
  newHit->SetParticleEnergy(aStep->GetPreStepPoint()->GetKineticEnergy() );

  HitID = scintillatorCollection->insert(newHit);

  if(posx == 0 && posy == 0 && posz == 0 && n>0){
  eki=ek;
  MomentumDiri.setX(momx);
  MomentumDiri.setY(momy); 
  MomentumDiri.setZ(momz);
  }

  if(particleName == "neutron" && Volume == "physWorld" && HasHitAr == 0){
  HasHit = HasHit+1;
  eke = ek;
  Position.setX(posx);
  Position.setY(posy); 
  Position.setZ(posz);
  MomentumDire.setX(momx);
  MomentumDire.setY(momy); 
  MomentumDire.setZ(momz);
  }
  
  if(particleName == "neutron" && Volume == "physAr"){
  HasHitAr = 1;
  }

  if(particleName == "neutron" && (Volume == "physS" || Volume == "physSap")){ //In case a neutron escapes the sapphire/poly and hits the poly/sapphire again
  HasHit = 0;
  eke = 0;
  Position.setX(0);
  Position.setY(0); 
  Position.setZ(0);
  MomentumDire.setX(0);
  MomentumDire.setY(0); 
  MomentumDire.setZ(0);
  }

  /*if(n1!=n){
  HasHit = 0;
  eke = 0;
  Position.setX(0);
  Position.setY(0); 
  Position.setZ(0);
  n1=n;}*/
  
  /*if(IntProcessName == "nCapture"){
  HasHit = HasHit+1;
  }

  if(IntProcessName == "nCapture" || IntProcessName == "hadElastic"){
  zCol = posz;
  }*/

  // Info << '\n' << n << "," << HasHitAr << "," << Volume;

  return true;
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::EndOfEvent(G4HCofThisEvent* HCE)
{

  G4String HCname = collectionName[0];
  static G4int HCID = -1;
  if(HCID<0)
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(HCname);
  HCE->AddHitsCollection(HCID,scintillatorCollection);



  //ROOT
  auto aMan = G4AnalysisManager::Instance();
  aMan->FillNtupleDColumn(0,HasHit);
  aMan->FillNtupleDColumn(1,Position.getX());
  aMan->FillNtupleDColumn(2,Position.getY());
  aMan->FillNtupleDColumn(3,Position.getZ());
  aMan->FillNtupleDColumn(4,eki);
  aMan->FillNtupleDColumn(5,eke);
  aMan->FillNtupleDColumn(6,MomentumDire.getX());
  aMan->FillNtupleDColumn(7,MomentumDire.getY());
  aMan->FillNtupleDColumn(8,MomentumDire.getZ());
  aMan->FillNtupleDColumn(9,HasHitAr);
  aMan->FillNtupleDColumn(10,MomentumDiri.getX());
  aMan->FillNtupleDColumn(11,MomentumDiri.getY());
  aMan->FillNtupleDColumn(12,MomentumDiri.getZ());
  aMan->AddNtupleRow();


  /*if(n==0){
  Info << n <<"," << HasHit <<"," << eki <<"," << eke;
  	if(HasHit > 0){
  	Info1 << n <<"," << eke <<"," << Position.getX()<<"," << Position.getY()<<"," << Position.getZ();}
  }
  if(n>0){
  Info << '\n' << n <<"," << HasHit <<"," << eki <<"," << eke;
  	if(HasHit > 0){
  	Info1 << '\n' << n <<"," << eke <<"," << Position.getX()<<"," << Position.getY()<<"," << Position.getZ();}
  }*/


  G4int nHits = scintillatorCollection->entries();
  if (verboseLevel>=1)
    G4cout << "     LXe collection: " <<  nHits << " hits" << G4endl;
    n++;
  if (verboseLevel>=2)
    scintillatorCollection->PrintAllHits();


  HasHit = 0;
  HasHitAr = 0;
  eke = 0;
  Position.setX(0);
  Position.setY(0); 
  Position.setZ(0);
  MomentumDire.setX(0);
  MomentumDire.setY(0);
  MomentumDire.setZ(0);
  MomentumDiri.setX(0);
  MomentumDiri.setY(0); 
  MomentumDiri.setZ(0);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::clear()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::DrawAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXScintSD::PrintAll()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

