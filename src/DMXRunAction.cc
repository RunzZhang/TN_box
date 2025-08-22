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
// History:
// 17 Jan 2002 Alex Howard Added Analysis
// 23 Oct 2009 Luciano Pandola Removed un-necessary calls from EndOfRun()
//
// RunAction program
// --------------------------------------------------------------

#include "DMXRunActionMessenger.hh"
#include "DMXRunAction.hh"

#include "G4Run.hh"
#include "G4ios.hh"
#include "G4AnalysisManager.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunAction::DMXRunAction()
{
  runMessenger = new DMXRunActionMessenger(this);

  // Logan ROOT
  /*auto aMan = G4AnalysisManager::Instance(); 
  aMan->SetVerboseLevel( 1 );
  //aMan>SetNtupleMerging( 1 );    
  aMan->CreateNtuple("tree", "tree");
  aMan->CreateNtupleDColumn("Hit"); //0
  aMan->CreateNtupleDColumn("x"); //1
  aMan->CreateNtupleDColumn("y"); //2
  aMan->CreateNtupleDColumn("z"); //3
  aMan->CreateNtupleDColumn("KEinitial"); //4
  aMan->CreateNtupleDColumn("KEescape"); //5
  aMan->CreateNtupleDColumn("px"); //6
  aMan->CreateNtupleDColumn("py"); //7
  aMan->CreateNtupleDColumn("pz"); //8
  aMan->CreateNtupleDColumn("ArHit"); //9
  aMan->CreateNtupleDColumn("pxi"); //10
  aMan->CreateNtupleDColumn("pyi"); //11
  aMan->CreateNtupleDColumn("pzi"); //12
  aMan->FinishNtuple();  
 */
  savehitsFile = "hits.out";
  savepmtFile  = "pmt.out";

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXRunAction::~DMXRunAction()
{
  delete runMessenger;
  runMessenger = 0;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunAction::BeginOfRunAction(const G4Run* aRun)
{
  
  // Logan ROOT
  /*auto aMan = G4AnalysisManager::Instance();
  aMan->OpenFile( "Data_ryan.root" );*/
  savehistFile = "/data/runzezhang/result/TN_box/dmx_Cfneutron_Ncry_1E6.root";
  G4int runID = aRun -> GetRunID();
  G4String filename1 = "/data/runzezhang/result/TN_box/dmx_Cfneutron_Ncry_1E6_" + std::to_string(runID) + ".root";
  savehistFile = filename1;
  //Master mode or sequential
  if (IsMaster())    
    G4cout << "### Run " << aRun->GetRunID() << " starts (master)." << G4endl;
  else
    G4cout << "### Run " << aRun->GetRunID() << " starts (worker)." << G4endl;
  
  // Book histograms and ntuples
  Book();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunAction::EndOfRunAction(const G4Run*)
{
  // Logan ROOT
  /*auto aMan = G4AnalysisManager::Instance();
  aMan->Write();
  aMan->CloseFile();*/
    auto man = G4AnalysisManager::Instance();
  man->Write();
  man->CloseFile();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXRunAction::Book()
{
// Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();

  man->SetVerboseLevel(2);
  man->SetDefaultFileType("root");

  // Open an output file
  man->OpenFile(savehistFile);
  man->SetFirstHistoId(1);
  man->SetFirstNtupleId(1);

   man->CreateNtuple("tree", "Scintillation Info");
  man->CreateNtupleDColumn("Event");
  man->CreateNtupleSColumn("name");
  man->CreateNtupleIColumn("Parent ID");
  man->CreateNtupleIColumn("Track ID");
  man->CreateNtupleIColumn("Step ID");
  man->CreateNtupleDColumn("X/mm");
  man->CreateNtupleDColumn("Y/mm");
  man->CreateNtupleDColumn("Z/mm");
  man->CreateNtupleDColumn("px/MeV");
  man->CreateNtupleDColumn("py/MeV");
  man->CreateNtupleDColumn("pz/MeV");
  man->CreateNtupleDColumn("Kinetic/keV");
  man->CreateNtupleDColumn("Recoiled/keV");
  man->CreateNtupleSColumn("Volume");
  man->CreateNtupleSColumn("Process");
  man->FinishNtuple();

  /*// Get/create analysis manager
  G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->SetDefaultFileType("root");
  
  // Open an output file
  man->OpenFile(savehistFile);
  man->SetFirstHistoId(1);
  man->SetFirstNtupleId(1);


  // ---- primary ntuple ------
  // id==1
  man->CreateNtuple("tree1", "Particle Source Energy");
  man->CreateNtupleDColumn("energy");
  man->FinishNtuple();

  // ---- secondary ntuple ------   
  //id==2
  man->CreateNtuple("tree2", "Scintillation Hits Info");
  man->CreateNtupleDColumn("Event");
  man->CreateNtupleDColumn("e_prim");
  man->CreateNtupleDColumn("tot_e");
  man->CreateNtupleDColumn("s_hits");
  man->CreateNtupleDColumn("xe_time");
  man->CreateNtupleDColumn("num_ph");
  man->CreateNtupleDColumn("avphtime");
  man->CreateNtupleDColumn("firstpart");
  man->CreateNtupleDColumn("firstparte");
  man->CreateNtupleDColumn("gamma");
  man->CreateNtupleDColumn("neutron");
  man->CreateNtupleDColumn("posi");
  man->CreateNtupleDColumn("elec");
  man->CreateNtupleDColumn("other");
  man->CreateNtupleDColumn("seed1");
  man->CreateNtupleDColumn("seed2");
  man->FinishNtuple();

  // ---- tertiary ntuple ------   
  //id==3
  man->CreateNtuple("tree3", "PMT Hits Info");
  man->CreateNtupleDColumn("event");
  man->CreateNtupleDColumn("hits");
  man->CreateNtupleDColumn("xpos");
  man->CreateNtupleDColumn("ypos");
  man->CreateNtupleDColumn("zpos");
  man->FinishNtuple();
 
  // Creating 1-dimensional histograms
  man->CreateH1("h1","Source Energy /keV",  1000,0.,10000.);
  man->CreateH1("h2","Energy Deposit /keV", 1000,0.,1000.);
  man->CreateH1("h3","Nuclear Recoil Edep /keV", 100,0.,100.);
  man->CreateH1("h4","Number of Photons - LowE", 200,0.,200.);
  man->CreateH1("h5","Number of Photons - HighE", 100,0.,10000.);
  man->CreateH1("h6","Average Photon Arrival/ns", 200,0.,200.);
  man->CreateH1("h7","1st event Photon Arrival", 200,0.,200.);
  man->CreateH1("h8","Gamma Energy Deposit/keV", 1000,0.,1000.);
  man->CreateH1("h9","Neutron Ener Deposit/keV", 1000,0.,1000.);
  man->CreateH1("h10","Electron Ener Deposit/keV",1000,0.,1000.);
  man->CreateH1("h11","Positron Ener Deposit/keV",1000,0.,1000.);
  man->CreateH1("h12","Other Ener Deposit/keV", 1000,0.,1000.);

  //Creating 2-dimensional histograms
  man->CreateH2("hh1","PMT Hit Pattern", 
		300 ,-30.,30.,300,-30.,30.);
  man->CreateH2("hh2","1st event PMT Hit Pattern", 
		300 ,-30.,30.,300,-30.,30.);*/

  return;

}


