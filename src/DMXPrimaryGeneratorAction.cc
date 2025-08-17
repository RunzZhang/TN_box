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
//
// PrimaryGeneratorAction program
// --------------------------------------------------------------

#include "DMXPrimaryGeneratorAction.hh"

#ifdef DMXENV_GPS_USE
#include "G4GeneralParticleSource.hh"
#else
#include "DMXParticleSource.hh"
#endif

#include "G4AnalysisManager.hh"
#include "G4Event.hh"

#include "Randomize.hh"

#include "globals.hh"

//	C/C++ includes
//
#include <vector>

//
//	GEANT4 includes
//
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4Poisson.hh"
G4int A = 252;
G4int Z = 98;
G4double multiplicity  = 3.75;
//------++++++------++++++------++++++------++++++------++++++------++++++------
//					BaccGeneratorCfFission()
//------++++++------++++++------++++++------++++++------++++++------++++++------
void DMXPrimaryGeneratorAction::BaccGeneratorCfFission()
{
    name = "CfFission";
    activityMultiplier = 1;
    neutronDef = G4Neutron::Definition();
    gammaDef = G4Gamma::Definition();

    averageGammaEnergy = -1.33 + 119.6 * pow(Z, 1. / 3) / A;

    static constexpr G4double gammaInfo[] = {
        0.141, 3.952,
        0.156, 5.091,
        0.171, 7.718,
        0.194, 8.609,
        0.251, 9.773,
        0.270, 8.303,
        0.301, 7.448,
        0.316, 7.722,
        0.339, 8.927,
        0.361, 10.700,
        0.377, 11.504,
        0.425, 9.035,
        0.448, 8.819,
        0.477, 8.736,
        0.499, 9.773,
        0.535, 8.799,
        0.556, 8.303,
        0.575, 8.455,
        0.594, 9.256,
        0.640, 7.722,
        0.667, 6.927,
        0.689, 6.927,
        0.697, 6.214,
        0.738, 5.998,
        0.776, 5.094,
        0.815, 4.171,
        0.853, 3.355,
        0.929, 2.798,
        0.968, 2.555,
        1.018, 2.305,
        1.079, 2.137,
        1.154, 2.026,
        1.235, 1.946,
        1.258, 1.757,
        1.350, 1.509,
        1.388, 1.378,
        1.426, 1.236,
        1.503, 1.149,
        1.579, 1.049,
        1.679, 0.950,
        1.743, 0.893,
        1.834, 0.855,
        1.880, 0.819,
        1.912, 0.761,
        1.970, 0.720,
        2.016, 0.689,
        2.133, 0.654,
        2.247, 0.620,
        2.267, 0.566,
        2.343, 0.566,
        2.400, 0.526,
        2.446, 0.474,
        2.611, 0.455,
        2.690, 0.383,
        2.840, 0.347,
        2.916, 0.317,
        3.018, 0.269,
        3.117, 0.230,
        3.222, 0.212,
        3.260, 0.190,
        3.356, 0.174,
        3.467, 0.159,
        3.520, 0.147,
        3.649, 0.123,
        3.785, 0.113,
        3.876, 0.103,
        3.966, 0.095,
        4.090, 0.085,
        4.178, 0.078,
        4.284, 0.070,
        4.420, 0.062,
        4.601, 0.050,
        4.738, 0.040,
        4.919, 0.034,
        5.100, 0.030,
        5.237, 0.024,
        5.373, 0.020,
        5.554, 0.017,
        5.736, 0.014,
        5.872, 0.011,
        6.053, 0.009,
        6.235, 0.007,
        6.507, 0.004,
        6.734, 0.003,
        7.051, 0.002,
        7.414, 4.522e-4
    };

    G4int numPoints = sizeof(gammaInfo) / sizeof(G4double) / 2;
    fCDFSize = numPoints;
    G4double* gammaPDF, totalArea = 0.;
    gammaPDF = new G4double[numPoints];
    for (G4int i = 0; i < numPoints; i++)
    {
        gammaEnergy[i] = gammaInfo[i * 2];
        gammaPDF[i] = gammaInfo[i * 2 + 1];
        totalArea += gammaPDF[i];
    }

    for (G4int i = 0; i < numPoints; i++)
    {
        gammaPDF[i] /= totalArea;
    }

    gammaCDF[0] = 0.;
    for (G4int i = 1; i < numPoints; i++)
    {
        gammaCDF[i] = gammaCDF[i - 1] + gammaPDF[i - 1];
    }
    gammaCDF[numPoints - 1] = 1.;

    delete[] gammaPDF;
}

DMXPrimaryGeneratorAction::DMXPrimaryGeneratorAction() {
  
#ifdef DMXENV_GPS_USE
  particleGun = new G4GeneralParticleSource();
#else
  particleGun = new DMXParticleSource();
#endif

  energy_pri=0;
  //  seeds=NULL;
  seeds[0] =-1;
  seeds[1] =-1;

  BaccGeneratorCfFission();

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

DMXPrimaryGeneratorAction::~DMXPrimaryGeneratorAction() {

  delete particleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo....

void DMXPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  energy_pri = 0.;

  // seeds
  seeds[0] = *G4Random::getTheSeeds();
  seeds[1] = *(G4Random::getTheSeeds()+1);

  //particleGun->GeneratePrimaryVertex(anEvent);

  //energy_pri = particleGun->GetParticleEnergy();

  //Fill ntuple #1
  /*G4AnalysisManager* man = G4AnalysisManager::Instance();
  man->FillNtupleDColumn(1,0,energy_pri);
  man->AddNtupleRow(1);*/
    /*
  //	Then create the gammas
  particleGun->SetParticleDefinition(gammaDef);

  //	The total energy of the gammas is 6.95 +/- 0.3 MeV. I take this to mean
  //	there is at least a 4% uncertainty in the total gamma energy. This
  //	generator therefore smears out the total energy by 4% so that we don't
  //	get a well-defined total gamma energy for any given multiplicity.
  G4double totalEnergy = (2.51 - 1.13e-5 * Z * Z * sqrt(A)) * multiplicity + 4.;
  totalEnergy = G4RandGauss::shoot(totalEnergy, totalEnergy * 0.3 / 6.95);

  G4double numGammas1 = totalEnergy / averageGammaEnergy;
  G4int numGammas2 = G4Poisson(numGammas1);
  G4double runningTotal = 0;
  std::vector<G4double> gammaEnergies;
  for (G4int i = 0; i < numGammas2; i++)
  {
      particleGun->GetCurrentSource()->GetAngDist()->SetParticleMomentumDirection(GetRandomDirection());

      G4double singleEnergy;
      if (i == numGammas2 - 1)
      {
          singleEnergy = totalEnergy;
      }
      else
      {
          singleEnergy = GetGammaEnergy(totalEnergy);
      }

      particleGun->GetCurrentSource()->GetEneDist()->SetMonoEnergy(
          singleEnergy);
      particleGun->GeneratePrimaryVertex(anEvent);

      //baccManager->AddPrimaryParticle(GetParticleInfo(particleGun));

      runningTotal += singleEnergy;
      totalEnergy -= singleEnergy;

      if (totalEnergy == 0)
      {
          break;
      }
  }*/
}


//------++++++------++++++------++++++------++++++------++++++------++++++------
//					GetGammaEnergy()
//------++++++------++++++------++++++------++++++------++++++------++++++------
G4double DMXPrimaryGeneratorAction::GetGammaEnergy(G4double limit)
{
    G4double energy = limit + 1;

    if (limit < gammaEnergy[0])
    {
        energy = limit;
    }
    else
    {
        while (energy > limit)
        {
            G4double prob = G4UniformRand();
            G4int indexLo = 0, indexHi = fCDFSize - 1;

            while ((gammaCDF[indexLo + 1] <= prob) || (gammaCDF[indexHi - 1] >= prob))
            {
                if (gammaCDF[(indexLo + indexHi) / 2] < prob)
                {
                    indexLo = (indexLo + indexHi) / 2;
                }
                else
                {
                    indexHi = (indexLo + indexHi) / 2;
                }
            }

            G4double split = (prob - gammaCDF[indexLo]) / (gammaCDF[indexHi] - gammaCDF[indexLo]);
            energy = gammaEnergy[indexLo] + split * (gammaEnergy[indexHi] - gammaEnergy[indexLo]);
        }
    }

    return (energy);
}

//------++++++------++++++------++++++------++++++------++++++------++++++------
//					GetGammaEnergy()
//------++++++------++++++------++++++------++++++------++++++------++++++------
G4ThreeVector DMXPrimaryGeneratorAction::GetRandomDirection() {
    // isotropic direction
    G4double cosTheta = 2.0*G4UniformRand() - 1.0; // -1 to 1
    G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
    G4double phi      = 2.0 * CLHEP::pi * G4UniformRand();
    G4double x = sinTheta * std::cos(phi);
    G4double y = sinTheta * std::sin(phi);
    G4double z = cosTheta;
    return G4ThreeVector(x, y, z);
}
