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
// $Id: B1SteppingAction.cc 74483 2013-10-09 13:37:06Z gcosmo $
//
/// \file B1SteppingAction.cc
/// \brief Implementation of the B1SteppingAction class

#include "B1SteppingAction.hh"
#include "B1EventAction.hh"
#include "B1DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"

#include <fstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::B1SteppingAction(B1EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1SteppingAction::~B1SteppingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void B1SteppingAction::UserSteppingAction(const G4Step* step)
{
  if (!fScoringVolume) {
    const B1DetectorConstruction* detectorConstruction
      = static_cast<const B1DetectorConstruction*>
        (G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    fScoringVolume = detectorConstruction->GetScoringVolume();
  }

  // output and kill photons on creation
  G4Track* track = step->GetTrack();
  std::ofstream outfile;
  if (track->GetDefinition()->GetParticleName() == "gamma") {
    outfile.open("photons_energy.txt", std::ios_base::app);
    outfile << track->GetTotalEnergy() << std::endl;
    outfile.close();
    outfile.open("photons_p.txt", std::ios_base::app);
    outfile << track->GetMomentum().getX() << " "
            << track->GetMomentum().getY() << " "
            << track->GetMomentum().getZ() << std::endl;
    outfile.close();
    track->SetTrackStatus(fStopAndKill);
  }

  // get volume of the current step
  G4LogicalVolume* volume
    = step->GetPreStepPoint()->GetTouchableHandle()
      ->GetVolume()->GetLogicalVolume();

  // check if we are in scoring volume
  if (volume != fScoringVolume) return;
  outfile.open("energy_deposit.txt");
  outfile << step->GetPostStepPoint()->GetPosition().getX() << " "
          << step->GetTotalEnergyDeposit() << std::endl;
  outfile.close();

  // collect energy deposited in this step
  G4double edepStep = step->GetTotalEnergyDeposit();
  fEventAction->AddEdep(edepStep);

  // electron
  if (track->GetDefinition()->GetParticleName() == "e-") {
    // only output if the particle escapes in the right direction
    if (step->GetPostStepPoint()->GetPosition().getX() > 0.00966) {
      outfile.open("electrons_energy.txt", std::ios_base::app);
      outfile << track->GetTotalEnergy() << std::endl;
      outfile.close();
      outfile.open("electrons_p.txt", std::ios_base::app);
      outfile << track->GetMomentum().getX() << " "
              << track->GetMomentum().getY() << " "
              << track->GetMomentum().getZ() << std::endl;
      outfile.close();
    }
    track->SetTrackStatus(fStopAndKill);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
