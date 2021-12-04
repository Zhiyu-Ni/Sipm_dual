#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "DR_PMTSD.hh"
#include "DetectorConstruction.hh"
#include "TString.h"
#include "TRandom3.h"
//#include "TCint.h"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4SteppingManager.hh"
#include <time.h>

#include "G4EventManager.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"
#include "G4UnitsTable.hh"

#include <iostream>
#include <fstream>
#include <vector>
#include "TTree.h"

//long int CreateSeed();

using namespace std;
using namespace CLHEP;


int to_int(string name)
{
  int Result; // int which will contain the result
  stringstream convert(name);
  string dummy;
  convert >> dummy;
  convert >> Result;
  return Result;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::SteppingAction(DetectorConstruction *detectorConstruction,
                               const G4int &scint, const G4int &cher) : fDetectorConstruction(detectorConstruction),
                                                                        propagateScintillation(scint),
                                                                        propagateCerenkov(cher)
{
  maxtracklength = 500000. * mm;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

SteppingAction::~SteppingAction()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void SteppingAction::UserSteppingAction(const G4Step *theStep)
{

  G4Track *theTrack = theStep->GetTrack();

  //  const G4ThreeVector& theTrackDirection = theTrack->GetMomentumDirection();
  //  const G4ThreeVector& theTrackVertexDirection = theTrack->GetVertexMomentumDirection();

  //  TrackInformation* theTrackInfo = (TrackInformation*)(theTrack->GetUserInformation());

  G4ParticleDefinition *particleType = theTrack->GetDefinition();
  //G4int trackID = theTrack->GetTrackID();

  G4StepPoint *thePrePoint = theStep->GetPreStepPoint();
  G4StepPoint *thePostPoint = theStep->GetPostStepPoint();
  const G4ThreeVector &thePrePosition = thePrePoint->GetPosition();
  G4VPhysicalVolume *thePrePV = thePrePoint->GetPhysicalVolume();
  G4VPhysicalVolume *thePostPV = thePostPoint->GetPhysicalVolume();
  G4String thePrePVName = "";
  if (thePrePV)
    thePrePVName = thePrePV->GetName();
  G4String thePostPVName = "";
  if (thePostPV)
    thePostPVName = thePostPV->GetName();

  //  G4VSolid* thePreS = thePrePV->GetLogicalVolume()->GetSolid();

  G4int nStep = theTrack->GetCurrentStepNumber();

  G4int TrPDGid = theTrack->GetDefinition()->GetPDGEncoding();

  //        cout << " step length = " << theStep->GetStepLength() << endl;
  //-------------

  // get position
  G4double global_x = thePrePosition.x() / mm;
  G4double global_y = thePrePosition.y() / mm;
  G4double global_z = thePrePosition.z() / mm;

  //if(thePrePVName.contains("ecalDetP_rr")) cout<<"z "<<global_z<<endl;

  G4double energy = theStep->GetTotalEnergyDeposit();
  G4double energyIon = 0;

  double weight = 1.;
  double birkCut = 0.1;
  double charge = thePrePoint->GetCharge();
  if (charge != 0. && theStep->GetStepLength() > 0)
  {
    const G4Material *matbirk = thePrePoint->GetMaterial();
    double density = matbirk->GetDensity();
    double dedx = theStep->GetTotalEnergyDeposit() / theStep->GetStepLength();
    double rkb = 0.03333 * g / (MeV * cm2) / density;
    if (dedx > 0)
    {
      weight = 1. - 0.253694 * log(rkb * dedx);
      if (weight < birkCut)
        weight = birkCut;
      else if (weight > 1.)
        weight = 1.;
    }
  }
  bool use_birk = false;

  if (use_birk)
  {
    energyIon = (energy - theStep->GetNonIonizingEnergyDeposit()) * weight;
  }
  else
  {
    energyIon = (energy - theStep->GetNonIonizingEnergyDeposit());
  }

  if (((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3.) > 75)
    energyIon = 0;
  //if(energy>0) cout << "time "<<((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3.)<<endl;
  //if(((thePrePoint->GetGlobalTime() / ns) - global_z / 300. - 19 / 3.)>75) cout << " step length = " << theStep->GetStepLength()  << " energy = " << energy/MeV<<" denominator = "<<(1+0.0126*mm*energy/MeV/theStep->GetStepLength())<<" energyIon = "<<energyIon<< endl;

  G4double energyElec = 0.;
  //total energy by particle types
  G4double energyPion_n = 0.;
  G4double energyPositron = 0.;
  G4double energyElectron = 0.;
  G4double energyPhoton = 0.;
  G4double energyPion_p = 0.;
  G4double energyKion = 0.;
  G4double energyNeutron = 0.;
  G4double energyProton = 0.;
  //ion energy by particle types
  G4double energyIonPion_n = 0.;
  G4double energyIonPositron = 0.;
  G4double energyIonElectron = 0.;
  G4double energyIonPhoton = 0.;
  G4double energyIonPion_p = 0.;
  G4double energyIonKion = 0.;
  G4double energyIonNeutron = 0.;
  G4double energyIonProton = 0.;

 if (TrPDGid == (11))
  {
    energyElectron = energy;
    energyIonElectron = energyIon;
  }

  energyElec = energyIonPositron + energyIonElectron;

  //std::cout<<"TrPDGid energy energyIon enegyElec are "<<TrPDGid<<" "<<energy<<" "<<energyIon<<" "<<energyElec<<std::endl;



  // optical photon

  if (particleType == G4OpticalPhoton::OpticalPhotonDefinition())
  {
    //if optics
    G4String processName = theTrack->GetCreatorProcess()->GetProcessName();
    float photWL = MyMaterials::fromEvToNm(theTrack->GetTotalEnergy() / eV);
    //only consider 300 to 1000nm
    if (photWL > 1000 || photWL < 300)
      theTrack->SetTrackStatus(fKillTrackAndSecondaries);
    else
    {
      //only consider Cerenkov and Scintillation
      if ( (processName == "Cerenkov") )
      {
        /*
        if(thePostPVName.contains("ecalGapP_f") ){
           CreateTree::Instance()->h_phot_lambda_ECAL_f_collect_Ceren0->Fill(photWL);
        }
        if(thePostPVName.contains("ecalDetWindowP_f") ){
           CreateTree::Instance()->h_phot_lambda_ECAL_f_collect_Ceren1->Fill(photWL);
        }*/
        G4OpBoundaryProcessStatus boundaryStatus = Undefined;
        static G4OpBoundaryProcess* boundary     = NULL;
        if( !boundary ){
          G4ProcessManager* pm = theStep->GetTrack()->GetDefinition()->GetProcessManager();
          G4int nprocesses     = pm->GetProcessListLength();
          G4ProcessVector* pv  = pm->GetProcessList();
          G4int i;
          for( i = 0; i < nprocesses; i++ ){
            if( ( *pv )[i]->GetProcessName() == "OpBoundary" ){
              boundary = (G4OpBoundaryProcess*)( *pv )[i];
              break;
            }
          }
        }
        boundaryStatus = boundary->GetStatus();
        switch( boundaryStatus ){
        case Detection:
        {
          G4EventManager::GetEventManager()->KeepTheCurrentEvent();
          G4SDManager* SDman = G4SDManager::GetSDMpointer();
          DR_PMTSD* pmtSD
            = (DR_PMTSD*)SDman->FindSensitiveDetector( "/DR_Det" );
          if( pmtSD ){
            if(thePostPVName.contains("ecalDetP_f") ){
              // std::cout << "detected pre"<<thePrePVName<<" post " << thePostPVName << std::endl;
              //cout<<"get one. yes!"<<endl;
              if (processName == "Cerenkov"){
                CreateTree::Instance()->h_phot_lambda_ECAL_collect_Ceren->Fill(photWL);
                //CreateTree::Instance()->tot_phot_cer_should_det3 += 1;
                pmtSD->ProcessHits_constStep( theStep, NULL );
                theTrack->SetTrackStatus(fKillTrackAndSecondaries);
              }
              if (processName == "Scintillation"){
                CreateTree::Instance()->h_phot_lambda_ECAL_collect_Scin->Fill(photWL);
                //CreateTree::Instance()->tot_phot_cer_should_det4 += 1;
                pmtSD->ProcessHits_constStep( theStep, NULL );
                theTrack->SetTrackStatus(fKillTrackAndSecondaries);
              }
            }
          }
          break;
        }
        default:  
          break;
        }
      /*
        if(thePostPVName.contains("ecalDetP_fr") || thePostPVName.contains("ecalDetP_ff") ){
              CreateTree::Instance()->h_phot_lambda_ECAL_collect_Scin->Fill(photWL);
              theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
      */
        if (thePrePVName.contains("world"))
        {
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
        else if (thePrePVName.contains("ecalGapP"))
        {
          //std::cout << "in Gap " << thePrePVName << std::endl;
        }
        else if (thePrePVName.contains("wrap"))
        {
          //std::cout << "Crystal " << thePrePVName << std::endl;
        }
        else
        {
          //std::cout << "weird PrePVName " << thePrePVName << std::endl;
        }
        if (thePostPVName.contains("world")){
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }
          //std::cout << "out of world " << thePrePVName << std::endl;
        G4double tracklength = theStep->GetTrack()->GetTrackLength();
        if (tracklength > maxtracklength)
        {
          theStep->GetTrack()->SetTrackStatus(fStopAndKill);
          std::cout << "maximum " << thePrePVName << std::endl;
        }
        if (!propagateCerenkov && (processName == "Cerenkov"))
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);

        if (!propagateScintillation && (processName == "Scintillation"))
          theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        /*
        if ((nStep == 1))
        {
          if (thePrePVName.contains("ecalCrystalP_f_fiber_cherenc"))
          {
            //CreateTree::Instance()->tot_phot_cer_ECAL_cheren_f_total += 1;
            //CreateTree::Instance()->h_photon_2D_produce_Ceren_r->Fill(sqrt(pow(global_x,2)+ pow(global_y,2)));
            //CreateTree::Instance()->h_photon_2D_produce_Ceren_z->Fill(global_z);
            //CreateTree::Instance()->h_phot_lambda_ECAL_f_produce_Ceren->Fill(photWL);
          }
          TrackInformation *theTrackInfo = (TrackInformation *)(theTrack->GetUserInformation());
          G4int aapdgid = theTrackInfo->GetParentPDGid();
//          CreateTree::Instance()->h_phot_cer_lambda_ECAL_f->Fill(photWL);
          //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
        }*/
      }
      else
      {

        theTrack->SetTrackStatus(fKillTrackAndSecondaries);
      }
    }
  }
  
   // non optical photon

  return;
}
