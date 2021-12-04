#include "CreateTree.hh"
#include <algorithm>

using namespace std;

CreateTree *CreateTree::fInstance = NULL;

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

CreateTree::CreateTree(TString name)
{
  if (fInstance)
  {
    return;
  }

  this->fInstance = this;
  this->fname = name;
  this->ftree = new TTree(name, name);

  this->GetTree()->Branch("Event", &this->Event, "Event/I");

  this->GetTree()->Branch("inputTrackerX0", &this->inputTrackerX0, "inputTrackerX0/F");
  this->GetTree()->Branch("inputServiceAlmm", &this->inputServiceAlmm, "inputServiceAlmm/F");
  this->GetTree()->Branch("inputTimingThick", &this->inputTimingThick, "inputTimingThick/F");
  this->GetTree()->Branch("inputE1Thick", &this->inputE1Thick, "inputE1Thick/F");
  this->GetTree()->Branch("inputE2Thick", &this->inputE2Thick, "inputE2Thick/F");
  this->GetTree()->Branch("inputE1Width", &this->inputE1Width, "inputE1Width/F");
  this->GetTree()->Branch("inputTimingECAL_dist", &this->inputTimingECAL_dist, "inputTimingECAL_dist/F");

  inputInitialPosition = new vector<float>(3, 0.);
  inputMomentum = new vector<float>(4, 0.);
  primaryPosT1 = new vector<float>(3, 0.);
  primaryMomT1 = new vector<float>(4, 0.);
  primaryPosE1 = new vector<float>(3, 0.);
  primaryMomE1 = new vector<float>(4, 0.);

  this->GetTree()->Branch("inputInitialPosition", "vector<float>", &inputInitialPosition);
  this->GetTree()->Branch("inputMomentum", "vector<float>", &inputMomentum);
  this->GetTree()->Branch("primaryPosT1", "vector<float>", &primaryPosT1);
  this->GetTree()->Branch("primaryMomT1", "vector<float>", &primaryMomT1);
  this->GetTree()->Branch("primaryPosE1", "vector<float>", &primaryPosE1);
  this->GetTree()->Branch("primaryMomE1", "vector<float>", &primaryMomE1);

  //h_parentID_Edep          = new vector<float>(11,0.);
  //h_parentID_Edep[11] = 0;
  //this -> GetTree() -> Branch("h_parentID_Edep",   "map<int,double>", &h_parentID_Edep);

  this->GetTree()->Branch("nTracksT1", &this->nTracksT1, "nTracksT1/I");
  this->GetTree()->Branch("nTracksT2", &this->nTracksT2, "nTracksT2/I");
  this->GetTree()->Branch("nTracksE1", &this->nTracksE1, "nTracksE1/I");
  this->GetTree()->Branch("nTracksE2", &this->nTracksE2, "nTracksE2/I");
  this->GetTree()->Branch("nTracksTRK", &this->nTracksTRK, "nTracksTRK[6]/F");

  //integrated per longitudinal layer
  this->GetTree()->Branch("depositedEnergyTotal", &this->depositedEnergyTotal, "depositedEnergyTotal/F");
  this->GetTree()->Branch("depositedEnergyEscapeWorld", &this->depositedEnergyEscapeWorld, "depositedEnergyEscapeWorld/F");
  this->GetTree()->Branch("depositedEnergyTiming_f", &this->depositedEnergyTiming_f, "depositedEnergyTiming_f/F");
  this->GetTree()->Branch("depositedEnergyTiming_r", &this->depositedEnergyTiming_r, "depositedEnergyTiming_r/F");
  this->GetTree()->Branch("depositedEnergyECAL_f", &this->depositedEnergyECAL_f, "depositedEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedEnergyECAL_r", &this->depositedEnergyECAL_r, "depositedEnergyECAL_r[3]/F");
  this->GetTree()->Branch("depositedEnergyHCALAct", &this->depositedEnergyHCALAct, "depositedEnergyHCALAct/F");
  this->GetTree()->Branch("depositedEnergyHCALPas", &this->depositedEnergyHCALPas, "depositedEnergyHCALPas/F");
  this->GetTree()->Branch("depositedEnergyWorld", &this->depositedEnergyWorld, "depositedEnergyWorld/F");
  this->GetTree()->Branch("depositedEnergyTimingGap", &this->depositedEnergyTimingGap, "depositedEnergyTimingGap/F");
  this->GetTree()->Branch("depositedEnergyServices", &this->depositedEnergyServices, "depositedEnergyServices/F");
  this->GetTree()->Branch("depositedEnergyEcalGap", &this->depositedEnergyEcalGap, "depositedEnergyEcalGap/F");
  this->GetTree()->Branch("depositedEnergyEcalDet", &this->depositedEnergyEcalDet, "depositedEnergyEcalDet/F");
  this->GetTree()->Branch("depositedEnergySolenoid", &this->depositedEnergySolenoid, "depositedEnergySolenoid/F");

  this->GetTree()->Branch("depositedEnergyECAL_absorb_f_particleID", &this->depositedEnergyECAL_absorb_f_particleID, "depositedEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_absorb_r_particleID", &this->depositedEnergyECAL_absorb_r_particleID, "depositedEnergyECAL_absorb_r_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_scinti_f_particleID", &this->depositedEnergyECAL_scinti_f_particleID, "depositedEnergyECAL_scinti_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_scinti_r_particleID", &this->depositedEnergyECAL_scinti_r_particleID, "depositedEnergyECAL_scinti_r_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_cheren_f_particleID", &this->depositedEnergyECAL_cheren_f_particleID, "depositedEnergyECAL_cheren_f_particleID[8]/F");
  this->GetTree()->Branch("depositedEnergyECAL_cheren_r_particleID", &this->depositedEnergyECAL_cheren_r_particleID, "depositedEnergyECAL_cheren_r_particleID[8]/F");

  this->GetTree()->Branch("depositedIonEnergyTotal", &this->depositedIonEnergyTotal, "depositedIonEnergyTotal/F");
  this->GetTree()->Branch("depositedIonEnergyTiming_f", &this->depositedIonEnergyTiming_f, "depositedIonEnergyTiming_f/F");
  this->GetTree()->Branch("depositedIonEnergyTiming_r", &this->depositedIonEnergyTiming_r, "depositedIonEnergyTiming_r/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_f", &this->depositedIonEnergyECAL_f, "depositedIonEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_r", &this->depositedIonEnergyECAL_r, "depositedIonEnergyECAL_r[3]/F");
  this->GetTree()->Branch("depositedIonEnergyHCALAct", &this->depositedIonEnergyHCALAct, "depositedIonEnergyHCALAct/F");
  this->GetTree()->Branch("depositedIonEnergyHCALPas", &this->depositedIonEnergyHCALPas, "depositedIonEnergyHCALPas/F");
  this->GetTree()->Branch("depositedIonEnergyWorld", &this->depositedIonEnergyWorld, "depositedIonEnergyWorld/F");
  this->GetTree()->Branch("depositedIonEnergyTimingGap", &this->depositedIonEnergyTimingGap, "depositedIonEnergyTimingGap/F");
  this->GetTree()->Branch("depositedIonEnergyServices", &this->depositedIonEnergyServices, "depositedIonEnergyServices/F");
  this->GetTree()->Branch("depositedIonEnergyEcalGap", &this->depositedIonEnergyEcalGap, "depositedIonEnergyEcalGap/F");
  this->GetTree()->Branch("depositedIonEnergyEcalDet", &this->depositedIonEnergyEcalDet, "depositedIonEnergyEcalDet/F");
  this->GetTree()->Branch("depositedIonEnergySolenoid", &this->depositedIonEnergySolenoid, "depositedIonEnergySolenoid/F");

  this->GetTree()->Branch("depositedIonEnergyECAL_absorb_f_particleID", &this->depositedIonEnergyECAL_absorb_f_particleID, "depositedIonEnergyECAL_absorb_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_absorb_r_particleID", &this->depositedIonEnergyECAL_absorb_r_particleID, "depositedIonEnergyECAL_absorb_r_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_scinti_f_particleID", &this->depositedIonEnergyECAL_scinti_f_particleID, "depositedIonEnergyECAL_scinti_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_scinti_r_particleID", &this->depositedIonEnergyECAL_scinti_r_particleID, "depositedIonEnergyECAL_scinti_r_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_cheren_f_particleID", &this->depositedIonEnergyECAL_cheren_f_particleID, "depositedIonEnergyECAL_cheren_f_particleID[8]/F");
  this->GetTree()->Branch("depositedIonEnergyECAL_cheren_r_particleID", &this->depositedIonEnergyECAL_cheren_r_particleID, "depositedIonEnergyECAL_cheren_r_particleID[8]/F");

  this->GetTree()->Branch("depositedElecEnergyTotal", &this->depositedElecEnergyTotal, "depositedElecEnergyTotal/F");
  this->GetTree()->Branch("depositedElecEnergyTiming_f", &this->depositedElecEnergyTiming_f, "depositedElecEnergyTiming_f/F");
  this->GetTree()->Branch("depositedElecEnergyTiming_r", &this->depositedElecEnergyTiming_r, "depositedElecEnergyTiming_r/F");
  this->GetTree()->Branch("depositedElecEnergyECAL_f", &this->depositedElecEnergyECAL_f, "depositedElecEnergyECAL_f[3]/F");
  this->GetTree()->Branch("depositedElecEnergyECAL_r", &this->depositedElecEnergyECAL_r, "depositedElecEnergyECAL_r[3]/F");
  this->GetTree()->Branch("depositedElecEnergyHCALAct", &this->depositedElecEnergyHCALAct, "depositedElecEnergyHCALAct/F");
  this->GetTree()->Branch("depositedElecEnergyHCALPas", &this->depositedElecEnergyHCALPas, "depositedElecEnergyHCALPas/F");
  this->GetTree()->Branch("depositedElecEnergyWorld", &this->depositedElecEnergyWorld, "depositedElecEnergyWorld/F");
  this->GetTree()->Branch("depositedElecEnergyTimingGap", &this->depositedElecEnergyTimingGap, "depositedElecEnergyTimingGap/F");
  this->GetTree()->Branch("depositedElecEnergyServices", &this->depositedElecEnergyServices, "depositedElecEnergyServices/F");
  this->GetTree()->Branch("depositedElecEnergyEcalGap", &this->depositedElecEnergyEcalGap, "depositedElecEnergyEcalGap/F");
  this->GetTree()->Branch("depositedElecEnergyEcalDet", &this->depositedElecEnergyEcalDet, "depositedElecEnergyEcalDet/F");
  this->GetTree()->Branch("depositedElecEnergySolenoid", &this->depositedElecEnergySolenoid, "depositedElecEnergySolenoid/F");

  //single channels
  this->GetTree()->Branch("Edep_Tracker_layer", &this->Edep_Tracker_layer, "Edep_Tracker_layer[6]/F");
  this->GetTree()->Branch("Edep_Timing_f_ch", &this->Edep_Timing_f_ch, "Edep_Timing_f_ch[18]/F");
  this->GetTree()->Branch("Edep_Timing_r_ch", &this->Edep_Timing_r_ch, "Edep_Timing_r_ch[18]/F");
  this->GetTree()->Branch("Edep_ECAL_f_ch", &this->Edep_ECAL_f_ch, "Edep_ECAL_f_ch[6400]/F");
  this->GetTree()->Branch("Edep_ECAL_r_ch", &this->Edep_ECAL_r_ch, "Edep_ECAL_r_ch[6400]/F");

  this->GetTree()->Branch("IonEdep_ECAL_f_ch", &this->IonEdep_ECAL_f_ch, "IonEdep_ECAL_f_ch[6400]/F");
  this->GetTree()->Branch("IonEdep_ECAL_r_ch", &this->IonEdep_ECAL_r_ch, "IonEdep_ECAL_r_ch[6400]/F");
//  this->GetTree()->Branch("ElecEdep_ECAL_f_ch", &this->ElecEdep_ECAL_f_ch, "ElecEdep_ECAL_f_ch[6400]/F");
//  this->GetTree()->Branch("ElecEdep_ECAL_r_ch", &this->ElecEdep_ECAL_r_ch, "ElecEdep_ECAL_r_ch[6400]/F");

//  this->GetTree()->Branch("Edep_ECAL_f_scinti_ch", &this->Edep_ECAL_f_scinti_ch, "Edep_ECAL_f_scinti_ch[6400]/F");
//  this->GetTree()->Branch("Edep_ECAL_r_scinti_ch", &this->Edep_ECAL_r_scinti_ch, "Edep_ECAL_r_scinti_ch[6400]/F");
//  this->GetTree()->Branch("IonEdep_ECAL_f_scinti_ch", &this->IonEdep_ECAL_f_scinti_ch, "IonEdep_ECAL_f_scinti_ch[6400]/F");
//  this->GetTree()->Branch("IonEdep_ECAL_r_scinti_ch", &this->IonEdep_ECAL_r_scinti_ch, "IonEdep_ECAL_r_scinti_ch[6400]/F");
//  this->GetTree()->Branch("ElecEdep_ECAL_f_scinti_ch", &this->ElecEdep_ECAL_f_scinti_ch, "ElecEdep_ECAL_f_scinti_ch[6400]/F");
//  this->GetTree()->Branch("ElecEdep_ECAL_r_scinti_ch", &this->ElecEdep_ECAL_r_scinti_ch, "ElecEdep_ECAL_r_scinti_ch[6400]/F");

//  this->GetTree()->Branch("Edep_ECAL_f_cheren_ch", &this->Edep_ECAL_f_cheren_ch, "Edep_ECAL_f_cheren_ch[6400]/F");
//  this->GetTree()->Branch("Edep_ECAL_r_cheren_ch", &this->Edep_ECAL_r_cheren_ch, "Edep_ECAL_r_cheren_ch[6400]/F");
//  this->GetTree()->Branch("IonEdep_ECAL_f_cheren_ch", &this->IonEdep_ECAL_f_cheren_ch, "IonEdep_ECAL_f_cheren_ch[6400]/F");
//  this->GetTree()->Branch("IonEdep_ECAL_r_cheren_ch", &this->IonEdep_ECAL_r_cheren_ch, "IonEdep_ECAL_r_cheren_ch[6400]/F");
//  this->GetTree()->Branch("ElecEdep_ECAL_f_cheren_ch", &this->ElecEdep_ECAL_f_cheren_ch, "ElecEdep_ECAL_f_cheren_ch[6400]/F");
//  this->GetTree()->Branch("ElecEdep_ECAL_r_cheren_ch", &this->ElecEdep_ECAL_r_cheren_ch, "ElecEdep_ECAL_r_cheren_ch[6400]/F");

  //Cerenkov photons
  //  this -> GetTree() -> Branch("tot_phot_sci_Timing",        &this->tot_phot_sci_Timing,               "tot_phot_sci_Timing/I");
  this->GetTree()->Branch("tot_phot_cer_Timing_f_total", &this->tot_phot_cer_Timing_f_total, "tot_phot_cer_Timing_f_total/I");
  this->GetTree()->Branch("tot_phot_cer_Timing_r_total", &this->tot_phot_cer_Timing_r_total, "tot_phot_cer_Timing_r_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_f_total", &this->tot_phot_cer_ECAL_scinti_f_total, "tot_phot_cer_ECAL_scinti_f_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_r_total", &this->tot_phot_cer_ECAL_scinti_r_total, "tot_phot_cer_ECAL_scinti_r_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_f_particleID", &this->tot_phot_cer_ECAL_scinti_f_particleID, "tot_phot_cer_ECAL_scinti_f_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_scinti_r_particleID", &this->tot_phot_cer_ECAL_scinti_r_particleID, "tot_phot_cer_ECAL_scinti_r_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_f_total", &this->tot_phot_cer_ECAL_cheren_f_total, "tot_phot_cer_ECAL_cheren_f_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_r_total", &this->tot_phot_cer_ECAL_cheren_r_total, "tot_phot_cer_ECAL_cheren_r_total/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_f_particleID", &this->tot_phot_cer_ECAL_cheren_f_particleID, "tot_phot_cer_ECAL_cheren_f_particleID[8]/I");
  this->GetTree()->Branch("tot_phot_cer_ECAL_cheren_r_particleID", &this->tot_phot_cer_ECAL_cheren_r_particleID, "tot_phot_cer_ECAL_cheren_r_particleID[8]/I");

  this->GetTree()->Branch("tot_phot_cer_HCAL", &this->tot_phot_cer_HCAL, "tot_phot_cer_HCAL/I");
  this->GetTree()->Branch("tot_phot_cer_SDdetected_ff", &this->tot_phot_cer_SDdetected_ff, "tot_phot_cer_SDdetected_ff/I");
  this->GetTree()->Branch("tot_phot_cer_SDdetected_fr", &this->tot_phot_cer_SDdetected_fr, "tot_phot_cer_SDdetected_fr/I");
  this->GetTree()->Branch("tot_phot_cer_SDdetected_rf", &this->tot_phot_cer_SDdetected_rf, "tot_phot_cer_SDdetected_rf/I");
  this->GetTree()->Branch("tot_phot_cer_SDdetected_rr", &this->tot_phot_cer_SDdetected_rr, "tot_phot_cer_SDdetected_rr/I");
  this->GetTree()->Branch("tot_phot_cer_should_det1", &this->tot_phot_cer_should_det1, "tot_phot_cer_should_det1/I");
  this->GetTree()->Branch("tot_phot_cer_should_det2", &this->tot_phot_cer_should_det2, "tot_phot_cer_should_det2/I");
  this->GetTree()->Branch("tot_phot_cer_should_det3", &this->tot_phot_cer_should_det3, "tot_phot_cer_should_det3/I");
  this->GetTree()->Branch("tot_phot_cer_should_det4", &this->tot_phot_cer_should_det4, "tot_phot_cer_should_det4/I");
  this->GetTree()->Branch("tot_phot_sin_should_det4", &this->tot_phot_sin_should_det4, "tot_phot_sin_should_det4/I");

  this->GetTree()->Branch("E_Zdep_0to5000mm_total", &this->E_Zdep_0to5000mm_total, "E_Zdep_0to5000mm_total[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Pion_n", &this->E_Zdep_0to5000mm_Pion_n, "E_Zdep_0to5000mm_Pion_n[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Positron", &this->E_Zdep_0to5000mm_Positron, "E_Zdep_0to5000mm_Positron[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Electron", &this->E_Zdep_0to5000mm_Electron, "E_Zdep_0to5000mm_Electron[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Photon", &this->E_Zdep_0to5000mm_Photon, "E_Zdep_0to5000mm_Photon[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Pion_p", &this->E_Zdep_0to5000mm_Pion_p, "E_Zdep_0to5000mm_Pion_p[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Kion", &this->E_Zdep_0to5000mm_Kion, "E_Zdep_0to5000mm_Kion[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Neutron", &this->E_Zdep_0to5000mm_Neutron, "E_Zdep_0to5000mm_Neutron[2500]/F");
  this->GetTree()->Branch("E_Zdep_0to5000mm_Proton", &this->E_Zdep_0to5000mm_Proton, "E_Zdep_0to5000mm_Proton[2500]/F");




  this->GetTree()->Branch("E_Tdep_0to5ns_total", &this->E_Tdep_0to5ns_total, "E_Tdep_0to5ns_total[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Pion_n", &this->E_Tdep_0to5ns_Pion_n, "E_Tdep_0to5ns_Pion_n[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Positron", &this->E_Tdep_0to5ns_Positron, "E_Tdep_0to5ns_Positron[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Electron", &this->E_Tdep_0to5ns_Electron, "E_Tdep_0to5ns_Electron[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Photon", &this->E_Tdep_0to5ns_Photon, "E_Tdep_0to5ns_Photon[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Pion_p", &this->E_Tdep_0to5ns_Pion_p, "E_Tdep_0to5ns_Pion_p[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Kion", &this->E_Tdep_0to5ns_Kion, "E_Tdep_0to5ns_Kion[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Neutron", &this->E_Tdep_0to5ns_Neutron, "E_Tdep_0to5ns_Neutron[2500]/F");
  this->GetTree()->Branch("E_Tdep_0to5ns_Proton", &this->E_Tdep_0to5ns_Proton, "E_Tdep_0to5ns_Proton[2500]/F");


/*
  h_time_z_st = new TH2F("h_time_z_st", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_time_z_mt = new TH2F("h_time_z_mt", "", 2501, 0., 5002., 5010, 0., 100.2);
  h_time_z_lt = new TH2F("h_time_z_lt", "", 2501, 0., 5002., 5010, 0., 5010.);

  h_localtime_z_st = new TH2F("h_localtime_z_st", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt = new TH2F("h_localtime_z_mt", "", 2501, 0., 5002., 5010, 0., 100.2);
  h_localtime_z_lt = new TH2F("h_localtime_z_lt", "", 2501, 0., 5002., 5010, 0., 5010.);

  h_localtime_z_st_en = new TH2F("h_localtime_z_st_en", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_en = new TH2F("h_localtime_z_mt_en", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_en = new TH2F("h_localtime_z_lt_en", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_localtime_z_st_ep = new TH2F("h_localtime_z_st_ep", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_ep = new TH2F("h_localtime_z_mt_ep", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_ep = new TH2F("h_localtime_z_lt_ep", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_localtime_z_st_gamma = new TH2F("h_localtime_z_st_gamma", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_gamma = new TH2F("h_localtime_z_mt_gamma", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_gamma = new TH2F("h_localtime_z_lt_gamma", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_localtime_z_st_n = new TH2F("h_localtime_z_st_n", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_n = new TH2F("h_localtime_z_mt_n", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_n = new TH2F("h_localtime_z_lt_n", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_localtime_z_st_p = new TH2F("h_localtime_z_st_p", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_p = new TH2F("h_localtime_z_mt_p", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_p = new TH2F("h_localtime_z_lt_p", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_localtime_z_st_pin = new TH2F("h_localtime_z_st_pin", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_pin = new TH2F("h_localtime_z_mt_pin", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_pin = new TH2F("h_localtime_z_lt_pin", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_localtime_z_st_pip = new TH2F("h_localtime_z_st_pip", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_mt_pip = new TH2F("h_localtime_z_mt_pip", "", 2501, 0., 5002., 5010, 0., 5.01);
  h_localtime_z_lt_pip = new TH2F("h_localtime_z_lt_pip", "", 2501, 0., 5002., 5010, 0., 5.01);

  h_parentID_Edep = new TH1F("h_parentID_Edep", "", 10000, -5000, 5000); //600,-300,300);

  h_phot_cer_lambda_Timing_f = new TH1F("h_phot_cer_lambda_Timing_f", "", 1250, 0., 1250.);
  h_phot_cer_lambda_Timing_r = new TH1F("h_phot_cer_lambda_Timing_r", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_collect_Ceren = new TH1F("h_phot_lambda_ECAL_f_collect_Ceren", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_produce_Ceren = new TH1F("h_phot_lambda_ECAL_f_produce_Ceren", "", 1250, 0., 1250.);
  h_phot_cer_lambda_HCAL = new TH1F("h_phot_cer_lambda_HCAL", "", 1250, 0., 1250.);
  h_phot_cer_parentID = new TH1F("h_phot_cer_parentID", "", 600, -300, 300);
*/

  h_phot_lambda_ECAL_f_collect_Ceren0 = new TH1F("h_phot_lambda_ECAL_f_collect_Ceren0", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_collect_Ceren = new TH1F("h_phot_lambda_ECAL_f_collect_Ceren", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_collect_Ceren1 = new TH1F("h_phot_lambda_ECAL_f_collect_Ceren1", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_collect_Ceren2 = new TH1F("h_phot_lambda_ECAL_f_collect_Ceren2", "", 1250, 0., 1250.);
  h_phot_lambda_ECAL_f_produce_Ceren = new TH1F("h_phot_lambda_ECAL_f_produce_Ceren", "", 1250, 0., 1250.);
  h_photon_2D_produce_Ceren_r = new TH1F("h_photon_2D_produce_Ceren_r", "", 500, 0., 15);
  h_photon_2D_receive_Ceren_r = new TH1F("h_photon_2D_receive_Ceren_r", "", 500, 0., 15);
  h_photon_2D_produce_Ceren_z = new TH1F("h_photon_2D_produce_Ceren_z", "", 500, 0., 200);
  h_photon_2D_receive_Ceren_z = new TH1F("h_photon_2D_receive_Ceren_z", "", 500, 0., 200);

 h_photon_2D_receive_Ceren = new TH2F("h_photon_2D_receive_Ceren", "", 500, -50, 50, 500, 0., 200);
  h_photon_2D_produce_Ceren = new TH2F("h_photon_2D_produce_Ceren", "", 500, -50, 50, 500, 0., 200);

    h_phot_lambda_ECAL_collect_Ceren = new TH1F("h_phot_lambda_ECAL_collect_Ceren", "", 1250, 0., 1250.);
    h_phot_lambda_ECAL_collect_Scin = new TH1F("h_phot_lambda_ECAL_collect_Scin", "", 1250, 0., 1250.);
    

  
  this->Clear();
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

CreateTree::~CreateTree()
{
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

int CreateTree::Fill()
{
//  this->GetTree()->Write(NULL, TObject::kOverwrite );
  return this->GetTree()->Fill();
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

bool CreateTree::Write(TFile *outfile)
{
  outfile->cd();
  ftree->Write();
  h_phot_lambda_ECAL_collect_Ceren->Write();
  h_phot_lambda_ECAL_collect_Scin->Write();
  
  h_phot_lambda_ECAL_f_collect_Ceren0->Write();
  h_phot_lambda_ECAL_f_collect_Ceren->Write();
  h_phot_lambda_ECAL_f_collect_Ceren1->Write();
  h_phot_lambda_ECAL_f_collect_Ceren2->Write();
  
  h_phot_lambda_ECAL_f_produce_Ceren->Write();
  h_photon_2D_produce_Ceren_r->Write();
  h_photon_2D_produce_Ceren_z->Write();
  h_photon_2D_receive_Ceren_r->Write();
  h_photon_2D_receive_Ceren_z->Write();
  
    h_photon_2D_produce_Ceren->Write();
  h_photon_2D_receive_Ceren->Write();
/*
  h_phot_cer_lambda_Timing_f->Write();
  h_phot_cer_lambda_Timing_r->Write();
  h_phot_lambda_ECAL_f_collect_Ceren->Write();
  h_phot_lambda_ECAL_f_produce_Ceren->Write();
  h_phot_cer_lambda_HCAL->Write();
  h_phot_cer_parentID->Write();

  h_time_z_st->Write();
  h_time_z_mt->Write();
  h_time_z_lt->Write();

  h_localtime_z_st->Write();
  h_localtime_z_mt->Write();
  h_localtime_z_lt->Write();

  h_localtime_z_st_en->Write();
  h_localtime_z_mt_en->Write();
  h_localtime_z_lt_en->Write();

  h_localtime_z_st_ep->Write();
  h_localtime_z_mt_ep->Write();
  h_localtime_z_lt_ep->Write();

  h_localtime_z_st_gamma->Write();
  h_localtime_z_mt_gamma->Write();
  h_localtime_z_lt_gamma->Write();

  h_localtime_z_st_n->Write();
  h_localtime_z_mt_n->Write();
  h_localtime_z_lt_n->Write();

  h_localtime_z_st_p->Write();
  h_localtime_z_mt_p->Write();
  h_localtime_z_lt_p->Write();

  h_localtime_z_st_pin->Write();
  h_localtime_z_mt_pin->Write();
  h_localtime_z_lt_pin->Write();

  h_localtime_z_st_pip->Write();
  h_localtime_z_mt_pip->Write();
  h_localtime_z_lt_pip->Write();

  h_parentID_Edep->Write();
*/
  return true;
}

// ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void CreateTree::Clear()
{
  Event = 0;

  nTracksT1 = 0;
  nTracksT2 = 0;
  nTracksE1 = 0;
  nTracksE2 = 0;

  for (int iLayer = 0; iLayer < 6; iLayer++)
  {
    nTracksTRK[iLayer] = 0;
  }

  depositedEnergyEscapeWorld = 0.;

  depositedEnergyTotal = 0.;
  depositedEnergyTiming_f = 0.;
  depositedEnergyTiming_r = 0.;
  for (int i = 0; i < 3; i++){
    depositedEnergyECAL_f[i] = 0.;
    depositedEnergyECAL_r[i] = 0.;
  }
  depositedEnergyHCALAct = 0.;
  depositedEnergyHCALPas = 0.;
  depositedEnergyWorld = 0.;
  depositedEnergyServices = 0.;
  depositedEnergyTimingGap = 0.;
  depositedEnergyEcalGap = 0.;
  depositedEnergyEcalDet = 0.;
  depositedEnergySolenoid = 0.;

  depositedIonEnergyTotal = 0.;
  depositedIonEnergyTiming_f = 0.;
  depositedIonEnergyTiming_r = 0.;
  for (int i = 0; i < 3; i++){
    depositedIonEnergyECAL_f[i] = 0.;
    depositedIonEnergyECAL_r[i] = 0.;
  }
  depositedIonEnergyHCALAct = 0.;
  depositedIonEnergyHCALPas = 0.;
  depositedIonEnergyWorld = 0.;
  depositedIonEnergyServices = 0.;
  depositedIonEnergyTimingGap = 0.;
  depositedIonEnergyEcalGap = 0.;
  depositedIonEnergyEcalDet = 0.;
  depositedIonEnergySolenoid = 0.;

  depositedElecEnergyTotal = 0.;
  depositedElecEnergyTiming_f = 0.;
  depositedElecEnergyTiming_r = 0.;
  for (int i = 0; i < 3; i++){
    depositedElecEnergyECAL_f[i] = 0.;
    depositedElecEnergyECAL_r[i] = 0.;
  }
  depositedElecEnergyHCALAct = 0.;
  depositedElecEnergyHCALPas = 0.;
  depositedElecEnergyWorld = 0.;
  depositedElecEnergyServices = 0.;
  depositedElecEnergyTimingGap = 0.;
  depositedElecEnergyEcalGap = 0.;
  depositedElecEnergyEcalDet = 0.;
  depositedElecEnergySolenoid = 0.;

  tot_phot_cer_Timing_f_total = 0.;
  tot_phot_cer_Timing_r_total = 0.;
  tot_phot_cer_ECAL_scinti_f_total = 0.;
  tot_phot_cer_ECAL_scinti_r_total = 0.;
  tot_phot_cer_ECAL_cheren_f_total = 0.;
  tot_phot_cer_ECAL_cheren_r_total = 0.;
  tot_phot_cer_HCAL = 0.;
  tot_phot_cer_SDdetected_ff = 0.;
  tot_phot_cer_SDdetected_fr = 0.;
  tot_phot_cer_SDdetected_rf = 0.;
  tot_phot_cer_SDdetected_rr = 0.;
  tot_phot_cer_should_det1 = 0.;
  tot_phot_cer_should_det2 = 0.;
  tot_phot_cer_should_det3 = 0.;
  tot_phot_cer_should_det4 = 0.;
  tot_phot_sin_should_det4 = 0.;
  for (int iparticle = 0; iparticle < 8; iparticle++)
  {
    depositedEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_absorb_r_particleID[iparticle] = 0.;
    depositedEnergyECAL_scinti_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_scinti_r_particleID[iparticle] = 0.;
    depositedEnergyECAL_cheren_f_particleID[iparticle] = 0.;
    depositedEnergyECAL_cheren_r_particleID[iparticle] = 0.;

    depositedIonEnergyECAL_absorb_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_absorb_r_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_scinti_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_scinti_r_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_cheren_f_particleID[iparticle] = 0.;
    depositedIonEnergyECAL_cheren_r_particleID[iparticle] = 0.;

    tot_phot_cer_ECAL_scinti_f_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_scinti_r_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_cheren_f_particleID[iparticle] = 0.;
    tot_phot_cer_ECAL_cheren_r_particleID[iparticle] = 0.;
  }

  for (int iLayer = 0; iLayer < 6; iLayer++)
  {
    Edep_Tracker_layer[iLayer] = 0.;
  }

  for (int iBar = 0; iBar < 18; iBar++)
  {
    Edep_Timing_f_ch[iBar] = 0.;
    Edep_Timing_r_ch[iBar] = 0.;
  }
  for (int iCh = 0; iCh < 6400; iCh++)
  {
    Edep_ECAL_f_ch[iCh] = 0.;
    Edep_ECAL_r_ch[iCh] = 0.;

    IonEdep_ECAL_f_ch[iCh] = 0.;
    IonEdep_ECAL_r_ch[iCh] = 0.;
/*
    ElecEdep_ECAL_f_ch[iCh] = 0.;
    ElecEdep_ECAL_r_ch[iCh] = 0.;

    Edep_ECAL_f_scinti_ch[iCh] = 0.;
    Edep_ECAL_r_scinti_ch[iCh] = 0.;

    IonEdep_ECAL_f_scinti_ch[iCh] = 0.;
    IonEdep_ECAL_r_scinti_ch[iCh] = 0.;

    ElecEdep_ECAL_f_scinti_ch[iCh] = 0.;
    ElecEdep_ECAL_r_scinti_ch[iCh] = 0.;

    Edep_ECAL_f_cheren_ch[iCh] = 0.;
    Edep_ECAL_r_cheren_ch[iCh] = 0.;

    IonEdep_ECAL_f_cheren_ch[iCh] = 0.;
    IonEdep_ECAL_r_cheren_ch[iCh] = 0.;

    ElecEdep_ECAL_f_cheren_ch[iCh] = 0.;
    ElecEdep_ECAL_r_cheren_ch[iCh] = 0.;
*/
  }

  for (int iZ = 0; iZ < 2500; iZ++)
  {
    E_Zdep_0to5000mm_total[iZ]=0.;
    E_Zdep_0to5000mm_Pion_n[iZ]=0.;
    E_Zdep_0to5000mm_Positron[iZ]=0.;
    E_Zdep_0to5000mm_Electron[iZ]=0.;
    E_Zdep_0to5000mm_Photon[iZ]=0.;
    E_Zdep_0to5000mm_Pion_p[iZ]=0.;
    E_Zdep_0to5000mm_Kion[iZ]=0.;
    E_Zdep_0to5000mm_Neutron[iZ]=0.;
    E_Zdep_0to5000mm_Proton[iZ]=0.;
  }
/*
  for (int iT = 0; iT < 2500; iT++)
  {
    E_Tdep_0to5ns_total[iT]=0.;
    E_Tdep_0to5ns_Pion_n[iT]=0.;
    E_Tdep_0to5ns_Positron[iT]=0.;
    E_Tdep_0to5ns_Electron[iT]=0.;
    E_Tdep_0to5ns_Photon[iT]=0.;
    E_Tdep_0to5ns_Pion_p[iT]=0.;
    E_Tdep_0to5ns_Kion[iT]=0.;
    E_Tdep_0to5ns_Neutron[iT]=0.;
    E_Tdep_0to5ns_Proton[iT]=0.;
  }
*/

  for (int i = 0; i < 3; ++i)
  {
    inputInitialPosition->at(i) = 0.;
    primaryPosT1->at(i) = 0.;
    primaryPosE1->at(i) = 0.;
  }
  for (int i = 0; i < 4; ++i)
  {
    inputMomentum->at(i) = 0.;
    primaryMomT1->at(i) = 0.;
    primaryMomE1->at(i) = 0.;
  }

  for (int i = 0; i < 3; ++i)
  {
    inputInitialPosition->at(i) = 0.;
    primaryPosT1->at(i) = 0.;
    primaryPosE1->at(i) = 0.;
  }
}
