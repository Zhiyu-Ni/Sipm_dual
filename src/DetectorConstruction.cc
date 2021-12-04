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
// * institutes, nor the agencies providing financial support for this *
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
// $Id: DetectorConstruction.cc, v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

#include "DetectorConstruction.hh"
#include "DR_PMTSD.hh"
#include "SurfaceProperty.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MagneticField.hh"
#include "G4UniformMagField.hh"
#include "G4FieldManager.hh"
#include "G4TransportationManager.hh"
#include "G4MagIntegratorStepper.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ClassicalRK4.hh"
#include "G4ExplicitEuler.hh"
#include "G4ChordFinder.hh"
#include "G4EqMagElectricField.hh"
#include "G4PropagatorInField.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4SubtractionSolid.hh"
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
#include "G4NistManager.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4SDManager.hh"
#include "G4MaterialPropertiesTable.hh"

#include "G4UserLimits.hh"

#include <G4TransportationManager.hh>
#include <G4MagneticField.hh>
#include <G4UniformMagField.hh>
#include <G4FieldManager.hh>
#include "CreateTree.hh"
#include <algorithm>
#include <string>
#include <sstream>
#include "G4MultiUnion.hh"

using namespace CLHEP;

DetectorConstruction::DetectorConstruction(const string &configFileName)
{
  //---------------------------------------
  //------------- Parameters --------------
  //---------------------------------------

  ConfigFile config(configFileName);

  config.readInto(checkOverlaps, "checkOverlaps");

  config.readInto(world_material, "world_material");
  config.readInto(bar_length, "bar_length");

  config.readInto(core_radius_x, "core_radius_x");
  config.readInto(core_radius_y, "core_radius_y");
  config.readInto(core_material, "core_material");
  config.readInto(core_rIndex, "core_rIndex");
  config.readInto(core_absLength, "core_absLength");

  config.readInto(gap_l, "gap_l");
  config.readInto(gap_size_x, "gap_size_x");
  config.readInto(gap_size_y, "gap_size_y");
  config.readInto(gap_material, "gap_material");

  config.readInto(det_l, "det_l");
  config.readInto(det_size_x, "det_size_x");
  config.readInto(det_size_y, "det_size_y");
  config.readInto(det_material, "det_material");
  config.readInto(ecal_det_size, "ecal_det_size");

  config.readInto(depth, "depth");
  config.readInto(cryst_dist, "cryst_dist");
  config.readInto(trackerX0, "trackerX0");
  config.readInto(services_thick, "services_thick");

  config.readInto(ecal_incline, "ecal_incline");
  config.readInto(ecal_xy_gap, "ecal_xy_gap");
  config.readInto(fiber_type, "fiber_type");
  config.readInto(ecal_material, "ecal_material");
  config.readInto(scinti_material, "scinti_material");
  config.readInto(Cherenc_material, "Cherenc_material");
  config.readInto(Cherenp_material, "Cherenp_material");
  config.readInto(wrap_material, "wrap_material");
  config.readInto(wrap_ref, "wrap_ref");
  config.readInto(rear_filter, "rear_filter");
  config.readInto(front_filter, "front_filter");
  config.readInto(ecal_front_length, "ecal_front_length");
  config.readInto(ecal_rear_length, "ecal_rear_length");
  config.readInto(ecal_front_face, "ecal_front_face");
  config.readInto(ecal_rear_face, "ecal_rear_face");
  config.readInto(hole_diameter, "hole_diameter");
  config.readInto(fiber_diameter, "fiber_diameter");
  config.readInto(wrapping_thick, "wrapping_thick");

  config.readInto(ecal_timing_distance, "ecal_timing_distance");

  config.readInto(hcal_width, "hcal_width");
  config.readInto(hcalTile_width, "hcalTile_width");
  config.readInto(hcalAbs_1_thick, "hcalAbs_1_thick");
  config.readInto(hcalAbs_2_thick, "hcalAbs_2_thick");
  config.readInto(solenoid_thick, "solenoid_thick");
  config.readInto(hcalTile_thick, "hcalTile_thick");

  B_field_intensity = config.read<double>("B_field_intensity") * tesla;

  expHall_x = 300. * cm;
  expHall_y = 300. * cm;
  expHall_z = 300. * cm;

  B_field_IsInitialized = false;

  initializeMaterials();
  initializeSurface();

  CreateTree::Instance()->inputTrackerX0 = trackerX0;
  CreateTree::Instance()->inputServiceAlmm = services_thick;
  CreateTree::Instance()->inputTimingThick = core_radius_x * 2;
  CreateTree::Instance()->inputE1Thick = ecal_front_length;
  CreateTree::Instance()->inputE2Thick = ecal_rear_length;
  CreateTree::Instance()->inputE1Width = ecal_front_face;
  CreateTree::Instance()->inputTimingECAL_dist = ecal_timing_distance;
}

//---- ---- ---- ---- ---- ---- ---- ---- ----  ---- ---- ---- ---- ---- ----

DetectorConstruction::~DetectorConstruction()
{
  delete stepLimit;
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

G4VPhysicalVolume *DetectorConstruction::Construct()
{
  G4cout << ">>>>>> DetectorConstruction::Construct ()::begin <<<<<<" << G4endl;

  //------------------------------------
  //------------- Geometry -------------
  //------------------------------------

  // The experimental Hall
  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  G4VSolid *worldS = new G4Box("worldS",  expHall_x,  expHall_y,  expHall_z);
  G4LogicalVolume *worldLV = new G4LogicalVolume(worldS, WoMaterial, "worldLV", 0, 0, 0);
  G4VPhysicalVolume *worldPV = new G4PVPlacement(0, G4ThreeVector(), worldLV, "worldPV", 0, false, 0, checkOverlaps);


  //SURFACE STATE IMPLEMENTATION
  /*
   //front face surface state (opposite to photodetector)
   G4Box* Front_skin 			= new G4Box("Front_skin", core_radius_x, core_radius_y, 0.5*depth);
   G4LogicalVolume* Front_skin_log 	= new G4LogicalVolume(Front_skin, MyMaterials::Air(), "Front_skin_log", 0, 0, 0);
//   G4VPhysicalVolume* Front_skin_phys 	= new G4PVPlacement(0, G4ThreeVector(0., 0., -(0.5*bar_length+0.5*depth)) , Front_skin_log, "Front_skin_phys", worldLV, false, 0);

   //side surface state
   G4VSolid* dummySk = new G4Box ("dummySk", core_radius_x+depth, core_radius_y+depth, 0.5*bar_length) ;
   G4VSolid* subSk   = new G4Box ("subSk", core_radius_x, core_radius_y, 0.51*bar_length);
   
   G4SubtractionSolid* Side_skin 	= new G4SubtractionSolid("Side_skin", dummySk, subSk);				
   G4LogicalVolume* Side_skin_log 	= new G4LogicalVolume(Side_skin, MyMaterials::Air(), "Side_skin_log", 0, 0, 0);
  // G4VPhysicalVolume* Side_skin_phys 	= new G4PVPlacement(0, G4ThreeVector(0., 0., 0.), Side_skin_log, "Side_skin_phys", worldLV, false, 0);

   G4LogicalBorderSurface *CrystalFrontSkin   	= NULL;
   G4LogicalBorderSurface *CrystalSideSkin   	= NULL;
   G4OpticalSurface *OpCrystalSurface		= NULL;
	
   ///-------CRYSTAL SURFACE-------
   OpCrystalSurface = new G4OpticalSurface("crystal");
   OpCrystalSurface->SetType(dielectric_dielectric);
   OpCrystalSurface->SetModel(unified);

*/

  //********************************************
  //  ELECTROMAGNETIC CALORIMETER
  //********************************************

  //  float pos_Z_E1 =  0*mm;
  //  float pos_Z_E2 =  pos_Z_T26*mm;
  //double pointingAngle = 0.001976598*M_PI;	//~0.36° -- arctan (29.31/4720)
  double pointingAngle = ecal_incline; //~0.36° -- arctan (29.31/4720)
                                       //  double alveola_thickness = 0.2 * mm;

  //  double pointingAngle = ecal_incline; //~0.36° -- arctan (29.31/4720)
  double alveola_thickness = ecal_xy_gap * mm;

  //  float total_ecal_length = ecal_front_length + ecal_rear_length;

  double window_l=0.8*mm;
  double Case_l=0.01*mm;
  //construct dream detector
  //

  //  double dchange_rotation =  sqrt(ecal_front_face* ecal_front_face + (ecal_front_length+2*det_l+2*gap_l) * (ecal_front_length+2*det_l+2*gap_l)) * abs(cos(pointingAngle * deg - atan(ecal_front_face / (ecal_front_length+2*det_l+2*gap_l)))) / 2 - (ecal_front_length+2*det_l+2*gap_l) * 0.5 + sqrt(ecal_front_face * ecal_front_face + (ecal_rear_length+2*det_l+2*gap_l) * (ecal_rear_length+2*det_l+2*gap_l)) * abs(cos(pointingAngle * deg - atan(ecal_front_face / (ecal_rear_length+2*det_l+2*gap_l)))) / 2 - (ecal_rear_length+2*det_l+2*gap_l) * 0.5;


  // ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
  //--------------------- build a dream detector here ------------------
  G4Box *ecalCrystalS_f_hole_s;
  G4Box *ecalCrystalS_f_fiber_s;
  G4Box *ecalCrystalS_f_outwrap_s;
  G4Box *ecalCrystalS_f_innerwrap_s;

  G4Tubs *ecalCrystalS_f_hole_c;
  G4Tubs *ecalCrystalS_f_fiber_c;
  G4Tubs *ecalCrystalS_f_outwrap_c;
  G4Tubs *ecalCrystalS_f_innerwrap_c;

  G4SubtractionSolid *ecalCrystalS_f_hollowwrap;
  G4SubtractionSolid *ecalCrystalS_f_absorb;

  G4SubtractionSolid *ecalCrystalS_f_fiberGap;


  G4MultiUnion *ecalCrystalS_f_fiber_scinti;
  G4MultiUnion *ecalCrystalS_f_fiber_cherenp;
  G4MultiUnion *ecalCrystalS_f_wrap;

  G4LogicalVolume *ecalCrystalL_f_absorb;
  G4LogicalVolume *ecalCrystalL_f_fiber_wrapping;
  G4LogicalVolume *ecalCrystalL_f_fiberGap;

  G4LogicalVolume *ecalCrystalL_f_fiber_scinti;
  G4LogicalVolume *ecalCrystalL_f_fiber_cherenc;
  G4LogicalVolume *ecalCrystalL_f_fiber_cherenp;


  G4Box *ecalGapS;
  G4Box *ecalDetS;
  G4Box *ecalDetWindowS;

  G4Box *ecalDetCaseOuterS;
  G4Box *ecalDetCaseInnerS;
  G4SubtractionSolid *ecalDetCase_fS;
  G4SubtractionSolid *ecalDetCase_rS;

  G4LogicalVolume *ecalGapL;
  G4LogicalVolume *ecalDetL;
  G4LogicalVolume *ecalDetWindowL;
  G4LogicalVolume *ecalDetCase_fL;
  G4LogicalVolume *ecalDetCase_rL;

    //--------------------hole & fiber-----------------------
    ecalCrystalS_f_hole_s = new G4Box("ecalCrystalS_f_hole", fiber_diameter * 0.5 + wrapping_thick + depth, fiber_diameter * 0.5 + wrapping_thick + depth, 0.5 * ecal_front_length);
    ecalCrystalS_f_fiber_s = new G4Box("ecalCrystalS_f_fiber", 0.5 * fiber_diameter, 0.5 * fiber_diameter, 0.5 * ecal_front_length);

    //---------------------- wrapping  --------------------------
    ecalCrystalS_f_outwrap_s = new G4Box("ecalCrystalS_f_outwrap_s", fiber_diameter * 0.5 + wrapping_thick + depth, fiber_diameter * 0.5 + wrapping_thick + depth, ecal_front_length * 0.5);
    ecalCrystalS_f_innerwrap_s = new G4Box("ecalCrystalS_f_innerwrap_s", fiber_diameter * 0.5 + depth, fiber_diameter * 0.5 + depth, ecal_front_length);
    ecalCrystalS_f_hollowwrap = new G4SubtractionSolid("ecalCrystalS_f_hollowwrap", ecalCrystalS_f_outwrap_s, ecalCrystalS_f_innerwrap_s);
  
    //----------------------- gap between fiber and wrap----------
    ecalCrystalS_f_fiberGap = new G4SubtractionSolid("ecalCrystalS_f_fiberGap", ecalCrystalS_f_hole_s, ecalCrystalS_f_hollowwrap);
    ecalCrystalS_f_fiberGap = new G4SubtractionSolid("ecalCrystalS_f_fiberGap", ecalCrystalS_f_fiberGap, ecalCrystalS_f_fiber_s);
    

    //----------------------detGap & Det & Detwindow-------------------
    ecalGapS = new G4Box("ecalGapS", (fiber_diameter-5) * 0.5, (fiber_diameter-5) * 0.5, 0.5 * (gap_l));
    ecalDetS = new G4Box("ecalDetS", (fiber_diameter-5) * 0.5, (fiber_diameter-5) * 0.5, 0.5 * (det_l));
    ecalDetWindowS = new G4Box("ecalDetWindowS", (fiber_diameter-5) * 0.5, (fiber_diameter-5) * 0.5, 0.5 * (window_l));
    //ecalGapS = new G4Box("ecalGapS", (fiber_diameter) * 0.5, (fiber_diameter) * 0.5, 0.5 * (gap_l));
    //ecalDetS = new G4Box("ecalDetS", (fiber_diameter) * 0.5, (fiber_diameter) * 0.5, 0.5 * (det_l));
    //ecalDetWindowS = new G4Box("ecalDetWindowS", (fiber_diameter) * 0.5, (fiber_diameter) * 0.5, 0.5 * (window_l));
    
    //---------------------det Case on two sides---------------------
    ecalDetCaseOuterS=new G4Box("ecalDetCaseOuterS", (fiber_diameter-5) * 0.5+Case_l, (fiber_diameter-5) * 0.5+Case_l, 0.5 * (det_l+window_l)+Case_l);
    ecalDetCaseInnerS=new G4Box("ecalDetCaseInnerS", (fiber_diameter-5) * 0.5, (fiber_diameter-5) * 0.5, 0.5 * (det_l+window_l)+Case_l);;
    //ecalDetCaseOuterS=new G4Box("ecalDetCaseOuterS", (fiber_diameter) * 0.5+Case_l, (fiber_diameter) * 0.5+Case_l, 0.5 * (det_l+window_l)+Case_l);
    //ecalDetCaseInnerS=new G4Box("ecalDetCaseInnerS", (fiber_diameter) * 0.5, (fiber_diameter) * 0.5, 0.5 * (det_l+window_l)+Case_l);;
    
    ecalDetCase_fS= new G4SubtractionSolid( "ecalDetCase_fS",ecalDetCaseOuterS,ecalDetCaseInnerS, 0, G4ThreeVector(0,0,Case_l) );
    ecalDetCase_rS= new G4SubtractionSolid( "ecalDetCase_rS",ecalDetCaseOuterS,ecalDetCaseInnerS, 0, G4ThreeVector(0,0,-1*Case_l) );

    //-----------------Fiber & (Det Gap & Det & window & Case Logical)--------
    ecalGapL = new G4LogicalVolume(ecalGapS, GaMaterial, "ecalGapL", 0, 0, 0);
    ecalDetL = new G4LogicalVolume(ecalDetS, DeMaterial, "ecalDetL", 0, 0, 0);
    ecalDetWindowL = new G4LogicalVolume(ecalDetWindowS, WindowMaterial, "ecalDetWindowL", 0, 0, 0);
    ecalDetCase_fL = new G4LogicalVolume(ecalDetCase_fS, WindowMaterial, "ecalDetCase_fL", 0, 0, 0);
    ecalDetCase_rL = new G4LogicalVolume(ecalDetCase_rS, WindowMaterial, "ecalDetCase_rL", 0, 0, 0);
    ecalCrystalL_f_fiber_cherenc = new G4LogicalVolume(ecalCrystalS_f_fiber_s, CherencMaterial, "ecalCrystalL_f_fiber_cherenc", 0, 0, 0);
    ecalCrystalL_f_fiberGap = new G4LogicalVolume(ecalCrystalS_f_fiberGap, GaMaterial, "ecalCrystalL_f_fiberGap", 0, 0, 0);
    ecalCrystalL_f_fiber_wrapping = new G4LogicalVolume(ecalCrystalS_f_hollowwrap, WrapMaterial, "ecalCrystalL_f_fiber_wrapping", 0, 0, 0);
    
  

  // ECAL physical placement
  const int NECAL_CRYST = 1; //6400; //2500;
  G4VPhysicalVolume *ecalCrystalP_f_wrap[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiber_scinti[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiber_cherenc[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiber_cherenp[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalGapP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetWindowP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_ff[NECAL_CRYST];
  G4VPhysicalVolume *ecalDetCaseP_fr[NECAL_CRYST];
  G4VPhysicalVolume *ecalCrystalP_f_fiberGap[NECAL_CRYST];

  G4LogicalSkinSurface* FiberWrap_f_Surface = new G4LogicalSkinSurface( "FiberWrap_f_Surface", ecalCrystalL_f_fiber_wrapping, fFiberWrapSurface );
  G4LogicalSkinSurface* PMTSurface = new G4LogicalSkinSurface( "PMTSurface",ecalDetL, fPMTSurface );
  G4LogicalSkinSurface* PMTCase_fSurface = new G4LogicalSkinSurface( "PMTCase_fSurface",ecalDetCase_fL, fPMTCaseSurface );
  G4LogicalSkinSurface* PMTCase_rSurface = new G4LogicalSkinSurface( "PMTCase_rSurface",ecalDetCase_rL, fPMTCaseSurface );
  
  G4LogicalBorderSurface* FilterSurface_ff[NECAL_CRYST];
  G4LogicalBorderSurface* FilterSurface_fr[NECAL_CRYST];
  G4LogicalBorderSurface* Fiber_cherenc_Surface[NECAL_CRYST];
//G4LogicalBorderSurface

  char name[60];
  G4double x_pos[NECAL_CRYST];
  G4double y_pos[NECAL_CRYST];
  int nArrayECAL = (int)sqrt(NECAL_CRYST);

  int iCrystal;
  for (int iX = 0; iX < nArrayECAL; iX++)
  {
    for (int iY = 0; iY < nArrayECAL; iY++)
    {

      G4RotationMatrix *piRotEcal = new G4RotationMatrix;
      piRotEcal->rotateX(pointingAngle * deg);
      G4RotationMatrix *piRotEcal_op = new G4RotationMatrix;
      piRotEcal_op->rotateX((pointingAngle+180) * deg);

      iCrystal = nArrayECAL * iX + iY;
      x_pos[iCrystal] = (iX - nArrayECAL / 2) * (ecal_front_face + alveola_thickness); // position the baricenter of crystals and then rotating them by
      y_pos[iCrystal] = (iY - nArrayECAL / 2) * (ecal_front_face / cos(pointingAngle * deg) + alveola_thickness);
      cout << " x_pos [" << iCrystal << "] = " << x_pos[iCrystal] << " :: y_pos[" << iCrystal << "] = " << y_pos[iCrystal] << " :: angle = [" << pointingAngle * iX << ", " << pointingAngle * iY << "] " << endl;

      //----------------------  wrap and fibergap--------------------------
      //add dream detector instead of the original one
      sprintf(name, "ecalCrystalP_f_wrap_%d", iCrystal);
      ecalCrystalP_f_wrap[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f_fiber_wrapping, name, worldLV, false, 0, checkOverlaps);

      sprintf(name, "ecalCrystalP_f_fiberGap_%d", iCrystal);
      ecalCrystalP_f_fiberGap[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f_fiberGap, name, worldLV, false, 0, checkOverlaps);

      //---------------------- fibers --------------------------
      sprintf(name, "ecalCrystalP_f_fiber_cherenc_%d", iCrystal);
      ecalCrystalP_f_fiber_cherenc[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal], ecal_timing_distance + ecal_front_length * 0.5), ecalCrystalL_f_fiber_cherenc, name, worldLV, false, 0, checkOverlaps);

      //-------------------surface between fiber and fibergap
      Fiber_cherenc_Surface[iCrystal] = new G4LogicalBorderSurface( "Fiber_cherenc_Surface",  ecalCrystalP_f_fiberGap[iCrystal], ecalCrystalP_f_fiber_cherenc[iCrystal], fIdealPolishedOpSurface );

      //---------------------- detectors at both ends --------------------------
      sprintf(name, "ecalDetP_ff_%d", iCrystal);
      ecalDetP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l+window_l + det_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg)), ecalDetL, name, worldLV, false, 0, checkOverlaps);

      sprintf(name, "ecalDetP_fr_%d", iCrystal);
      ecalDetP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l+window_l + det_l * 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg)), ecalDetL, name, worldLV, false, 0, checkOverlaps);
      cout << " z_pos [" << iCrystal << "] = " << ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg) << " :: z_pos[" << iCrystal << "] = " << ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+window_l + det_l * 0.5) * cos(pointingAngle * deg) << " " << endl;

      //---------------------- Det PMT windows --------------------------
      sprintf(name, "ecalDetWindowP_ff_%d", iCrystal);
      ecalDetWindowP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l+window_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+window_l * 0.5) * cos(pointingAngle * deg)), ecalDetWindowL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetWindowP_fr_%d", iCrystal);
      ecalDetWindowP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l+window_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+window_l* 0.5) * cos(pointingAngle * deg)), ecalDetWindowL, name, worldLV, false, 0, checkOverlaps);

      //----------------------Det PMT Case --------------------------
      sprintf(name, "ecalDetCaseP_ff_%d", iCrystal);
      ecalDetCaseP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l+0.5 * (det_l+window_l)+Case_l) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l+0.5 * (det_l+window_l)+Case_l) * cos(pointingAngle * deg)), ecalDetCase_fL, name, worldLV, false, 0, checkOverlaps);
      sprintf(name, "ecalDetCaseP_fr_%d", iCrystal);
      ecalDetCaseP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l+0.5 * (det_l+window_l)+Case_l) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l+0.5 * (det_l+window_l)+Case_l) * cos(pointingAngle * deg)), ecalDetCase_rL, name, worldLV, false, 0, checkOverlaps);

      //---------------------- gap between detectors and fibers --------------------------
      sprintf(name, "ecalGapP_ff_%d", iCrystal);
      ecalGapP_ff[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] - (0.5 * ecal_front_length + gap_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 - (ecal_front_length * 0.5 + gap_l* 0.5) * cos(pointingAngle * deg)), ecalGapL, name, worldLV, false, iCrystal, checkOverlaps);
      sprintf(name, "ecalGapP_fr_%d", iCrystal);
      ecalGapP_fr[iCrystal] = new G4PVPlacement(piRotEcal, G4ThreeVector(x_pos[iCrystal], y_pos[iCrystal] + (0.5 * ecal_front_length + gap_l* 0.5) * sin(pointingAngle * deg), ecal_timing_distance + ecal_front_length * 0.5 + (ecal_front_length * 0.5 + gap_l* 0.5) * cos(pointingAngle * deg)), ecalGapL, name, worldLV, false, iCrystal, checkOverlaps);
      
      //--------------------u330/ug5 surface between detgap and detwindow
      FilterSurface_ff[iCrystal] = new G4LogicalBorderSurface( "FilterSurface_ff", ecalGapP_ff[iCrystal], ecalDetWindowP_ff[iCrystal], fFilterSurface_ff );
      FilterSurface_fr[iCrystal] = new G4LogicalBorderSurface( "FilterSurface_fr", ecalGapP_fr[iCrystal], ecalDetWindowP_fr[iCrystal], fFilterSurface_fr );

    }
  }

  //-----------------------------------------------------
  //------------- sensitive detector & surface --------------
  //-----------------------------------------------------

  auto sdman = G4SDManager::GetSDMpointer();
  auto fDRDetSD = new DR_PMTSD("/DR_Det");
  sdman->AddNewDetector(fDRDetSD);
  ecalDetL->SetSensitiveDetector(fDRDetSD);


  if (B_field_intensity > 0.1 * tesla)
  ConstructField();

  G4cout << ">>>>>> DetectorConstruction::Construct ()::end <<< " << G4endl;
  return worldPV;
}

void DetectorConstruction::initializeMaterials()
{
  //-----------------
  // define materials

  WoMaterial = NULL;
  if (world_material == 0)
    WoMaterial = MyMaterials::Vacuum();
  else if (world_material == 1)
    WoMaterial = MyMaterials::Air();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre world material specifier " << world_material << G4endl;
    exit(-1);
  }
  G4cout << "Wo. material: " << WoMaterial << G4endl;

  CoMaterial = NULL;
  if (core_material == 1)
    CoMaterial = MyMaterials::Quartz();
  else if (core_material == 2)
    CoMaterial = MyMaterials::SiO2();
  else if (core_material == 3)
    CoMaterial = MyMaterials::SiO2_Ce();
  else if (core_material == 4)
    CoMaterial = MyMaterials::LuAG_Ce();
  else if (core_material == 5)
    CoMaterial = MyMaterials::YAG_Ce();
  else if (core_material == 6)
    CoMaterial = MyMaterials::LSO();
  else if (core_material == 7)
    CoMaterial = MyMaterials::LYSO();
  else if (core_material == 8)
    CoMaterial = MyMaterials::LuAG_undoped();
  else if (core_material == 9)
    CoMaterial = MyMaterials::GAGG_Ce();
  else if (core_material == 11)
    CoMaterial = MyMaterials::LuAG_Pr();
  else if (core_material == 12)
    CoMaterial = MyMaterials::PbF2();
  else if (core_material == 13)
    CoMaterial = MyMaterials::PlasticBC408();
  else if (core_material == 14)
    CoMaterial = MyMaterials::PlasticBC418();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << core_material << G4endl;
    exit(-1);
  }
  G4cout << "Co. material: " << CoMaterial << G4endl;

  EcalMaterial = NULL;
  if (ecal_material == 1)
    EcalMaterial = MyMaterials::Quartz();
  else if (ecal_material == 2)
    EcalMaterial = MyMaterials::SiO2();
  else if (ecal_material == 3)
    EcalMaterial = MyMaterials::SiO2_Ce();
  else if (ecal_material == 4)
    EcalMaterial = MyMaterials::LuAG_Ce();
  else if (ecal_material == 5)
    EcalMaterial = MyMaterials::YAG_Ce();
  else if (ecal_material == 6)
    EcalMaterial = MyMaterials::LSO();
  else if (ecal_material == 7)
    EcalMaterial = MyMaterials::LYSO();
  else if (ecal_material == 8)
    EcalMaterial = MyMaterials::LuAG_undoped();
  else if (ecal_material == 9)
    EcalMaterial = MyMaterials::GAGG_Ce();
  else if (ecal_material == 10)
    EcalMaterial = MyMaterials::LuAG_Pr();
  else if (ecal_material == 11)
    EcalMaterial = MyMaterials::PbF2();
  else if (ecal_material == 12)
    EcalMaterial = MyMaterials::PlasticBC408();
  else if (ecal_material == 13)
    EcalMaterial = MyMaterials::PlasticBC418();
  else if (ecal_material == 14)
    EcalMaterial = MyMaterials::PWO();
  else if (ecal_material == 15)
    EcalMaterial = MyMaterials::Acrylic();
  else if (ecal_material == 16)
    EcalMaterial = MyMaterials::copper();
  else if (ecal_material == 17)
    EcalMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << ecal_material << G4endl;
    exit(-1);
  }
  G4cout << "ECAL material: " << EcalMaterial << G4endl;

  /************************************************************************************/
  ScintiMaterial = NULL;
  if (scinti_material == 1)
    ScintiMaterial = MyMaterials::Quartz();
  else if (scinti_material == 2)
    ScintiMaterial = MyMaterials::SiO2();
  else if (scinti_material == 3)
    ScintiMaterial = MyMaterials::SiO2_Ce();
  else if (scinti_material == 4)
    ScintiMaterial = MyMaterials::LuAG_Ce();
  else if (scinti_material == 5)
    ScintiMaterial = MyMaterials::YAG_Ce();
  else if (scinti_material == 6)
    ScintiMaterial = MyMaterials::LSO();
  else if (scinti_material == 7)
    ScintiMaterial = MyMaterials::LYSO();
  else if (scinti_material == 8)
    ScintiMaterial = MyMaterials::LuAG_undoped();
  else if (scinti_material == 9)
    ScintiMaterial = MyMaterials::GAGG_Ce();
  else if (scinti_material == 10)
    ScintiMaterial = MyMaterials::LuAG_Pr();
  else if (scinti_material == 11)
    ScintiMaterial = MyMaterials::PbF2();
  else if (scinti_material == 12)
    ScintiMaterial = MyMaterials::PlasticBC408();
  else if (scinti_material == 13)
    ScintiMaterial = MyMaterials::PlasticBC418();
  else if (scinti_material == 14)
    ScintiMaterial = MyMaterials::PWO();
  else if (scinti_material == 15)
    ScintiMaterial = MyMaterials::Acrylic();
  else if (scinti_material == 16)
    ScintiMaterial = MyMaterials::copper();
  else if (scinti_material == 17)
    ScintiMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << scinti_material << G4endl;
    exit(-1);
  }

  CherencMaterial = NULL;
  if (Cherenc_material == 1)
    CherencMaterial = MyMaterials::Quartz();
  else if (Cherenc_material == 2)
    CherencMaterial = MyMaterials::SiO2();
  else if (Cherenc_material == 3)
    CherencMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenc_material == 4)
    CherencMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenc_material == 5)
    CherencMaterial = MyMaterials::YAG_Ce();
  else if (Cherenc_material == 6)
    CherencMaterial = MyMaterials::LSO();
  else if (Cherenc_material == 7)
    CherencMaterial = MyMaterials::LYSO();
  else if (Cherenc_material == 8)
    CherencMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenc_material == 9)
    CherencMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenc_material == 10)
    CherencMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenc_material == 11)
    CherencMaterial = MyMaterials::PbF2();
  else if (Cherenc_material == 12)
    CherencMaterial = MyMaterials::PlasticBC408();
  else if (Cherenc_material == 13)
    CherencMaterial = MyMaterials::PlasticBC418();
  else if (Cherenc_material == 14)
    CherencMaterial = MyMaterials::PWO();
  else if (Cherenc_material == 15)
    CherencMaterial = MyMaterials::Acrylic();
  else if (Cherenc_material == 16)
    CherencMaterial = MyMaterials::copper();
  else if (Cherenc_material == 17)
    CherencMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenc_material << G4endl;
    exit(-1);
  }

  CherenpMaterial = NULL;
  if (Cherenp_material == 1)
    CherenpMaterial = MyMaterials::Quartz();
  else if (Cherenp_material == 2)
    CherenpMaterial = MyMaterials::SiO2();
  else if (Cherenp_material == 3)
    CherenpMaterial = MyMaterials::SiO2_Ce();
  else if (Cherenp_material == 4)
    CherenpMaterial = MyMaterials::LuAG_Ce();
  else if (Cherenp_material == 5)
    CherenpMaterial = MyMaterials::YAG_Ce();
  else if (Cherenp_material == 6)
    CherenpMaterial = MyMaterials::LSO();
  else if (Cherenp_material == 7)
    CherenpMaterial = MyMaterials::LYSO();
  else if (Cherenp_material == 8)
    CherenpMaterial = MyMaterials::LuAG_undoped();
  else if (Cherenp_material == 9)
    CherenpMaterial = MyMaterials::GAGG_Ce();
  else if (Cherenp_material == 10)
    CherenpMaterial = MyMaterials::LuAG_Pr();
  else if (Cherenp_material == 11)
    CherenpMaterial = MyMaterials::PbF2();
  else if (Cherenp_material == 12)
    CherenpMaterial = MyMaterials::PlasticBC408();
  else if (Cherenp_material == 13)
    CherenpMaterial = MyMaterials::PlasticBC418();
  else if (Cherenp_material == 14)
    CherenpMaterial = MyMaterials::PWO();
  else if (Cherenp_material == 15)
    CherenpMaterial = MyMaterials::Acrylic();
  else if (Cherenp_material == 16)
    CherenpMaterial = MyMaterials::copper();
  else if (Cherenp_material == 17)
    CherenpMaterial = MyMaterials::EJ200();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << Cherenp_material << G4endl;
    exit(-1);
  }

  WrapMaterial = NULL;
  if (wrap_material == 1)
    WrapMaterial = MyMaterials::Quartz();
  else if (wrap_material == 2)
    WrapMaterial = MyMaterials::SiO2();
  else if (wrap_material == 3)
    WrapMaterial = MyMaterials::SiO2_Ce();
  else if (wrap_material == 4)
    WrapMaterial = MyMaterials::LuAG_Ce();
  else if (wrap_material == 5)
    WrapMaterial = MyMaterials::YAG_Ce();
  else if (wrap_material == 6)
    WrapMaterial = MyMaterials::LSO();
  else if (wrap_material == 7)
    WrapMaterial = MyMaterials::LYSO();
  else if (wrap_material == 8)
    WrapMaterial = MyMaterials::LuAG_undoped();
  else if (wrap_material == 9)
    WrapMaterial = MyMaterials::GAGG_Ce();
  else if (wrap_material == 10)
    WrapMaterial = MyMaterials::LuAG_Pr();
  else if (wrap_material == 11)
    WrapMaterial = MyMaterials::PbF2();
  else if (wrap_material == 12)
    WrapMaterial = MyMaterials::PlasticBC408();
  else if (wrap_material == 13)
    WrapMaterial = MyMaterials::PlasticBC418();
  else if (wrap_material == 14)
    WrapMaterial = MyMaterials::PWO();
  else if (wrap_material == 15)
    WrapMaterial = MyMaterials::Acrylic();
  else if (wrap_material == 16)
    WrapMaterial = MyMaterials::copper();
  else if (wrap_material == 17)
    WrapMaterial = MyMaterials::Epoxy();
  else if (wrap_material == 18)
    WrapMaterial = G4NistManager::Instance()->FindOrBuildMaterial( "G4_Al" );
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid fibre clad material specifier " << wrap_material << G4endl;
    exit(-1);
  }

  /************************************************************************************/
  WindowMaterial = MyMaterials::PyrexGlass();

  GaMaterial = NULL;
  if (gap_material == 1)
    GaMaterial = MyMaterials::Air();
  else if (gap_material == 2)
    GaMaterial = MyMaterials::OpticalGrease();
  else if (gap_material == 3)
    GaMaterial = MyMaterials::MeltMount168();
  else if (gap_material == 4)
    GaMaterial = MyMaterials::OpticalGrease155();
  else if (gap_material == 5)
    GaMaterial = MyMaterials::silicone();
  else if (gap_material == 6)
    GaMaterial = MyMaterials::PyrexGlass();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid gap material specifier " << gap_material << G4endl;
    exit(-1);
  }
  G4cout << "Gap material: " << gap_material << G4endl;

  DeMaterial = NULL;
  if (det_material == 1)
    DeMaterial = MyMaterials::Silicon();
  else if (det_material == 2)
    DeMaterial = MyMaterials::Quartz();
  else if (det_material == 3)
    DeMaterial = MyMaterials::Air();
  else if (det_material == 4)
    DeMaterial = MyMaterials::Bialkali();
  else
  {
    G4cerr << "<DetectorConstructioninitializeMaterials>: Invalid detector material specifier " << det_material << G4endl;
    exit(-1);
  }
  G4cout << "Detector material: " << det_material << G4endl;

  //------------------
  // change properties

  if (core_absLength > 0)
  {
    const G4int nEntries_ABS = 2;
    G4double PhotonEnergy_ABS[nEntries_ABS] = {1. * eV, 10. * eV};
    G4double Absorption[nEntries_ABS] = {core_absLength * mm, core_absLength * mm};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("ABSLENGTH");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("ABSLENGTH", PhotonEnergy_ABS, Absorption, nEntries_ABS);
  }
  if (core_rIndex > 0)
  {
    const G4int nEntries_RI = 2;
    G4double PhotonEnergy_RI[nEntries_RI] = {1. * eV, 10. * eV};
    G4double RefractiveIndex[nEntries_RI] = {core_rIndex, core_rIndex};

    CoMaterial->GetMaterialPropertiesTable()->RemoveProperty("RINDEX");
    CoMaterial->GetMaterialPropertiesTable()->AddProperty("RINDEX", PhotonEnergy_RI, RefractiveIndex, nEntries_RI);
  }
}

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----

void DetectorConstruction::ConstructField()
{
  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::begin <<<<<<" << G4endl;

  static G4TransportationManager *trMgr = G4TransportationManager::GetTransportationManager();

  // A field object is held by a field manager
  // Find the global Field Manager
  G4FieldManager *globalFieldMgr = trMgr->GetFieldManager();

  if (!B_field_IsInitialized)
  {
    // magnetic field parallel to the beam direction (w/ tilt)
    G4ThreeVector fieldVector(0.0522 * B_field_intensity, 0.0522 * B_field_intensity, 0.9973 * B_field_intensity);

    B_field = new G4UniformMagField(fieldVector);
    globalFieldMgr->SetDetectorField(B_field);
    globalFieldMgr->CreateChordFinder(B_field);
    globalFieldMgr->GetChordFinder()->SetDeltaChord(0.005 * mm);
    B_field_IsInitialized = true;
  }

  G4cout << ">>>>>> DetectorConstruction::ConstructField ()::end <<< " << G4endl;
  return;
}

void DetectorConstruction::SetMaxStep(G4double maxStep)
{
  if ((stepLimit) && (maxStep > 0.))
    stepLimit->SetMaxAllowedStep(maxStep);
}

// //-----------------------------------------------
// //------------- Surface properties --------------
// //-----------------------------------------------

void DetectorConstruction::initializeSurface()
{
  fFiberWrapSurface = MakeS_wrap();
  fPMTSurface = MakeS_PMT();
  fPMTCaseSurface = MakeS_IdealWhiteSurface();

  static const unsigned nentries = 2;
  static double phoE[nentries]   = {0.01*eV, 10.0*eV};
  double reflectivity[nentries]  = {wrap_ref, wrap_ref};
  G4MaterialPropertiesTable* table = fFiberWrapSurface->GetMaterialPropertiesTable();
  if( table ){
    table->RemoveProperty( "REFLECTIVITY" );
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
  } else {
    table = new G4MaterialPropertiesTable();
    table->AddProperty( "REFLECTIVITY", phoE, reflectivity, nentries );
    fFiberWrapSurface->SetMaterialPropertiesTable( table );
  }
  fIdealPolishedOpSurface = MakeS_IdealPolished();
  fFilterSurface_ff = MakeS_IdealPolished();
  fFilterSurface_fr = MakeS_IdealPolished();
/*
  const G4int nEntries_tran_u330 = 18;
  G4double PhotonEnergy_tran_u330[nEntries_tran_u330] = {4.103235582*eV, 4.048912761*eV,3.944470896    *eV,3.81331407      *eV,3.750947295     *eV,3.690598697     *eV,3.603625797     *eV,3.534216059     *eV,3.456978466     *eV,3.380554073     *eV,3.309819592     *eV,3.230572067     *eV,3.185706353     *eV,3.131340814     *eV,3.087086539     *eV,3.050146549     *eV,2.992445212     *eV,2.933127681     *eV};
  G4double transIndex_u330[nEntries_tran_u330] = {0.201372        ,0.202705        ,0.211043        ,0.227125        ,0.234102        ,0.233987        ,0.235942        ,0.235798        ,0.229958        ,0.219856        ,0.199831        ,0.155664        ,0.115833        ,0.068881        ,0.037554        ,0.017592        ,0.00466 ,0.000935        };

  const G4int nEntries_tran_ug5 = 25;
  G4double PhotonEnergy_tran_ug5[nEntries_tran_ug5] = {4.092143963	*eV,4.034407045	*eV,3.943981544	*eV,3.825267479	*eV,3.743879546	*eV,3.62234533	*eV,3.530100421	*eV,3.414187953	*eV,3.300875562	*eV,3.233225815	*eV,3.168293273	*eV,3.137871012	*eV,3.099604675	*eV,3.074608111	*eV,3.041899835	*eV,3.001980276	*eV,2.947821354	*eV,2.873756177	*eV,2.764363413	*eV,2.697530944	*eV,2.642988727	*eV,2.564470256	*eV,2.529030177	*eV,2.498643446	*eV,2.482374634	*eV};
  G4double transIndex_ug5[nEntries_tran_ug5] = {0.197435	,0.199462	,0.209962	,0.224666	,0.230184	,0.23489	,0.238226	,0.232947	,0.219825	,0.204682	,0.176017	,0.153858	,0.125979	,0.106681	,0.089491	,0.06658	,0.04575	,0.031245	,0.029397	,0.029832	,0.02319	,0.016426	,0.012699	,0.009703	,0.008199	};
*/
  const G4int nEntries_tran_u330 = 21;
  G4double PhotonEnergy_tran_u330[nEntries_tran_u330] = {	4.13280661444001*eV,	3.99949027203872*eV,	3.87450620103751*eV,	3.75709692221819*eV,	3.64659407156471*eV,	3.54240566952001*eV,	3.44400551203334*eV,	3.35092428197839*eV,	3.26274206403159*eV,	3.17908201110770*eV,	3.09960496083001*eV,	3.02400483983415*eV,	2.95200472460001*eV,	2.88335345193489*eV,	2.81782269166364*eV,	2.75520440962667*eV,	2.69530866159131*eV,	2.63796166879150*eV,	2.58300413402501*eV,	2.53028976394286*eV,	2.47968396866401*eV};
  G4double transIndex_u330[nEntries_tran_u330] = {0.89, 0.89, 0.89, 0.9, 0.89, 0.88,0.87,0.82,0.7,0.45,0.18,0.05,0.02,0.007,0.005,0,0,0,0,0,0};
  //{0.291390000000000,	0.316220000000000	,0.331920000000000	,0.342950000000000	,0.343670000000000	,0.339970000000000	,0.332720000000000	,0.312390000000000	,0.266490000000000	,0.169380000000000	,0.0668400000000000	,0.0181100000000000	,0.00703000000000000	,0.00238000000000000	,0.00163000000000000	,0	,0	,0	,0	,0	,0	};

  const G4int nEntries_tran_ug5 = 21;
  G4double PhotonEnergy_tran_ug5[nEntries_tran_ug5] = {	4.13280661444001*eV,	3.99949027203872*eV,	3.87450620103751*eV,	3.75709692221819*eV,	3.64659407156471*eV,	3.54240566952001*eV,	3.44400551203334*eV,	3.35092428197839*eV,	3.26274206403159*eV,	3.17908201110770*eV,	3.09960496083001*eV,	3.02400483983415*eV,	2.95200472460001*eV,	2.88335345193489*eV,	2.81782269166364*eV,	2.75520440962667*eV,	2.69530866159131*eV,	2.63796166879150*eV,	2.58300413402501*eV,	2.53028976394286*eV,	2.47968396866401*eV};
  G4double transIndex_ug5[nEntries_tran_ug5] = {	0.284840000000000,	0.312660000000000,	0.331920000000000	,0.339130000000000	,0.343670000000000	,0.343830000000000,	0.340370000000000,	0.331440000000000	,0.315980000000000	,0.271010000000000,	0.193110000000000,	0.119550000000000	,0.0738300000000000,	0.0545000000000000	,0.0488800000000000	,0.0490900000000000	,0.0483600000000000	,0.0395000000000000	,0.0321600000000000,	0.0257100000000000,	0.0196200000000000};


    G4MaterialPropertiesTable* Filtertable_u330 = new G4MaterialPropertiesTable();
    Filtertable_u330->AddProperty("REFLECTIVITY", PhotonEnergy_tran_u330, transIndex_u330, nEntries_tran_u330);
    G4MaterialPropertiesTable* Filtertable_ug5 = new G4MaterialPropertiesTable();
    Filtertable_ug5->AddProperty("REFLECTIVITY", PhotonEnergy_tran_ug5, transIndex_ug5, nEntries_tran_ug5);

    if(front_filter== -1){
      std::cout<<"No front Filter!"<<std::endl;
    }else if(front_filter==0){
      fFilterSurface_ff->SetMaterialPropertiesTable( Filtertable_u330 );
    }else if(front_filter==1){
      fFilterSurface_ff->SetMaterialPropertiesTable( Filtertable_ug5 );
    }
    if(rear_filter== -1){
      std::cout<<"No rear Filter!"<<std::endl;
    }else if(rear_filter==0){
      fFilterSurface_fr->SetMaterialPropertiesTable( Filtertable_u330 );
    }else if(rear_filter==1){
      fFilterSurface_fr->SetMaterialPropertiesTable( Filtertable_ug5 );
    }

  return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ---- ----
