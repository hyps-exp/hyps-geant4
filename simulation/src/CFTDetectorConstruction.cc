// ********************************************************************
// * DISCLAIMER                                                       *
// *                                                                  *
// * The following disclaimer summarizes all the specific disclaimers *
// * of contributors to this software. The specific disclaimers,which *
// * govern, are listed with their locations in:                      *
// *   http://cern.ch/geant4/license                                  *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.                                                             *
// *                                                                  *
// * This  code  implementation is the  intellectual property  of the *
// * GEANT4 collaboration.                                            *
// * By copying,  distributing  or modifying the Program (or any work *
// * based  on  the Program)  you indicate  your  acceptance of  this *
// * statement, and all its terms.                                    *
// ********************************************************************
//
//
// $Id: ExN03DetectorConstruction.cc,v 1.19 2003/11/25 14:23:44 gcosmo Exp $
// GEANT4 tag $Name: geant4-06-00-patch-01 $
//
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CFTDetectorConstruction.hh"
#include "MaterialList.hh"
#include "RadDeg.hh"

#include "ConfMan.hh"
#include "DetectorID.hh"
#include "DCGeomMan.hh"

#include "CFTFiberSD.hh"
#include "CrystalSD.hh"
#include "SKSChamberSD.hh"
#include "SKSCounterSD.hh"
#include "SKSVirtualPlaneSD.hh"

#include "SKSField.hh"

#include "common.h"

#include "G4FieldManager.hh"
#include "G4ChordFinder.hh"
#include "G4TransportationManager.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Sphere.hh"
#include "G4Torus.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4UniformMagField.hh"
#include "G4UserLimits.hh"

#include "G4ThreeVector.hh"
#include "G4RotationMatrix.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SystemOfUnits.hh"

#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTDetectorConstruction::CFTDetectorConstruction()
  : mList_(0), solidWorld(0),logicWorld(0),physWorld(0), chamberSD(0), counterSD(0), virtualPlaneSD(0)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

CFTDetectorConstruction::~CFTDetectorConstruction()
{
 }

MaterialList *CFTDetectorConstruction::DefineMaterials()
{
  // not safe for exceotion
  if(mList_) delete mList_;
  return new MaterialList();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* CFTDetectorConstruction::Construct()
{
  // Clean old geometry, if any
  //
  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  // Define materials
  mList_ = DefineMaterials();


  // Magnetic field
  ConfMan *confMan = ConfMan::GetConfManager();
  if (confMan->FieldMapName() == "") {
    G4cout << "Field Map is not set. Now there is no magnetic field Ok?" << G4endl;
    G4cout << "Press return " ;
    getchar();
  }

  static G4bool fieldIsInitialized = false;
  if(!fieldIsInitialized)
  {
    SKSField* myField = new SKSField(confMan->FieldMapName());
    double MapScale = confMan->MapScale();
    if (MapScale != 0.)
      myField->SetMapScale(MapScale);

    G4FieldManager* fieldMgr;

    fieldMgr
      = G4TransportationManager::GetTransportationManager()
        ->GetFieldManager();

    fieldMgr->SetDetectorField(myField);
    fieldMgr->CreateChordFinder(myField);
    fieldMgr->GetChordFinder()->SetDeltaChord( 1.0E-3*mm );

    fieldIsInitialized = true;

    for (double z=-400.; z<400.; z+=10.) {
      double pos[3] = {0.*mm, 0.*mm, z*mm};
      double value[3];
      myField->GetFieldValue(pos, value);
      printf("z = %.1lf Bx = %lf, By = %lf, Bz = %lf\n", z, value[0]/tesla, value[1]/tesla, value[2]/tesla);
    }

  }


  //
  // World
  //
  G4double WorldSizeX=10.0*m;
  G4double WorldSizeYZ=10.0*m;

  solidWorld = new G4Box("World",				//its name
                   WorldSizeX/2,WorldSizeYZ/2,WorldSizeYZ/2);	//its size

  logicWorld = new G4LogicalVolume(solidWorld,		//its solid
                                   mList_->Air,        	//its material
				   //mList_->Vacuum,        	//its material
                                   "World");		//its name

  physWorld = new G4PVPlacement(0,			//no rotation
  				 G4ThreeVector(),	//at (0,0,0)
                                 logicWorld,		//its logical volume
                                 "World",		//its name
                                 0,			//its mother  volume
                                 false,			//no boolean operation
                                 0);			//copy number


  // Construct CFT
  ConstructCFT(physWorld);
  // Construct LH2 target
  ConstructTarget(physWorld);

  // Construct KURAMA magnet
  ConstructKURAMA(physWorld);
  // Construct SdcIn
  ConstructSdcIn(physWorld);
  // Construct SdcOut
  ConstructSdcOut(physWorld);
  // Construct CH
  //ConstructCH(physWorld);
  //ConstructCHwithBeamHole(physWorld);
  // Construct SH
  //ConstructSH(physWorld);
  // Construct TOF
  ConstructTOF(physWorld);
  // Construct T0
  ConstructT0(physWorld);
  // Construct AC
  ConstructAC(physWorld);
  ConstructeVeto(physWorld);
  // Construct BcOut
  // ConstructBcOut(physWorld);
  // Construct Virtual Plane
  ConstructVP(physWorld);
  // Construct BeamDump
  //ConstructBeamDump(physWorld);
  // Construct S-2S
  //ConstructS2S(physWorld);


  //
  // Visualization attributes
  //
  auto invisibility = new G4VisAttributes();
  invisibility->SetVisibility(false);
  logicWorld->SetVisAttributes (invisibility);
  // logicWorld->SetVisAttributes (G4VisAttributes::invisibility);
  // logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
  //G4VisAttributes* simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  //simpleBoxVisAtt->SetVisibility(true);
  //logicCalor->SetVisAttributes(simpleBoxVisAtt);

  //
  //always return the physical World
  //
  return physWorld;

}

void  CFTDetectorConstruction::ConstructCFT(G4VPhysicalVolume *pMother)
{
  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);


  char name[100];

  G4double r1ContCFT = 46.*mm;
  //G4double r2ContCFT = 96.*mm;
  G4double r2ContCFT = 170.*mm;
  G4double zContCFT  = 560.*mm;
  //G4double zContCFT  = 800.*mm;

  G4Tubs *solidContCFT =
    new G4Tubs( "solidContCFT",
		r1ContCFT,
		r2ContCFT,
		zContCFT/2,
		0.*degree, 360.*degree );

  G4LogicalVolume *logicContCFT = new G4LogicalVolume(solidContCFT,//its solid
						       mList_->Vacuum,        	//its material
						       "logicContCFT");		//its name
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum = 61; // CFT_U1
  G4ThreeVector gloPosContCFT = geomMan.GetGlobalPosition(lnum)*mm;
  G4double RA1 = geomMan.GetRotAngle1(lnum);
  G4double RA2 = geomMan.GetRotAngle2(lnum);

  G4RotationMatrix rotContCFT;
  rotContCFT.rotateX(RA1*deg);
  rotContCFT.rotateZ(RA2*deg);

  G4VPhysicalVolume* physContCFT =
    new G4PVPlacement( G4Transform3D(rotContCFT, gloPosContCFT),
		       "physContCFT",
		       logicContCFT,
		       pMother,
		       false,
		       0 );


  G4double r1SupportCFT = 46.*mm;
  //G4double r2SupportCFT = 96.*mm;
  G4double r2SupportCFT = 85.9*mm;
  G4double zSupportCFT  = 10.*mm;

  G4Tubs *solidSupportCFT =
    new G4Tubs( "solidSupportCFT",
		r1SupportCFT,
		r2SupportCFT,
		zSupportCFT/2,
		0.*degree, 360.*degree );


  G4LogicalVolume *logicSupportCFT = new G4LogicalVolume(solidSupportCFT,//its solid
							 //mList_->Acrylic,        	//its material
							 mList_->Al,        	//its material
							 "logicSupportCFT");		//its name

  G4VPhysicalVolume* physSupportCFT[2];

  for (int i=0; i<2; i++) {
    sprintf(name, "physSupportCFT-%d", i);
    double z;
    if (i==0)
      z = (zContCFT/2 - zSupportCFT/2.);
    else
      z = -(zContCFT/2 - zSupportCFT/2.);

    G4ThreeVector pos(0*mm,0*mm, z);

    physSupportCFT[i] = new G4PVPlacement( 0,
					   pos,
					   name,
					   logicSupportCFT,
					   physContCFT,
					   false,
					   i );

  }

  double zCFTFiber = 400.*mm;

  G4Tubs *solidSupportCFT2 =
    new G4Tubs( "solidSupportCFT",
		r1SupportCFT,
		r2SupportCFT,
		zSupportCFT/2,
		0.*degree, 360.*degree );


  G4LogicalVolume *logicSupportCFT2 = new G4LogicalVolume(solidSupportCFT,//its solid
							  //mList_->Al,        	//its material
							  mList_->Air,        	//its material
							 "logicSupportCFT2");		//its name

  G4VPhysicalVolume* physSupportCFT2[2];

  for (int i=0; i<2; i++) {
    sprintf(name, "physSupportCFT2-%d", i);
    double z;
    if (i==0)
      z = (zCFTFiber/2 + 5*mm + zSupportCFT/2.);
    else
      z = -(zCFTFiber/2 + 5*mm + zSupportCFT/2.);

    G4ThreeVector pos(0*mm,0*mm, z);

    physSupportCFT2[i] = new G4PVPlacement( 0,
					    pos,
					    name,
					    logicSupportCFT2,
					    physContCFT,
					    false,
					    i );

  }


  G4Tubs *solidCFTFiberPHI[4];
  G4LogicalVolume *logicCFTFiberPHI[4];
  G4VPhysicalVolume* physCFTFiberPHI[4][NumOfSegCFT_PHI4];

  for (int l=0; l<NumOfPlaneCFT; l++) {
    int layerPhi=0;
    if (l==CFT_PHI1)
      layerPhi=0;
    else if (l==CFT_PHI2)
      layerPhi=1;
    else if (l==CFT_PHI3)
      layerPhi=2;
    else if (l==CFT_PHI4)
      layerPhi=3;
    else
      continue;


    sprintf(name, "solidCFTFiberPHI%d", layerPhi+1);
    solidCFTFiberPHI[layerPhi] = new G4Tubs( name,
					     0,
					     DiameterOfCFTFiber[l]/2.*mm,
					     zCFTFiber/2,
					     0.*degree, 360.*degree );

    sprintf(name, "logicCFTFiberPHI%d", layerPhi+1);
    logicCFTFiberPHI[layerPhi] = new G4LogicalVolume(solidCFTFiberPHI[layerPhi],          //its solid
						     mList_->Scin,          //its material
						     name);                 //its name


    for (int i=0; i<NumOfSegCFT[l]; i++) {
      sprintf(name, "physCFTFiberPHI%d_%d", layerPhi+1, i);
      double R=RadiusOfCFT[l];

      //if (DiameterOfCFTFiber[l] == 0.75) {
      if (DiameterOfCFTFiber[l] <= 0.85) {
	if (i%2 == 0)
	  R -= 0.4;
	else
	  R += 0.4;
      } else if ((DiameterOfCFTFiber[l] == 1.0)) {
	if (i%2 == 0)
	  R -= 0.5;
	else
	  R += 0.5;
      } else {
	fprintf(stderr, "CFTDetectorConstruction::ConstructCFT No such FiberDiameter %lf\n", DiameterOfCFTFiber[l]);
	exit(-1);
      }

      double theta = (360./(double)NumOfSegCFT[l])*i;

      G4ThreeVector pos(R*cos(theta*Deg2Rad)*mm,R*sin(theta*Deg2Rad)*mm, 0*mm);

      physCFTFiberPHI[layerPhi][i] = new G4PVPlacement( 0,
							pos,
							name,
							logicCFTFiberPHI[layerPhi],
							physContCFT,
							false,
							i );

    }
  }


  G4Tubs *solidCFTFiberUV[4];
  G4LogicalVolume *logicCFTFiberUV[4];
  G4VPhysicalVolume* physCFTFiberUV[4];

  for (int l=0; l<NumOfPlaneCFT; l++) {
    int layerUV=0;
    if (l==CFT_U1)
      layerUV=0;
    else if (l==CFT_V2)
      layerUV=1;
    else if (l==CFT_U3)
      layerUV=2;
    else if (l==CFT_V4)
      layerUV=3;
    else
      continue;


    sprintf(name, "solidCFTFiberUV%d", layerUV+1);
    solidCFTFiberUV[layerUV] = new G4Tubs( name,
					   (RadiusOfCFT[l]-DiameterOfCFTFiber[l]/2.)*mm,
					   (RadiusOfCFT[l]+DiameterOfCFTFiber[l]/2.)*mm,
					   zCFTFiber/2,
					   0.*degree, 360.*degree );

    sprintf(name, "logicCFTFiberUV%d", layerUV+1);
    logicCFTFiberUV[layerUV] = new G4LogicalVolume(solidCFTFiberUV[layerUV],          //its solid
						   mList_->Scin,          //its material
						     name);                 //its name

    sprintf(name, "physCFTFiberUV%d", layerUV+1);
    physCFTFiberUV[layerUV] = new G4PVPlacement( 0,
						  G4ThreeVector(0.*mm, 0*mm, 0*mm),
						  name,
						  logicCFTFiberUV[layerUV],
						  physContCFT,
						  false,
						  0 );


  }

  G4Tubs *solidBGO_VP = new G4Tubs( "solidBGO_VP",
				    (RadiusOfBGOSurface-2.)*mm,
				    (RadiusOfBGOSurface-1.)*mm,
				      zContCFT/2,
				      0.*degree, 360.*degree );

  G4LogicalVolume *logicBGO_VP = new G4LogicalVolume(solidBGO_VP,          //its solid
						     mList_->Air,          //its material
						     "logicBGO_VP");                 //its name

  G4VPhysicalVolume* physBGO_VP = new G4PVPlacement( 0,
						  G4ThreeVector(0.*mm, 0*mm, 0*mm),
						  "physBGO_VP",
						  logicBGO_VP,
						  physContCFT,
						  false,
						  0 );




  G4Box *solidBGO =     new G4Box("solidBGO", //its name
				  BGO_X/2,
				  BGO_Y/2,
				  BGO_Z/2);	//its size


  G4LogicalVolume *logicBGO
    =  new G4LogicalVolume(solidBGO, //its solid
			   mList_->BGO,   //its material
			   //mList_->Air,   //its material
			   "logicBGO");  //its name

  G4Box *solidBGO_PMT =     new G4Box("solidBGO_PMT", //its name
				      BGO_PMT_X/2,
				      BGO_PMT_Y/2,
				      BGO_PMT_Z/2);	//its size

  G4LogicalVolume *logicBGO_PMT
    =  new G4LogicalVolume(solidBGO_PMT, //its solid
			   mList_->Fe,   //its material
			   //mList_->Air,   //its material
			   "logicBGO_PMT");  //its name


#if 0
  G4VPhysicalVolume* physBGO[NumOfBGOUnit*NumOfBGOInOneUnit];

  G4RotationMatrix rotBGOUnitCont[NumOfBGOUnit];

  for (int i=0; i<NumOfBGOUnit; i++) {
    double theta = 22.5+(double)i*45.;

    rotBGOUnitCont[i].rotateZ(-90.*deg);
    rotBGOUnitCont[i].rotateZ(theta*deg);

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-1)*BGO_X;

      G4ThreeVector pos((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
			(x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
			0*mm);

      sprintf(name, "physBGO%d", j+3*i);
      G4cout << "i " << i << ", j" << j << ", " << j+3*i << ", " << name << G4endl;

      /*
      physBGO[j+3*i] = new G4PVPlacement( 0,
					  G4ThreeVector((double)(1-j)*BGO_X*mm, 0*mm, 0*mm),
					  name,
					  logicBGO,
					  physBGOUnitCont[i],
					  false,
					  0 );
      */

      physBGO[j+3*i] = new G4PVPlacement( G4Transform3D(rotBGOUnitCont[i], pos),
					    name,
					    logicBGO,
					    physContCFT,
					    false,
					    0 );
    }


  }
#endif

  const G4double zOffsetBGO=60*mm;

  G4VPhysicalVolume* physBGO[NumOfBGOUnit*(NumOfBGOInOneUnit+NumOfBGOInOneUnit2)];
  G4VPhysicalVolume* physBGO_PMT[NumOfBGOUnit*(NumOfBGOInOneUnit+NumOfBGOInOneUnit2)];

  G4RotationMatrix rotBGOUnitCont[NumOfBGOUnit];
  G4RotationMatrix rotBGOUnitCont2[NumOfBGOUnit];

  for (int i=0; i<NumOfBGOUnit; i++) {

    double theta = (double)i*45.;

    rotBGOUnitCont[i].rotateZ(-90.*deg);
    rotBGOUnitCont[i].rotateZ(theta*deg);

    for (int j=0; j<NumOfBGOInOneUnit; j++) {
      double x0 = RadiusOfBGOSurface+BGO_Y/2;
      double y0 = (double)(j-0.5)*BGO_X;

      G4ThreeVector pos((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
			(x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
			//0*mm);
			//60*mm);
			zOffsetBGO);

      sprintf(name, "physBGO%d", j+3*i);
      G4cout << "i " << i << ", j" << j << ", " << j+3*i << ", " << name << G4endl;

      /*
      physBGO[j+3*i] = new G4PVPlacement( 0,
					  G4ThreeVector((double)(1-j)*BGO_X*mm, 0*mm, 0*mm),
					  name,
					  logicBGO,
					  physBGOUnitCont[i],
					  false,
					  0 );
      */


      ///////////
      //if (i!=0) {
      if (1) {
	physBGO[j+3*i] = new G4PVPlacement( G4Transform3D(rotBGOUnitCont[i], pos),
					    name,
					    logicBGO,
					    physContCFT,
					    false,
					    0 );



	G4ThreeVector pos_PMT((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
			      (x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
			      zOffsetBGO-(BGO_Z+BGO_PMT_Z)/2.*mm);

	sprintf(name, "physBGO_PMT%d", j+3*i);
	physBGO_PMT[j+3*i] = new G4PVPlacement( G4Transform3D(rotBGOUnitCont[i], pos_PMT),
						name,
						logicBGO_PMT,
						physContCFT,
						false,
						0 );


      }
    }


    theta = 22.5+(double)i*45.;

    rotBGOUnitCont2[i].rotateZ(-90.*deg);
    rotBGOUnitCont2[i].rotateZ(theta*deg);

    for (int j=0; j<NumOfBGOInOneUnit2; j++) {
      double x0 = RadiusOfBGOSurface2+BGO_Y/2;
      double y0 = (double)(j)*BGO_X;

      G4ThreeVector pos((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
			(x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
			//0*mm);
			//60*mm);
			zOffsetBGO);

      sprintf(name, "physBGO%d", j+3*i+NumOfBGOInOneUnit);
      G4cout << "i " << i << ", j" << j << ", " << j+3*i << ", " << name << G4endl;

      /*
      physBGO[j+3*i] = new G4PVPlacement( 0,
					  G4ThreeVector((double)(1-j)*BGO_X*mm, 0*mm, 0*mm),
					  name,
					  logicBGO,
					  physBGOUnitCont[i],
					  false,
					  0 );
      */

      physBGO[j+NumOfBGOInOneUnit+3*i] = new G4PVPlacement( G4Transform3D(rotBGOUnitCont2[i], pos),
							    name,
							    logicBGO,
							    physContCFT,
							    false,
							    0 );

      G4ThreeVector pos_PMT((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
			    (x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
			    //0*mm);
			    //60*mm);
			    zOffsetBGO-(BGO_Z+BGO_PMT_Z)/2.*mm);

      sprintf(name, "physBGO_PMT%d", j+3*i+NumOfBGOInOneUnit);

      physBGO_PMT[j+NumOfBGOInOneUnit+3*i] = new G4PVPlacement( G4Transform3D(rotBGOUnitCont2[i], pos_PMT),
								name,
								logicBGO_PMT,
								physContCFT,
								false,
								0 );

    }
  }

  G4Tubs *solidBGO_MagShld = new G4Tubs( "solidBGO_MagShld",
					 BGO_MagShld_R*mm,
					 (BGO_MagShld_R+BGO_MagShld_Thick)*mm,
					 BGO_MagShld_Z/2*mm,
					 0.*degree, 360.*degree );

  G4LogicalVolume *logicBGO_MagShld
    =  new G4LogicalVolume(solidBGO_MagShld, //its solid
			   mList_->Fe,   //its material
			   //mList_->Air,   //its material
			   "logicBGO_MagShld");  //its name
  logicBGO_MagShld->SetVisAttributes(counter_att);

  G4VPhysicalVolume* physBGO_MagShld = new G4PVPlacement( 0,
							  G4ThreeVector(0.*mm, 0*mm, zOffsetBGO-(BGO_Z-BGO_MagShld_Z)/2.*mm),
							  "physBGO_MagShld",
							  logicBGO_MagShld,
							  physContCFT,
							  false,
							  0 );





  /* E10 BGO rate study*/
#if 0
  /*
  G4Box *solidBGO2 =     new G4Box("solidBGO2", //its name
				   400./2*mm,
				   25./2*mm,
				   32./2*mm);	//its size


  G4LogicalVolume *logicBGO2
    =  new G4LogicalVolume(solidBGO2, //its solid
			   mList_->BGO,   //its material
			   "logicBGO2");  //its name

  G4RotationMatrix rotBGOUnitCont2;
  rotBGOUnitCont2.rotateZ(-8.*deg);

  G4ThreeVector localpos(-200.*mm,
			 -30*mm,
			 -200*mm);
  */
  int lnumTgt = geomMan.GetDetectorId("Target");
  /*
  G4ThreeVector gloPosBGO2 = geomMan.Local2GlobalPos(lnumTgt, localpos);
  G4VPhysicalVolume* physBGO2 = new G4PVPlacement( G4Transform3D(rotBGOUnitCont2, gloPosBGO2),
						   "physBGO24",
						   logicBGO2,
						   pMother,
						   false,
						   0 );
  */
  G4Box *solidBH2 =     new G4Box("solidBH2", //its name
				   5./2*mm,
				   250./2*mm,
				   100./2*mm);	//its size

  G4LogicalVolume *logicBH2
    =  new G4LogicalVolume(solidBH2, //its solid
			   mList_->Scin,   //its material
			   "logicBH2");  //its name

  G4RotationMatrix rotBH2;
  rotBH2.rotateZ(-8.*deg);

  G4ThreeVector localpos2(0.*mm,
			 0*mm,
			 -420*mm);

  G4ThreeVector gloPosBH2 = geomMan.Local2GlobalPos(lnumTgt, localpos2);
  G4VPhysicalVolume* physBH2 = new G4PVPlacement( G4Transform3D(rotBH2, gloPosBH2),
						   "physBH2",
						   logicBH2,
						   pMother,
						   false,
						   0 );

#endif



  G4Box *solidPiV =     new G4Box("solidPiV", //its name
				  PiV_X/2,
				  PiV_Y/2,
				  PiV_Z/2);	//its size


  G4LogicalVolume *logicPiV
    =  new G4LogicalVolume(solidPiV, //its solid
			   mList_->Scin,   //its material
			   "logicPiV");  //its name

  G4VPhysicalVolume* physPiV[NumOfPiVUnit*NumOfPiVInOneUnit];

  G4RotationMatrix rotPiVUnitCont[NumOfPiVUnit];


  for (int i=0; i<NumOfPiVUnit; i++) {
    //double theta = 22.5+(double)i*45.;
    double theta = (double)i*45.;

    rotPiVUnitCont[i].rotateZ(-90.*deg);
    rotPiVUnitCont[i].rotateZ(theta*deg);

    for (int j=0; j<NumOfPiVInOneUnit; j++) {
      double x0 = RadiusOfPiVSurface+PiV_Y/2;
      double y0 = (double)(j-1)*PiV_X;

      G4ThreeVector pos((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
			(x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
			//0*mm);
			//60*mm);
			zOffsetBGO);

      sprintf(name, "physPiV%d", j+3*i);
      G4cout << "i " << i << ", j" << j << ", " << j+3*i << ", " << name << G4endl;

      /*
      physBGO[j+3*i] = new G4PVPlacement( 0,
					  G4ThreeVector((double)(1-j)*BGO_X*mm, 0*mm, 0*mm),
					  name,
					  logicBGO,
					  physBGOUnitCont[i],
					  false,
					  0 );
      */

      physPiV[j+3*i] = new G4PVPlacement( G4Transform3D(rotPiVUnitCont[i], pos),
					    name,
					    logicPiV,
					    physContCFT,
					    false,
					    0 );
    }


  }

  G4Box *solidPiV2 =     new G4Box("solidPiV2", //its name
				  PiV2_X/2,
				  PiV2_Y/2,
				  PiV2_Z/2);	//its size



  G4LogicalVolume *logicPiV2
    =  new G4LogicalVolume(solidPiV2, //its solid
			   mList_->Scin,   //its material
			   "logicPiV2");  //its name

  G4VPhysicalVolume* physPiV2[NumOfPiVUnit];

  G4RotationMatrix rotPiV2UnitCont[NumOfPiVUnit];

  for (int i=0; i<NumOfPiVUnit; i++) {
    //double theta = (double)i*45.;
    double theta = 22.5+(double)i*45.;

    rotPiV2UnitCont[i].rotateZ(-90.*deg);
    rotPiV2UnitCont[i].rotateZ(theta*deg);


    double x0 = RadiusOfPiV2Surface+PiV2_Y/2;
    double y0 = 0.;

    G4ThreeVector pos((x0*cos(theta*Deg2Rad) - y0*sin(theta*Deg2Rad))*mm,
		      (x0*sin(theta*Deg2Rad) + y0*cos(theta*Deg2Rad))*mm,
		      //0*mm);
		      //60*mm);
		      zOffsetBGO);

    sprintf(name, "physPiV%d", NumOfPiVUnit*NumOfPiVInOneUnit+i);

    physPiV2[i] = new G4PVPlacement( G4Transform3D(rotPiV2UnitCont[i], pos),
				     name,
				     logicPiV2,
				     physContCFT,
				     false,
				     0 );

  }


  //////// Sensitive Detectors ////////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  CFTFiberSD *cftFiberSD = new CFTFiberSD( "cftFiberSD" );
  SDMan->AddNewDetector( cftFiberSD );

  for (int i=0; i<4; i++)
    if (logicCFTFiberPHI[i])
      logicCFTFiberPHI[i]->SetSensitiveDetector( cftFiberSD );

  for (int i=0; i<4; i++)
    logicCFTFiberUV[i]->SetSensitiveDetector( cftFiberSD );


  CrystalSD *crystalSD = new CrystalSD( "crystalSD" );
  SDMan->AddNewDetector( crystalSD );

  logicBGO->SetSensitiveDetector( crystalSD );

#if 0
  /* E10 BGO rate study*/
  logicBGO2->SetSensitiveDetector( crystalSD );
#endif

  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  logicPiV->SetSensitiveDetector( counterSD );
  logicPiV2->SetSensitiveDetector( counterSD );

  logicBGO_VP->SetSensitiveDetector( counterSD );

}

void  CFTDetectorConstruction::ConstructTarget(G4VPhysicalVolume *pMother)
{
  G4Colour aqua(0.247, 0.8, 1.0);
  G4Colour magenta(1.0, 0.0, 1.0);
  G4Colour green(0.0, 1.0, 0.0);

  G4VisAttributes *LH2_att = new G4VisAttributes(aqua);
  LH2_att->SetForceWireframe(true);

  //-----------------LH2 Target
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum = geomMan.GetDetectorId("Target");

  G4RotationMatrix rotLH2;
  rotLH2.rotateX(90.*deg);
  rotLH2.rotateZ(geomMan.GetRotAngle2(lnum)*deg);

  G4ThreeVector localPosLH2( 0., 0., 0.);
  G4ThreeVector gloPosLH2 = geomMan.Local2GlobalPos(lnum, localPosLH2);

  G4double TargetRmin = 0.0 * mm;
  G4double TargetRmax = LH2TgtR/2.*mm;
  G4double TargetZ    = Target_Length/2.*mm;
  G4double TargetPhi1 = 0.0 * deg;
  G4double TargetPhi2 = 360.0 * deg;

  G4Tubs* LH2_Tubs = new G4Tubs("LH2_Tubs",     // name
				TargetRmin,   // RMin
				TargetRmax,   // RMax
				TargetZ,      // Dz
				TargetPhi1,   // Start Angle
				TargetPhi2);  // Spanning Angle


  G4LogicalVolume* LH2_log = new G4LogicalVolume(LH2_Tubs,
						 mList_->LH,
						 "LH2_log",0,0,0);
  LH2_log->SetVisAttributes(LH2_att);
  G4double maxStep=2.0*mm;
  LH2_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* LH2_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2),
		       "LH2_phys",
		       LH2_log,
		       pMother,
		       false,
		       0 );


  //=======================================Target vessel
  G4VisAttributes *TargetVess_att = new G4VisAttributes(green);
  TargetVess_att->SetForceWireframe(true);

  G4double TargetVessRmin = LH2TgtR/2.*mm;
  G4double TargetVessRmax = (LH2TgtR/2. + TargetVessThickness)*mm;
  G4double TargetVessZ    = Target_Length/2.*mm;
  G4double TargetVessPhi1 = 0.0 * deg;
  G4double TargetVessPhi2 = 360.0 * deg;

  G4Tubs* TargetVess_Tubs = new G4Tubs("targetVess",     // name
				       TargetVessRmin,   // RMin
				       TargetVessRmax,   // RMax
				       TargetVessZ,      // Dz
				       TargetVessPhi1,   // Start Angle
				       TargetVessPhi2);  // Spanning Angle

  G4LogicalVolume* TargetVess_log = new G4LogicalVolume(TargetVess_Tubs,   // Solid
							mList_->Mylar,      // Material
							"TargetVess_log",   // name
							0,              // FieldMgr
							0,              // SDetector
							0);             // ULimits

  TargetVess_log->SetVisAttributes(TargetVess_att);

  G4VPhysicalVolume* TargetVess_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2),
		       "TargetVess_phys",
		       TargetVess_log,
		       pMother,
		       false,
		       0 );


  //=========== target edge
  G4double TargetUp_Z    = 35/2.*mm;

  G4ThreeVector localPosLH2_UpEdge( 0., 0., -Target_Length/2-TargetUp_Z);
  G4ThreeVector gloPosLH2_UpEdge = geomMan.Local2GlobalPos(lnum, localPosLH2_UpEdge);


  G4Tubs* LH2_UpEdge_Tubs = new G4Tubs("LH2_UpEdge_Tubs",     // name
				TargetRmin,   // RMin
				TargetRmax,   // RMax
				TargetUp_Z,      // Dz
				TargetPhi1,   // Start Angle
				TargetPhi2);  // Spanning Angle


  G4LogicalVolume* LH2_UpEdge_log = new G4LogicalVolume(LH2_UpEdge_Tubs,
							mList_->LH,
							"LH2_UpEdge_log",0,0,0);
  LH2_UpEdge_log->SetVisAttributes(LH2_att);
  LH2_UpEdge_log->SetUserLimits(new G4UserLimits(maxStep));


  G4VPhysicalVolume* LH2_UpEdge_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_UpEdge),
		       "LH2_UpEdge_phys",
		       LH2_UpEdge_log,
		       pMother,
		       false,
		       0 );


  G4Tubs* LH2_UpEdgeVess_Tubs = new G4Tubs("LH2_UpEdgeVess_Tubs",     // name
					   TargetVessRmin,   // RMin
					   TargetVessRmax,   // RMax
					   TargetUp_Z,      // Dz
					   TargetPhi1,   // Start Angle
					   TargetPhi2);  // Spanning Angle


  G4LogicalVolume* LH2_UpEdgeVess_log = new G4LogicalVolume(LH2_UpEdgeVess_Tubs,
							     mList_->Mylar,
							     "LH2_UpEdgeVess_log",0,0,0);
  LH2_UpEdgeVess_log->SetVisAttributes(TargetVess_att);


  G4VPhysicalVolume* LH2_UpEdgeVess_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_UpEdge),
		       "LH2_UpEdgeVess_phys",
		       LH2_UpEdgeVess_log,
		       pMother,
		       false,
		       0 );


  G4ThreeVector localPosLH2_UpSphere( 0., 0., -Target_Length/2-TargetUp_Z*2);
  G4ThreeVector gloPosLH2_UpSphere = geomMan.Local2GlobalPos(lnum, localPosLH2_UpSphere);

  G4Sphere* LH2_UpEdge_Sphere = new G4Sphere("LH2_UpEdge_Spere",     // name
					     TargetRmin,   // RMin
					     TargetRmax,   // RMax
					     0*degree,      // StartPhi
					     360*degree,    // EndPhi
					     90*degree,     // Start Theta
					     180*degree);   // End Theta

  G4LogicalVolume* LH2_UpEdge_Sphere_log = new G4LogicalVolume(LH2_UpEdge_Sphere,
							       mList_->LH,
							       "LH2_UpEdge_Sphere_log",0,0,0);
  LH2_UpEdge_Sphere_log->SetVisAttributes(LH2_att);
  LH2_UpEdge_Sphere_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* LH2_UpEdge_Sphere_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_UpSphere),
		       "LH2_UpEdge_Sphere_phys",
		       LH2_UpEdge_Sphere_log,
		       pMother,
		       false,
		       0 );


  G4Sphere* LH2_UpEdgeVess_Sphere = new G4Sphere("LH2_UpEdgeVess_Spere",     // name
					     TargetVessRmin,   // RMin
					     TargetVessRmax,   // RMax
					     0*degree,      // StartPhi
					     360*degree,    // EndPhi
					     90*degree,     // Start Theta
					     180*degree);   // End Theta

  G4LogicalVolume* LH2_UpEdgeVess_Sphere_log = new G4LogicalVolume(LH2_UpEdgeVess_Sphere,
								   mList_->Mylar,
								   "LH2_UpEdgeVess_Sphere_log",0,0,0);
  LH2_UpEdgeVess_Sphere_log->SetVisAttributes(TargetVess_att);

  G4VPhysicalVolume* LH2_UpEdgeVess_Sphere_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_UpSphere),
		       "LH2_UpEdgeVess_Sphere_phys",
		       LH2_UpEdgeVess_Sphere_log,
		       pMother,
		       false,
		       0 );


  /*
  // Upstream SUS cover case
  G4double SUS_Cap_thickness = 0.05*mm;
  G4ThreeVector localPosLH2_UpSUSCap( 0., 0., -Target_Length/2-TargetUp_Z*2-SUS_Cap_thickness/2);
  G4ThreeVector gloPosLH2_UpSUSCap = geomMan.Local2GlobalPos(lnum, localPosLH2_UpSUSCap);


  G4Tubs* LH2_UpSUSCap_Tubs = new G4Tubs("LH2_UpSUSCap_Tubs",     // name
				TargetRmin,   // RMin
				TargetRmax,   // RMax
				SUS_Cap_thickness/2,      // Dz
				TargetPhi1,   // Start Angle
				TargetPhi2);  // Spanning Angle


  G4LogicalVolume* LH2_UpSUSCap_log = new G4LogicalVolume(LH2_UpSUSCap_Tubs,
							mList_->SUS316L,
							"LH2_UpSUSCap_log",0,0,0);

  G4VPhysicalVolume* LH2_UpSUSCap_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_UpSUSCap),
		       "LH2_UpSUSCap_phys",
		       LH2_UpSUSCap_log,
		       pMother,
		       false,
		       0 );

  */



  G4ThreeVector localPosLH2_DownSphere( 0., 0., Target_Length/2);
  G4ThreeVector gloPosLH2_DownSphere = geomMan.Local2GlobalPos(lnum, localPosLH2_DownSphere);

  G4Sphere* LH2_DownEdge_Sphere = new G4Sphere("LH2_DownEdge_Spere",     // name
					     TargetRmin,   // RMin
					     TargetRmax,   // RMax
					     0*degree,      // StartPhi
					     360*degree,    // EndPhi
					     0*degree,     // Start Theta
					     90*degree);   // End Theta

  G4LogicalVolume* LH2_DownEdge_Sphere_log = new G4LogicalVolume(LH2_DownEdge_Sphere,
								 mList_->LH,
								 "LH2_DownEdge_Sphere_log",0,0,0);
  LH2_DownEdge_Sphere_log->SetVisAttributes(LH2_att);
  LH2_DownEdge_Sphere_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* LH2_DownEdge_Sphere_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_DownSphere),
		       "LH2_DownEdge_Sphere_phys",
		       LH2_DownEdge_Sphere_log,
		       pMother,
		       false,
		       0 );


  G4Sphere* LH2_DownEdgeVess_Sphere = new G4Sphere("LH2_DownEdgeVess_Spere",     // name
					     TargetVessRmin,   // RMin
					     TargetVessRmax,   // RMax
					     0*degree,      // StartPhi
					     360*degree,    // EndPhi
					     0*degree,     // Start Theta
					     90*degree);   // End Theta

  G4LogicalVolume* LH2_DownEdgeVess_Sphere_log = new G4LogicalVolume(LH2_DownEdgeVess_Sphere,
							       mList_->Mylar,
							       "LH2_DownEdgeVess_Sphere_log",0,0,0);
  LH2_DownEdgeVess_Sphere_log->SetVisAttributes(TargetVess_att);

  G4VPhysicalVolume* LH2_DownEdgeVess_Sphere_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2_DownSphere),
		       "LH2_DownEdgeVess_Sphere_phys",
		       LH2_DownEdgeVess_Sphere_log,
		       pMother,
		       false,
		       0 );


#if 1
  //=========== target holder

  // G10
  G4double TgtHolder_Z    = TgtG10HolderZ/2.*mm;
 G4double TgtHolder_Thicnkess = TgtG10HolderThickness;
  G4Material *HolderMater = mList_->G10;

  /*
  // SUS
  G4double TgtHolder_Z    = TgtSUSHolderZ/2.*mm;
  G4double TgtHolder_Thicnkess = TgtSUSHolderThickness;
  G4Material *HolderMater = mList_->SUS316L;
  */

  G4double TgtHolderRmin = (LH2TgtR/2. + TargetVessThickness)*mm;
  G4double TgtHolderRmax = (LH2TgtR/2. + TargetVessThickness + TgtHolder_Thicnkess)*mm;

  G4ThreeVector localPosTgtHolder( 0., 0., -Target_Length/2-TgtHolder_Z);
  G4ThreeVector gloPosTgtHolder = geomMan.Local2GlobalPos(lnum, localPosTgtHolder);

  G4Tubs* TgtHolder_Tubs = new G4Tubs("TgtHolder_Tubs",     // name
				TgtHolderRmin,   // RMin
				TgtHolderRmax,   // RMax
				TgtHolder_Z,      // Dz
				TargetPhi1,   // Start Angle
				TargetPhi2);  // Spanning Angle

  G4LogicalVolume* TgtHolder_log = new G4LogicalVolume(TgtHolder_Tubs,
						       HolderMater,
						       "TgtHolder_log",0,0,0);
  TgtHolder_log->SetVisAttributes(TargetVess_att);

  G4VPhysicalVolume* TgtHolder_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosTgtHolder),
		       "TgtHolder_phys",
		       TgtHolder_log,
		       pMother,
		       false,
		       0 );
#endif

#if 1
  //=======================================Vaccum chamber window
  G4VisAttributes *VaccumWin_att = new G4VisAttributes(magenta);
  VaccumWin_att->SetForceWireframe(true);

  G4double VaccumWinRmin = VaccumChamWinR/2.*mm ;
  G4double VaccumWinRmax = (VaccumChamWinR/2. + VaccumChamThickness)*mm;
  G4double VaccumWinZ    = VaccumChamWinZ/2.*mm;
  G4double VaccumWinPhi1 = 0.0 * deg;
  G4double VaccumWinPhi2 = 360.0 * deg;

  G4Tubs* VaccumWin_Tubs = new G4Tubs("VaccumWin",     // name
				       VaccumWinRmin,   // RMin
				       VaccumWinRmax,   // RMax
				       VaccumWinZ,      // Dz
				       VaccumWinPhi1,   // Start Angle
				       VaccumWinPhi2);  // Spanning Angle

  G4LogicalVolume* VaccumWin_log = new G4LogicalVolume(VaccumWin_Tubs,    // Solid
						       mList_->CFRP,             // Material
				       "VaccumWin_log",   // name
				       0,              // FieldMgr
				       0,              // SDetector
				       0);             // ULimits

  VaccumWin_log->SetVisAttributes(VaccumWin_att);

  G4VPhysicalVolume* VaccumWin_phys =
    new G4PVPlacement( G4Transform3D(rotLH2, gloPosLH2),
		       "VaccumWin_phys",
		       VaccumWin_log,
		       pMother,
		       false,
		       0 );

#endif






}

void  CFTDetectorConstruction::ConstructKURAMA(G4VPhysicalVolume *pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  // Vis Attribure
  G4Colour blue(0.0, 0.0, 1.0);
  G4Colour red(1.0, 0.0, 0.0);

  G4VisAttributes *magnet_att = new G4VisAttributes(blue);
  magnet_att->SetForceWireframe(true);

  // define magnet size and position

  G4double size_MFIELD[XYZ];
  G4double size_YOKE_UD[XYZ];
  G4double size_YOKE_LR[XYZ];
  G4double size_YOKE_LR_GapSpace[XYZ];
  G4double size_UGUARD_UD[XYZ];
  G4double size_UGUARD_LR[XYZ];
  G4double size_UGUARD_R1[XYZ];
  G4double size_UGUARD_R2[XYZ];
  G4double size_DGUARD_UD[XYZ];
  G4double size_DGUARD_LR[XYZ];

  // KURAMA
  //size_MFIELD[XCOORD] = 400.0*mm;
  //size_MFIELD[YCOORD] = 500.0*mm;
  //size_MFIELD[ZCOORD] = 400.0*mm;
  // LEPS
  size_MFIELD[XCOORD] = 400.0*mm;
  size_MFIELD[YCOORD] = 675.0*mm;
  size_MFIELD[ZCOORD] = 275.0*mm;

  G4double size_COIL1[3];
  size_COIL1[XCOORD] = 280./2*mm;
  size_COIL1[YCOORD] = 900./2*mm;
  size_COIL1[ZCOORD] = 193./2*mm;

  G4double size_COIL2[3];
  size_COIL2[XCOORD] = 280.0/2.*mm;
  size_COIL2[YCOORD] = 193./2.*mm;
  size_COIL2[ZCOORD] = 58.5*mm;


  G4double size_COIL3[3];
  size_COIL3[XCOORD] = 280.0/2.*mm;
  size_COIL3[YCOORD] = 193.0/2.*mm;
  size_COIL3[ZCOORD] = 68.5*mm;

  G4double size_COIL4[3];
  size_COIL4[XCOORD] = 740.0/2*mm;
  size_COIL4[YCOORD] = 193./2*mm;
  size_COIL4[ZCOORD] = 280./2*mm;

  G4double size_COIL5[3];
  size_COIL5[XCOORD] = 280./2*mm;
  size_COIL5[YCOORD] = 900.0/2*mm;
  size_COIL5[ZCOORD] = 214./2*mm;


  size_YOKE_UD[XCOORD] = 400.0*mm;
  //size_YOKE_UD[YCOORD] = 1000.0*mm;
  size_YOKE_UD[YCOORD] = 1100.0*mm;
  //size_YOKE_UD[ZCOORD] = 250.0*mm;
  size_YOKE_UD[ZCOORD] = 370.0/2.*mm;

  size_YOKE_LR[XCOORD] = 400.0*mm;
  //size_YOKE_LR[YCOORD] = 250.0*mm;
  size_YOKE_LR[YCOORD] = 200.0*mm;
  //size_YOKE_LR[ZCOORD] = 250.0*mm;
  size_YOKE_LR[ZCOORD] = 400.0*mm;

  size_YOKE_LR_GapSpace[XCOORD] = 400.0*mm;
  size_YOKE_LR_GapSpace[YCOORD] = 220.0/2.*mm;
  size_YOKE_LR_GapSpace[ZCOORD] = 120.0/2.*mm;

  size_UGUARD_UD[XCOORD] =  50.0*mm;
  size_UGUARD_UD[YCOORD] = 950.0*mm;
  size_UGUARD_UD[ZCOORD] = 310.0*mm;

  size_UGUARD_LR[XCOORD] =  50.0*mm;
  //size_UGUARD_LR[YCOORD] = 400.0*mm;
  size_UGUARD_LR[YCOORD] = 325.0*mm;
  size_UGUARD_LR[ZCOORD] = 150.0*mm;

  size_UGUARD_R1[XCOORD] =  50.0*mm;
  size_UGUARD_R1[YCOORD] = 590.0/2*mm;
  size_UGUARD_R1[ZCOORD] = 30.0*mm;

  size_UGUARD_R2[XCOORD] =  50.0*mm;
  size_UGUARD_R2[YCOORD] = 325.0*mm;
  size_UGUARD_R2[ZCOORD] = 60.0*mm;

  size_DGUARD_UD[XCOORD] =  50.0*mm;
  size_DGUARD_UD[YCOORD] = 800.0*mm;
  size_DGUARD_UD[ZCOORD] = 210.0*mm;

  size_DGUARD_LR[XCOORD] =  50.0*mm;
  size_DGUARD_LR[YCOORD] = 125.0*mm;
  //size_DGUARD_LR[ZCOORD] = 350.0*mm;
  size_DGUARD_LR[ZCOORD] = 550.0*mm;

  G4double pos_MFIELD[XYZ];
  G4double pos_YOKE_U[XYZ];
  G4double pos_YOKE_D[XYZ];
  G4double pos_YOKE_L[XYZ];
  G4double pos_YOKE_R[XYZ];
  G4double pos_YOKE_L_GapSpace[XYZ];
  G4double pos_YOKE_R_GapSpace[XYZ];
  G4double pos_UGUARD_U[XYZ];
  G4double pos_UGUARD_D[XYZ];
  G4double pos_UGUARD_L[XYZ];
  G4double pos_UGUARD_R[XYZ];
  G4double pos_UGUARD_R1[XYZ];
  G4double pos_UGUARD_R2U[XYZ];
  G4double pos_UGUARD_R2D[XYZ];
  G4double pos_DGUARD_U[XYZ];
  G4double pos_DGUARD_D[XYZ];
  G4double pos_DGUARD_L[XYZ];
  G4double pos_DGUARD_R[XYZ];

  pos_MFIELD[XCOORD] = 0.0*mm;
  pos_MFIELD[YCOORD] = 0.0*mm;
  pos_MFIELD[ZCOORD] = 0.0*mm;

  pos_YOKE_U[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_U[YCOORD] = 0.0*mm;
  pos_YOKE_U[ZCOORD] = size_MFIELD[ZCOORD]+size_YOKE_UD[ZCOORD];

  pos_YOKE_D[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_D[YCOORD] = 0.0*mm;
  pos_YOKE_D[ZCOORD] = -(size_MFIELD[ZCOORD]+size_YOKE_UD[ZCOORD]);

  pos_YOKE_L[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_L[YCOORD] = size_MFIELD[YCOORD]+size_YOKE_LR[YCOORD]+200.*mm;
  pos_YOKE_L[ZCOORD] = 0.0*mm;

  pos_YOKE_R[XCOORD] = pos_MFIELD[XCOORD];
  pos_YOKE_R[YCOORD] = -(size_MFIELD[YCOORD]+size_YOKE_LR[YCOORD]+200.*mm);
  pos_YOKE_R[ZCOORD] = 0.0*mm;

  pos_YOKE_L_GapSpace[XCOORD] = pos_YOKE_L[XCOORD];
  pos_YOKE_L_GapSpace[YCOORD] = pos_YOKE_L[YCOORD]-size_YOKE_LR[YCOORD]-size_YOKE_LR_GapSpace[YCOORD];
  pos_YOKE_L_GapSpace[ZCOORD] = 0.0*mm;

  pos_YOKE_R_GapSpace[XCOORD] = pos_YOKE_R[XCOORD];
  pos_YOKE_R_GapSpace[YCOORD] = pos_YOKE_R[YCOORD]+size_YOKE_LR[YCOORD]+size_YOKE_LR_GapSpace[YCOORD];
  pos_YOKE_R_GapSpace[ZCOORD] = 0.0*mm;

  pos_UGUARD_U[XCOORD] = pos_MFIELD[XCOORD]-820.0*mm + size_UGUARD_UD[XCOORD];
  //pos_UGUARD_U[YCOORD] = pos_MFIELD[YCOORD] - 150.0*mm;
  pos_UGUARD_U[YCOORD] = pos_MFIELD[YCOORD]; // E07
  pos_UGUARD_U[ZCOORD] = size_UGUARD_LR[ZCOORD]+size_UGUARD_UD[ZCOORD];

  pos_UGUARD_D[XCOORD] = pos_MFIELD[XCOORD]- 820.0*mm + size_UGUARD_UD[XCOORD];
  //pos_UGUARD_D[YCOORD] = pos_MFIELD[YCOORD] - 150.0*mm;
  pos_UGUARD_D[YCOORD] = pos_MFIELD[YCOORD]; // E07
  pos_UGUARD_D[ZCOORD] = -(size_UGUARD_LR[ZCOORD]+size_UGUARD_UD[ZCOORD]);

  pos_UGUARD_L[XCOORD] = pos_MFIELD[XCOORD]- 820.0*mm + size_UGUARD_LR[XCOORD];
  //pos_UGUARD_L[YCOORD] = pos_MFIELD[YCOORD]- 150.0*mm + (150.0*mm + size_UGUARD_LR[YCOORD]);
  pos_UGUARD_L[YCOORD] =  pos_UGUARD_U[YCOORD]+size_UGUARD_UD[YCOORD] - size_UGUARD_LR[YCOORD];
  pos_UGUARD_L[ZCOORD] = 0.0*mm;

  pos_UGUARD_R[XCOORD] = pos_MFIELD[XCOORD]- 820.0*mm + size_UGUARD_LR[XCOORD];
  //pos_UGUARD_R[YCOORD] = pos_MFIELD[YCOORD] - 150.0*mm - (150.0*mm + size_UGUARD_LR[YCOORD]);
  pos_UGUARD_R[YCOORD] = pos_UGUARD_U[YCOORD] -(size_UGUARD_UD[YCOORD] - size_UGUARD_LR[YCOORD]);
  pos_UGUARD_R[ZCOORD] = 0.0*mm;

  pos_UGUARD_R1[XCOORD] = pos_MFIELD[XCOORD]- 820.0*mm + size_UGUARD_LR[XCOORD];
  pos_UGUARD_R1[YCOORD] = pos_UGUARD_U[YCOORD] -(size_UGUARD_UD[YCOORD] - size_UGUARD_R1[YCOORD]);
  pos_UGUARD_R1[ZCOORD] = 0.0*mm;

  pos_UGUARD_R2U[XCOORD] = pos_MFIELD[XCOORD]- 820.0*mm + size_UGUARD_LR[XCOORD];
  pos_UGUARD_R2U[YCOORD] = pos_UGUARD_U[YCOORD] -(size_UGUARD_UD[YCOORD] - size_UGUARD_R2[YCOORD]);
  pos_UGUARD_R2U[ZCOORD] = size_UGUARD_R1[ZCOORD]+size_UGUARD_R2[ZCOORD];

  pos_UGUARD_R2D[XCOORD] = pos_MFIELD[XCOORD]- 820.0*mm + size_UGUARD_LR[XCOORD];
  pos_UGUARD_R2D[YCOORD] = pos_UGUARD_U[YCOORD] -(size_UGUARD_UD[YCOORD] - size_UGUARD_R2[YCOORD]);
  pos_UGUARD_R2D[ZCOORD] = -(size_UGUARD_R1[ZCOORD]+size_UGUARD_R2[ZCOORD]);

  pos_DGUARD_U[XCOORD] = pos_MFIELD[XCOORD]+ 820.0*mm - size_DGUARD_UD[XCOORD];
  pos_DGUARD_U[YCOORD] = pos_MFIELD[YCOORD];
  pos_DGUARD_U[ZCOORD] = size_DGUARD_LR[ZCOORD]+size_DGUARD_UD[ZCOORD];

  pos_DGUARD_D[XCOORD] = pos_MFIELD[XCOORD]+ 820.0*mm - size_DGUARD_UD[XCOORD];
  pos_DGUARD_D[YCOORD] = pos_MFIELD[YCOORD];
  pos_DGUARD_D[ZCOORD] = -(size_DGUARD_LR[ZCOORD]+size_DGUARD_UD[ZCOORD]);

  pos_DGUARD_L[XCOORD] = pos_MFIELD[XCOORD]+ 820.0*mm - size_DGUARD_LR[XCOORD];
  pos_DGUARD_L[YCOORD] = pos_MFIELD[YCOORD] + (550.0*mm + size_DGUARD_LR[YCOORD]);
  pos_DGUARD_L[ZCOORD] = 0.0*mm;

  pos_DGUARD_R[XCOORD] = pos_MFIELD[XCOORD]+ 820.0*mm - size_DGUARD_LR[XCOORD];
  pos_DGUARD_R[YCOORD] = pos_MFIELD[YCOORD]- (550.0*mm + size_DGUARD_LR[YCOORD]);
  pos_DGUARD_R[ZCOORD] = 0.0*mm;


  G4double pos_COIL1[3];
  pos_COIL1[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-(size_COIL1[XCOORD]+20.*mm);
  pos_COIL1[YCOORD] = pos_MFIELD[YCOORD];
  pos_COIL1[ZCOORD] = size_MFIELD[ZCOORD]+size_YOKE_UD[ZCOORD]*2. -(size_COIL1[ZCOORD]+20.*mm);

  G4double pos_COIL4L[3];
  pos_COIL4L[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL4L[YCOORD] = pos_MFIELD[YCOORD]+size_MFIELD[YCOORD]+size_COIL4[YCOORD];
  pos_COIL4L[ZCOORD] = size_MFIELD[ZCOORD]/2;

  G4double pos_COIL4R[3];
  pos_COIL4R[XCOORD] = pos_MFIELD[XCOORD];
  pos_COIL4R[YCOORD] = pos_MFIELD[YCOORD]-size_MFIELD[YCOORD]-size_COIL4[YCOORD];
  pos_COIL4R[ZCOORD] = size_MFIELD[ZCOORD]/2;

  G4double pos_COIL2L[3];
  pos_COIL2L[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-(size_COIL1[XCOORD]+20.*mm);
  pos_COIL2L[YCOORD] = pos_MFIELD[YCOORD]+size_MFIELD[YCOORD]+size_COIL4[YCOORD];
  pos_COIL2L[ZCOORD] = (pos_COIL4L[ZCOORD]+pos_COIL1[ZCOORD])/2+(size_COIL4[ZCOORD]-size_COIL1[ZCOORD])/2;


  G4double pos_COIL2R[3];
  pos_COIL2R[XCOORD] = pos_MFIELD[XCOORD]-size_MFIELD[XCOORD]-(size_COIL1[XCOORD]+20.*mm);
  pos_COIL2R[YCOORD] = pos_MFIELD[YCOORD]-size_MFIELD[YCOORD]-size_COIL4[YCOORD];
  pos_COIL2R[ZCOORD] = (pos_COIL4R[ZCOORD]+pos_COIL1[ZCOORD])/2+(size_COIL4[ZCOORD]-size_COIL1[ZCOORD])/2;


  G4double pos_COIL5[3];
  pos_COIL5[XCOORD] = pos_MFIELD[XCOORD] + size_MFIELD[XCOORD]+(size_COIL5[XCOORD]+21.*mm);
  pos_COIL5[YCOORD] = pos_MFIELD[YCOORD];
  pos_COIL5[ZCOORD] = size_MFIELD[ZCOORD]+size_YOKE_UD[ZCOORD]*2. -(size_COIL5[ZCOORD]);



  G4double pos_COIL3L[3];
  pos_COIL3L[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+(size_COIL5[XCOORD]+21.*mm);
  pos_COIL3L[YCOORD] = pos_MFIELD[YCOORD]+size_MFIELD[YCOORD]+size_COIL4[YCOORD];
  pos_COIL3L[ZCOORD] = (pos_COIL4L[ZCOORD]+pos_COIL5[ZCOORD])/2+(size_COIL4[ZCOORD]-size_COIL5[ZCOORD])/2;


  G4double pos_COIL3R[3];
  pos_COIL3R[XCOORD] = pos_MFIELD[XCOORD]+size_MFIELD[XCOORD]+(size_COIL5[XCOORD]+21.*mm);
  pos_COIL3R[YCOORD] = pos_MFIELD[YCOORD]-size_MFIELD[YCOORD]-size_COIL4[YCOORD];
  pos_COIL3R[ZCOORD] = (pos_COIL4R[ZCOORD]+pos_COIL5[ZCOORD])/2+(size_COIL4[ZCOORD]-size_COIL5[ZCOORD])/2;


  // Construct KURAMA Magnet
#if 0
  //-------------------- Upstraam End Guard
  G4Box* upGuard_UD_box = new G4Box("upGuard_UD_box",
				    size_UGUARD_UD[XCOORD],size_UGUARD_UD[YCOORD],size_UGUARD_UD[ZCOORD]);
  G4LogicalVolume* upGuard_U_log = new G4LogicalVolume(upGuard_UD_box,
						    mList_->Fe,
						    "upGuard_U_log",0,0,0);
  upGuard_U_log->SetVisAttributes(magnet_att);

  maxStep=0.00001*mm;
  //upGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));

  G4VPhysicalVolume* upGuard_U_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_U[XCOORD],pos_UGUARD_U[YCOORD],pos_UGUARD_U[ZCOORD]),
		      "upGuard_U_phys",
		      upGuard_U_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* upGuard_D_log = new G4LogicalVolume(upGuard_UD_box,
						       mList_->Fe,
						       "upGuard_D_log",0,0,0);
  upGuard_D_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //upGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* upGuard_D_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_D[XCOORD],pos_UGUARD_D[YCOORD],pos_UGUARD_D[ZCOORD]),
		      "upGuard_D_phys",
		      upGuard_D_log,
		      pMother,
		      false,
		      0);

  G4Box* upGuard_LR_box =
    new G4Box("upGuard_LR_box",
	      size_UGUARD_LR[XCOORD],size_UGUARD_LR[YCOORD],size_UGUARD_LR[ZCOORD]);

  G4LogicalVolume* upGuard_L_log =  new G4LogicalVolume(upGuard_LR_box,
						      mList_->Fe,
						      "upGuard_L_log",0,0,0);
  upGuard_L_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* upGuard_L_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_L[XCOORD],pos_UGUARD_L[YCOORD],pos_UGUARD_L[ZCOORD]),
		      "upGuard_L_phys",
		      upGuard_L_log,
		      pMother,
		      false,
		      0);


  /*
  G4LogicalVolume* upGuard_R_log = new G4LogicalVolume(upGuard_LR_box,
						     mList_->Fe,
						     "upGuard_R_log",0,0,0);
  upGuard_R_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //upGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* upGuard_R_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_R[XCOORD],pos_UGUARD_R[YCOORD],pos_UGUARD_R[ZCOORD]),
		      "upGuard_R_phys",
		      upGuard_R_log,
		      pMother,
		      false,
		      0);
  */

  G4Box* upGuard_R1_box =
    new G4Box("upGuard_R1_box",
	      size_UGUARD_R1[XCOORD],size_UGUARD_R1[YCOORD],size_UGUARD_R1[ZCOORD]);

  G4Box* upGuard_R2_box =
    new G4Box("upGuard_R2_box",
	      size_UGUARD_R2[XCOORD],size_UGUARD_R2[YCOORD],size_UGUARD_R2[ZCOORD]);


  G4LogicalVolume* upGuard_R1_log = new G4LogicalVolume(upGuard_R1_box,
							mList_->Fe,
							"upGuard_R1_log",0,0,0);
  upGuard_R1_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //upGuard_R1_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* upGuard_R1_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_R1[XCOORD],pos_UGUARD_R1[YCOORD],pos_UGUARD_R1[ZCOORD]),
		      "upGuard_R1_phys",
		      upGuard_R1_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* upGuard_R2U_log = new G4LogicalVolume(upGuard_R2_box,
							mList_->Fe,
							"upGuard_R2U_log",0,0,0);
  upGuard_R2U_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //upGuard_R2U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* upGuard_R2U_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_R2U[XCOORD],pos_UGUARD_R2U[YCOORD],pos_UGUARD_R2U[ZCOORD]),
		      "upGuard_R2U_phys",
		      upGuard_R2U_log,
		      pMother,
		      false,
		      0);


  G4LogicalVolume* upGuard_R2D_log = new G4LogicalVolume(upGuard_R2_box,
							mList_->Fe,
							"upGuard_R2D_log",0,0,0);
  upGuard_R2D_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //upGuard_R2U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* upGuard_R2D_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_UGUARD_R2D[XCOORD],pos_UGUARD_R2D[YCOORD],pos_UGUARD_R2D[ZCOORD]),
		      "upGuard_R2D_phys",
		      upGuard_R2D_log,
		      pMother,
		      false,
		      0);


#endif

  //////////////coil1U
  G4Box* Coil1_box = new G4Box("Coil1_box",
			       size_COIL1[XCOORD],size_COIL1[YCOORD],size_COIL1[ZCOORD]);
  G4LogicalVolume*  Coil1_log = new G4LogicalVolume(Coil1_box, mList_->Cu, "Coil1_log",0,0,0);
  Coil1_log->SetVisAttributes(red);


  G4PVPlacement* Coil1U_phys = new G4PVPlacement(0,
						 G4ThreeVector(pos_COIL1[XCOORD],pos_COIL1[YCOORD],pos_COIL1[ZCOORD]),
						 "Coil1U_phys",
						 Coil1_log,
						 pMother,
						 false,
						 0);
  G4PVPlacement* Coil1D_phys = new G4PVPlacement(0,
						 G4ThreeVector(pos_COIL1[XCOORD],pos_COIL1[YCOORD], -pos_COIL1[ZCOORD]),
						 "Coil1D_phys",
						 Coil1_log,
						 pMother,
						 false,
						 0);



  //////////////coil4RLUD
  G4Box* Coil4_box = new G4Box("Coil4_box",
				    size_COIL4[XCOORD],size_COIL4[YCOORD],size_COIL4[ZCOORD]);
  G4LogicalVolume*  Coil4_log = new G4LogicalVolume(Coil4_box, mList_->Cu, "Coil4_log",0,0,0);
  Coil4_log->SetVisAttributes(red);


  G4PVPlacement* Coil4UR_phys = new G4PVPlacement(0,
						  G4ThreeVector(pos_COIL4L[XCOORD],pos_COIL4L[YCOORD],pos_COIL4L[ZCOORD]),
						  "Coil4UR_phys",
						  Coil4_log,
						  pMother,
						  false,
						  0);
  G4PVPlacement* Coil4UL_phys = new G4PVPlacement(0,
						  G4ThreeVector(pos_COIL4R[XCOORD],pos_COIL4R[YCOORD],pos_COIL4R[ZCOORD]),
						  "Coil4UL_phys",
						  Coil4_log,
						  pMother,
						  false,
						  0);
  G4PVPlacement* Coil4DR_phys = new G4PVPlacement(0,
						  G4ThreeVector(pos_COIL4L[XCOORD],pos_COIL4L[YCOORD],-pos_COIL4L[ZCOORD]),
						  "Coil4DR_phys",
						  Coil4_log,
						  pMother,
						  false,
						  0);
  G4PVPlacement* Coil4DL_phys = new G4PVPlacement(0,
						  G4ThreeVector(pos_COIL4R[XCOORD],pos_COIL4R[YCOORD],-pos_COIL4R[ZCOORD]),
						  "Coil4UL_phys",
						  Coil4_log,
						  pMother,
						  false,
						  0);


  //////////////coil5UD
  G4Box* Coil5_box = new G4Box("Coil5_box",
				    size_COIL5[XCOORD],size_COIL5[YCOORD],size_COIL5[ZCOORD]);
  G4LogicalVolume*  Coil5_log = new G4LogicalVolume(Coil5_box, mList_->Cu, "Coil5_log",0,0,0);
  Coil5_log->SetVisAttributes(red);


  G4PVPlacement* Coil5U_phys = new G4PVPlacement(0,
						 G4ThreeVector(pos_COIL5[XCOORD],pos_COIL5[YCOORD],pos_COIL5[ZCOORD]),
						 "Coil5U_phys",
						 Coil5_log,
						 pMother,
						 false,
						 0);
  G4PVPlacement* Coil5D_phys = new G4PVPlacement(0,
						 G4ThreeVector(pos_COIL5[XCOORD],pos_COIL5[YCOORD],-pos_COIL5[ZCOORD]),
						 "Coil5D_phys",
						 Coil5_log,
						 pMother,
						 false,
						 0);

  //////////////coil6RLUD
  G4double size_COIL6[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL6[0] = 50.0*mm;
  size_COIL6[1] = size_COIL1[ZCOORD]*2+size_COIL6[0];
  //  size_COIL6[1] = 330.*mm;
  size_COIL6[2] = 280./2*mm;
  size_COIL6[3] = 90.*deg;

  G4double pos_COIL6LU[3];
  G4double pos_COIL6RU[3];
  G4double pos_COIL6LD[3];
  G4double pos_COIL6RD[3];
  //LU
  pos_COIL6LU[XCOORD] = pos_MFIELD[XCOORD] - size_MFIELD[XCOORD]-(size_COIL6[2]+21.*mm);
  pos_COIL6LU[YCOORD] = pos_MFIELD[YCOORD] +size_MFIELD[YCOORD]-size_COIL6[0];
  pos_COIL6LU[ZCOORD] = pos_COIL1[ZCOORD]  -(size_COIL6[0]+size_COIL1[ZCOORD]);

  //RU
  pos_COIL6RU[XCOORD] = pos_MFIELD[XCOORD] - size_MFIELD[XCOORD]-(size_COIL6[2]+21.*mm);
  pos_COIL6RU[YCOORD] = pos_MFIELD[YCOORD] -size_MFIELD[YCOORD]+size_COIL6[0];
  pos_COIL6RU[ZCOORD] = pos_COIL1[ZCOORD]  -(size_COIL6[0]+size_COIL1[ZCOORD]);

  //LD
  pos_COIL6LD[XCOORD] = pos_MFIELD[XCOORD]  - size_MFIELD[XCOORD]-(size_COIL6[2]+21.*mm);
  pos_COIL6LD[YCOORD] = pos_MFIELD[YCOORD]  +size_MFIELD[YCOORD]-size_COIL6[0];
  pos_COIL6LD[ZCOORD] = -pos_COIL1[ZCOORD] +(size_COIL6[0]+size_COIL1[ZCOORD]);

  //RD
  pos_COIL6RD[XCOORD] = pos_MFIELD[XCOORD]  - size_MFIELD[XCOORD]-(size_COIL6[2]+21.*mm);
  pos_COIL6RD[YCOORD] = pos_MFIELD[YCOORD]  -size_MFIELD[YCOORD]+size_COIL6[0];
  pos_COIL6RD[ZCOORD] = -pos_COIL1[ZCOORD] +(size_COIL6[0]+size_COIL1[ZCOORD]);



  G4Tubs* Coil6_tub = new G4Tubs("Coil6_tubs",
				size_COIL6[0],size_COIL6[1],size_COIL6[2],0.,size_COIL6[3]);
  G4LogicalVolume*  Coil6_log = new G4LogicalVolume(Coil6_tub, mList_->Cu, "Coil6_log",0,0,0);
  Coil6_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix *rotcoil6lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil6rd = new G4RotationMatrix();
  /*
  rotcoil6lu->rotateY(-fSpectrometerAngle);
  rotcoil6ru->rotateY(-fSpectrometerAngle);
  rotcoil6ld->rotateY(-fSpectrometerAngle);
  rotcoil6rd->rotateY(-fSpectrometerAngle);
  */

  rotcoil6lu->rotateY(90.*deg);
  rotcoil6lu->rotateX(0.*deg);

  G4PVPlacement* Coil6LU_phys
    = new G4PVPlacement(rotcoil6lu,
			G4ThreeVector(pos_COIL6LU[XCOORD],
				      pos_COIL6LU[YCOORD],
				      pos_COIL6LU[ZCOORD]),
			"Coil6LU_phys",
			Coil6_log,
			pMother,
			false,
			0);

  rotcoil6ru->rotateY(90.*deg);
  rotcoil6ru->rotateZ(90.*deg);
  G4PVPlacement* Coil6RU_phys
    = new G4PVPlacement(rotcoil6ru,
			G4ThreeVector(pos_COIL6RU[XCOORD],
				      pos_COIL6RU[YCOORD],
				      pos_COIL6RU[ZCOORD]),
			"Coil6RU_phys",
			Coil6_log,
			pMother,
			false,
			0);

  rotcoil6ld->rotateY(90.*deg);
  rotcoil6ld->rotateZ(-180.*deg);
  G4PVPlacement* Coil6RD_phys
    = new G4PVPlacement(rotcoil6ld,
			G4ThreeVector(pos_COIL6RD[XCOORD],
				      pos_COIL6RD[YCOORD],
				      pos_COIL6RD[ZCOORD]),
			"Coil6RD_phys",
			Coil6_log,
			pMother,
			false,
			0);

  rotcoil6rd->rotateY(90.*deg);
  rotcoil6rd->rotateZ(270.*deg);
  G4PVPlacement* Coil6LD_phys
    = new G4PVPlacement(rotcoil6rd,
			G4ThreeVector(pos_COIL6LD[XCOORD],
				      pos_COIL6LD[YCOORD],
				      pos_COIL6LD[ZCOORD]),
			"Coil6LD_phys",
			Coil6_log,
			pMother,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);


  //////////////coil8RLUD
  G4double size_COIL8[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL8[0] = 50.0*mm;
  size_COIL8[1] = size_COIL5[ZCOORD]*2.+size_COIL8[0];
  size_COIL8[2] = 280./2*mm;
  size_COIL8[3] = 90.*deg;

  G4double pos_COIL8LU[3];
  G4double pos_COIL8RU[3];
  G4double pos_COIL8LD[3];
  G4double pos_COIL8RD[3];
  //LU
  pos_COIL8LU[XCOORD] = pos_MFIELD[XCOORD] + size_MFIELD[XCOORD]+(size_COIL8[2]+21.*mm);
  pos_COIL8LU[YCOORD] = pos_MFIELD[YCOORD] +size_MFIELD[YCOORD]-size_COIL8[0];
  pos_COIL8LU[ZCOORD] = pos_COIL5[ZCOORD]  -(size_COIL8[0]+size_COIL5[ZCOORD]);

  //RU
  pos_COIL8RU[XCOORD] = pos_MFIELD[XCOORD] + size_MFIELD[XCOORD]+(size_COIL8[2]+21.*mm);
  pos_COIL8RU[YCOORD] = pos_MFIELD[YCOORD] -size_MFIELD[YCOORD]+size_COIL8[0];
  pos_COIL8RU[ZCOORD] = pos_COIL5[ZCOORD]  -(size_COIL8[0]+size_COIL5[ZCOORD]);

  //LD
  pos_COIL8LD[XCOORD] = pos_MFIELD[XCOORD]  + size_MFIELD[XCOORD]+(size_COIL8[2]+21.*mm);
  pos_COIL8LD[YCOORD] = pos_MFIELD[YCOORD]  +size_MFIELD[YCOORD]-size_COIL8[0];
  pos_COIL8LD[ZCOORD] = -pos_COIL5[ZCOORD] +(size_COIL8[0]+size_COIL5[ZCOORD]);

  //RD
  pos_COIL8RD[XCOORD] = pos_MFIELD[XCOORD]  + size_MFIELD[XCOORD]+(size_COIL8[2]+21.*mm);
  pos_COIL8RD[YCOORD] = pos_MFIELD[YCOORD]  -size_MFIELD[YCOORD]+size_COIL8[0];
  pos_COIL8RD[ZCOORD] = -pos_COIL5[ZCOORD] +(size_COIL8[0]+size_COIL5[ZCOORD]);


  G4Tubs* Coil8_tub = new G4Tubs("Coil8_tubs",
				size_COIL8[0],size_COIL8[1],size_COIL8[2],0.,size_COIL8[3]);
  G4LogicalVolume*  Coil8_log = new G4LogicalVolume(Coil8_tub, mList_->Cu, "Coil8_log",0,0,0);
  Coil8_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix *rotcoil8lu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8ld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil8rd = new G4RotationMatrix();
  /*
  rotcoil8lu->rotateY(-fSpectrometerAngle);
  rotcoil8ru->rotateY(-fSpectrometerAngle);
  rotcoil8ld->rotateY(-fSpectrometerAngle);
  rotcoil8rd->rotateY(-fSpectrometerAngle);
  */
  rotcoil8lu->rotateY(90.*deg);
  rotcoil8lu->rotateZ(0.*deg);
  G4PVPlacement* Coil8LU_phys
    = new G4PVPlacement(rotcoil8lu,
			G4ThreeVector(pos_COIL8LU[XCOORD],
				      pos_COIL8LU[YCOORD],
				      pos_COIL8LU[ZCOORD]),
			"Coil8LU_phys",
			Coil8_log,
			pMother,
			false,
			0);
  rotcoil8ru->rotateY(90.*deg);
  rotcoil8ru->rotateZ(90.*deg);
  G4PVPlacement* Coil8RU_phys
    = new G4PVPlacement(rotcoil8ru,
			G4ThreeVector(pos_COIL8RU[XCOORD],
				      pos_COIL8RU[YCOORD],
				      pos_COIL8RU[ZCOORD]),
			"Coil8RU_phys",
			Coil8_log,
			pMother,
			false,
			0);

  rotcoil8ld->rotateY(90.*deg);
  rotcoil8ld->rotateZ(-180.*deg);
  G4PVPlacement* Coil8RD_phys
    = new G4PVPlacement(rotcoil8ld,
			G4ThreeVector(pos_COIL8RD[XCOORD],
				      pos_COIL8RD[YCOORD],
				      pos_COIL8RD[ZCOORD]),
			"Coil8RD_phys",
			Coil8_log,
			pMother,
			false,
			0);
  rotcoil8rd->rotateY(90.*deg);
  rotcoil8rd->rotateZ(270.*deg);
  G4PVPlacement* Coil8LD_phys
    = new G4PVPlacement(rotcoil8rd,
			G4ThreeVector(pos_COIL8LD[XCOORD],
				      pos_COIL8LD[YCOORD],
				      pos_COIL8LD[ZCOORD]),
			"Coil8LD_phys",
			Coil8_log,
			pMother,
			false,
			0);
  //  rotForwardSp->rotateZ(90.*deg);
  //  rotForwardSp->rotateZ(90.*deg);


  //////////////coil7RLUD
  G4double size_COIL7[4];
  //0:in
  //1:out
  //2:z
  //3:angle
  size_COIL7[0] = 50.0*mm;
  size_COIL7[1] = size_COIL4[ZCOORD]*2.+size_COIL7[0];
  size_COIL7[2] = size_COIL4[YCOORD];
  size_COIL7[3] = 90.*deg;

  G4double pos_COIL7ULU[3];
  G4double pos_COIL7URU[3];
  G4double pos_COIL7ULD[3];
  G4double pos_COIL7URD[3];
  G4double pos_COIL7DLU[3];
  G4double pos_COIL7DRU[3];
  G4double pos_COIL7DLD[3];
  G4double pos_COIL7DRD[3];
  //ULU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] - size_COIL4[XCOORD];
  pos_COIL7ULU[YCOORD] = pos_COIL4L[YCOORD];
  pos_COIL7ULU[ZCOORD] = pos_COIL4L[ZCOORD] +(size_COIL7[0]+size_COIL4[ZCOORD]);

  //URU
  pos_COIL7URU[XCOORD] = pos_MFIELD[XCOORD] - size_COIL4[XCOORD];
  pos_COIL7URU[YCOORD] = pos_COIL4R[YCOORD];
  pos_COIL7URU[ZCOORD] = pos_COIL4R[ZCOORD] +(size_COIL7[0]+size_COIL4[ZCOORD]);

  //ULD
  pos_COIL7ULD[XCOORD] = pos_MFIELD[XCOORD]  - size_COIL4[XCOORD];
  pos_COIL7ULD[YCOORD] = pos_COIL4L[YCOORD];
  pos_COIL7ULD[ZCOORD] = -pos_COIL4L[ZCOORD] -(size_COIL7[0]+size_COIL4[ZCOORD]);

  //URD
  pos_COIL7URD[XCOORD] = pos_MFIELD[XCOORD] - size_COIL4[XCOORD];
  pos_COIL7URD[YCOORD] = pos_COIL4R[YCOORD];
  pos_COIL7URD[ZCOORD] = -pos_COIL4R[ZCOORD] -(size_COIL7[0]+size_COIL4[ZCOORD]);



  //DLU
  //  pos_COIL7ULU[XCOORD] = pos_MFIELD[XCOORD] +size_MFIELD[XCOORD]+size_COIL7[0];
  pos_COIL7DLU[XCOORD] = pos_MFIELD[XCOORD] + size_COIL4[XCOORD];
  pos_COIL7DLU[YCOORD] = pos_COIL4L[YCOORD];
  pos_COIL7DLU[ZCOORD] = pos_COIL4L[ZCOORD] +(size_COIL7[0]+size_COIL4[ZCOORD]);

  //DRU
  pos_COIL7DRU[XCOORD] = pos_MFIELD[XCOORD] + size_COIL4[XCOORD];
  pos_COIL7DRU[YCOORD] = pos_COIL4R[YCOORD];
  pos_COIL7DRU[ZCOORD] = pos_COIL4R[ZCOORD] +(size_COIL7[0]+size_COIL4[ZCOORD]);

  //DLD
  pos_COIL7DLD[XCOORD] = pos_MFIELD[XCOORD]  + size_COIL4[XCOORD];
  pos_COIL7DLD[YCOORD] = pos_COIL4L[YCOORD];
  pos_COIL7DLD[ZCOORD] = -pos_COIL4L[ZCOORD] -(size_COIL7[0]+size_COIL4[ZCOORD]);

  //DRD
  pos_COIL7DRD[XCOORD] = pos_MFIELD[XCOORD] + size_COIL4[XCOORD];
  pos_COIL7DRD[YCOORD] = pos_COIL4R[YCOORD];
  pos_COIL7DRD[ZCOORD] = -pos_COIL4R[ZCOORD] -(size_COIL7[0]+size_COIL4[ZCOORD]);



  G4Tubs* Coil7_tub = new G4Tubs("Coil7_tubs",
				size_COIL7[0],size_COIL7[1],size_COIL7[2],0.,size_COIL7[3]);
  G4LogicalVolume*  Coil7_log = new G4LogicalVolume(Coil7_tub, mList_->Cu, "Coil7_log",0,0,0);
  Coil7_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4RotationMatrix *rotcoil7ulu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7uld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7urd = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dlu = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dru = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7dld = new G4RotationMatrix();
  G4RotationMatrix *rotcoil7drd = new G4RotationMatrix();
  /*
  rotcoil7ulu->rotateY(-fSpectrometerAngle);
  rotcoil7uru->rotateY(-fSpectrometerAngle);
  rotcoil7uld->rotateY(-fSpectrometerAngle);
  rotcoil7urd->rotateY(-fSpectrometerAngle);
  rotcoil7dlu->rotateY(-fSpectrometerAngle);
  rotcoil7dru->rotateY(-fSpectrometerAngle);
  rotcoil7dld->rotateY(-fSpectrometerAngle);
  rotcoil7drd->rotateY(-fSpectrometerAngle);
  */

  rotcoil7ulu->rotateX(90.*deg);
  rotcoil7ulu->rotateZ(-90.*deg);
  /*
  rotcoil7ulu->rotateY(90.*deg);
  rotcoil7ulu->rotateX(0.*deg);
  rotcoil7ulu->rotateY(180.*deg);
  rotcoil7ulu->rotateX(90.*deg);
  */
  G4PVPlacement* Coil7ULU_phys
    = new G4PVPlacement(rotcoil7ulu,
			G4ThreeVector(pos_COIL7ULU[XCOORD],
				      pos_COIL7ULU[YCOORD],
				      pos_COIL7ULU[ZCOORD]),
			"Coil7ULU_phys",
			Coil7_log,
			pMother,
			false,
			0);



  rotcoil7uru->rotateX(90.*deg);
  rotcoil7uru->rotateZ(-90.*deg);

  /*
  rotcoil7ulu->rotateY(90.*deg);
  rotcoil7uru->rotateX(0.*deg);
  rotcoil7uru->rotateY(180.*deg);
  rotcoil7uru->rotateX(90.*deg);
  */
  G4PVPlacement* Coil7URU_phys
    = new G4PVPlacement(rotcoil7uru,
			G4ThreeVector(pos_COIL7URU[XCOORD],
				      pos_COIL7URU[YCOORD],
				      pos_COIL7URU[ZCOORD]),
			"Coil7URU_phys",
			Coil7_log,
			pMother,
			false,
			0);

  rotcoil7urd->rotateX(-90.*deg);
  rotcoil7urd->rotateZ(-90.*deg);

  /*
  rotcoil7urd->rotateZ(0.*deg);
  rotcoil7urd->rotateX(90.*deg);
  rotcoil7urd->rotateY(90.*deg);
  */
  G4PVPlacement* Coil7URD_phys
    = new G4PVPlacement(rotcoil7urd,
			G4ThreeVector(pos_COIL7URD[XCOORD],
				      pos_COIL7URD[YCOORD],
				      pos_COIL7URD[ZCOORD]),
			"Coil7URD_phys",
			Coil7_log,
			pMother,
			false,
			0);
  rotcoil7uld->rotateX(-90.*deg);
  rotcoil7uld->rotateZ(-90.*deg);
  /*
  rotcoil7uld->rotateZ(0.*deg);
  rotcoil7uld->rotateX(90.*deg);
  rotcoil7uld->rotateY(90.*deg);
  */
  G4PVPlacement* Coil7ULD_phys
    = new G4PVPlacement(rotcoil7uld,
			G4ThreeVector(pos_COIL7ULD[XCOORD],
				      pos_COIL7ULD[YCOORD],
				      pos_COIL7ULD[ZCOORD]),
			"Coil7ULD_phys",
			Coil7_log,
			pMother,
			false,
			0);


  ///down

  /*
  rotcoil7dlu->rotateZ(0.*deg);
  rotcoil7dlu->rotateX(-90.*deg);
  rotcoil7dlu->rotateY(90.*deg);
  */
  rotcoil7dlu->rotateX(90.*deg);
  G4PVPlacement* Coil7DLU_phys
    = new G4PVPlacement(rotcoil7dlu,
			G4ThreeVector(pos_COIL7DLU[XCOORD],
				      pos_COIL7DLU[YCOORD],
				      pos_COIL7DLU[ZCOORD]),
			"Coil7DLU_phys",
			Coil7_log,
			pMother,
			false,
			0);

  /*
  rotcoil7dru->rotateZ(0.*deg);
  rotcoil7dru->rotateX(-90.*deg);
  rotcoil7dru->rotateY(90.*deg);
  */
  rotcoil7dru->rotateX(90.*deg);
  G4PVPlacement* Coil7DRU_phys
    = new G4PVPlacement(rotcoil7dru,
			G4ThreeVector(pos_COIL7DRU[XCOORD],
				      pos_COIL7DRU[YCOORD],
				      pos_COIL7DRU[ZCOORD]),
			"Coil7DRU_phys",
			Coil7_log,
			pMother,
			false,
			0);

  /*
  rotcoil7drd->rotateZ(0.*deg);
  rotcoil7drd->rotateX(0.*deg);
  rotcoil7drd->rotateY(90.*deg);
  */

  rotcoil7drd->rotateX(-90.*deg);
  G4PVPlacement* Coil7DRD_phys
    = new G4PVPlacement(rotcoil7drd,
			G4ThreeVector(pos_COIL7DRD[XCOORD],
				      pos_COIL7DRD[YCOORD],
				      pos_COIL7DRD[ZCOORD]),
			"Coil7DRD_phys",
			Coil7_log,
			pMother,
			false,
			0);
  /*
  rotcoil7dld->rotateZ(0.*deg);
  rotcoil7dld->rotateX(0.*deg);
  rotcoil7dld->rotateY(90.*deg);
  */

  rotcoil7dld->rotateX(-90.*deg);
  G4PVPlacement* Coil7DLD_phys
    = new G4PVPlacement(rotcoil7dld,
			G4ThreeVector(pos_COIL7DLD[XCOORD],
				      pos_COIL7DLD[YCOORD],
				      pos_COIL7DLD[ZCOORD]),
			"Coil7DLD_phys",
			Coil7_log,
			pMother,
			false,
			0);


  ///coil2
  //////////////coil2
  G4Box* Coil2_box = new G4Box("Coil2_box",
				    size_COIL2[XCOORD],size_COIL2[YCOORD],size_COIL2[ZCOORD]);
  G4LogicalVolume*  Coil2_log = new G4LogicalVolume(Coil2_box, mList_->Cu, "Coil2_log",0,0,0);
  Coil2_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4PVPlacement* Coil2UL_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL2L[XCOORD],pos_COIL2L[YCOORD],pos_COIL2L[ZCOORD]),
		      "Coil2UL_phys",
		      Coil2_log,
		      pMother,
		      false,
		      0);

  G4PVPlacement* Coil2UR_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL2R[XCOORD],pos_COIL2R[YCOORD],pos_COIL2R[ZCOORD]),
		      "Coil2UR_phys",
		      Coil2_log,
		      pMother,
		      false,
		      0);


  G4PVPlacement* Coil2DL_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL2L[XCOORD],pos_COIL2L[YCOORD],-pos_COIL2L[ZCOORD]),
		      "Coil2DL_phys",
		      Coil2_log,
		      pMother,
		      false,
		      0);

  G4PVPlacement* Coil2DR_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL2R[XCOORD],pos_COIL2R[YCOORD],-pos_COIL2R[ZCOORD]),
		      "Coil2DR_phys",
		      Coil2_log,
		      pMother,
		      false,
		      0);

  ///coil3
  //////////////coil3
  G4Box* Coil3_box = new G4Box("Coil3_box",
				    size_COIL3[XCOORD],size_COIL3[YCOORD],size_COIL3[ZCOORD]);
  G4LogicalVolume*  Coil3_log = new G4LogicalVolume(Coil3_box, mList_->Cu, "Coil3_log",0,0,0);
  Coil3_log->SetVisAttributes(red);
  maxStep=0.00001*mm;
  //upGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));

  G4PVPlacement* Coil3UL_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL3L[XCOORD],pos_COIL3L[YCOORD],pos_COIL3L[ZCOORD]),
		      "Coil3UL_phys",
		      Coil3_log,
		      pMother,
		      false,
		      0);

  G4PVPlacement* Coil3UR_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL3R[XCOORD],pos_COIL3R[YCOORD],pos_COIL3R[ZCOORD]),
		      "Coil3UR_phys",
		      Coil3_log,
		      pMother,
		      false,
		      0);


  G4PVPlacement* Coil3DL_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL3L[XCOORD],pos_COIL3L[YCOORD],-pos_COIL3L[ZCOORD]),
		      "Coil3DL_phys",
		      Coil3_log,
		      pMother,
		      false,
		      0);

  G4PVPlacement* Coil3DR_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_COIL3R[XCOORD],pos_COIL3R[YCOORD],-pos_COIL3R[ZCOORD]),
		      "Coil3DR_phys",
		      Coil3_log,
		      pMother,
		      false,
		      0);


  //-----------------Magetic Field
  /*
  G4Box* MField_box = new G4Box("MField_box",
				size_MFIELD[XCOORD],size_MFIELD[YCOORD],size_MFIELD[ZCOORD]);
  MField_log = new G4LogicalVolume(MField_box, Air, "MField_log",0,0,0);
  maxStep = 1.0*mm;
  //MField_log->SetUserLimits(new G4UserLimits(maxStep));
  MField_log->SetVisAttributes(magnet_att);
  MField_phys = new G4PVPlacement(0,
				  G4ThreeVector(pos_MFIELD[XCOORD],pos_MFIELD[YCOORD],pos_MFIELD[ZCOORD]),
				  MField_log,
				  "MField_phys",
				  experimentalHall_log,
				  false,
				  0);
  */
  //-------------------- Yoke
  G4Box* Yoke_UD_box = new G4Box("Yoke_UD_box",
				 size_YOKE_UD[XCOORD],size_YOKE_UD[YCOORD],size_YOKE_UD[ZCOORD]);
  G4LogicalVolume* Yoke_U_log = new G4LogicalVolume(Yoke_UD_box,
						    mList_->Fe,
						    "Yoke_U_log",0,0,0);
  Yoke_U_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //Yoke_U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume *Yoke_U_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_YOKE_U[XCOORD],pos_YOKE_U[YCOORD],pos_YOKE_U[ZCOORD]),
		      "Yoke_U_phys",
		      Yoke_U_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* Yoke_D_log = new G4LogicalVolume(Yoke_UD_box,
						    mList_->Fe,
						    "Yoke_D_log",0,0,0);
  Yoke_D_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //Yoke_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume *Yoke_D_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_YOKE_D[XCOORD],pos_YOKE_D[YCOORD],pos_YOKE_D[ZCOORD]),
		      "Yoke_D_phys",
		      Yoke_D_log,
		      pMother,
		      false,
		      0);

  G4Box* Yoke_LR_box = new G4Box("Yoke_LR_box",
				 size_YOKE_LR[XCOORD],size_YOKE_LR[YCOORD],size_YOKE_LR[ZCOORD]);
  G4LogicalVolume* Yoke_L_log = new G4LogicalVolume(Yoke_LR_box,
						    mList_->Fe,
						    "Yoke_L_log",0,0,0);
  Yoke_L_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_L_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_YOKE_L[XCOORD],pos_YOKE_L[YCOORD],pos_YOKE_L[ZCOORD]),
		      "Yoke_L_phys",
		      Yoke_L_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* Yoke_R_log = new G4LogicalVolume(Yoke_LR_box,
						    mList_->Fe,
						    "Yoke_R_log",0,0,0);
  Yoke_R_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //Yoke_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_R_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_YOKE_R[XCOORD],pos_YOKE_R[YCOORD],pos_YOKE_R[ZCOORD]),
		      "Yoke_R_phys",
		      Yoke_R_log,
		      pMother,
		      false,
		      0);

  G4Box* Yoke_LR_GapSpace_box = new G4Box("Yoke_LR_GapSpace_box",
					  size_YOKE_LR_GapSpace[XCOORD],
					  size_YOKE_LR_GapSpace[YCOORD],
					  size_YOKE_LR_GapSpace[ZCOORD]);
  G4LogicalVolume* Yoke_L_GapSpace_log = new G4LogicalVolume(Yoke_LR_GapSpace_box,
							     mList_->Fe,
							     "Yoke_L_GapSpace_log",0,0,0);
  Yoke_L_GapSpace_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //Yoke_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* Yoke_L_GapSpace_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_YOKE_L_GapSpace[XCOORD],pos_YOKE_L_GapSpace[YCOORD],pos_YOKE_L_GapSpace[ZCOORD]),
		      "Yoke_L_GapSpace_phys",
		      Yoke_L_GapSpace_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* Yoke_R_GapSpace_log = new G4LogicalVolume(Yoke_LR_GapSpace_box,
							     mList_->Fe,
							     "Yoke_R_GapSpace_log",0,0,0);
  Yoke_R_GapSpace_log->SetVisAttributes(magnet_att);
  G4VPhysicalVolume* Yoke_R_GapSpace_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_YOKE_R_GapSpace[XCOORD],pos_YOKE_R_GapSpace[YCOORD],pos_YOKE_R_GapSpace[ZCOORD]),
		      "Yoke_R_GapSpace_phys",
		      Yoke_R_GapSpace_log,
		      pMother,
		      false,
		      0);


#if 0
  //-------------------- Downstraam End Guard
  G4Box* downGuard_UD_box = new G4Box("downGuard_UD_box",
				      size_DGUARD_UD[XCOORD],size_DGUARD_UD[YCOORD],size_DGUARD_UD[ZCOORD]);
  G4LogicalVolume* downGuard_U_log = new G4LogicalVolume(downGuard_UD_box,
							 mList_->Fe,
							 "downGuard_U_log",0,0,0);
  downGuard_U_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //downGuard_U_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_U_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_DGUARD_U[XCOORD],pos_DGUARD_U[YCOORD],pos_DGUARD_U[ZCOORD]),
		      "downGuard_U_phys",
		      downGuard_U_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* downGuard_D_log = new G4LogicalVolume(downGuard_UD_box,
							 mList_->Fe,
							 "downGuard_D_log",0,0,0);
  downGuard_D_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //downGuard_D_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_D_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_DGUARD_D[XCOORD],pos_DGUARD_D[YCOORD],pos_DGUARD_D[ZCOORD]),
		      "downGuard_D_phys",
		      downGuard_D_log,
		      pMother,
		      false,
		      0);

  G4Box* downGuard_LR_box = new G4Box("downGuard_LR_box",
				      size_DGUARD_LR[XCOORD],size_DGUARD_LR[YCOORD],size_DGUARD_LR[ZCOORD]);
  G4LogicalVolume* downGuard_L_log = new G4LogicalVolume(downGuard_LR_box,
							 mList_->Fe,
							 "downGuard_L_log",0,0,0);
  downGuard_L_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //downGuard_L_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_L_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_DGUARD_L[XCOORD],pos_DGUARD_L[YCOORD],pos_DGUARD_L[ZCOORD]),
		      "downGuard_L_phys",
		      downGuard_L_log,
		      pMother,
		      false,
		      0);

  G4LogicalVolume* downGuard_R_log = new G4LogicalVolume(downGuard_LR_box,
						       mList_->Fe,
						       "downGuard_R_log",0,0,0);
  downGuard_R_log->SetVisAttributes(magnet_att);
  maxStep=0.00001*mm;
  //downGuard_R_log->SetUserLimits(new G4UserLimits(maxStep));
  G4VPhysicalVolume* downGuard_R_phys =
    new G4PVPlacement(0,
		      G4ThreeVector(pos_DGUARD_R[XCOORD],pos_DGUARD_R[YCOORD],pos_DGUARD_R[ZCOORD]),
		      "downGuard_R_phys",
		      downGuard_R_log,
		      pMother,
		      false,
		      0);
#endif
}

void  CFTDetectorConstruction::ConstructSdcIn(G4VPhysicalVolume *pMother)
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour green(0.0, 1.0, 0.0);
  G4Colour blue(0.0, 0.0, 1.0);
  G4VisAttributes* chamber_att = new G4VisAttributes(green);
  chamber_att->SetForceWireframe(true);
  G4VisAttributes* magnet_att = new G4VisAttributes(blue);
  magnet_att->SetForceWireframe(true);


  double size_SDC1[XYZ];

  // for LEPS
  // size_SDC1[XCOORD] = 600.0/2.*mm;
  // size_SDC1[YCOORD] = 59.7/2.*mm;
  // size_SDC1[ZCOORD] = 200.0/2.*mm;

  //for HYPS
  size_SDC1[XCOORD] = 600.0/2.*mm;
  size_SDC1[YCOORD] = 59.7/2.*mm;
  size_SDC1[ZCOORD] = 200.0/2.*mm;

  double size_SDC1Plane[XYZ];
  // original
  // for LEPS
  // size_SDC1Plane[XCOORD] = 580./2*mm;
  // size_SDC1Plane[YCOORD] = 0.0001/2.*mm;
  // size_SDC1Plane[ZCOORD] = 200./2*mm;
  size_SDC1Plane[XCOORD] = 560./2*mm;
  size_SDC1Plane[YCOORD] = 0.0001/2.*mm;
  size_SDC1Plane[ZCOORD] = 150./2*mm;

  lnum1 = geomMan.GetDetectorId("SDC1-x-1");
  G4ThreeVector posSdc1X1 = geomMan.GetGlobalPosition(lnum1);
  lnum2 = geomMan.GetDetectorId("SDC1-u-2");
  G4ThreeVector posSdc1U2 = geomMan.GetGlobalPosition(lnum2);

  G4ThreeVector offsetSDC1(0., geomMan.GetOffset(lnum1), 0.);
  posSdc1X1 += offsetSDC1;
  posSdc1U2 += offsetSDC1;

  G4RotationMatrix* rotSdc1 = new G4RotationMatrix();
  rotSdc1->rotateZ(-geomMan.GetRotAngle2(lnum1)*deg);

  G4ThreeVector SDC1Pos = G4ThreeVector((posSdc1X1.x()+posSdc1U2.x())/2.0 * mm,
					(posSdc1X1.y()+posSdc1U2.y())/2.0 * mm,
					(posSdc1X1.z()+posSdc1U2.z())/2.0 * mm);

  G4cout << "SDC1 CenterPos = ( " << SDC1Pos.x() << ", " << SDC1Pos.y()
	 << ", " << SDC1Pos.z() << ")" << G4endl;


  G4Box* SDC1_box = new G4Box("SDC1_box",size_SDC1[XCOORD],size_SDC1[YCOORD],size_SDC1[ZCOORD]);
  G4LogicalVolume* SDC1_log =
    new G4LogicalVolume(SDC1_box, mList_->ArGas, "SDC1_log",0,0,0);
  SDC1_log->SetVisAttributes(chamber_att);
  G4VPhysicalVolume* SDC1_phys = new G4PVPlacement(rotSdc1,
						   SDC1Pos,
						   "SDC1_phys",
						   SDC1_log,
						   pMother,
						   false,
						   0);

  G4LogicalVolume* SDC1Plane_log[NumOfPlaneSDC1];
  G4VPhysicalVolume* SDC1Plane_phys[NumOfPlaneSDC1];

  //---------SDC1x1
  G4Box* SDC1Plane_box = new G4Box("SDC1Plane_box",
				   size_SDC1Plane[XCOORD],
				   size_SDC1Plane[YCOORD],
				   size_SDC1Plane[ZCOORD]);
  for (int i=lnum1; i<=lnum2; i++) {
    G4ThreeVector posSdc1plane = geomMan.GetGlobalPosition(i);
    G4ThreeVector localPos;

    G4cout << "SDC " << i << " Pos = ( " << posSdc1plane.x() << ", "
	   << posSdc1plane.y() << ", " << posSdc1plane.z() << ")" << G4endl;


    switch (i-lnum1) {
    case 0:
      sprintf(name_log, "SDC1X1_log");
      sprintf(name_phys, "SDC1-x-1");
      //localPos = G4ThreeVector(0.0*mm, GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      localPos = G4ThreeVector(0.0*mm, (SDC1Pos.x()-posSdc1plane.x()*mm), 0.0*mm);
     break;
    case 1:
      sprintf(name_log, "SDC1X1_log");
      sprintf(name_phys, "SDC1-x-2");
      localPos = G4ThreeVector(0.0*mm, (SDC1Pos.x()-posSdc1plane.x()*mm), 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    case 2:
      sprintf(name_log, "SDC1U1_log");
      sprintf(name_phys, "SDC1-u-1");
      localPos = G4ThreeVector(0.0*mm, (SDC1Pos.x()-posSdc1plane.x()*mm), 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    case 3:
      sprintf(name_log, "SDC1U2_log");
      sprintf(name_phys, "SDC1-u-2");
      localPos = G4ThreeVector(0.0*mm, (SDC1Pos.x()-posSdc1plane.x()*mm), 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, -GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    default:
      fprintf(stderr, "DetectorConstruction::ConstructForwardSpectrometer()\n");
      fprintf(stderr, "No such plane %d in SDC1\n", i);
      exit(-1);
    }


    G4cout << "SDC(local) " << i << " Pos = ( " << localPos.x() << ", "
	   << localPos.y() << ", " << localPos.z() << ")" << G4endl;

    SDC1Plane_log[i-lnum1]  = new G4LogicalVolume(SDC1Plane_box,
						  mList_->ArGas,
						  name_log,
						  0,0,0);

    if ((i-lnum1)==0)
      SDC1Plane_log[i-lnum1]->SetVisAttributes(magnet_att);

    SDC1Plane_phys[i-lnum1] = new G4PVPlacement(0,
						localPos,
						SDC1Plane_log[i-lnum1],
						name_phys,
						SDC1_log,
						false,
						0);
  }


  //--------------SDC2
  /*
  double size_SDC2[XYZ];
  size_SDC2[XCOORD] = 310.0*mm;
  size_SDC2[YCOORD] = 70.0/2.*mm;
  size_SDC2[ZCOORD] = 210.0*mm;

  double size_SDC2Plane[XYZ];
  size_SDC2Plane[XCOORD] = 10.0*48.0*0.5*mm;
  size_SDC2Plane[YCOORD] = 0.0001/2.*mm;
  //size_SDC2Plane[ZCOORD] = 10.0*32.0*0.5*mm;
  size_SDC2Plane[ZCOORD] = 10.0*40.0*0.5*mm;
  */
  double size_SDC2[XYZ];
  //original
  /*
  size_SDC2[XCOORD] = 289.0*mm;
  size_SDC2[YCOORD] = 73.*mm;
  size_SDC2[ZCOORD] = 215.0*mm;
  */
  // for LEPS
  // size_SDC2[XCOORD] = 600.0/2.*mm;
  // size_SDC2[YCOORD] = 278./2.*mm;
  // size_SDC2[ZCOORD] = 400.0/2.*mm;
  size_SDC2[XCOORD] = 614.0/2.*mm;
  size_SDC2[YCOORD] = 276./2.*mm;
  size_SDC2[ZCOORD] = 310.0/2.*mm;

  double size_SDC2Plane[XYZ];
  // original
  /*
  size_SDC2Plane[XCOORD] = 6.0*67.0*0.5*mm;
  size_SDC2Plane[YCOORD] = 0.0001/2.*mm;
  //size_SDC2Plane[ZCOORD] = 10.0*32.0*0.5*mm;
  size_SDC2Plane[ZCOORD] = 6.0*42.0*0.5*mm;
  */
  // for LEPS
  size_SDC2Plane[XCOORD] = 576./2*mm;
  size_SDC2Plane[YCOORD] = 0.0001/2.*mm;
  //size_SDC2Plane[ZCOORD] = 10.0*32.0*0.5*mm;
  // size_SDC2Plane[ZCOORD] = 400./2*mm;
  size_SDC2Plane[ZCOORD] = 310./2*mm;

  lnum1 = geomMan.GetDetectorId("SDC2-x-1");
  G4ThreeVector posSdc2U1 = geomMan.GetGlobalPosition(lnum1);
  lnum2 = geomMan.GetDetectorId("SDC2-x-3");
  G4ThreeVector posSdc2X2 = geomMan.GetGlobalPosition(lnum2);

  G4ThreeVector offsetSDC2(0., geomMan.GetOffset(lnum2), 0.);
  posSdc2U1 += offsetSDC2;
  posSdc2X2 += offsetSDC2;

  G4RotationMatrix* rotSdc2 = new G4RotationMatrix();
  rotSdc2->rotateZ(-geomMan.GetRotAngle2(lnum1)*deg);

  G4ThreeVector SDC2Pos = G4ThreeVector((posSdc2U1.x()+posSdc2X2.x())/2.0 * mm,
					(posSdc2U1.y()+posSdc2X2.y())/2.0 * mm,
					(posSdc2U1.z()+posSdc2X2.z())/2.0 * mm);

  G4cout << "SDC2 CenterPos = ( " << SDC2Pos.x() << ", " << SDC2Pos.y()
	 << ", " << SDC2Pos.z() << ")" << G4endl;


  G4Box* SDC2_box = new G4Box("SDC2_box",size_SDC2[XCOORD],size_SDC2[YCOORD],size_SDC2[ZCOORD]);
  G4LogicalVolume* SDC2_log =
    new G4LogicalVolume(SDC2_box, mList_->ArGas, "SDC2_log",0,0,0);
  SDC2_log->SetVisAttributes(chamber_att);
  G4VPhysicalVolume* SDC2_phys = new G4PVPlacement(rotSdc2,
						   SDC2Pos,
						   "SDC2_phys",
						   SDC2_log,
						   pMother,
						   false,
						   0);

  G4LogicalVolume* SDC2Plane_log[NumOfPlaneSDC2];
  G4VPhysicalVolume* SDC2Plane_phys[NumOfPlaneSDC2];

  //---------SDC2u1
  G4Box* SDC2Plane_box = new G4Box("SDC2Plane_box",
				   size_SDC2Plane[XCOORD],
				   size_SDC2Plane[YCOORD],
				   size_SDC2Plane[ZCOORD]);
  for (int i=lnum1; i<=lnum2; i++) {
    G4ThreeVector posSdc2plane = geomMan.GetGlobalPosition(i);
    G4ThreeVector localPos;
    // G4ThreeVector offsetSDC2plane(0., geomMan.GetOffset(i), 0.);
    // posSdc2plane += offsetSDC2plane;

    G4cout << "SDC " << i << " Pos = ( " << posSdc2plane.x() << ", "
	   << posSdc2plane.y() << ", " << posSdc2plane.z() << ")" << G4endl;


    switch (i-lnum1) {
    case 0:
      sprintf(name_log, "SDC2X1_log");
      sprintf(name_phys, "SDC2-x-1");
      //localPos = G4ThreeVector(0.0*mm, GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      localPos = G4ThreeVector(0.0*mm, SDC2Pos.x()-posSdc2plane.x()*mm, 0.0*mm);
     break;
    case 1:
      sprintf(name_log, "SDC2V1_log");
      sprintf(name_phys, "SDC2-v-1");
      localPos = G4ThreeVector(0.0*mm, SDC2Pos.x()-posSdc2plane.x()*mm, 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    case 2:
      sprintf(name_log, "SDC2U1_log");
      sprintf(name_phys, "SDC2-u-1");
      localPos = G4ThreeVector(0.0*mm, SDC2Pos.x()-posSdc2plane.x()*mm, 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    case 3:
      sprintf(name_log, "SDC2U2_log");
      sprintf(name_phys, "SDC2-u-2");
      localPos = G4ThreeVector(0.0*mm, SDC2Pos.x()-posSdc2plane.x()*mm, 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, -GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    case 4:
      sprintf(name_log, "SDC2X2_log");
      sprintf(name_phys, "SDC2-x-2");
      localPos = G4ThreeVector(0.0*mm, SDC2Pos.x()-posSdc2plane.x()*mm, 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, -GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    case 5:
      sprintf(name_log, "SDC2X3_log");
      sprintf(name_phys, "SDC2-x-3");
      localPos = G4ThreeVector(0.0*mm, SDC2Pos.x()-posSdc2plane.x()*mm, 0.0*mm);
      //localPos = G4ThreeVector(0.0*mm, -GetDistance(SDC2Pos, posSdc2plane)*mm, 0.0*mm);
      break;
    default:
      fprintf(stderr, "DetectorConstruction::ConstructForwardSpectrometer()\n");
      fprintf(stderr, "No such plane %d in SDC2\n", i);
      exit(-1);
    }


    G4cout << "SDC(local) " << i << " Pos = ( " << localPos.x() << ", "
	   << localPos.y() << ", " << localPos.z() << ")" << G4endl;

    SDC2Plane_log[i-lnum1]  = new G4LogicalVolume(SDC2Plane_box,
						  mList_->ArGas,
						  name_log,
						  0,0,0);

    if ((i-lnum1)==0)
      SDC2Plane_log[i-lnum1]->SetVisAttributes(magnet_att);

    SDC2Plane_phys[i-lnum1] = new G4PVPlacement(0,
						localPos,
						SDC2Plane_log[i-lnum1],
						name_phys,
						SDC2_log,
						false,
						0);
  }

  // VP2
  G4Material *VPMater = mList_->Air;

  G4Box *solidVP2 =
    new G4Box( "VP2_Box",   size_SDC2[XCOORD], 0.03/2.0*mm, size_SDC2[ZCOORD] );
  G4LogicalVolume *logVP2 =
   new G4LogicalVolume( solidVP2,    VPMater,   "logVP2",   0, 0, 0 );

  /*
  G4int lnum_VP2;
  lnum_VP2 = geomMan.GetDetectorId("VP2");

  G4ThreeVector posVP2 = geomMan.GetGlobalPosition( lnum_VP2 );
  //G4ThreeVector offsetVP2(0., geomMan.GetOffset(lnum_VP2), 0.);
  //posVP2 += offsetVP2;
  G4ThreeVector localPos = G4ThreeVector(0.0*mm, SDC2Pos.x() - posVP2.x()*mm, 0.0*mm);

  G4VPhysicalVolume *physVP2 =
    new G4PVPlacement( 0,
		       localPos,
		       logVP2,
		       "VP2",
		       SDC2_log,
		       false, 0 );
  */

  //////// Sensitive Detectors ////////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!chamberSD) {
    chamberSD = new SKSChamberSD( "chamberSD");
    SDMan->AddNewDetector( chamberSD );
  }

  for (int i=0; i<NumOfPlaneSDC1; i++)
    SDC1Plane_log[i]->SetSensitiveDetector(chamberSD);

  for (int i=0; i<NumOfPlaneSDC2; i++)
    SDC2Plane_log[i]->SetSensitiveDetector(chamberSD);

  /*
  if (!virtualPlaneSD) {
    virtualPlaneSD = new SKSVirtualPlaneSD( "virtualPlaneSD");
    SDMan->AddNewDetector( virtualPlaneSD );
  }
  logVP2->SetSensitiveDetector( virtualPlaneSD );
  */
}

void  CFTDetectorConstruction::ConstructSdcOut(G4VPhysicalVolume* pMother)
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour green(0.0, 1.0, 0.0);
  G4VisAttributes* chamber_att = new G4VisAttributes(green);
  chamber_att->SetForceWireframe(true);

  double size_SDC3[XYZ], size_SDC4[XYZ];
  double size_SDC3Plane[XYZ], size_SDC4Plane[XYZ];

  // for HYPS
  size_SDC3[XCOORD] = 2560.0/2.0*mm;
  size_SDC3[YCOORD] = 1120.0/2.0*mm;
  size_SDC3[ZCOORD] = 270.0/2.0*mm;

  size_SDC3Plane[XCOORD] = 2410.0/2.0*mm;
  size_SDC3Plane[YCOORD] = 910.0/2.0*mm;
  size_SDC3Plane[ZCOORD] = 0.0001/2.*mm;

  size_SDC4[XCOORD] = 2560.0/2.0*mm;
  size_SDC4[YCOORD] = 1120.0/2.0*mm;
  size_SDC4[ZCOORD] = 270.0/2.0*mm;

  size_SDC4Plane[XCOORD] = 2410.0/2.0*mm;
  size_SDC4Plane[YCOORD] = 910.0/2.0*mm;
  size_SDC4Plane[ZCOORD] = 0.0001/2.*mm;


  //--------------SDC3

  lnum1 = geomMan.GetDetectorId("SDC3-v-1");
  G4ThreeVector posSDC3X1 = geomMan.GetGlobalPosition(lnum1);
  G4ThreeVector offsetSDC3X1(0., geomMan.GetOffset(lnum1), 0.);
  posSDC3X1 += offsetSDC3X1;

  lnum2 = geomMan.GetDetectorId("SDC3-x-2");
  G4ThreeVector posSDC3X2 = geomMan.GetGlobalPosition(lnum2);
  G4ThreeVector offsetSDC3X2(0., geomMan.GetOffset(lnum2), 0.);
  posSDC3X2 += offsetSDC3X2;

  G4RotationMatrix* rotSDC3 = new G4RotationMatrix();
  rotSDC3->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotSDC3->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4ThreeVector SDC3Pos = G4ThreeVector((posSDC3X2.x()+posSDC3X1.x())/2.0 * mm,
				       (posSDC3X2.y()+posSDC3X1.y())/2.0 * mm,
				       (posSDC3X2.z()+posSDC3X1.z())/2.0 * mm);

  G4Box* SDC3_box = new G4Box("SDC3_box",size_SDC3[XCOORD],size_SDC3[YCOORD],size_SDC3[ZCOORD]);
  G4LogicalVolume* SDC3_log = new G4LogicalVolume(SDC3_box,
						  mList_->ArGas,
						  "SDC3_log",0,0,0);
  SDC3_log->SetVisAttributes(chamber_att);
  G4VPhysicalVolume* SDC3_phys = new G4PVPlacement(rotSDC3,
						   SDC3Pos,
						   "SDC3_phys",
						   SDC3_log,
						   pMother,
						   false,
						   0);

  //---------SDC3 planes
  G4Box* SDC3Plane_box = new G4Box("SDC3Plane_box",
				   size_SDC3Plane[XCOORD],
				   size_SDC3Plane[YCOORD],
				   size_SDC3Plane[ZCOORD]);

  G4LogicalVolume* SDC3Plane_log[NumOfPlaneSDC3];
  G4VPhysicalVolume* SDC3Plane_phys[NumOfPlaneSDC3];

  lnum1 = geomMan.GetDetectorId("SDC3-v-1");
  lnum2 = geomMan.GetDetectorId("SDC3-x-2");

  for (int i=lnum1; i<=lnum2; i++) {
    G4ThreeVector posSDC3plane = geomMan.GetGlobalPosition(i);
    G4ThreeVector localPos;

    switch (i-lnum1) {
    case 0:
      sprintf(name_log, "SDC3V1_log");
      sprintf(name_phys, "SDC3-v-1");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC3Pos.x()-posSDC3plane.x()*mm);
      break;
    case 1:
      sprintf(name_log, "SDC3U1_log");
      sprintf(name_phys, "SDC3-u-1");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC3Pos.x()-posSDC3plane.x()*mm);
      break;
    case 2:
      sprintf(name_log, "SDC3U2_log");
      sprintf(name_phys, "SDC3-u-2");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC3Pos.x()-posSDC3plane.x()*mm);
      break;
    case 3:
      sprintf(name_log, "SDC3X1_log");
      sprintf(name_phys, "SDC3-x-1");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC3Pos.x()-posSDC3plane.x()*mm);
      break;
    case 4:
      sprintf(name_log, "SDC3X2_log");
      sprintf(name_phys, "SDC3-x-2");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC3Pos.x()-posSDC3plane.x()*mm);
      break;
    default:
      fprintf(stderr, "DetectorConstruction::ConstructForwardSpectrometer()\n");
      fprintf(stderr, "No such plane %d in SDC3\n", i);
      exit(-1);
    }
    SDC3Plane_log[i-lnum1]  = new G4LogicalVolume(SDC3Plane_box,
						  mList_->ArGas,
						  name_log,
						  0,0,0);

    SDC3Plane_phys[i-lnum1] = new G4PVPlacement(0,
						localPos,
						SDC3Plane_log[i-lnum1],
						name_phys,
						SDC3_log,
						false,
						0);
  }



  //--------------SDC4

  lnum1 = geomMan.GetDetectorId("SDC4-v-1");
  G4ThreeVector posSDC4X1 = geomMan.GetGlobalPosition(lnum1);
  G4ThreeVector offsetSDC4X1(0., geomMan.GetOffset(lnum1), 0.);
  posSDC4X1 += offsetSDC4X1;

  lnum2 = geomMan.GetDetectorId("SDC4-x-2");
  G4ThreeVector posSDC4X2 = geomMan.GetGlobalPosition(lnum2);
  G4ThreeVector offsetSDC4X2(0., geomMan.GetOffset(lnum2), 0.);
  posSDC4X2 += offsetSDC4X2;

  G4RotationMatrix* rotSDC4 = new G4RotationMatrix();
  rotSDC4->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotSDC4->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4ThreeVector SDC4Pos = G4ThreeVector((posSDC4X1.x()+posSDC4X2.x())/2.0 * mm,
				       (posSDC4X1.y()+posSDC4X2.y())/2.0 * mm,
				       (posSDC4X1.z()+posSDC4X2.z())/2.0 * mm);

  G4Box* SDC4_box = new G4Box("SDC4_box",size_SDC4[XCOORD],size_SDC4[YCOORD],size_SDC4[ZCOORD]);
  G4LogicalVolume* SDC4_log = new G4LogicalVolume(SDC4_box,
						  mList_->ArGas,
						  "SDC4_log",0,0,0);
  SDC4_log->SetVisAttributes(chamber_att);
  G4VPhysicalVolume* SDC4_phys = new G4PVPlacement(rotSDC4,
						   SDC4Pos,
						   "SDC4_phys",
						   SDC4_log,
						   pMother,
						   false,
						   0);

  //---------SDC4 planes
  G4Box* SDC4Plane_box = new G4Box("SDC4Plane_box",
				  size_SDC4Plane[XCOORD],
				  size_SDC4Plane[YCOORD],
				  size_SDC4Plane[ZCOORD]);

  G4LogicalVolume* SDC4Plane_log[NumOfPlaneSDC4];
  G4VPhysicalVolume* SDC4Plane_phys[NumOfPlaneSDC4];

  lnum1 = geomMan.GetDetectorId("SDC4-v-1");
  lnum2 = geomMan.GetDetectorId("SDC4-x-2");

  for (int i=lnum1; i<=lnum2; i++) {
    G4ThreeVector posSDC4plane = geomMan.GetGlobalPosition(i);
    G4ThreeVector localPos;


    switch (i-lnum1) {
    case 0:
      sprintf(name_log, "SDC4V1_log");
      sprintf(name_phys, "SDC4-v-1");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC4Pos.x()-posSDC4plane.x()*mm);
      break;
    case 1:
      sprintf(name_log, "SDC4U1_log");
      sprintf(name_phys, "SDC4-u-1");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC4Pos.x()-posSDC4plane.x()*mm);
      break;
    case 2:
      sprintf(name_log, "SDC4U2_log");
      sprintf(name_phys, "SDC4-u-2");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC4Pos.x()-posSDC4plane.x()*mm);
      break;
    case 3:
      sprintf(name_log, "SDC4X1_log");
      sprintf(name_phys, "SDC4-x-1");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC4Pos.x()-posSDC4plane.x()*mm);
      break;
    case 4:
      sprintf(name_log, "SDC4X2_log");
      sprintf(name_phys, "SDC4-x-2");
      localPos = G4ThreeVector(0.0*mm, 0.0*mm, SDC4Pos.x()-posSDC4plane.x()*mm);
      break;
    default:
      fprintf(stderr, "DetectorConstruction::ConstructForwardSpectrometer()\n");
      fprintf(stderr, "No such plane %d in BD2\n", i);
      exit(-1);
    }
    SDC4Plane_log[i-lnum1]  = new G4LogicalVolume(SDC4Plane_box,
						  mList_->ArGas,
						  name_log,
						  0,0,0);

    SDC4Plane_phys[i-lnum1] = new G4PVPlacement(0,
						localPos,
						SDC4Plane_log[i-lnum1],
						name_phys,
						SDC4_log,
						false,
						0);
  }

  //////// Sensitive Detectors ////////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!chamberSD) {
    chamberSD = new SKSChamberSD( "chamberSD");
    SDMan->AddNewDetector( chamberSD );
  }

  for (int i=0; i<NumOfPlaneSDC3; i++)
    SDC3Plane_log[i]->SetSensitiveDetector(chamberSD);
  for (int i=0; i<NumOfPlaneSDC4; i++)
    SDC4Plane_log[i]->SetSensitiveDetector(chamberSD);


}

void  CFTDetectorConstruction::ConstructBcOut(G4VPhysicalVolume *pMother)
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour green(0.0, 1.0, 0.0);
  G4Colour blue(0.0, 0.0, 1.0);
  G4VisAttributes* chamber_att = new G4VisAttributes(green);
  chamber_att->SetForceWireframe(true);
  G4VisAttributes* magnet_att = new G4VisAttributes(blue);
  magnet_att->SetForceWireframe(true);

  //--------------BC3
  double size_BC3[XYZ];
  size_BC3[XCOORD] = 220.0/2.*mm;// tekitou
  size_BC3[YCOORD] = 60.0/2.*mm; // tekitou
  size_BC3[ZCOORD] = 170.0/2.*mm;// tekitou

  double size_BC3Plane[XYZ];
  size_BC3Plane[XCOORD] = 192./2.*mm;
  size_BC3Plane[YCOORD] = 0.0001/2.*mm;
  size_BC3Plane[ZCOORD] = 150./2.*mm;


  lnum1 = geomMan.GetDetectorId("BC3-x-1");
  G4ThreeVector posBc3X1 = geomMan.GetGlobalPosition(lnum1);
  lnum2 = geomMan.GetDetectorId("BC3-u-2");
  G4ThreeVector posBc3U2 = geomMan.GetGlobalPosition(lnum2);

  G4RotationMatrix* rotBc3 = new G4RotationMatrix();
  rotBc3->rotateZ(-geomMan.GetRotAngle2(lnum1)*deg);

  G4ThreeVector BC3Pos = G4ThreeVector((posBc3X1.x()+posBc3U2.x())/2.0 * mm,
				       (posBc3X1.y()+posBc3U2.y())/2.0 * mm,
				       (posBc3X1.z()+posBc3U2.z())/2.0 * mm);

  G4cout << "BC3 CenterPos = ( " << BC3Pos.x() << ", " << BC3Pos.y()
	 << ", " << BC3Pos.z() << ")" << G4endl;


  G4Box* BC3_box = new G4Box("BC3_box",size_BC3[XCOORD],size_BC3[YCOORD],size_BC3[ZCOORD]);
  G4LogicalVolume* BC3_log =
    new G4LogicalVolume(BC3_box, mList_->ArGas, "BC3_log",0,0,0);
  BC3_log->SetVisAttributes(chamber_att);
  G4VPhysicalVolume* BC3_phys = new G4PVPlacement(rotBc3,
						  BC3Pos,
						  "BC3_phys",
						  BC3_log,
						  pMother,
						  false,
						  0);

  G4LogicalVolume* BC3Plane_log[NumOfPlaneBC3];
  G4VPhysicalVolume* BC3Plane_phys[NumOfPlaneBC3];

  //---------SDC2u1
  G4Box* BC3Plane_box = new G4Box("BC3Plane_box",
				  size_BC3Plane[XCOORD],
				  size_BC3Plane[YCOORD],
				  size_BC3Plane[ZCOORD]);
  for (int i=lnum1; i<=lnum2; i++) {
    G4ThreeVector posBc3plane = geomMan.GetGlobalPosition(i);
    G4ThreeVector localPos;

    G4cout << "BC " << i << " Pos = ( " << posBc3plane.x() << ", "
	   << posBc3plane.y() << ", " << posBc3plane.z() << ")" << G4endl;


    switch (i-lnum1) {
    case 0:
      sprintf(name_log, "BC3X1_log");
      sprintf(name_phys, "BC3-x-1");
      localPos = G4ThreeVector(0.0*mm, GetDistance(BC3Pos, posBc3plane)*mm, 0.0*mm);
     break;
    case 1:
      sprintf(name_log, "BC3X2_log");
      sprintf(name_phys, "BC3-x-2");
      localPos = G4ThreeVector(0.0*mm, GetDistance(BC3Pos, posBc3plane)*mm, 0.0*mm);
      break;
    case 2:
      sprintf(name_log, "BC3V1_log");
      sprintf(name_phys, "BC3-v-1");
      localPos = G4ThreeVector(0.0*mm, GetDistance(BC3Pos, posBc3plane)*mm, 0.0*mm);
      break;
    case 3:
      sprintf(name_log, "BC3V2_log");
      sprintf(name_phys, "BC3-v-2");
      localPos = G4ThreeVector(0.0*mm, -GetDistance(BC3Pos, posBc3plane)*mm, 0.0*mm);
      break;
    case 4:
      sprintf(name_log, "BC3U1_log");
      sprintf(name_phys, "BC3-u-1");
      localPos = G4ThreeVector(0.0*mm, -GetDistance(BC3Pos, posBc3plane)*mm, 0.0*mm);
      break;
    case 5:
      sprintf(name_log, "BC3U2_log");
      sprintf(name_phys, "BC3-u-2");
      localPos = G4ThreeVector(0.0*mm, -GetDistance(BC3Pos, posBc3plane)*mm, 0.0*mm);
      break;
    default:
      fprintf(stderr, "DetectorConstruction::ConstructForwardSpectrometer()\n");
      fprintf(stderr, "No such plane %d in BC3\n", i);
      exit(-1);
    }

    G4cout << "BC(local) " << i << " Pos = ( " << localPos.x() << ", "
	   << localPos.y() << ", " << localPos.z() << ")" << G4endl;

    BC3Plane_log[i-lnum1]  = new G4LogicalVolume(BC3Plane_box,
						 mList_->ArGas,
						 name_log,
						 0,0,0);

    if ((i-lnum1)==0)
      BC3Plane_log[i-lnum1]->SetVisAttributes(magnet_att);

    BC3Plane_phys[i-lnum1] = new G4PVPlacement(0,
					       localPos,
					       BC3Plane_log[i-lnum1],
					       name_phys,
					       BC3_log,
					       false,
					       0);
  }


  //--------------BC4
  double size_BC4[XYZ];
  size_BC4[XCOORD] = 220.0/2.*mm;// tekitou
  size_BC4[YCOORD] = 60.0/2.*mm; // tekitou
  size_BC4[ZCOORD] = 170.0/2.*mm;// tekitou

  double size_BC4Plane[XYZ];
  size_BC4Plane[XCOORD] = 192./2.*mm;
  size_BC4Plane[YCOORD] = 0.0001/2.*mm;
  size_BC4Plane[ZCOORD] = 150./2.*mm;


  lnum1 = geomMan.GetDetectorId("BC4-u-1");
  G4ThreeVector posBc4U1 = geomMan.GetGlobalPosition(lnum1);
  lnum2 = geomMan.GetDetectorId("BC4-x-2");
  G4ThreeVector posBc4X2 = geomMan.GetGlobalPosition(lnum2);

  G4RotationMatrix* rotBc4 = new G4RotationMatrix();
  rotBc4->rotateZ(-geomMan.GetRotAngle2(lnum1)*deg);

  G4ThreeVector BC4Pos = G4ThreeVector((posBc4U1.x()+posBc4X2.x())/2.0 * mm,
				       (posBc4U1.y()+posBc4X2.y())/2.0 * mm,
				       (posBc4U1.z()+posBc4X2.z())/2.0 * mm);

  G4cout << "BC4 CenterPos = ( " << BC4Pos.x() << ", " << BC4Pos.y()
	 << ", " << BC4Pos.z() << ")" << G4endl;


  G4Box* BC4_box = new G4Box("BC4_box",size_BC4[XCOORD],size_BC4[YCOORD],size_BC4[ZCOORD]);
  G4LogicalVolume* BC4_log =
    new G4LogicalVolume(BC4_box, mList_->ArGas, "BC4_log",0,0,0);
  BC4_log->SetVisAttributes(chamber_att);
  G4VPhysicalVolume* BC4_phys = new G4PVPlacement(rotBc4,
						  BC4Pos,
						  "BC4_phys",
						  BC4_log,
						  pMother,
						  false,
						  0);

  G4LogicalVolume* BC4Plane_log[NumOfPlaneBC4];
  G4VPhysicalVolume* BC4Plane_phys[NumOfPlaneBC4];

  //---------SDC2u1
  G4Box* BC4Plane_box = new G4Box("BC4Plane_box",
				  size_BC4Plane[XCOORD],
				  size_BC4Plane[YCOORD],
				  size_BC4Plane[ZCOORD]);
  for (int i=lnum1; i<=lnum2; i++) {
    G4ThreeVector posBc4plane = geomMan.GetGlobalPosition(i);
    G4ThreeVector localPos;

    G4cout << "BC " << i << " Pos = ( " << posBc4plane.x() << ", "
	   << posBc4plane.y() << ", " << posBc4plane.z() << ")" << G4endl;


    switch (i-lnum1) {
    case 0:
      sprintf(name_log, "BC4U1_log");
      sprintf(name_phys, "BC4-u-1");
      localPos = G4ThreeVector(0.0*mm, GetDistance(BC4Pos, posBc4plane)*mm, 0.0*mm);
     break;
    case 1:
      sprintf(name_log, "BC4U2_log");
      sprintf(name_phys, "BC4-u-2");
      localPos = G4ThreeVector(0.0*mm, GetDistance(BC4Pos, posBc4plane)*mm, 0.0*mm);
      break;
    case 2:
      sprintf(name_log, "BC4V1_log");
      sprintf(name_phys, "BC4-v-1");
      localPos = G4ThreeVector(0.0*mm, GetDistance(BC4Pos, posBc4plane)*mm, 0.0*mm);
      break;
    case 3:
      sprintf(name_log, "BC4V2_log");
      sprintf(name_phys, "BC4-v-2");
      localPos = G4ThreeVector(0.0*mm, -GetDistance(BC4Pos, posBc4plane)*mm, 0.0*mm);
      break;
    case 4:
      sprintf(name_log, "BC4X1_log");
      sprintf(name_phys, "BC4-x-1");
      localPos = G4ThreeVector(0.0*mm, -GetDistance(BC4Pos, posBc4plane)*mm, 0.0*mm);
      break;
    case 5:
      sprintf(name_log, "BC4X2_log");
      sprintf(name_phys, "BC4-x-2");
      localPos = G4ThreeVector(0.0*mm, -GetDistance(BC4Pos, posBc4plane)*mm, 0.0*mm);
      break;
    default:
      fprintf(stderr, "DetectorConstruction::ConstructForwardSpectrometer()\n");
      fprintf(stderr, "No such plane %d in BC4\n", i);
      exit(-1);
    }

    G4cout << "BC(local) " << i << " Pos = ( " << localPos.x() << ", "
	   << localPos.y() << ", " << localPos.z() << ")" << G4endl;

    BC4Plane_log[i-lnum1]  = new G4LogicalVolume(BC4Plane_box,
						 mList_->ArGas,
						 name_log,
						 0,0,0);

    if ((i-lnum1)==0)
      BC4Plane_log[i-lnum1]->SetVisAttributes(magnet_att);

    BC4Plane_phys[i-lnum1] = new G4PVPlacement(0,
					       localPos,
					       BC4Plane_log[i-lnum1],
					       name_phys,
					       BC4_log,
					       false,
					       0);
  }


  //////// Sensitive Detectors ////////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!chamberSD) {
    chamberSD = new SKSChamberSD( "chamberSD");
    SDMan->AddNewDetector( chamberSD );
  }


  for (int i=0; i<NumOfPlaneBC3; i++)
    BC3Plane_log[i]->SetSensitiveDetector(chamberSD);

  for (int i=0; i<NumOfPlaneBC4; i++)
    BC4Plane_log[i]->SetSensitiveDetector(chamberSD);

}

void  CFTDetectorConstruction::ConstructVP(G4VPhysicalVolume *pMother)
{
  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour green(0.0, 1.0, 0.0);
  G4Colour blue(0.0, 0.0, 1.0);
  G4VisAttributes* chamber_att = new G4VisAttributes(green);
  chamber_att->SetForceWireframe(true);
  G4VisAttributes* magnet_att = new G4VisAttributes(blue);
  magnet_att->SetForceWireframe(true);


  G4Material *VPMater = mList_->Air;

  // VP1
  G4Box *solidVP1 =
    new G4Box( "VP1_Box",   600.0/2.0*mm, 0.03/2.0*mm, 600.0/2.0*mm );
  G4LogicalVolume *logVP1 =
   new G4LogicalVolume( solidVP1,    VPMater,   "logVP1",   0, 0, 0 );

  G4int lnum_VP1;
  lnum_VP1 = geomMan.GetDetectorId("VP1");
  G4double RotateAngleVP1;
  RotateAngleVP1 = geomMan.GetRotAngle2( lnum_VP1 );

  G4RotationMatrix RM_VP1;
  RM_VP1.rotateZ( RotateAngleVP1*degree );

  G4ThreeVector posVP1 = geomMan.GetGlobalPosition( lnum_VP1 );
  G4ThreeVector offsetVP1(0., geomMan.GetOffset(lnum_VP1), 0.);
  posVP1 += offsetVP1;

  G4VPhysicalVolume *physVP1 =
    new G4PVPlacement( G4Transform3D(RM_VP1,
                                     G4ThreeVector( posVP1.x()*mm,
                                                    posVP1.y()*mm,
                                                    posVP1.z()*mm )),
		       "VP1", logVP1, pMother, false, 0 );


  // VP2
  G4Box *solidVP2 =
    //new G4Box( "VP2_Box",   600.0/2.0*mm, 0.03/2.0*mm, 600.0/2.0*mm );
    new G4Box( "VP2_Box",   1350.0/2.0*mm, 0.03/2.0*mm, 600.0/2.0*mm );
  G4LogicalVolume *logVP2 =
   new G4LogicalVolume( solidVP2,    VPMater,   "logVP2",   0, 0, 0 );

  G4int lnum_VP2;
  lnum_VP2 = geomMan.GetDetectorId("VP2");

  G4ThreeVector posVP2 = geomMan.GetGlobalPosition( lnum_VP2 );
  G4ThreeVector offsetVP2(0., geomMan.GetOffset(lnum_VP2), 0.);
  posVP2 += offsetVP2;

  G4VPhysicalVolume *physVP2 =
    new G4PVPlacement( G4Transform3D(RM_VP1,
                                     G4ThreeVector( posVP2.x()*mm,
                                                    posVP2.y()*mm,
                                                    posVP2.z()*mm )),
		       "VP2", logVP2, pMother, false, 0 );


  // VP3
  G4Box *solidVP3 =
    //new G4Box( "VP3_Box",   600.0/2.0*mm, 0.03/2.0*mm, 600.0/2.0*mm );
    new G4Box( "VP3_Box",   1350.0/2.0*mm, 0.03/2.0*mm, 600.0/2.0*mm );
  G4LogicalVolume *logVP3 =
   new G4LogicalVolume( solidVP3,    VPMater,   "logVP3",   0, 0, 0 );

  G4int lnum_VP3;
  lnum_VP3 = geomMan.GetDetectorId("VP3");

  G4ThreeVector posVP3 = geomMan.GetGlobalPosition( lnum_VP3 );
  G4ThreeVector offsetVP3(0., geomMan.GetOffset(lnum_VP3), 0.);
  posVP3 += offsetVP3;

  G4VPhysicalVolume *physVP3 =
    new G4PVPlacement( G4Transform3D(RM_VP1,
                                     G4ThreeVector( posVP3.x()*mm,
                                                    posVP3.y()*mm,
                                                    posVP3.z()*mm )),
		       "VP3", logVP3, pMother, false, 0 );


  // VP4
  G4Box *solidVP4 =
    new G4Box( "VP4_Box",   3000.0/2.0*mm, 0.03/2.0*mm, 2000.0/2.0*mm );
  G4LogicalVolume *logVP4 =
   new G4LogicalVolume( solidVP4,    VPMater,   "logVP4",   0, 0, 0 );

  G4int lnum_VP4;
  lnum_VP4 = geomMan.GetDetectorId("VP4");

  G4ThreeVector posVP4 = geomMan.GetGlobalPosition( lnum_VP4 );
  G4ThreeVector offsetVP4(0., geomMan.GetOffset(lnum_VP4), 0.);
  posVP4 += offsetVP4;

  G4VPhysicalVolume *physVP4 =
    new G4PVPlacement( G4Transform3D(RM_VP1,
                                     G4ThreeVector( posVP4.x()*mm,
                                                    posVP4.y()*mm,
                                                    posVP4.z()*mm )),
		       "VP4", logVP4, pMother, false, 0 );

  // VP5
  G4Box *solidVP5 =
    new G4Box( "VP5_Box",   3000.0/2.0*mm, 0.03/2.0*mm, 2000.0/2.0*mm );
  G4LogicalVolume *logVP5 =
   new G4LogicalVolume( solidVP5,    VPMater,   "logVP5",   0, 0, 0 );

  G4int lnum_VP5;
  lnum_VP5 = geomMan.GetDetectorId("VP5");

  G4ThreeVector posVP5 = geomMan.GetGlobalPosition( lnum_VP5 );
  G4ThreeVector offsetVP5(0., geomMan.GetOffset(lnum_VP5), 0.);
  posVP5 += offsetVP5;

  G4VPhysicalVolume *physVP5 =
    new G4PVPlacement( G4Transform3D(RM_VP1,
                                     G4ThreeVector( posVP5.x()*mm,
                                                    posVP5.y()*mm,
                                                    posVP5.z()*mm )),
		       "VP5", logVP5, pMother, false, 0 );




  //////// Sensitive Detectors ////////
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!virtualPlaneSD) {
    virtualPlaneSD = new SKSVirtualPlaneSD( "virtualPlaneSD");
    SDMan->AddNewDetector( virtualPlaneSD );
  }

  logVP1->SetSensitiveDetector( virtualPlaneSD );
  logVP2->SetSensitiveDetector( virtualPlaneSD );
  logVP3->SetSensitiveDetector( virtualPlaneSD );
  logVP4->SetSensitiveDetector( virtualPlaneSD );
  logVP5->SetSensitiveDetector( virtualPlaneSD );
}


void CFTDetectorConstruction::ConstructCH(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  double size_CHWall[XYZ], size_CH[XYZ];
  size_CHWall[XCOORD] = 11.5*NumOfSegmentCH/2.0*mm;
  size_CHWall[YCOORD] = 240.0*mm;
  size_CHWall[ZCOORD] = 4.0/2.0*mm;

  //size_CH[XCOORD] = 5.0*mm;
  size_CH[XCOORD] = 11.5/2*mm;
  size_CH[YCOORD] = 200.0*mm;
  size_CH[ZCOORD] = 1.0*mm;

  //--------------CH
  lnum1 = geomMan.GetDetectorId("CH");
  G4ThreeVector CHWallPos = geomMan.GetGlobalPosition(lnum1);
  G4ThreeVector offsetCH(0., geomMan.GetOffset(lnum1), 0.);
  CHWallPos += offsetCH;

  G4RotationMatrix* rotCH = new G4RotationMatrix();
  rotCH->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotCH->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* CH_box = new G4Box("CH_box",
			    size_CHWall[XCOORD],
			    size_CHWall[YCOORD],
			    size_CHWall[ZCOORD]);

  G4LogicalVolume* CH_log = new G4LogicalVolume(CH_box,         // Solid
						mList_->Air,               // Material
						"CH_log",       // name
						0,                // FieldMgr
						0,                // SDetector
						0);               // ULimits

  G4VPhysicalVolume* CH_phys =
    new G4PVPlacement(rotCH,         // Rotaion
		      CHWallPos,     // location
		      "CH",                                    // name
		      CH_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo


  //-------------------------------each CH
  G4Box* sCH_box = new G4Box("sCH_box",
			       size_CH[XCOORD],
			       size_CH[YCOORD],
			       size_CH[ZCOORD]);

  G4LogicalVolume* sCH_log[NumOfSegmentCH];
  G4VPhysicalVolume* sCH_phys[NumOfSegmentCH];
  G4VPhysicalVolume* sCH_phys_D[NumOfSegmentCH];

  for (int i=1; i<=NumOfSegmentCH; i++) {
    sprintf(name_log, "sCH_log%d", i);
    sCH_log[i-1] = new G4LogicalVolume(sCH_box,   // Solid
				       mList_->Scin,     // Material
				       name_log,   // name
				       0,          // FieldMgr
				       0,          // SDetector
				       0);         // ULimits
    sCH_log[i-1]->SetVisAttributes(counter_att);
    //// Span is 17.5mm, stagging 2.5mm
    // Span is 10.5mm, stagging 1.0mm
    double localYpos;
    if (i%2 == 0)
      localYpos = -1.0*mm;
    else
      localYpos = +1.0*mm;

    /* Normal Setup */
    G4ThreeVector chSegLocalPos =
      //G4ThreeVector((double)(i-32.5)*9.0*mm, localYpos, 0.0);
      G4ThreeVector((double)(i-32.5)*10.5*mm, 0, localYpos);

    sprintf(name_phys, "sCH_phys%d", i);

    sCH_phys[i-1] =
      new G4PVPlacement(0,                        // Rotaion
			chSegLocalPos,           // location
			sCH_log[i-1],            // Logical
			name_phys,                // name
			CH_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo

  }

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  for (int i=0; i<NumOfSegmentCH; i++) {
    sCH_log[i]->SetSensitiveDetector( counterSD );
  }

}

#if 1

void CFTDetectorConstruction::ConstructSH(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4Colour red(1.0, 0.0, 0.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  G4VisAttributes* mark_att = new G4VisAttributes(red);

  counter_att->SetForceWireframe(true);

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  double size_SH1Wall[XYZ], size_SH[XYZ];
  size_SH1Wall[XCOORD] = 1100./2.0*mm;
  size_SH1Wall[YCOORD] = 1100.0/2.0*mm;
  size_SH1Wall[ZCOORD] = 60.0/2.0*mm;

  size_SH[XCOORD] = 6.0/2*mm;
  size_SH[YCOORD] = 1000.0/2.*mm;
  size_SH[ZCOORD] = 2.0/2.*mm;

  //--------------SH1
  lnum1 = geomMan.GetDetectorId("SH1-0");
  G4ThreeVector posSH1_0 = geomMan.GetGlobalPosition(lnum1);

  lnum2 = geomMan.GetDetectorId("SH1-3");
  G4ThreeVector posSH1_3 = geomMan.GetGlobalPosition(lnum2);

  G4ThreeVector SH1Pos = G4ThreeVector((posSH1_0.x()+posSH1_3.x())/2.*mm,
				       0.*mm,
				       (posSH1_0.z()+posSH1_3.z())/2.*mm);

  G4RotationMatrix* rotSH1 = new G4RotationMatrix();
  rotSH1->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotSH1->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* SH1_box = new G4Box("SH1_box",
			    size_SH1Wall[XCOORD],
			    size_SH1Wall[YCOORD],
			    size_SH1Wall[ZCOORD]);

  G4LogicalVolume* SH1_log = new G4LogicalVolume(SH1_box,         // Solid
						 mList_->Air,               // Material
						 "SH1_log",       // name
						 0,                // FieldMgr
						 0,                // SDetector
						 0);               // ULimits

  G4VPhysicalVolume* SH1_phys =
    new G4PVPlacement(rotSH1,         // Rotaion
		      SH1Pos,     // location
		      "SH1",                                    // name
		      SH1_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo


  //-------------------------------each SH1
  G4Box* sSH1_box = new G4Box("sSH1_box",
			       size_SH[XCOORD],
			       size_SH[YCOORD],
			       size_SH[ZCOORD]);

  G4LogicalVolume* sSH1_log[NumOfPlaneSH1][NumOfSegmentSH1];
  G4VPhysicalVolume* sSH1_phys[NumOfPlaneSH1][NumOfSegmentSH1];
  G4VPhysicalVolume* sSH1_phys_D[NumOfPlaneSH1][NumOfSegmentSH1];

  G4RotationMatrix* rotSH_Z = new G4RotationMatrix();
  rotSH_Z->rotateZ(90*deg);

  for (int i=1; i<=NumOfPlaneSH1; i++) {
    double w0 = 23.;
    if (i==2 || i==4)
      w0 = 23.5;

    for (int j=0; j<NumOfSegmentSH1; j++) {
      sprintf(name_log, "sSH1_%d_log%d", i, j);
      sSH1_log[i-1][j] = new G4LogicalVolume(sSH1_box,   // Solid
					     //mList_->Scin,     // Material
					     mList_->Air,     // Material
					     name_log,   // name
					     0,          // FieldMgr
					     0,          // SDetector
					     0);         // ULimits
      if (j==0)
	sSH1_log[i-1][j]->SetVisAttributes(mark_att);
      else
	sSH1_log[i-1][j]->SetVisAttributes(counter_att);
      // Span is 4.0mm, stagging 2.0mm
      double localZpos;
      if (j%2 == 0)
	localZpos = -1.0*mm;
      else
	localZpos = +1.0*mm;

      G4RotationMatrix* rotMatrix;

      double localXpos, localYpos;
      if (i==1 || i==2) {
	localXpos = ((double)j-w0)*4.0*mm + posSH1_0.y() + geomMan.GetOffset(lnum1);
	localYpos = 0;
	rotMatrix = 0;
      } else {
	localXpos = 0;
	localYpos = -((double)j-w0)*4.0*mm; // - is determined to match actural geometry
	rotMatrix = rotSH_Z;
      }

      if (i==1)
	localZpos += 22.5;
      else if (i==2)
	localZpos += 7.5;
      else if (i==3)
	localZpos += -7.5;
      else if (i==4)
	localZpos += -22.5;

      /* Normal Setup */
      G4ThreeVector sh1SegLocalPos =
	G4ThreeVector(localXpos, localYpos, localZpos);



      sprintf(name_phys, "sSH1_%d_phys%d", i, j);

      sSH1_phys[i-1][j] =
	new G4PVPlacement(rotMatrix,                        // Rotaion
			  sh1SegLocalPos,           // location
			  sSH1_log[i-1][j],            // Logical
			  name_phys,                // name
			  SH1_log,                  // Mother
			  false,                    // Many
			  0);                       // CopyNo

    }
  }

  double size_SH2Wall[XYZ];
  size_SH2Wall[XCOORD] = 2240./2.0*mm;
  size_SH2Wall[YCOORD] = 1240.0/2.0*mm;
  size_SH2Wall[ZCOORD] = 60.0/2.0*mm;

  double size_SH2[XYZ];
  size_SH2[XCOORD] = 6.0/2*mm;
  size_SH2[YCOORD] = 2240.0/2.*mm;
  size_SH2[ZCOORD] = 2.0/2.*mm;

  //--------------SH2
  lnum2 = geomMan.GetDetectorId("SH2-3");
  G4ThreeVector posSH2_3 = geomMan.GetGlobalPosition(lnum2);

  lnum1 = geomMan.GetDetectorId("SH2-0");
  G4ThreeVector posSH2_0 = geomMan.GetGlobalPosition(lnum1);


  G4ThreeVector SH2Pos = G4ThreeVector((posSH2_0.x()+posSH2_3.x())/2.*mm,
				       0.*mm,
				       (posSH2_0.z()+posSH2_3.z())/2.*mm);

  G4RotationMatrix* rotSH2 = new G4RotationMatrix();
  rotSH2->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotSH2->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* SH2_box = new G4Box("SH2_box",
			    size_SH2Wall[XCOORD],
			    size_SH2Wall[YCOORD],
			    size_SH2Wall[ZCOORD]);

  G4LogicalVolume* SH2_log = new G4LogicalVolume(SH2_box,         // Solid
						 mList_->Air,               // Material
						 "SH2_log",       // name
						 0,                // FieldMgr
						 0,                // SDetector
						 0);               // ULimits

  G4VPhysicalVolume* SH2_phys =
    new G4PVPlacement(rotSH2,         // Rotaion
		      SH2Pos,     // location
		      "SH2",                                    // name
		      SH2_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo


  //-------------------------------each SH2
  G4Box* sSH2_box = new G4Box("sSH2_box",
			       size_SH2[XCOORD],
			       size_SH2[YCOORD],
			       size_SH2[ZCOORD]);

  G4LogicalVolume* sSH2_log[NumOfPlaneSH2][NumOfSegmentSH2];
  G4VPhysicalVolume* sSH2_phys[NumOfPlaneSH2][NumOfSegmentSH2];
  G4VPhysicalVolume* sSH2_phys_D[NumOfPlaneSH2][NumOfSegmentSH2];

  for (int i=1; i<=NumOfPlaneSH2; i++) {
    double w0 = 23.;
    if (i==2 || i==4)
      w0 = 23.5;

    for (int j=0; j<NumOfSegmentSH2; j++) {
      sprintf(name_log, "sSH2_%d_log%d", i, j);
      sSH2_log[i-1][j] = new G4LogicalVolume(sSH2_box,   // Solid
					     //mList_->Scin,     // Material
					     mList_->Air,     // Material
					     name_log,   // name
					     0,          // FieldMgr
					     0,          // SDetector
					     0);         // ULimits
      sSH2_log[i-1][j]->SetVisAttributes(counter_att);
      // Span is 4.0mm, stagging 2.0mm
      double localZpos;
      if (j%2 == 0)
	localZpos = -1.0*mm;
      else
	localZpos = +1.0*mm;

      G4RotationMatrix* rotMatrix;
      double localXpos, localYpos;
      if (i==1 || i==2) {
	localXpos = 0;
	localYpos = -((double)j-w0)*4.0*mm;
	rotMatrix = rotSH_Z;
      } else {
	localXpos = ((double)j-w0)*4.0*mm + posSH2_0.y() + geomMan.GetOffset(lnum2);
	localYpos = 0;
	rotMatrix = 0;
      }

      if (i==1)
	localZpos += 22.5;
      else if (i==2)
	localZpos += 7.5;
      else if (i==3)
	localZpos += -7.5;
      else if (i==4)
	localZpos += -22.5;

      /* Normal Setup */
      G4ThreeVector sh2SegLocalPos =
	G4ThreeVector(localXpos, localYpos, localZpos);

      sprintf(name_phys, "sSH2_%d_phys%d", i, j);

      sSH2_phys[i-1][j] =
	new G4PVPlacement(rotMatrix,                        // Rotaion
			  sh2SegLocalPos,           // location
			  sSH2_log[i-1][j],            // Logical
			  name_phys,                // name
			  SH2_log,                  // Mother
			  false,                    // Many
			  0);                       // CopyNo

    }
  }
  /*
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!chamberSD) {
    chamberSD = new SKSChamberSD( "chamberSD");
    SDMan->AddNewDetector( chamberSD );
  }

  for (int i=0; i<NumOfPlaneSH1; i++) {
    for (int j=0; j<NumOfSegmentSH1; j++) {
      sSH1_log[i][j]->SetSensitiveDetector( chamberSD );
    }
  }

  for (int i=0; i<NumOfPlaneSH2; i++) {
    for (int j=0; j<NumOfSegmentSH2; j++) {
      sSH2_log[i][j]->SetSensitiveDetector( chamberSD );
    }
  }
  */
}

#endif

void CFTDetectorConstruction::ConstructCHwithBeamHole(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  double size_CHWall[XYZ], size_CH[XYZ];
  size_CHWall[XCOORD] = 20.0*24./2.0*mm;
  size_CHWall[YCOORD] = 165.0*mm;
  size_CHWall[ZCOORD] = 4.0/2.0*mm;

  size_CH[XCOORD] = 10.0*mm;
  size_CH[YCOORD] = 165.0*mm;
  size_CH[ZCOORD] = 1.0*mm;

  //--------------CH
  lnum1 = geomMan.GetDetectorId("CH");
  G4ThreeVector CHWallPos = geomMan.GetGlobalPosition(lnum1);

  G4RotationMatrix* rotCH = new G4RotationMatrix();
  rotCH->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotCH->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* CH_box = new G4Box("CH_box",
			    size_CHWall[XCOORD],
			    size_CHWall[YCOORD],
			    size_CHWall[ZCOORD]);

  G4LogicalVolume* CH_log = new G4LogicalVolume(CH_box,         // Solid
						mList_->Air,               // Material
						"CH_log",       // name
						0,                // FieldMgr
						0,                // SDetector
						0);               // ULimits

  G4VPhysicalVolume* CH_phys =
    new G4PVPlacement(rotCH,         // Rotaion
		      CHWallPos,     // location
		      "CH",                                    // name
		      CH_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo


  //-------------------------------each CH
  G4Box* sCH_box = new G4Box("sCH_box",
			       size_CH[XCOORD],
			       size_CH[YCOORD],
			       size_CH[ZCOORD]);

  G4LogicalVolume* sCH_log[NumOfSegmentCH];
  G4VPhysicalVolume* sCH_phys[NumOfSegmentCH];
  G4VPhysicalVolume* sCH_physD[NumOfSegmentCH];

  int MinHoleSeg=3;
  //int MaxHoleSeg=6;
  int MaxHoleSeg=7;
  double HoleSizeY=15.*mm;

  for (int i=1; i<=NumOfSegmentCH; i++) {
    if (!(i>=MinHoleSeg && i<=MaxHoleSeg)) {
      sprintf(name_log, "sCH_log%d", i);
      sCH_log[i-1] = new G4LogicalVolume(sCH_box,   // Solid
					 mList_->Scin,     // Material
					 name_log,   // name
					 0,          // FieldMgr
					 0,          // SDetector
					 0);         // ULimits
      sCH_log[i-1]->SetVisAttributes(counter_att);
      // Span is 17.5mm, stagging 2.5mm
      double localZpos;
      if (i%2 == 0)
	localZpos = -1.0*mm;
      else
	localZpos = +1.0*mm;


      /* Normal Setup */
      G4ThreeVector chSegLocalPos =
	G4ThreeVector((double)(i-12.5)*17.5*mm, 0*mm, localZpos);

      sprintf(name_phys, "sCH_phys%d", i);

      sCH_phys[i-1] =
	new G4PVPlacement(0,                        // Rotaion
			  chSegLocalPos,           // location
			  sCH_log[i-1],            // Logical
			  name_phys,                // name
			  CH_log,                  // Mother
			  false,                    // Many
			  0);                       // CopyNo

    }
  }

  G4Box* sCH2_box = new G4Box("sCH2_box",
			      size_CH[XCOORD],
			      (size_CH[YCOORD]-HoleSizeY)/2.,
			      size_CH[ZCOORD]);

  for (int i=MinHoleSeg; i<=MaxHoleSeg; i++) {
    sprintf(name_log, "sCH_log%d", i);
    sCH_log[i-1] = new G4LogicalVolume(sCH2_box,   // Solid
				       mList_->Scin,     // Material
				       name_log,   // name
				       0,          // FieldMgr
				       0,          // SDetector
				       0);         // ULimits
    sCH_log[i-1]->SetVisAttributes(counter_att);
    // Span is 17.5mm, stagging 2.5mm
    double localZpos;
    if (i%2 == 0)
      localZpos = -1.0*mm;
    else
      localZpos = +1.0*mm;

    double localYpos = HoleSizeY + (size_CH[YCOORD]-HoleSizeY)/2.;

    /* Normal Setup */
    G4ThreeVector chSegLocalPos =
      G4ThreeVector((double)(i-12.5)*17.5*mm, localYpos, localZpos);

    sprintf(name_phys, "sCH_phys%d", i);

    sCH_phys[i-1] =
      new G4PVPlacement(0,                        // Rotaion
			chSegLocalPos,           // location
			sCH_log[i-1],            // Logical
			name_phys,                // name
			CH_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo

    G4ThreeVector chSegLocalPosD =
      G4ThreeVector((double)(i-12.5)*17.5*mm, -localYpos, localZpos);
    sprintf(name_phys, "sCH_physD%d", i);
    sCH_physD[i-1] =
      new G4PVPlacement(0,                        // Rotaion
			chSegLocalPosD,           // location
			sCH_log[i-1],            // Logical
			name_phys,                // name
			CH_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo


  }



  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  for (int i=0; i<NumOfSegmentCH; i++) {
    sCH_log[i]->SetSensitiveDetector( counterSD );
  }

}


void CFTDetectorConstruction::ConstructT0(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  double size_T0[XYZ];

  //for HYPS
  size_T0[XCOORD] = 600.0/2.*mm;
  size_T0[YCOORD] = 172.0/2.*mm;
  size_T0[ZCOORD] = 10.0/2.*mm;


  //--------------T0
  lnum1 = geomMan.GetDetectorId("T0");
  G4ThreeVector T0Pos = geomMan.GetGlobalPosition(lnum1);
  G4ThreeVector offsetT0(0., geomMan.GetOffset(lnum1), 0.);
  T0Pos += offsetT0;

  G4RotationMatrix* rotT0 = new G4RotationMatrix();
  rotT0->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotT0->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* T0_box = new G4Box("T0_box",
			    size_T0[XCOORD],
			    size_T0[YCOORD],
			    size_T0[ZCOORD]);

  G4LogicalVolume* T0_log = new G4LogicalVolume(T0_box,         // Solid
						mList_->Air, // Material
						"T0_log",       // name
						0,                // FieldMgr
						0,                // SDetector
						0);               // ULimits

  G4VPhysicalVolume* T0_phys =
    new G4PVPlacement(rotT0,         // Rotaion
		      T0Pos,     // location
		      "T0",                                    // name
		      T0_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo

  //-------------------------------each T0
  double size_T0_each[XYZ];
  size_T0_each[XCOORD] = size_T0[XCOORD]/NumOfSegmentT0;
  size_T0_each[YCOORD] = size_T0[YCOORD];
  size_T0_each[ZCOORD] = size_T0[ZCOORD];

  G4Box* sT0_box = new G4Box("sT0_box",
			     size_T0_each[XCOORD],
			     size_T0_each[YCOORD],
			     size_T0_each[ZCOORD]);

  G4LogicalVolume* sT0_log = new G4LogicalVolume(sT0_box,         // Solid
						 mList_->Scin, // Material
						 "sT0_log",       // name
						 0,                // FieldMgr
						 0,                // SDetector
						 0);               // ULimits


  double size_T0_Hole[XYZ];
  /*
  double xHoleWidth = 80*mm;
  double yHoleWidth = 40*mm;
  */
  // double xHoleWidth = 40./2.*mm;
  // double yHoleWidth = 30./2.*mm;

  /*
  //for HYPS
  double xHoleWidth = 0./2.*mm;
  double yHoleWidth = 0./2.*mm;
  size_T0_Hole[XCOORD] = xHoleWidth;
  size_T0_Hole[YCOORD] = yHoleWidth;
  size_T0_Hole[ZCOORD] = size_T0_each[ZCOORD];

  G4Box* sT0_Hole_box = new G4Box("sT0_Hole_box",
				  size_T0_Hole[XCOORD],
				  size_T0_Hole[YCOORD],
				  size_T0_Hole[ZCOORD]);

  G4LogicalVolume* sT0_Hole_log = new G4LogicalVolume(sT0_Hole_box,         // Solid
						      mList_->Air, // Material
						      "sT0_Hole_log",       // name
						      0,                // FieldMgr
						      0,                // SDetector
						      0);               // ULimits

  */

  G4VPhysicalVolume* sT0_phys[NumOfSegmentT0];
  // G4VPhysicalVolume* sT0_Hole_phys;


  // NumOfSegmentT0 should be 1
  for (int i=0; i<NumOfSegmentT0; i++) {
      /* Normal Setup */
    G4ThreeVector LocalPos =
      G4ThreeVector(((double)i-((double)(NumOfSegmentT0-1))/2.)*size_T0[XCOORD]*2, 0.*mm, 0.*mm);

    sprintf(name_phys, "sT0_phys%d", i);

    sT0_phys[i] =
      new G4PVPlacement(0,                        // Rotaion
			LocalPos,           // location
			sT0_log,            // Logical
			name_phys,                // name
			T0_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo
    /*
    double xHoleCenter = 0.*mm;
    double yHoleCenter = 0.*mm;

    if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7) {
      //xHoleCenter = 75.*mm - 210.*mm;

      xHoleCenter = -size_T0[XCOORD] + size_T0_Hole[XCOORD];
      yHoleCenter = 0.*mm;
    } else {
      //xHoleCenter = -380.*mm - 210.*mm;
      xHoleCenter = -size_T0[XCOORD] + size_T0_Hole[XCOORD];
      yHoleCenter = 0.*mm;
    }
    G4ThreeVector LocalPos_Hole = G4ThreeVector(xHoleCenter, yHoleCenter, 0);

    sT0_Hole_phys =
      new G4PVPlacement(0,                        // Rotaion
			LocalPos_Hole,           // location
			sT0_Hole_log,            // Logical
			"sT0_Hole_phys",                // name
			sT0_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo


    */

  }
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  T0_log->SetSensitiveDetector( counterSD );

}


void CFTDetectorConstruction::ConstructAC(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  double size_ACWall[XYZ];
  /*
  size_ACWall[XCOORD] = 2740.0/2.*mm;
  size_ACWall[YCOORD] = 1740.0/2.*mm;
  size_ACWall[ZCOORD] = 80.0/2.*mm;
  */

  // size_ACWall[XCOORD] = 450.0/2.*mm;
  // size_ACWall[YCOORD] = 350.0/2.*mm;
  // size_ACWall[ZCOORD] = 80.0/2.*mm;

  //for HYPS
  size_ACWall[XCOORD] = 120.0/2.*mm;
  size_ACWall[YCOORD] = 110.0/2.*mm;
  size_ACWall[ZCOORD] = 60.0/2.*mm;


  //--------------AC
  lnum1 = geomMan.GetDetectorId("AC");
  G4ThreeVector ACWallPos = geomMan.GetGlobalPosition(lnum1);
  G4ThreeVector offsetAC(0., geomMan.GetOffset(lnum1), 0.);
  ACWallPos += offsetAC;

  G4RotationMatrix* rotAC = new G4RotationMatrix();
  rotAC->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotAC->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* AC_box = new G4Box("AC_box",
			    size_ACWall[XCOORD],
			    size_ACWall[YCOORD],
			    size_ACWall[ZCOORD]);

  G4LogicalVolume* AC_log = new G4LogicalVolume(AC_box,         // Solid
						mList_->Air, // Material
						"AC_log",       // name
						0,                // FieldMgr
						0,                // SDetector
						0);               // ULimits

  G4VPhysicalVolume* AC_phys =
    new G4PVPlacement(rotAC,         // Rotaion
		      ACWallPos,     // location
		      "AC",                                    // name
		      AC_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo



  //-------------------------------each AC
  double size_AC[XYZ];
  size_AC[XCOORD] = size_ACWall[XCOORD]/NumOfSegmentAC;
  size_AC[YCOORD] = size_ACWall[YCOORD];
  size_AC[ZCOORD] = size_ACWall[ZCOORD];

  G4Box* sAC_box = new G4Box("sAC_box",
			     size_AC[XCOORD],
			     size_AC[YCOORD],
			     size_AC[ZCOORD]);

  G4LogicalVolume* sAC_log = new G4LogicalVolume(sAC_box,         // Solid
						 mList_->Aerogel, // Material
						 "sAC_log",       // name
						 0,                // FieldMgr
						 0,                // SDetector
						 0);               // ULimits


  double size_AC_Hole[XYZ];
  /*
  double xHoleWidth = 80*mm;
  double yHoleWidth = 40*mm;
  */
  // double xHoleWidth = 40./2.*mm;
  // double yHoleWidth = 30./2.*mm;

  /*
  //for HYPS
  double xHoleWidth = 0./2.*mm;
  double yHoleWidth = 0./2.*mm;
  size_AC_Hole[XCOORD] = xHoleWidth;
  size_AC_Hole[YCOORD] = yHoleWidth;
  size_AC_Hole[ZCOORD] = size_AC[ZCOORD];

  G4Box* sAC_Hole_box = new G4Box("sAC_Hole_box",
				  size_AC_Hole[XCOORD],
				  size_AC_Hole[YCOORD],
				  size_AC_Hole[ZCOORD]);

  G4LogicalVolume* sAC_Hole_log = new G4LogicalVolume(sAC_Hole_box,         // Solid
						      mList_->Air, // Material
						      "sAC_Hole_log",       // name
						      0,                // FieldMgr
						      0,                // SDetector
						      0);               // ULimits


  */

  G4VPhysicalVolume* sAC_phys[NumOfSegmentAC];
  // G4VPhysicalVolume* sAC_Hole_phys;

  // NumOfSegmentAC should be 1
  for (int i=0; i<NumOfSegmentAC; i++) {
      /* Normal Setup */
    G4ThreeVector LocalPos =
      G4ThreeVector(((double)i-((double)(NumOfSegmentAC-1))/2.)*size_AC[XCOORD]*2, 0.*mm, 0.*mm);

    sprintf(name_phys, "sAC_phys%d", i);

    sAC_phys[i] =
      new G4PVPlacement(0,                        // Rotaion
			LocalPos,           // location
			sAC_log,            // Logical
			name_phys,                // name
			AC_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo
    /*
    double xHoleCenter = 0.*mm;
    double yHoleCenter = 0.*mm;

    if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7) {
      //xHoleCenter = 75.*mm - 210.*mm;

      xHoleCenter = -size_ACWall[XCOORD] + size_AC_Hole[XCOORD];
      yHoleCenter = 0.*mm;
    } else {
      //xHoleCenter = -380.*mm - 210.*mm;
      xHoleCenter = -size_ACWall[XCOORD] + size_AC_Hole[XCOORD];
      yHoleCenter = 0.*mm;
    }
    G4ThreeVector LocalPos_Hole = G4ThreeVector(xHoleCenter, yHoleCenter, 0);

    sAC_Hole_phys =
      new G4PVPlacement(0,                        // Rotaion
			LocalPos_Hole,           // location
			sAC_Hole_log,            // Logical
			"sAC_Hole_phys",                // name
			sAC_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo

    */


  }

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  sAC_log->SetSensitiveDetector( counterSD );

}

void CFTDetectorConstruction::ConstructeVeto(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  double size_eVeto[XYZ];

  //for HYPS
  size_eVeto[XCOORD] = 1850.0/2.*mm;
  size_eVeto[YCOORD] = 40.0/2.*mm;
  size_eVeto[ZCOORD] = 20.0/2.*mm;


  //--------------eVeto
  lnum1 = geomMan.GetDetectorId("eVeto");
  G4ThreeVector eVetoPos = geomMan.GetGlobalPosition(lnum1);
  G4ThreeVector offseteVeto(0., geomMan.GetOffset(lnum1)+3.3, 0.);
  eVetoPos += offseteVeto;

  G4RotationMatrix* roteVeto = new G4RotationMatrix();
  roteVeto->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  roteVeto->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* eVeto_box = new G4Box("eVeto_box",
			    size_eVeto[XCOORD],
			    size_eVeto[YCOORD],
			    size_eVeto[ZCOORD]);

  G4LogicalVolume* eVeto_log = new G4LogicalVolume(eVeto_box,         // Solid
						mList_->Air, // Material
						"eVeto_log",       // name
						0,                // FieldMgr
						0,                // SDetector
						0);               // ULimits

  G4VPhysicalVolume* eVeto_phys =
    new G4PVPlacement(roteVeto,         // Rotaion
		      eVetoPos,     // location
		      "eVeto",                                    // name
		      eVeto_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo

  //-------------------------------each eVeto
  double size_eVeto_each[XYZ];
  size_eVeto_each[XCOORD] = size_eVeto[XCOORD]/NumOfSegmenteVeto;
  size_eVeto_each[YCOORD] = size_eVeto[YCOORD];
  size_eVeto_each[ZCOORD] = size_eVeto[ZCOORD];

  G4Box* seVeto_box = new G4Box("seVeto_box",
			     size_eVeto_each[XCOORD],
			     size_eVeto_each[YCOORD],
			     size_eVeto_each[ZCOORD]);

  G4LogicalVolume* seVeto_log = new G4LogicalVolume(seVeto_box,         // Solid
						 mList_->Scin, // Material
						 "seVeto_log",       // name
						 0,                // FieldMgr
						 0,                // SDetector
						 0);               // ULimits


  double size_eVeto_Hole[XYZ];
  /*
  double xHoleWidth = 80*mm;
  double yHoleWidth = 40*mm;
  */
  // double xHoleWidth = 40./2.*mm;
  // double yHoleWidth = 30./2.*mm;

  //for HYPS
  double xHoleWidth = 50./2.*mm;
  double yHoleWidth = 20./2.*mm;
  size_eVeto_Hole[XCOORD] = xHoleWidth;
  size_eVeto_Hole[YCOORD] = yHoleWidth;
  size_eVeto_Hole[ZCOORD] = size_eVeto_each[ZCOORD];

  G4Box* seVeto_Hole_box = new G4Box("seVeto_Hole_box",
				  size_eVeto_Hole[XCOORD],
				  size_eVeto_Hole[YCOORD],
				  size_eVeto_Hole[ZCOORD]);

  G4LogicalVolume* seVeto_Hole_log = new G4LogicalVolume(seVeto_Hole_box,         // Solid
						      mList_->Air, // Material
						      "seVeto_Hole_log",       // name
						      0,                // FieldMgr
						      0,                // SDetector
						      0);               // ULimits


  G4VPhysicalVolume* seVeto_phys[NumOfSegmenteVeto];
  G4VPhysicalVolume* seVeto_Hole_phys;


  // NumOfSegmenteVeto should be 1
  for (int i=0; i<NumOfSegmenteVeto; i++) {
      /* Normal Setup */
    G4ThreeVector LocalPos =
      G4ThreeVector(((double)i-((double)(NumOfSegmenteVeto-1))/2.)*size_eVeto[XCOORD]*2, 0.*mm, 0.*mm);

    sprintf(name_phys, "seVeto_phys%d", i);

    seVeto_phys[i] =
      new G4PVPlacement(0,                        // Rotaion
			LocalPos,           // location
			seVeto_log,            // Logical
			name_phys,                // name
			eVeto_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo
    double xHoleCenter = 0.*mm;
    double yHoleCenter = 0.*mm;

    // if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7) {
    //   //xHoleCenter = 75.*mm - 210.*mm;

    //   xHoleCenter = -size_eVeto[XCOORD] + size_eVeto_Hole[XCOORD];
    //   yHoleCenter = 0.*mm;
    // } else {
    //   //xHoleCenter = -380.*mm - 210.*mm;
    //   xHoleCenter = -size_eVeto[XCOORD] + size_eVeto_Hole[XCOORD];
    //   yHoleCenter = 0.*mm;
    // }
    G4ThreeVector LocalPos_Hole = G4ThreeVector(xHoleCenter, yHoleCenter, 0);

    seVeto_Hole_phys =
      new G4PVPlacement(0,                        // Rotaion
			LocalPos_Hole,           // location
			seVeto_Hole_log,            // Logical
			"seVeto_Hole_phys",                // name
			seVeto_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo


  }
  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  eVeto_log->SetSensitiveDetector( counterSD );

}

void  CFTDetectorConstruction::ConstructTOF(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  double size_TOFWall[XYZ], size_TOF[XYZ];
  //size_TOFWall[XCOORD] = 80.0*(double)(NumOfSegmentTOF+1)/2.0*mm;
  // size_TOFWall[XCOORD] = 120.0*(double)(NumOfSegmentTOF+1)/2.0*mm;
  size_TOFWall[XCOORD] = 93.175*(double)(NumOfSegmentTOF)/2.0*mm; // tmp.
  size_TOFWall[YCOORD] = 1600.0/2.0*mm;
  size_TOFWall[ZCOORD] = 4.0/2.0*mm;

  //size_TOF[XCOORD] = 40.0*mm;
  // size_TOF[XCOORD] = 120.0/2.0*mm;
  size_TOF[XCOORD] = 93.175/2.0*mm;
  size_TOF[YCOORD] = 1600.0/2.0*mm;
  size_TOF[ZCOORD] = 4.0/2*mm;

  double size_TOF_Hole[XYZ];
  size_TOF_Hole[XCOORD] = 20.0*mm;
  size_TOF_Hole[YCOORD] = size_TOF[YCOORD];
  size_TOF_Hole[ZCOORD] = size_TOF[ZCOORD];


  //--------------TOF
  lnum1 = geomMan.GetDetectorId("FTOF");
  G4ThreeVector TOFLocalPos(geomMan.GetOffset(lnum1), 0, 0.);
  G4ThreeVector FTofWallPos = geomMan.Local2GlobalPos(lnum1, TOFLocalPos);

  /*
  G4double ofs = geomMan.GetOffset(lnum1);
  G4ThreeVector FTofWallOffset =
    G4ThreeVector( ofs*cos(geomMan.GetRotAngle2(lnum1)*Deg2Rad),
		   ofs*sin(geomMan.GetRotAngle2(lnum1)*Deg2Rad),
		   0.);
  FTofWallPos += FTofWallOffset;
  */

  G4RotationMatrix* rotTOF = new G4RotationMatrix();
  rotTOF->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotTOF->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* TOF_box = new G4Box("TOF_box",
			      size_TOFWall[XCOORD],
			      size_TOFWall[YCOORD],
			      size_TOFWall[ZCOORD]);

  G4LogicalVolume* TOF_log = new G4LogicalVolume(TOF_box,         // Solid
						  mList_->Air,               // Material
						  "TOF_log",       // name
						  0,                // FieldMgr
						  0,                // SDetector
						  0);               // ULimits

  G4VPhysicalVolume* TOF_phys =
    new G4PVPlacement(rotTOF,         // Rotaion
		      FTofWallPos,     // location
		      "TOF",                                    // name
		      TOF_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo


  //-------------------------------each TOF
  G4Box* sTOF_box = new G4Box("sTOF_box",
			       size_TOF[XCOORD],
			       size_TOF[YCOORD],
			       size_TOF[ZCOORD]);

  G4LogicalVolume* sTOF_log[NumOfSegmentTOF];
  G4VPhysicalVolume* sTOF_phys[NumOfSegmentTOF];

  G4Box* sTOF_Hole_box = new G4Box("sTOF_Hole_box",
				    size_TOF_Hole[XCOORD],
				    size_TOF_Hole[YCOORD],
				    size_TOF_Hole[ZCOORD]);

  G4LogicalVolume* sTOF_Hole_log[NumOfSegmentTOF];
  G4VPhysicalVolume* sTOF_Hole_phys[NumOfSegmentTOF];


  for (int i=1; i<=NumOfSegmentTOF; i++) {
    sprintf(name_log, "sTOF_log%d", i);
    sTOF_log[i-1] = new G4LogicalVolume(sTOF_box,   // Solid
					 mList_->Scin,     // Material
					 name_log,   // name
					 0,          // FieldMgr
					 0,          // SDetector
					 0);         // ULimits
    sTOF_log[i-1]->SetVisAttributes(counter_att);

    //original

    // Span is 75mm, stagging 5mm
    double localYpos = 0*mm;
    /*
    if (i%2 == 0)
      localYpos = -17.5*mm;
    else
      localYpos = +17.5*mm;
    */

    G4ThreeVector tofSegLocalPos =
      //G4ThreeVector((double)(i-12.5)*75.0*mm, 0.0, localYpos);
      // G4ThreeVector((double)(i-20.5)*110.0*mm, 0.0, localYpos);
      G4ThreeVector((double)(i-24-0.5)*2*size_TOF[XCOORD]*mm, 0.0, localYpos);

    /*
    G4ThreeVector tofSegLocalPos =
      G4ThreeVector((double)(i-(double)NumOfSegmentTOF/2.-0.5)*80.0*mm, 0.0, 0.0);
    */

    sprintf(name_phys, "sFTOF_phys%d", i);

    sTOF_phys[i-1] =
      new G4PVPlacement(0,                        // Rotaion
			tofSegLocalPos,           // location
			sTOF_log[i-1],            // Logical
			name_phys,                // name
			TOF_log,                  // Mother
			false,                    // Many
			0);                       // CopyNo

    /*
    if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7 || ReactionMode == 14 || (ReactionMode >= 16 && ReactionMode <= 22)) {
      // origin start from 1
      //if (i>=11&&i<=13) {
      if (i>=14&&i<=16) { // 2018/1/5 E40 detector setup and E07 magnetic field
	sprintf(name_log, "sTOF_Hole_log%d", i);
	sTOF_Hole_log[i-1] = new G4LogicalVolume(sTOF_Hole_box,   // Solid
						  mList_->Air,     // Material
						  //mList_->Scin,     // Material
						  name_log,   // name
						  0,          // FieldMgr
						  0,          // SDetector
						  0);         // ULimits

	G4ThreeVector tofSegHoleLocalPos = G4ThreeVector(0.0, 0.0, 0.0);

	sprintf(name_phys, "sTOF_Hole_phys%d", i);

	sTOF_Hole_phys[i-1] =
	  new G4PVPlacement(0,                        // Rotaion
			    tofSegHoleLocalPos,           // location
			    sTOF_Hole_log[i-1],            // Logical
			    name_phys,                // name
			    sTOF_log[i-1],                  // Mother
			    false,                    // Many
			    0);                       // CopyNo


      }

    } else {
      // origin start from 1
      if (i>=5&&i<=7) {
	sprintf(name_log, "sTOF_Hole_log%d", i);
	sTOF_Hole_log[i-1] = new G4LogicalVolume(sTOF_Hole_box,   // Solid
						  mList_->Air,     // Material
						  //mList_->Scin,     // Material
						  name_log,   // name
						  0,          // FieldMgr
						  0,          // SDetector
						  0);         // ULimits

	G4ThreeVector tofSegHoleLocalPos = G4ThreeVector(0.0, 0.0, 0.0);

	sprintf(name_phys, "sTOF_Hole_phys%d", i);

	sTOF_Hole_phys[i-1] =
	  new G4PVPlacement(0,                        // Rotaion
			    tofSegHoleLocalPos,           // location
			    sTOF_Hole_log[i-1],            // Logical
			    name_phys,                // name
			    sTOF_log[i-1],                  // Mother
			    false,                    // Many
			    0);                       // CopyNo


      }
    }
    */

  }

  G4SDManager *SDMan = G4SDManager::GetSDMpointer();
  if (!counterSD) {
    counterSD = new SKSCounterSD( "counterSD");
    SDMan->AddNewDetector( counterSD );
  }

  for (int i=0; i<NumOfSegmentTOF; i++) {
    sTOF_log[i]->SetSensitiveDetector( counterSD );
  }

}


void  CFTDetectorConstruction::ConstructBeamDump(G4VPhysicalVolume* pMother)
{
  char name_log[30], name_phys[30];
  G4double maxStep;

  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  const DCGeomMan & geomMan=DCGeomMan::GetInstance();
  int lnum1, lnum2;

  G4Colour aqua(0.247, 0.8, 1.0);
  G4VisAttributes* counter_att = new G4VisAttributes(aqua);
  counter_att->SetForceWireframe(true);

  double size_Dump[XYZ];
  size_Dump[XCOORD] = 2000./2.0*mm;
  size_Dump[YCOORD] = 1800.0/2.*mm;
  size_Dump[ZCOORD] = 500.0/2.0*mm;


  //--------------TOF
  lnum1 = geomMan.GetDetectorId("FTOF");
  G4ThreeVector DumpLocalPos(geomMan.GetOffset(lnum1), 0, size_Dump[ZCOORD]+100.);
  G4ThreeVector DumpPos = geomMan.Local2GlobalPos(lnum1, DumpLocalPos);


  G4RotationMatrix* rotDump = new G4RotationMatrix();
  rotDump->rotateX(geomMan.GetRotAngle1(lnum1)*deg);
  rotDump->rotateY(geomMan.GetRotAngle2(lnum1)*deg);

  G4Box* Dump_box = new G4Box("Dump_box",
			      size_Dump[XCOORD],
			      size_Dump[YCOORD],
			      size_Dump[ZCOORD]);

  G4LogicalVolume* Dump_log = new G4LogicalVolume(Dump_box,         // Solid
						  mList_->Fe,               // Material
						  "Dump_log",       // name
						  0,                // FieldMgr
						  0,                // SDetector
						  0);               // ULimits

  G4VPhysicalVolume* Dump_phys =
    new G4PVPlacement(rotDump,         // Rotaion
		      DumpPos,     // location
		      "Dump",                                    // name
		      Dump_log,             // Logical
		      pMother,                     // Mother
		      false,                                    // Many
		      0);                                       // CopyNo


}


void  CFTDetectorConstruction::ConstructS2S(G4VPhysicalVolume* pMother)
{
  // Vis Attribure
  G4Colour blue(0.0, 0.0, 1.0);
  G4Colour red(1.0, 0.0, 0.0);

  G4VisAttributes *magnet_att = new G4VisAttributes(blue);
  magnet_att->SetForceWireframe(true);

  ConfMan *confMan = ConfMan::GetConfManager();
  int ReactionMode = confMan->ReactionMode();

  double S2Spos[XYZ] = {1870.*mm, 2676.*mm, 0.*mm};

  double GapSize = 500.*mm;
  double r1 = 1400.*mm;
  double r2 = 2086.*mm;

  // Inner yoke
  G4Tubs* InnerYoke_tub = new G4Tubs("InnerYoke_tubs", r1, r2, GapSize/2, 270.*deg, 70.*deg);
  G4LogicalVolume*  InnerYoke_log = new G4LogicalVolume(InnerYoke_tub, mList_->Fe, "InnerYoke_log",0,0,0);
  InnerYoke_log->SetVisAttributes(magnet_att);

  G4PVPlacement* InnerYoke_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]),
			"InnerYoke_phys",
			InnerYoke_log,
			pMother,
			false,
			0);

  double r3 = 2498.*mm;


  G4Tubs* InnerCoil_tub = new G4Tubs("InnerCoil_tubs", r2, r3, GapSize/4, 270.*deg, 70.*deg);

  // Inner Up Coil
  G4LogicalVolume*  InnerUpCoil_log = new G4LogicalVolume(InnerCoil_tub, mList_->Cu, "InnerUpCoil_log",0,0,0);
  InnerUpCoil_log->SetVisAttributes(red);

  G4PVPlacement* InnerUpCoil_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]+GapSize/4),
			"InnerUpCoil_phys",
			InnerUpCoil_log,
			pMother,
			false,
			0);


  // Inner Down Coil
  G4LogicalVolume*  InnerDownCoil_log = new G4LogicalVolume(InnerCoil_tub, mList_->Cu, "InnerDownCoil_log",0,0,0);
  InnerDownCoil_log->SetVisAttributes(red);

  G4PVPlacement* InnerDownCoil_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]-GapSize/4),
			"InnerDownCoil_phys",
			InnerDownCoil_log,
			pMother,
			false,
			0);



  double r4 = 3498.*mm;
  double r5 = 3910.*mm;

  G4Tubs* OuterCoil_tub = new G4Tubs("OuterCoil_tubs", r4, r5, GapSize/4, 270.*deg, 70.*deg);

  // Outer Up Coil
  G4LogicalVolume*  OuterUpCoil_log = new G4LogicalVolume(OuterCoil_tub, mList_->Cu, "OuterUpCoil_log",0,0,0);
  OuterUpCoil_log->SetVisAttributes(red);

  G4PVPlacement* OuterUpCoil_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]+GapSize/4),
			"OuterUpCoil_phys",
			OuterUpCoil_log,
			pMother,
			false,
			0);

  // Outer Down Coil
  G4LogicalVolume*  OuterDownCoil_log = new G4LogicalVolume(OuterCoil_tub, mList_->Cu, "OuterDownCoil_log",0,0,0);
  OuterDownCoil_log->SetVisAttributes(red);

  G4PVPlacement* OuterDownCoil_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]-GapSize/4),
			"OuterDownCoil_phys",
			OuterDownCoil_log,
			pMother,
			false,
			0);

  // Outer yoke
  double r6 = 4598.*mm;


  G4Tubs* OuterYoke_tub = new G4Tubs("OuterYoke_tubs", r5, r6, GapSize/2, 270.*deg, 70.*deg);
  G4LogicalVolume*  OuterYoke_log = new G4LogicalVolume(OuterYoke_tub, mList_->Fe, "OuterYoke_log",0,0,0);
  OuterYoke_log->SetVisAttributes(magnet_att);

  G4PVPlacement* OuterYoke_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]),
			"OuterYoke_phys",
			OuterYoke_log,
			pMother,
			false,
			0);


  // Upper Yoke
  double Z_UpperYoke = 800.*mm;

  G4Tubs* UpperYoke_tub = new G4Tubs("UpperYoke_tubs", r1, r6, Z_UpperYoke/2, 270.*deg, 70.*deg);
  G4LogicalVolume*  UpperYoke_log = new G4LogicalVolume(UpperYoke_tub, mList_->Fe, "UpperYoke_log",0,0,0);
  UpperYoke_log->SetVisAttributes(magnet_att);

  G4PVPlacement* UpperYoke_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]+GapSize/2.+Z_UpperYoke/2),
			"UpperYoke_phys",
			UpperYoke_log,
			pMother,
			false,
			0);

  // Bottom Yoke
  double Z_BottomYoke = Z_UpperYoke;

  G4Tubs* BottomYoke_tub = new G4Tubs("BottomYoke_tubs", r1, r6, Z_BottomYoke/2, 270.*deg, 70.*deg);
  G4LogicalVolume*  BottomYoke_log = new G4LogicalVolume(BottomYoke_tub, mList_->Fe, "BottomYoke_log",0,0,0);
  BottomYoke_log->SetVisAttributes(magnet_att);

  G4PVPlacement* BottomYoke_phys
    = new G4PVPlacement(0,
			G4ThreeVector(S2Spos[XCOORD],
				      S2Spos[YCOORD],
				      S2Spos[ZCOORD]-GapSize/2.-Z_UpperYoke/2),
			"BottomYoke_phys",
			BottomYoke_log,
			pMother,
			false,
			0);

  double WallThickness = 100.*mm;
  double Wall1Size[XYZ] = {WallThickness/2, 500./2*mm, GapSize/2};
  G4Box* Wall1_box = new G4Box("Wall1_box",
			       Wall1Size[XCOORD],
			       Wall1Size[YCOORD],
			       Wall1Size[ZCOORD]);

  double Wall1Pos[XYZ] = {S2Spos[XCOORD]+Wall1Size[XCOORD], S2Spos[YCOORD]-r3-Wall1Size[YCOORD], S2Spos[ZCOORD]};

  if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7 || ReactionMode == 14 || (ReactionMode >= 16 && ReactionMode <= 22)) {
    Wall1Pos[YCOORD] = S2Spos[YCOORD]-r4+Wall1Size[YCOORD];
  }

  G4LogicalVolume*  Wall1_log = new G4LogicalVolume(Wall1_box, mList_->Fe, "Wall1_log",0,0,0);
  Wall1_log->SetVisAttributes(magnet_att);
  G4PVPlacement* Wall1_phys
    = new G4PVPlacement(0,
			G4ThreeVector(Wall1Pos[XCOORD],
				      Wall1Pos[YCOORD],
				      Wall1Pos[ZCOORD]),
			"Wall1_phys",
			Wall1_log,
			pMother,
			false,
			0);

  double WallGapSize = 300.*mm;
  double Wall2Size[XYZ] = {WallThickness/2, 500./2*mm, (GapSize/2-WallGapSize/2)/2};

  G4Box* Wall2_box = new G4Box("Wall2_box",
			       Wall2Size[XCOORD],
			       Wall2Size[YCOORD],
			       Wall2Size[ZCOORD]);


  double Wall2UPos[XYZ] = {S2Spos[XCOORD]+Wall1Size[XCOORD], S2Spos[YCOORD]-r3-Wall1Size[YCOORD]*2-Wall2Size[YCOORD], S2Spos[ZCOORD]+GapSize/2-Wall2Size[ZCOORD]};

  if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7 || ReactionMode == 14 || (ReactionMode >= 16 && ReactionMode <= 22)) {
    Wall2UPos[YCOORD] = S2Spos[YCOORD]-r4+Wall1Size[YCOORD]*2+Wall2Size[YCOORD];
  }



  G4LogicalVolume*  Wall2U_log = new G4LogicalVolume(Wall2_box, mList_->Fe, "Wall2U_log",0,0,0);
  Wall2U_log->SetVisAttributes(magnet_att);
  G4PVPlacement* Wall2U_phys
    = new G4PVPlacement(0,
			G4ThreeVector(Wall2UPos[XCOORD],
				      Wall2UPos[YCOORD],
				      Wall2UPos[ZCOORD]),
			"Wall2U_phys",
			Wall2U_log,
			pMother,
			false,
			0);

  double Wall2DPos[XYZ] = {S2Spos[XCOORD]+Wall1Size[XCOORD], S2Spos[YCOORD]-r3-Wall1Size[YCOORD]*2-Wall2Size[YCOORD], S2Spos[ZCOORD]-GapSize/2+Wall2Size[ZCOORD]};

  if (ReactionMode == 4 || ReactionMode == 6 || ReactionMode == 7 || ReactionMode == 14 || (ReactionMode >= 16 && ReactionMode <= 22)) {
    Wall2DPos[YCOORD] = S2Spos[YCOORD]-r4+Wall1Size[YCOORD]*2+Wall2Size[YCOORD];
  }


  G4LogicalVolume*  Wall2D_log = new G4LogicalVolume(Wall2_box, mList_->Fe, "Wall2D_log",0,0,0);
  Wall2D_log->SetVisAttributes(magnet_att);
  G4PVPlacement* Wall2D_phys
    = new G4PVPlacement(0,
			G4ThreeVector(Wall2DPos[XCOORD],
				      Wall2DPos[YCOORD],
				      Wall2DPos[ZCOORD]),
			"Wall2D_phys",
			Wall2D_log,
			pMother,
			false,
			0);



}


G4double CFTDetectorConstruction::GetDistance(G4ThreeVector pos1, G4ThreeVector pos2)
{
  G4double distance;

  distance = sqrt (pow(pos1.x()-pos2.x(),2.0)
		   + pow(pos1.y()-pos2.y(),2.0));

  return distance;
}
