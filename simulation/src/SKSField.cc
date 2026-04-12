//
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

#include "G4ThreeVector.hh"

#include "SKSField.hh"
#include "SksFieldMap.hh"
#include "FieldElements.hh"

#include "G4SystemOfUnits.hh"

SKSField::SKSField(const std::string &FieldMapName)
  : FieldMapName_(FieldMapName), MapScale_(2.0)
{
  const std::string funcName="[SKSField::SKSField]";
  if (!Initialize()) {
    std::cerr << funcName << " Error Initialize" << std::endl;
    exit(1);
  }
}

SKSField::~SKSField()
{
  delete Sksmap_;
}

bool SKSField::Initialize( void )
{
  SksFieldMap *pSks = new SksFieldMap( FieldMapName_.c_str() );
  if( pSks ){
    /*
      delete Sksmap_;
    */
    Sksmap_ = pSks;
  }
  if( Sksmap_ )
    return Sksmap_->Initialize();
  else
    return false;
}

void SKSField::GetFieldValue(const double Point[3],double *Bfield) const
{
  /*
  if (Point[0]<-400) {
    Bfield[0] = 0.*tesla;
    Bfield[1] = 0.*tesla;
    Bfield[2] = 0.*tesla;
  return;
  }
  */

    
  G4ThreeVector field( 0.0, 0.0, 0.0 );
  if( Sksmap_ ){
    double p[3], b[3];
    p[0]=Point[0]/mm*0.1; p[1]=Point[1]/mm*0.1; p[2]=Point[2]/mm*0.1;

    if (p[1]<-40.)
      p[1] = -40.;
    else if (p[1]>40)
      p[1] = 40.;
    
    if( Sksmap_->GetFieldValue( p, b ) ){
      field.setX(b[0]); field.setY(b[1]); field.setZ(b[2]);
    }
  }


#if 1
  G4ThreeVector position = G4ThreeVector(Point[0], Point[1], Point[2]);
  FEIterator end = elemList_.end();
  for( FEIterator itr=elemList_.begin(); itr!=end; ++itr ){
    if( (*itr)->ExistField( position ) )
      field += (*itr)->GetField( position );
  }
#endif
  /*
  Bfield[0] = -0.825/0.55*field.x()*tesla;
  Bfield[1] = -0.825/0.55*field.y()*tesla;
  Bfield[2] = -0.825/0.55*field.z()*tesla;
  */


  Bfield[0] = MapScale_*field.x()*tesla;
  Bfield[1] = MapScale_*field.y()*tesla;
  Bfield[2] = MapScale_*field.z()*tesla;

  /*
  double factor=0.7;
  Bfield[0] = factor*field.x()*tesla;
  Bfield[1] = factor*field.y()*tesla;
  Bfield[2] = factor*field.z()*tesla;
  */

  /*
  double S2SPos_X = 1870.*mm, S2SPos_Y = 2676.*mm;
  double GapSize = 500.*mm;
  double tmp_val = (Point[0]-S2SPos_X)*(Point[0]-S2SPos_X)+(Point[1]-S2SPos_Y)*(Point[1]-S2SPos_Y);
  double r1 = 2498.*mm, r2 = 3498.*mm;

  if (Point[0]>S2SPos_X &&  tmp_val >= r1*r1 && tmp_val <= r2*r2 && fabs(Point[2])<GapSize/2) {
    Bfield[0] = 0*tesla;
    Bfield[1] = 0*tesla;
    Bfield[2] = -2*tesla;
  }
  */


}


