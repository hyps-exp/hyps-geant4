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

#ifndef SKSField_H
#define SKSField_H 1

#include "globals.hh"
#include "G4MagneticField.hh"
#include "common.h"

#include <string>
#include <vector>

class SksFieldMap;
class FieldElements;

class SKSField : public G4MagneticField
{
  public:
    SKSField(const std::string &FieldMapName);
    ~SKSField();

    void SetFileName( const char *filename ) { FieldMapName_=filename; }
    void SetFileName( const std::string &filename ) { FieldMapName_=filename; }

    bool Initialize( void );

    void GetFieldValue( const  double Point[3],
                               double *Bfield ) const;
    void SetMapScale(double scale) {MapScale_ = scale;}
  private:
    std::string FieldMapName_;
    SksFieldMap *Sksmap_;
    typedef std::vector <FieldElements *> FEContainer;
    typedef std::vector <FieldElements *>::const_iterator FEIterator;

    FEContainer elemList_;
  double MapScale_;

};

#endif

