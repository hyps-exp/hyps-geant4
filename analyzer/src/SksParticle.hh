/*
  SksParticle.hh

*/

#ifndef SksParticle_h

#define SksParticle_h 1

#include "ThreeVector.hh"


class SksTrack;
class HodoCluster;
class DCLocalTrack;

class SksParticle
{
public:
  SksParticle( SksTrack *track, HodoCluster *Tof, DCLocalTrack *Out/*, HodoCluster *Lc*/ )
    : Track_(track), Tof_(Tof), Out_(Out)/*, Lc_(Lc)*/
  {}
  ~SksParticle() {}
  SksParticle( SksTrack *track, HodoCluster *Tof/*, HodoCluster *Lc*/ )
    : Track_(track), Tof_(Tof)/*, Lc_(Lc)*/
  {}
private:
  SksTrack *Track_;
  HodoCluster *Tof_/*, *Lc_*/;
  DCLocalTrack *Out_;

public:
  SksTrack * GetTrack( void ) { return Track_; }
  HodoCluster * GetTof( void ) { return Tof_; }
  DCLocalTrack * GetOut( void ) { return Out_; }
  //HodoCluster * GetLc( void ) { return Lc_; }

  ThreeVector Momentum( void ) const;
  double Polarity( void ) const;
  ThreeVector Position( void ) const;
  double PathLengthToTOF( void ) const;
  //double PathLengthTotal( void ) const;
  double MassSquare( void ) const; 
  double yTof( void ) const; 
  double xTof( void ) const; 
};  

#endif
