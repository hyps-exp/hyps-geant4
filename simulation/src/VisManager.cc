/*
  VisManager.cc
  2007/4  K.Shirotori
*/

#include "VisManager.hh"
#include "G4OpenGL.hh"
#include "G4ASCIITree.hh"
#include "G4DAWNFILE.hh"
//#include "G4GAGTree.hh"
// #include "G4HepRepFile.hh"
// #include "G4HepRep.hh"
#include "G4RayTracer.hh"
// #include "G4VRML1File.hh"
#include "G4VRML2File.hh"

#include "G4TrajectoryFilterFactories.hh"
#include "G4TrajectoryModelFactories.hh"

#ifdef G4VIS_USE_DAWN
#include "G4FukuiRenderer.hh"
#endif

#ifdef G4VIS_USE_DAWNFILE
#include "G4DAWNFILE.hh"
#endif

#ifdef G4VIS_USE_OPACS
#include "G4Wo.hh"
#include "G4Xo.hh"
#endif

#ifdef G4VIS_USE_OPENGLX
#include "G4OpenGLImmediateX.hh"
#include "G4OpenGLStoredX.hh"
#endif

#ifdef G4VIS_USE_OPENGLXM
#include "G4OpenGLImmediateXm.hh"
#include "G4OpenGLStoredXm.hh"
#endif

#ifdef G4VIS_USE_OIX
#include "G4OpenInventorX.hh"
#endif

#ifdef G4VIS_USE_RAYX
#include "G4RayX.hh"
#endif

#ifdef G4VIS_USE_VRML
#include "G4VRML1.hh"
#include "G4VRML2.hh"
#endif

#ifdef G4VIS_USE_VRMLFILE
// #include "G4VRML1File.hh"
#include "G4VRML2File.hh"
#endif

#ifdef G4VIS_USE_RAYTRACER
#include "G4RayTracer.hh"
#endif

VisManager::VisManager() 
  : G4VisManager()
{}

VisManager::~VisManager() 
{}

void VisManager::RegisterGraphicsSystems() 
{

#ifdef G4VIS_USE_DAWN
  //G4cout << "*  Visualizaton Driver: Fukui Renderer" << G4endl;
  //RegisterGraphicsSystem(new G4FukuiRenderer);
#endif

#ifdef G4VIS_USE_DAWNFILE
  G4cout << "*  Visualizaton Driver: DAWN" << G4endl;
  RegisterGraphicsSystem(new G4DAWNFILE);
#endif

#ifdef G4VIS_USE_OPACS
  RegisterGraphicsSystem(new G4Wo);
  RegisterGraphicsSystem(new G4Xo);
#endif

#ifdef G4VIS_USE_OPENGLX
  RegisterGraphicsSystem(new G4OpenGLImmediateX);
  RegisterGraphicsSystem(new G4OpenGLStoredX);
#endif

#ifdef G4VIS_USE_OPENGLXM
  RegisterGraphicsSystem(new G4OpenGLImmediateXm);
  RegisterGraphicsSystem(new G4OpenGLStoredXm);
#endif

#ifdef G4VIS_USE_OIX
  RegisterGraphicsSystem(new G4OpenInventorX);
#endif

#ifdef G4VIS_USE_RAYX
  RegisterGraphicsSystem(new G4RayX);
#endif

#ifdef G4VIS_USE_VRML
  G4cout << "*  Visualizaton Driver: VRML1/2" << G4endl;
  RegisterGraphicsSystem(new G4VRML1);
  RegisterGraphicsSystem(new G4VRML2);
#endif

#ifdef G4VIS_USE_VRMLFILE
  G4cout << "*  Visualizaton Driver: VRML1/2 File" << G4endl;
  // RegisterGraphicsSystem(new G4VRML1File);
  RegisterGraphicsSystem(new G4VRML2File);
#endif
}

void VisManager::RegisterModelFactories()
{
   // Trajectory draw models
   RegisterModelFactory(new G4TrajectoryGenericDrawerFactory());       
   RegisterModelFactory(new G4TrajectoryDrawByChargeFactory());
   RegisterModelFactory(new G4TrajectoryDrawByParticleIDFactory());
   RegisterModelFactory(new G4TrajectoryDrawByOriginVolumeFactory());  

   // Trajectory filter models
   RegisterModelFactory(new G4TrajectoryChargeFilterFactory());
   RegisterModelFactory(new G4TrajectoryParticleFilterFactory());
   RegisterModelFactory(new G4TrajectoryOriginVolumeFilterFactory());
}
