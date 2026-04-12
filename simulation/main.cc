/*
  Main.cc
  2007/4  K.Shirotori
*/

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4UIterminal.hh"
#include "G4UItcsh.hh"

#include "Randomize.hh"
#include "GetNumberFromKernelEntropyPool.hh"
#include "ConfMan.hh"

#include "G4VisExecutive.hh"

#include "G4UIExecutive.hh"

// #include "RunAction.hh"
#include "CFTDetectorConstruction.hh"
#include "PhysicsList.hh"
#include "CFTPrimaryGeneratorAction.hh"

#include "EventAction.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "Analysis.hh"

#define RANDOM_SEED 1

// enum EArgv { kProcess, kConfFile, kOutFile, kG4Macro, kArgc };
enum EArgv { kProcess, kConfFile, kRootFile, kDatFile, kG4Macro, kArgc };

//_____________________________________________________________________________
int main( int argc, char **argv )
{

  if (argc != kArgc-1 && argc != kArgc) {
    G4cout << "Usage: " << argv[kProcess]
           // << " [ConfFile] [OutputName] (G4Macro)" << G4endl;
           << " [ConfFile] [OutFile(.root)] [OutFile(.dat)] (G4Macro <default: init_vis.mac>)" << G4endl;
    return EXIT_SUCCESS;
  }

  ConfMan *confManager = new ConfMan(argv[kConfFile]);
  if (!confManager->Initialize()) {
    return EXIT_FAILURE;
  }

  // Initialization on Random Number Generator
  int32_t initSeed = 1;
  G4cout << "Initial Seed for Geant4 Randomize" << G4endl;
  G4cout << "Randomized or Fixed : ";
#if RANDOM_SEED
  initSeed = GetIntFromKernelEntropyPool() & 0x7FFFFFFF;
  CLHEP::HepRandom::setTheSeed(initSeed);
  G4cout << "Randomized" << G4endl;
#else
  CLHEP::HepRandom::setTheSeed(initSeed);
  G4cout << "Fixed" << G4endl;
#endif

  G4cout << "Initial Seed Value : " << CLHEP::HepRandom::getTheSeed() << G4endl;

  CLHEP::HepRandom::showEngineStatus();

  //----- Analysis -----
  Analysis *anaManager = new Analysis(argv[kRootFile]);
  anaManager->SetDataFile(argv[kDatFile]);

  //----- Run Manager -----
  G4RunManager *runManager = new G4RunManager();
  runManager->SetUserInitialization(new CFTDetectorConstruction);
  runManager->SetUserInitialization(new PhysicsList);

  G4VisManager* visManager = new G4VisExecutive;
  visManager->SetVerboseLevel(0);
  visManager->Initialize();

  //----- Primary Generator -----
  CFTPrimaryGeneratorAction *priGen = new CFTPrimaryGeneratorAction(anaManager);
  // RunAction *runAction = new RunAction(anaManager);
  EventAction *eventAction = new EventAction(anaManager);
  SteppingAction* stepAction = new SteppingAction(anaManager);
  TrackingAction* trackAction = new TrackingAction(anaManager);

#if 0
  char *effstudydata;
  if (argc >= 6
      && (confManager->ReactionMode() == 95
	  || confManager->ReactionMode() == 99
	  || confManager->ReactionMode() == 100
	  || confManager->ReactionMode() == 101)) {
    effstudydata = argv[5];
    char outputname[100];
    sprintf(outputname, "data/OutputEffStudy");
    // char *p = strchr(effstudydata, '.');
    // strcat(outputname, p);
    // anaManager->OpenEffStudyFile( outputname );
    priGen->OpenEffStudyFile( effstudydata );
  }
#endif

  runManager->SetUserAction( priGen );
  // runManager->SetUserAction( runAction );
  runManager->SetUserAction( eventAction );
  runManager->SetUserAction( stepAction );
  runManager->SetUserAction( trackAction );
  runManager->Initialize();

  G4UImanager *uiManager = G4UImanager::GetUIpointer();
  if ( argc == kArgc - 1) {
    auto ui = new G4UIExecutive(argc, argv);
    uiManager->ApplyCommand("/control/execute init_vis.mac");
    ui->SessionStart();
    delete ui;
  } else {
    //--- Batch Mode ---
    G4String command = "/control/execute ";
    G4String fileName = argv[kG4Macro];
    uiManager->ApplyCommand(command + fileName);
  }

  delete anaManager;
  delete visManager;
  delete runManager;
  return EXIT_SUCCESS;
}
