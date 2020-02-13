#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4Sphere.hh>
#include <G4NistManager.hh>
#include <G4ParticleGun.hh>
#include <G4ParticleTable.hh>
#include <G4ParticleDefinition.hh>
#include <G4PVPlacement.hh>
#include <G4RunManager.hh>
#include <G4Event.hh>
#include <G4SDManager.hh>
#include <G4VisExecutive.hh>
#include <G4VisManager.hh>
#include <G4VSensitiveDetector.hh>
#include <G4VUserPhysicsList.hh>
#include <G4VUserDetectorConstruction.hh>
#include <G4VUserPrimaryGeneratorAction.hh>
#include <G4UIExecutive.hh>
#include <G4UImanager.hh>
#include <G4UIsession.hh>
#include <G4UserEventAction.hh>
#include <QGSP_BERT_HP.hh>
#include <G4RandomDirection.hh>
#include <G4VPhysicsConstructor.hh>
#include <G4UserSteppingAction.hh>
#include <G4OpticalPhoton.hh>
#include <G4Scintillation.hh>
#include <G4ProcessManager.hh>
#include <G4OpAbsorption.hh>
#include <G4OpRayleigh.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4MaterialPropertiesTable.hh>
#include <G4MaterialPropertyVector.hh>
#include <G4Material.hh>
#include <G4UserRunAction.hh>
#include <G4Run.hh>
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4BosonConstructor.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4ProcessManager.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4UserEventAction.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"
#include "globals.hh"
#include "G4UImessenger.hh"
#include <QGSP_BERT_HP.hh>
#include <G4ProcessManager.hh>
#include <G4Cerenkov.hh>
#include <G4Scintillation.hh>
#include <G4OpAbsorption.hh>
#include <G4OpRayleigh.hh>
#include <G4OpMieHG.hh>
#include <G4OpBoundaryProcess.hh>
#include <G4LossTableManager.hh>
#include <G4EmSaturation.hh>
#include <G4IonTable.hh>


using namespace std;
using CLHEP::eV;
using CLHEP::gram;
using CLHEP::cm3;
using CLHEP::MeV;
using CLHEP::cm;
using CLHEP::mm;
using CLHEP::degree;
using CLHEP::ns;
using CLHEP::m;

int nphotons; 
std::vector<int> histogram;
const int nbins=1000000;
const int bin_size = 40;
const int MassNumber = 107;
const double Energy = 0.66166;
const int ChargeNumber = 47;
char file_output[] = "/home/golodovka/Рабочий стол/661.66.txt";


class my_event_action : public G4UserEventAction
{
public:
  my_event_action()
  {
    histogram.resize(nbins);
  }
  void BeginOfEventAction (const G4Event*)
  {
    // int ID = G4RunManager::GetRunManager()->GetCurrentEvent ()->GetEventID ();
    // G4cout << ">>> event has started " << ID << "\n" << std::flush;
    // fprintf(stderr,">>> event %d has started\n",ID);
    nphotons=0;
  }
  
  void EndOfEventAction (const G4Event*)
  {
    int ID = G4RunManager::GetRunManager()->GetCurrentEvent ()->GetEventID ();
    // G4cout << ">>> event has ended " << ID << "\n" << std::flush;
    fprintf(stderr,">>> event %d has finished ",ID);
    fprintf(stderr,"with number of opticalphotons %d \n",nphotons);
    int bin = nphotons/bin_size;
    if (bin<nbins)
    {
      histogram[bin]++;
    }
  }
};

class my_run_action : public G4UserRunAction
{
  void BegionOfRunAction (const G4Run*)
  {
    CLHEP::HepRandom::setTheSeed(1);
    histogram.clear();
  }
  void EndOfRunAction (const G4Run*)
  {
    ofstream fout(file_output);
    ifstream fin;
    std::vector<int>::iterator it;
    for (int i=0;i<histogram.size(); i++)
    {
      fout << histogram[i]<<"\n";
    }
    fout.close();
    histogram.clear();

  }

};

class my_stepping_action : public G4UserSteppingAction
{
  public: 
    void UserSteppingAction (const G4Step *step)
    { 

      const std::vector<const G4Track*> *sv = step ->GetSecondaryInCurrentStep();
      int nsecondary = step->GetNumberOfSecondariesInCurrentStep();

      for (int j=0;j<nsecondary;j++)
        {
          const G4Track *t = sv->at (j);
          G4ParticleDefinition *def=t->GetDefinition();

          if (def==G4OpticalPhoton::Definition())
          { 
            nphotons ++;
            // double E = t->GetKineticEnergy() / eV;
            // printf ("optical photon with Energy = %f eV \n", E );
          }
        }


    }
};

class OpticalPhysics : public QGSP_BERT_HP
{
public:
  void ConstructProcess()
  {
    QGSP_BERT_HP::ConstructProcess();
    ConstructOp();
  }
  
private:

  void ConstructOp();
  static G4ThreadLocal G4int fVerboseLevel;
  // static G4ThreadLocal G4int fMaxNumPhotonStep;

  // static G4ThreadLocal G4Cerenkov* fCerenkovProcess;
  static G4ThreadLocal G4Scintillation* fScintillationProcess;
  static G4ThreadLocal G4OpAbsorption* fAbsorptionProcess;
  static G4ThreadLocal G4OpRayleigh* fRayleighScatteringProcess;
  static G4ThreadLocal G4OpMieHG* fMieHGScatteringProcess;
  // static G4ThreadLocal G4OpBoundaryProcess* fBoundaryProcess;
};

G4ThreadLocal G4int OpticalPhysics::fVerboseLevel = 1;
// G4ThreadLocal G4int OpticalPhysics::fMaxNumPhotonStep = 20;
// G4ThreadLocal G4Cerenkov* OpticalPhysics::fCerenkovProcess = 0;
G4ThreadLocal G4Scintillation* OpticalPhysics::fScintillationProcess = 0;
G4ThreadLocal G4OpAbsorption* OpticalPhysics::fAbsorptionProcess = 0;
G4ThreadLocal G4OpRayleigh* OpticalPhysics::fRayleighScatteringProcess = 0;
G4ThreadLocal G4OpMieHG* OpticalPhysics::fMieHGScatteringProcess = 0;
// G4ThreadLocal G4OpBoundaryProcess* OpticalPhysics::fBoundaryProcess = 0;

void OpticalPhysics::ConstructOp()
{
  // fCerenkovProcess = new G4Cerenkov("Cerenkov");
  // fCerenkovProcess->SetMaxNumPhotonsPerStep(fMaxNumPhotonStep);
  // fCerenkovProcess->SetMaxBetaChangePerStep(10.0);
  // fCerenkovProcess->SetTrackSecondariesFirst(true);
  fScintillationProcess = new G4Scintillation("Scintillation");
  fScintillationProcess->SetScintillationYieldFactor(1.);
  fScintillationProcess->SetTrackSecondariesFirst(true);
  fAbsorptionProcess = new G4OpAbsorption();
  fRayleighScatteringProcess = new G4OpRayleigh();
  fMieHGScatteringProcess = new G4OpMieHG();
  // fBoundaryProcess = new G4OpBoundaryProcess();

  // fCerenkovProcess->SetVerboseLevel(fVerboseLevel);
  fScintillationProcess->SetVerboseLevel(fVerboseLevel);
  fAbsorptionProcess->SetVerboseLevel(fVerboseLevel);
  fRayleighScatteringProcess->SetVerboseLevel(fVerboseLevel);
  fMieHGScatteringProcess->SetVerboseLevel(fVerboseLevel);
  // fBoundaryProcess->SetVerboseLevel(fVerboseLevel);
  
  // // Use Birks Correction in the Scintillation process
  // if(G4Threading::IsMasterThread())
  // {
  //   G4EmSaturation* emSaturation =
  //             G4LossTableManager::Instance()->EmSaturation();
  //     fScintillationProcess->AddSaturation(emSaturation);
  // }

  auto particleIterator=GetParticleIterator();
  particleIterator->reset();
  while( (*particleIterator)() ){
    G4ParticleDefinition* particle = particleIterator->value();
    G4ProcessManager* pmanager = particle->GetProcessManager();
    G4String particleName = particle->GetParticleName();
    // if (fCerenkovProcess->IsApplicable(*particle)) {
    //   pmanager->AddProcess(fCerenkovProcess);
    //   pmanager->SetProcessOrdering(fCerenkovProcess,idxPostStep);
    // }
    if (fScintillationProcess->IsApplicable(*particle)) {
      pmanager->AddProcess(fScintillationProcess);
      pmanager->SetProcessOrderingToLast(fScintillationProcess, idxAtRest);
      pmanager->SetProcessOrderingToLast(fScintillationProcess, idxPostStep);
    }
    if (particleName == "opticalphoton") {
      G4cout << " AddDiscreteProcess to OpticalPhoton " << G4endl;
      pmanager->AddDiscreteProcess(fAbsorptionProcess);
      pmanager->AddDiscreteProcess(fRayleighScatteringProcess);
      pmanager->AddDiscreteProcess(fMieHGScatteringProcess);
      // pmanager->AddDiscreteProcess(fBoundaryProcess);
    }
  }
}

class sensitive_detector : public G4VSensitiveDetector
{
private:
  double edep;
public:
  sensitive_detector (G4String name) : G4VSensitiveDetector (name)
  {
  }

  virtual void Initialize (G4HCofThisEvent*)
  {
    edep = 0;
  }

  G4bool ProcessHits (G4Step *step, G4TouchableHistory*)
  {
    edep += step->GetTotalEnergyDeposit ();
    return true;
  }
  
  virtual void EndOfEvent (G4HCofThisEvent*)
  { 

    int e = G4RunManager::GetRunManager ()->GetCurrentEvent ()->GetEventID ();
    if (edep > 0)
      {
        G4cout << ">>> event " << e
               << " energy deposition " << edep / MeV << " MeV"
               << endl;
      }
  }
};

class detector_construction : public G4VUserDetectorConstruction
{
private:
  G4VPhysicalVolume *Construct ()
  {
    const G4int NUMENTRIES = 1;
    G4double Scnt_PP[NUMENTRIES] = 
    {0.329*eV};
    G4double Scnt_FAST[NUMENTRIES] = 
    {1.*ns};
    G4double Scnt_SLOW[NUMENTRIES] = 
    {10.*ns};
    G4double Scnt_ABSL[NUMENTRIES] = 
    {0.5*m};
    G4double Scnt_RINDEX[NUMENTRIES] = 
    {2.6};
    std::vector<G4String> ee = {"Zn","Se"};
    std::vector<G4int> aa= {1,1}; 
    G4Material* Scnt = G4NistManager::Instance()->ConstructNewMaterial("ZINC_SELENIDE",ee,aa,5.27*gram/cm3);
    G4MaterialPropertiesTable* Scnt_MPT = new G4MaterialPropertiesTable();
    Scnt_MPT->AddProperty("FASTCOMPONENT", Scnt_PP, Scnt_FAST, NUMENTRIES);
    Scnt_MPT->AddProperty("ABSLENGHT", Scnt_PP,Scnt_ABSL,NUMENTRIES)->SetSpline(true);
    Scnt_MPT->AddProperty("RINDEX", Scnt_PP,Scnt_RINDEX,NUMENTRIES)->SetSpline(true);
    Scnt_MPT->AddProperty("SLOWCOMPONENT", Scnt_PP, Scnt_SLOW, NUMENTRIES);
    Scnt_MPT->AddConstProperty("SCINTILLATIONYIELD", 28300./MeV);
    Scnt_MPT->AddConstProperty("RESOLUTIONSCALE", 1700.);
    Scnt_MPT->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
    Scnt_MPT->AddConstProperty("SLOWTIMECONSTANT", 10.*ns);
    Scnt_MPT->AddConstProperty("YIELDRATIO", 0.8);
    Scnt->SetMaterialPropertiesTable(Scnt_MPT);

    G4Material *au = G4NistManager::Instance()->FindOrBuildMaterial ("G4_Au");
    G4Material *air =G4NistManager::Instance()->FindOrBuildMaterial ("G4_AIR");
    // G4Material *vacum =G4NistManager::Instance()->FindOrBuildMaterial ("G4_VACUM");
    G4Material *ge = G4NistManager::Instance()->FindOrBuildMaterial ("G4_Ge");
    G4Material *zn = G4NistManager::Instance()->FindOrBuildMaterial ("G4_Zn");

    G4Box *world_box =     
      new G4Box ("world_box", 50*cm, 50*cm, 50*cm);
    G4LogicalVolume *world_logical_volume =
      new G4LogicalVolume (world_box, air, "world_logical_vol");
    G4VPhysicalVolume *world_physical_volume = 
      new G4PVPlacement (0, G4ThreeVector (), world_logical_volume,
                         "world_physical_vol", 0, false, 0);

    G4Box *scintillator = 
      new G4Box ("scintillator", 1*cm, 1*cm, 1*mm);
    G4LogicalVolume *scintillator_logical_volume =
      new G4LogicalVolume (scintillator, Scnt, "scintillator_logical_vol");
    G4VPhysicalVolume *scintillator_physical_volume =
      new G4PVPlacement (0, G4ThreeVector (0, 0, 2.85*cm), scintillator_logical_volume,
                         "scintillator_physical_vol", world_logical_volume, false, 0);

    G4Sphere *target_box =  
      new G4Sphere ("target_box", 0, 0.5*mm, 0, 360*degree, 0, 180*degree);
    G4LogicalVolume *target_logical_volume =
      new G4LogicalVolume (target_box, air, "target_logical_vol");
    G4VPhysicalVolume *target_physical_volume =
      new G4PVPlacement (0, G4ThreeVector (0, 0, 3*cm), target_logical_volume,
                         "target_physical_vol", world_logical_volume, false, 0);

    G4Tubs *detector_tubs =
      new G4Tubs ("detector_tubs", 0, 3.5*cm, 4.5*mm, 0, 360*degree);
    G4LogicalVolume *detector_logical_volume =
      new G4LogicalVolume (detector_tubs, ge, "detector_logical_vol");
    G4VPhysicalVolume *detector_physical_volume =
      new G4PVPlacement (0, G4ThreeVector (0, 0, 0.5*cm), detector_logical_volume,
                         "detector_physical_vol", world_logical_volume, false, 0);





    // sensitive_detector *detector = new sensitive_detector ("detector");
    // G4SDManager::GetSDMpointer ()->AddNewDetector (detector);
    // detector_logical_volume->SetSensitiveDetector (detector);





    // G4OpticalSurface* Scint2World = new G4OpticalSurface("ScintWrap");
    // new G4LogicalBorderSurface("Scint2World", scintillator_physical_volume,world_physical_volume,Scint2World);
    // Scint2World->SetType(dielectric_metal);
    // Scint2World->SetFinish(polishedair);
    // Scint2World->SetModel(glisur);
    // const G4int num = 1;
    // G4double pp[num] = {0.329*eV};
    // G4double reflectivity[num] = {1.3};
    // G4double efficiency[num] = {1.1};
    // G4MaterialPropertiesTable* WrapProperty
    // = new G4MaterialPropertiesTable();
    // WrapProperty->AddProperty("REFLECTIVITY",Scnt_PP,reflectivity,num);
    // WrapProperty->AddProperty("EFFICIENCY",Scnt_PP,efficiency,num);
    // Scint2World->SetMaterialPropertiesTable(WrapProperty);
    // new G4LogicalSkinSurface("MirrorSurface",scintillator_logical_volume,Scint2World);

    return world_physical_volume;
  }
};


class primary_generator_action : public G4VUserPrimaryGeneratorAction
{
public:
  primary_generator_action ()
  {
    particle_gun = new G4ParticleGun (1);

    G4ParticleTable* particle_table = G4ParticleTable::GetParticleTable();
    G4ParticleDefinition* particle = particle_table->FindParticle("gamma");
  
    particle_gun->SetParticleDefinition (particle);
    particle_gun->SetParticleEnergy ( Energy * MeV);
  }

  ~primary_generator_action ()
  {
    delete particle_gun;
  }

  void GeneratePrimaries (G4Event *event)
  {

    G4ParticleDefinition* particle = particle_gun->GetParticleDefinition();
  if (particle == G4ChargedGeantino::ChargedGeantino()) {
    //fluorine 
    int Z = ChargeNumber, A = MassNumber;
    G4double ionCharge   = ChargeNumber;
    G4double excitEnergy = 0.*MeV;
    
    G4ParticleDefinition* ion
       = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    particle_gun->SetParticleDefinition(ion);
    particle_gun->SetParticleCharge(ionCharge);
  }
    double x = 0;
    double y = 0;
    double z = drand48 () * 0.1 * mm + 3 * cm;
    particle_gun->SetParticlePosition (G4ThreeVector (x, y, z));
    particle_gun->SetParticleMomentumDirection (G4RandomDirection ());
    
    particle_gun->GeneratePrimaryVertex (event);
  }

private:
  G4ParticleGun *particle_gun;
};


int main (int argc, char *argv [])
{
  G4RunManager *rm = new G4RunManager;
  setlinebuf (stdout);

  rm->SetUserInitialization (new detector_construction);
  rm->SetUserInitialization (new OpticalPhysics);
  rm->SetUserAction (new my_event_action);
  rm->SetUserAction (new my_stepping_action);
  rm->SetUserAction (new my_run_action);
  rm->SetUserAction (new primary_generator_action);
  rm->Initialize ();

  G4UImanager *ui = G4UImanager::GetUIpointer ();

  G4VisManager *vm = new G4VisExecutive;
  vm->Initialize ();

  G4UIExecutive *uie = new G4UIExecutive (argc, argv);

  ui->ApplyCommand ("/control/execute vis1.mac");
  uie->SessionStart ();

  delete rm;
  return 0;
}