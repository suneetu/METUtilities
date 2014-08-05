// $Id: METRebuilderTester.cxx
// Based on CPToolTester.cxx 300804 2014-06-04 16:49:29Z krasznaa $

// System includes:
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <sstream>

// ROOT includes:
#include <TFile.h>

// Infrastructure includes:
#ifdef ROOTCORE
#   include "xAODRootAccess/Init.h"
#   include "xAODRootAccess/TEvent.h"
#   include "xAODRootAccess/TStore.h"
#endif // ROOTCORE

// EDM includes:
#include "xAODEventInfo/EventInfo.h"

#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETComposition.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETComponentMap.h"

#include "xAODMuon/MuonContainer.h"
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODEgamma/ElectronxAODHelpers.h"

#include "xAODJet/JetContainer.h"
#include "xAODTau/TauJetContainer.h"

#include "xAODTracking/VertexFwd.h"
#include "xAODTracking/TrackParticleFwd.h"
#include "xAODTracking/TrackParticle.h"
#include "EventPrimitives/EventPrimitivesHelpers.h"

// Local includes:
#include "METUtilities/METRebuilder.h"

using namespace std;

enum METterm_t {
	rebuild_ele,
	rebuild_gamma,
	rebuild_tau,
	rebuild_jet,
	rebuild_muon,
	n_rebuild_METterms
};

string collections[] = {
	"ElectronCollection",
	"PhotonCollection",
	"TauRecContainer",
	"AntiKt4LCTopoJets",
	"Muons",
	"_COLLECTIONS_" };

string terms[] = {
	"RefEle",
	"RefGamma",
	"RefTau",
	"RefJet",
	"Muons",
	"_TERMS_" };

const string inputMap = "METMap_RefFinal";
const string softTerm = "PVSoftTrk";

const string outputFilename = "AOD.METRebuild.root";

bool print_all_debug = false;

int main(int argc, char* argv[]) {

	// The application's name:
	const char* APP_NAME = argv[0];

	string fileName = "";
	int userEvtLimit = 0;
	string rebuildTerm = "";

	for (int i = 1; i < argc; i++) {
		if (strcmp(argv[i], "-n") == 0)
			userEvtLimit = atoi(argv[++i]);
		else if (strcmp(argv[i], "-f") == 0)
			fileName = argv[++i];
	}

	// Check if we received strings:
	string usageExample = string("  Usage: ") + APP_NAME + " -f [xAOD file name] -n [opitonal: number of events to process]";
	if (fileName.compare("") == 0) {
		cout << APP_NAME << " No file name received!" << endl;
		cout << usageExample << endl;
		return 1;
	}

	// Initialise the application:
	const bool chk_init = xAOD::Init(APP_NAME);

	if (chk_init) cout << APP_NAME << " initialized." << endl;
	else cout << APP_NAME << " failed." << endl;

	// Open the input file:
	cout << APP_NAME << " opening file: " << fileName.data() << endl;
	TFile* ifile = new TFile(fileName.data(), "READ");

	if (ifile->IsZombie()) {
		cout << APP_NAME << " Error opening file " << fileName.data() << endl;
		return 1;
	}

	// Create a TEvent object:
	xAOD::TEvent* eventStore = new xAOD::TEvent(xAOD::TEvent::kClassAccess);//!
	const bool evt_read = eventStore->readFrom(ifile);
	if (!evt_read) {
		cout << APP_NAME << " Error reading file " << fileName.data() << endl;
		return 1;
	}
	cout << APP_NAME << " Number of events in the file: " << eventStore->getEntries() << endl;

	// Create a transient object store. Needed for the tools.
	xAOD::TStore store;

	// Decide how many events to run over:
	Long64_t entries = eventStore->getEntries();
	if (userEvtLimit > 0 && userEvtLimit < entries) {
		entries = userEvtLimit;
	}

	// Create the tools to test:
	met::METRebuilder* metTool = new met::METRebuilder("RebuildMET");
	metTool->initialize();

	// Create a new container and its auxiliary store.
	xAOD::MissingETContainer* outContainer = new xAOD::MissingETContainer();
	xAOD::MissingETAuxContainer* outContainerAux = new xAOD::MissingETAuxContainer();
	outContainer->setStore(outContainerAux); //< Connect the two

	// Create an output file
	TFile* outputFile = new TFile(outputFilename.data(), "RECREATE");
	if (!eventStore->writeTo(outputFile)) cout << "ERROR";

	// Loop over the events:
	for (Long64_t entry = 0; entry < entries; entry++) {
		// Tell the object which entry to look at:
		eventStore->getEntry(entry);

		const xAOD::EventInfo* evtInfo = 0;
		if (!eventStore->retrieve(evtInfo, "EventInfo")) {// retrieve
			cout << APP_NAME << " Error retrieving EventInfo" << endl;
			return 1;
		}
		if (print_all_debug || (entry % 100) == 0) {
			cout << "===>>>  start processing event " << static_cast<int>(evtInfo->eventNumber());
			cout << " run #" << evtInfo->runNumber() << " " << entry << " of " << entries << "  <<<===" << endl;
		}

		// Copy EventInfo into the new file without modifications.
		if (!eventStore->copy("EventInfo")) cout << "ERROR";

		xAOD::MissingET* metFinal = new xAOD::MissingET(0., 0., 0.);
		metFinal->setName("Final");

		const xAOD::MissingETComponentMap* metMap = 0;
		if (eventStore->retrieve(metMap, inputMap).isFailure()) {
			cout << APP_NAME << "Unable to retrieve MissingETComponentMap: " << inputMap << endl;
			return 1;
		}

		const xAOD::MissingET* oldSoft = xAOD::MissingETComposition::getMissingET(metMap, softTerm);

		xAOD::MissingET* metSoft = new  xAOD::MissingET(-1 * oldSoft->mpx(),
			-1 * oldSoft->mpy(),
			-1 * oldSoft->sumet(),
			oldSoft->name(), oldSoft->source());

		bool doTracks = MissingETBase::Source::isTrackTerm(metSoft->source());

		const xAOD::MissingETComponent* metComp = 0;

		for (int currentTerm = 0; currentTerm < n_rebuild_METterms; currentTerm++) {
			if (currentTerm == rebuild_ele) {
				xAOD::MissingET* metEle = new xAOD::MissingET(0., 0., 0.);
				if (print_all_debug) cout << APP_NAME << " Rebuild electron MET term" << endl;
				outContainer->push_back(metEle);

				const xAOD::ElectronContainer* elec = 0;
				if (eventStore->retrieve(elec, collections[rebuild_ele]).isFailure()) {
					cout << APP_NAME << " Unable to retrieve ElectronContainer : " << collections[rebuild_ele] << endl;
					return 1;
				}

				metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_ele].data());
				if (!metComp) {
					cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
					return 1;
				}
				metEle->setName(terms[rebuild_ele]);
				if (!metTool->rebuildMET(metEle, elec, metComp, doTracks)) {
					cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_ele] << endl;
					return 1;
				}
				(*metFinal) += *metEle;
			}
			else if (currentTerm == rebuild_gamma) {
				xAOD::MissingET* metGamma = new xAOD::MissingET(0., 0., 0.);
				if (print_all_debug) cout << APP_NAME << " Rebuild photon MET term" << endl;
				outContainer->push_back(metGamma);
				//
				const xAOD::PhotonContainer* gamma = 0;
				if (eventStore->retrieve(gamma, collections[rebuild_gamma]).isFailure()) {
					cout << APP_NAME << " Unable to retrieve GammaContainer: " << collections[rebuild_gamma] << endl;
					return 1;
				}

				metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_gamma]);
				if (!metComp) {
					cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
					return 1;
				}
				metGamma->setName(terms[rebuild_gamma]);
				if (!metTool->rebuildMET(metGamma, gamma, metComp, doTracks)) {
					cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_gamma] << endl;
					return 1;
				}
				(*metFinal) += *metGamma;
			}
			else if (currentTerm == rebuild_tau) {
				xAOD::MissingET* metTau = new xAOD::MissingET(0., 0., 0.);
				if (print_all_debug) cout << APP_NAME << " Rebuild tau MET term" << endl;
				outContainer->push_back(metTau);
				//
				const xAOD::TauJetContainer* taujet = 0;
				if (eventStore->retrieve(taujet, collections[rebuild_tau]).isFailure()) {
					cout << APP_NAME << " Unable to retrieve TauJetContainer: " << collections[rebuild_tau] << endl;
					return 1;
				}

				metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_tau]);
				if (!metComp) {
					cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
					return 1;
				}
				metTau->setName(terms[rebuild_tau]);
				if (!metTool->rebuildMET(metTau, taujet, metComp, doTracks)) {
					cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_tau] << endl;
					return 1;
				}
				(*metFinal) += *metTau;
			}
			else if (currentTerm == rebuild_muon) {
				// May need implementation of Eloss correction
				// Place in separate tool (?)
				xAOD::MissingET* metMuon = new xAOD::MissingET(0., 0., 0.);
				if (print_all_debug) cout << APP_NAME << " Rebuild muon MET term" << endl;
				outContainer->push_back(metMuon);
				//
				const xAOD::MuonContainer* muon = 0;
				if (eventStore->retrieve(muon, collections[rebuild_muon]).isFailure()) {
					cout << APP_NAME << " Unable to retrieve MuonContainer: " << collections[rebuild_muon] << endl;
					return 1;
				}

				metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_muon]);
				if (!metComp) {
					cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
					return 1;
				}
				metMuon->setName(terms[rebuild_muon]);
				if (!metTool->rebuildMET(metMuon, muon, metComp, doTracks)) {
					cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_muon] << endl;
					return 1;
				}
				(*metFinal) += *metMuon;
			}
			else if (currentTerm == rebuild_jet) {
				// Needs implementation of the jet/soft term rebuilding too.
				// Place in separate tool (?)
				xAOD::MissingET* metJet = new xAOD::MissingET(0., 0., 0.);
				if (print_all_debug) cout << APP_NAME << " Rebuild jet and soft MET terms" << endl;
				outContainer->push_back(metJet);
				//
				const xAOD::JetContainer* jet = 0;
				if (eventStore->retrieve(jet, collections[rebuild_jet]).isFailure()) {
					cout << APP_NAME << " Unable to retrieve JetContainer: " << collections[rebuild_jet] << endl;
					return 1;
				}

				metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_jet]);
				if (!metComp) {
					cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
					return 1;
				}
				metJet->setName(terms[rebuild_jet]);
				if (!metTool->rebuildJetMET(metJet, metSoft, jet, metComp, doTracks)) {
					cout << APP_NAME << " rebuildJetMET() failed on term " << terms[rebuild_jet] << endl;
					return 1;
				}
				(*metFinal) += *metJet;
			}
			else {
				cout << APP_NAME << " Error finding MET term selection: (index)" << currentTerm << endl;
			}
		}

		outContainer->push_back(metFinal);
		if (!eventStore->record(outContainer, "MET_MyRefFinal")) return 1;// Don't delete containers that have been
		if (!eventStore->record(outContainerAux, "MET_MyRefFinalAux")) return 1;// recorded + owned by the TEvent.

		// Save the event:
		eventStore->fill();
	}

	// finalize and close our output xAOD file:
	if (!eventStore->finishWritingTo(outputFile)) cout << "ERROR writing output file.";
	outputFile->Save();

	cout << "Wrote " << outputFilename << endl;

	// Finalizers
	delete outputFile;
	delete metTool;
	delete eventStore;
	delete ifile;

	return 0;
}