#include <EventLoop/Job.h>
#include <EventLoop/StatusCode.h>
#include <EventLoop/Worker.h>
#include <METUtilities/METRebuilderEventLoop.h>

// System includes:
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <string>
#include <sstream>

// Infrastructure include(s):
#include "xAODRootAccess/TStore.h"
#include "EventLoop/OutputStream.h"
#include "EventLoop/StatusCode.h"

// ROOT includes:
#include "TFile.h"

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

// this is needed to distribute the algorithm to the workers
ClassImp(METRebuilderEventLoop)

METRebuilderEventLoop::METRebuilderEventLoop()
{
	// Here you put any code for the base initialization of variables,
	// e.g. initialize all pointers to 0.  Note that you should only put
	// the most basic initialization here, since this method will be
	// called on both the submission and the worker node.  Most of your
	// initialization code will go into histInitialize() and
	// initialize().
}



EL::StatusCode METRebuilderEventLoop::setupJob(EL::Job& job)
{
	// Here you put code that sets up the job on the submission object
	// so that it is ready to work with your algorithm, e.g. you can
	// request the D3PDReader service or add output files.  Any code you
	// put here could instead also go into the submission script.  The
	// sole advantage of putting it here is that it gets automatically
	// activated/deactivated when you add/remove the algorithm from your
	// job, which may or may not be of value to you.
	job.useXAOD();

	// let's initialize the algorithm to use the xAODRootAccess package
	xAOD::Init("METRebuilderEventLoop").ignore(); // call before opening first file

	m_metTool = new met::METRebuilder("RebuildMET");
	m_metTool->initialize();

	// tell EventLoop about our output xAOD:
	EL::OutputStream out("outputMyMET");
	job.outputAdd(out);

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::histInitialize()
{
	// Here you do everything that needs to be done at the very
	// beginning on each worker node, e.g. create histograms and output
	// trees.  This method gets called before any input files are
	// connected.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::fileExecute()
{
	// Here you do everything that needs to be done exactly once for every
	// single file, e.g. collect a list of all lumi-blocks processed
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::changeInput(bool firstFile)
{
	// Here you do everything you need to do when we change input files,
	// e.g. resetting branch addresses on trees.  If you are using
	// D3PDReader or a similar service this method is not needed.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::initialize()
{
	// Here you do everything that you need to do after the first input
	// file has been connected and before the first event is processed,
	// e.g. create additional histograms based on which variables are
	// available in the input files.  You can also create all of your
	// histograms and trees in here, but be aware that this method
	// doesn't get called if no events are processed.  So any objects
	// you create here won't be available in the output if you have no
	// input events.
	m_eventCounter = 0;

	m_event = wk()->xaodEvent();

	// as a check, let's see the number of events in our xAOD
	Info("initialize()", "Number of events = %lli", m_event->getEntries()); // print long long int

	// output xAOD
	TFile* file = wk()->getOutputFile("outputMyMET");
	m_event->writeTo(file);

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::execute()
{
	// Here you do everything that needs to be done on every single
	// events, e.g. read input variables, apply cuts, and fill
	// histograms and trees.  This is where most of your actual analysis
	// code will go.

	string APP_NAME = "METRebuilderEventLoop";

	// copy full container(s) to new xAOD
	// without modifying the contents of it: 
	if (!m_event->copy("EventInfo")) cout << "ERROR";

	// Create a new container and its auxiliary store.
	xAOD::MissingETContainer* outContainer = new xAOD::MissingETContainer();
	xAOD::MissingETAuxContainer* outContainerAux = new xAOD::MissingETAuxContainer();
	outContainer->setStore(outContainerAux); //< Connect the two

	xAOD::MissingET* metFinal = new xAOD::MissingET(0., 0., 0.);
	metFinal->setName("Final");

	const xAOD::MissingETComponentMap* metMap = 0;
	if (m_event->retrieve(metMap, inputMap).isFailure()) {
		cout << APP_NAME << "Unable to retrieve MissingETComponentMap: " << inputMap << endl;
		return EL::StatusCode::FAILURE;
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
			if (m_event->retrieve(elec, collections[rebuild_ele]).isFailure()) {
				cout << APP_NAME << " Unable to retrieve ElectronContainer : " << collections[rebuild_ele] << endl;
				return EL::StatusCode::FAILURE;
			}

			metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_ele].data());
			if (!metComp) {
				cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
				return EL::StatusCode::FAILURE;
			}
			metEle->setName(terms[rebuild_ele]);
			//rebuild// if (!m_metTool->rebuildMET(metEle, elec, metComp, doTracks)) {
			//rebuild// 	cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_ele] << endl;
			//rebuild// 	return EL::StatusCode::FAILURE;
			//rebuild// }
			(*metFinal) += *metEle;
		}
		else if (currentTerm == rebuild_gamma) {
			xAOD::MissingET* metGamma = new xAOD::MissingET(0., 0., 0.);
			if (print_all_debug) cout << APP_NAME << " Rebuild photon MET term" << endl;
			outContainer->push_back(metGamma);
			//
			const xAOD::PhotonContainer* gamma = 0;
			//rebuild// if (m_event->retrieve(gamma, collections[rebuild_gamma]).isFailure()) {
			//rebuild// 	cout << APP_NAME << " Unable to retrieve GammaContainer: " << collections[rebuild_gamma] << endl;
			//rebuild// 	return EL::StatusCode::FAILURE;
			//rebuild// }

			metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_gamma]);
			if (!metComp) {
				cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
				return EL::StatusCode::FAILURE;
			}
			metGamma->setName(terms[rebuild_gamma]);
			//rebuild// if (!m_metTool->rebuildMET(metGamma, gamma, metComp, doTracks)) {
			//rebuild// 	cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_gamma] << endl;
			//rebuild// 	return EL::StatusCode::FAILURE;
			//rebuild// }
			(*metFinal) += *metGamma;
		}
		else if (currentTerm == rebuild_tau) {
			xAOD::MissingET* metTau = new xAOD::MissingET(0., 0., 0.);
			if (print_all_debug) cout << APP_NAME << " Rebuild tau MET term" << endl;
			outContainer->push_back(metTau);
			//
			const xAOD::TauJetContainer* taujet = 0;
			if (m_event->retrieve(taujet, collections[rebuild_tau]).isFailure()) {
				cout << APP_NAME << " Unable to retrieve TauJetContainer: " << collections[rebuild_tau] << endl;
				return EL::StatusCode::FAILURE;
			}

			metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_tau]);
			if (!metComp) {
				cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
				return EL::StatusCode::FAILURE;
			}
			metTau->setName(terms[rebuild_tau]);
			//rebuild// if (!m_metTool->rebuildMET(metTau, taujet, metComp, doTracks)) {
			//rebuild// 	cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_tau] << endl;
			//rebuild// 	return EL::StatusCode::FAILURE;
			//rebuild// }
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
			if (m_event->retrieve(muon, collections[rebuild_muon]).isFailure()) {
				cout << APP_NAME << " Unable to retrieve MuonContainer: " << collections[rebuild_muon] << endl;
				return EL::StatusCode::FAILURE;
			}

			metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_muon]);
			if (!metComp) {
				cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
				return EL::StatusCode::FAILURE;
			}
			metMuon->setName(terms[rebuild_muon]);
			//rebuild// if (!m_metTool->rebuildMET(metMuon, muon, metComp, doTracks)) {
			//rebuild// 	cout << APP_NAME << " rebuildMET() failed on term " << terms[rebuild_muon] << endl;
			//rebuild// 	return EL::StatusCode::FAILURE;
			//rebuild// }
			(*metFinal) += *metMuon;
		}
		////else if (currentTerm == rebuild_jet) {
		////	// Needs implementation of the jet/soft term rebuilding too.
		////	// Place in separate tool (?)
		////	xAOD::MissingET* metJet = new xAOD::MissingET(0., 0., 0.);
		////	if (print_all_debug) cout << APP_NAME << " Rebuild jet and soft MET terms" << endl;
		////	outContainer->push_back(metJet);
		////	//
		////	const xAOD::JetContainer* jet = 0;
		////	if (m_event->retrieve(jet, collections[rebuild_jet]).isFailure()) {
		////		cout << APP_NAME << " Unable to retrieve JetContainer: " << collections[rebuild_jet] << endl;
		////		return EL::StatusCode::FAILURE;
		////	}
		////
		////	metComp = xAOD::MissingETComposition::getComponent(metMap, terms[rebuild_jet]);
		////	if (!metComp) {
		////		cout << APP_NAME << " Could not find current METComponent in MET Map!" << endl;
		////		return EL::StatusCode::FAILURE;
		////	}
		////	metJet->setName(terms[rebuild_jet]);
		////	//rebuild// if (!m_metTool->rebuildJetMET(metJet, metSoft, jet, metComp, doTracks)) {
		////	//rebuild// 	cout << APP_NAME << " rebuildJetMET() failed on term " << terms[rebuild_jet] << endl;
		////	//rebuild// 	return EL::StatusCode::FAILURE;
		////	//rebuild// }
		////	(*metFinal) += *metJet;
		////}
		else {
			if (print_all_debug) cout << APP_NAME << " Error finding MET term selection (index): " << currentTerm << endl;
		}
	}

	outContainer->push_back(metFinal);
	if (!m_event->record(outContainer, "MET_MyRefFinal")) return EL::StatusCode::FAILURE;;// Don't delete containers that have been
	if (!m_event->record(outContainerAux, "MET_MyRefFinalAux")) return EL::StatusCode::FAILURE;;// recorded + owned by the TEvent.

	// Save the event:
	m_event->fill();

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::postExecute()
{
	// Here you do everything that needs to be done after the main event
	// processing.  This is typically very rare, particularly in user
	// code.  It is mainly used in implementing the NTupleSvc.
	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::finalize()
{
	// This method is the mirror image of initialize(), meaning it gets
	// called after the last event has been processed on the worker node
	// and allows you to finish up any objects you created in
	// initialize() before they are written to disk.  This is actually
	// fairly rare, since this happens separately for each worker node.
	// Most of the time you want to do your post-processing on the
	// submission node after all your histogram outputs have been
	// merged.  This is different from histFinalize() in that it only
	// gets called on worker nodes that processed input events.

	// finalize and close our output xAOD file:
	TFile *file = wk()->getOutputFile("outputMyMET");
	if (!m_event->finishWritingTo(file)) cout << "ERROR: Failed to write output file!" << endl;

	if (m_metTool) {
		delete m_metTool;
		m_metTool = 0;
	}

	return EL::StatusCode::SUCCESS;
}



EL::StatusCode METRebuilderEventLoop::histFinalize()
{
	// This method is the mirror image of histInitialize(), meaning it
	// gets called after the last event has been processed on the worker
	// node and allows you to finish up any objects you created in
	// histInitialize() before they are written to disk.  This is
	// actually fairly rare, since this happens separately for each
	// worker node.  Most of the time you want to do your
	// post-processing on the submission node after all your histogram
	// outputs have been merged.  This is different from finalize() in
	// that it gets called on all worker nodes regardless of whether
	// they processed input events.
	return EL::StatusCode::SUCCESS;
}
