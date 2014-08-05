#ifndef METUtilities_METRebuilderEventLoop_H
#define METUtilities_METRebuilderEventLoop_H

// Infrastructure include(s):
#include "xAODRootAccess/Init.h"
#include "xAODRootAccess/TEvent.h"
#include "EventLoop/Algorithm.h"

// METRebuilder
namespace met {
	class METRebuilder;
}

class METRebuilderEventLoop : public EL::Algorithm
{
	// put your configuration variables here as public variables.
	// that way they can be set directly from CINT and python.
public:
	xAOD::TEvent *m_event; //!
	int m_eventCounter; //!

#ifndef __CINT__
	met::METRebuilder* m_metTool; //!
#endif // not __CINT__

	// variables that don't get filled at submission time should be
	// protected from being send from the submission node to the worker
	// node (done by the //!)
public:
	// Tree *myTree; //!
	// TH1 *myHist; //!

	// this is a standard constructor
	METRebuilderEventLoop();

	// these are the functions inherited from Algorithm
	virtual EL::StatusCode setupJob(EL::Job& job);
	virtual EL::StatusCode fileExecute();
	virtual EL::StatusCode histInitialize();
	virtual EL::StatusCode changeInput(bool firstFile);
	virtual EL::StatusCode initialize();
	virtual EL::StatusCode execute();
	virtual EL::StatusCode postExecute();
	virtual EL::StatusCode finalize();
	virtual EL::StatusCode histFinalize();

	// this is needed to distribute the algorithm to the workers
	ClassDef(METRebuilderEventLoop, 1);
};

#endif
