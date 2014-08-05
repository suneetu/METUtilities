///////////////////////// -*- C++ -*- /////////////////////////////
// METRebuilder.cxx 
// Implementation file for class METRebuilder
// Author: S.Binet<binet@cern.ch>
/////////////////////////////////////////////////////////////////// 

// METUtilities includes
#include "METUtilities/METRebuilder.h"

// MET EDM
#include "xAODMissingET/MissingETContainer.h"
#include "xAODMissingET/MissingETComposition.h"
#include "xAODMissingET/MissingETAuxContainer.h"
#include "xAODMissingET/MissingETComponentMap.h"

// EDM includes
#include "xAODEgamma/ElectronContainer.h"
#include "xAODEgamma/PhotonContainer.h"
#include "xAODTau/TauJetContainer.h"
#include "xAODMuon/MuonContainer.h"

// Egamma EDM
#include "xAODEgamma/ElectronxAODHelpers.h"

// Tracking EDM
#include "xAODTracking/TrackParticle.h"

// Track errors
#include "EventPrimitives/EventPrimitivesHelpers.h"

namespace met {

  using std::vector;

  using xAOD::MissingET;
  using xAOD::MissingETContainer;
  using xAOD::MissingETComponent;
  using xAOD::MissingETComponentMap;
  using xAOD::MissingETAuxContainer;
  //
  using xAOD::IParticle;
  using xAOD::IParticleContainer;
  //
  using xAOD::JetContainer;
  using xAOD::JetConstituentVector;
  //
  using xAOD::TrackParticle;

  /////////////////////////////////////////////////////////////////// 
  // Public methods: 
  /////////////////////////////////////////////////////////////////// 

  // Constructors
  ////////////////
  METRebuilder::METRebuilder(const std::string& name) : 
    AsgTool(name),
    m_doEle(false),
    m_doGamma(false),
    m_doTau(false),
    m_doJet(false),
    m_doMuon(false),
    m_objLinkAcc("originalObjectLink")
  {
    //
    // Property declaration
    // 
    declareProperty( "EleColl",         m_eleColl    = "ElectronCollection" );
    declareProperty( "GammaColl",       m_gammaColl  = "PhotonCollection"   );
    declareProperty( "TauColl",         m_tauColl    = "TauRecContainer"    );
    declareProperty( "JetColl",         m_jetColl    = "AntiKt4LCTopoJets"  );
    declareProperty( "MuonColl",        m_muonColl   = "Muons"              );
    //
    declareProperty( "EleTerm",         m_eleTerm    = "RefEle"             );
    declareProperty( "GammaTerm",       m_gammaTerm  = "RefGamma"           );
    declareProperty( "TauTerm",         m_tauTerm    = "RefTau"             );
    declareProperty( "JetTerm",         m_jetTerm    = "RefJet"             );
    declareProperty( "MuonTerm",        m_muonTerm   = "Muons"              );
    declareProperty( "SoftTerm",        m_softTerm   = "PVSoftTrk"          );
    //
    declareProperty( "InputMap",        m_inputMap   = "METMap_RefFinal"    );
    declareProperty( "OutputContainer", m_outMETCont = "MET_MyRefFinal"     );
    declareProperty( "OutputTotal",     m_outMETTerm = "Final"              );

    // migrate to new tool at some point
    declareProperty( "CalibJetPtCut",   m_jetPtCut   = 20e3                 );
    // should be in a track tool -- also should use track isolation decorations
    declareProperty( "DoTrackPVSel",    m_trk_doPVsel = true                );
    declareProperty( "VertexColl",      m_vtxColl    = "PrimaryVertices"    );
    declareProperty( "ClusterColl",     m_clusColl   = "CaloCalTopoCluster" );
    declareProperty( "TrackD0Max",      m_trk_d0Max = 1.5                   );
    declareProperty( "TrackZ0Max",      m_trk_z0Max = 1.5                   );

  }

  // Destructor
  ///////////////
  METRebuilder::~METRebuilder()
  {}

  // Athena algtool's Hooks
  ////////////////////////////
  StatusCode METRebuilder::initialize()
  {
    ATH_MSG_INFO ("Initializing " << name() << "...");

    if( m_inputMap.size()==0 ) {
      ATH_MSG_FATAL("Input MissingETComponentMap name must be provided.");
      return StatusCode::FAILURE;
    }
    ATH_MSG_INFO ("Input MET Map: " << m_inputMap);

    if( m_outMETCont.size()==0 ) {
      ATH_MSG_FATAL("Output MissingETContainer name must be provided.");
      return StatusCode::FAILURE;
    }
    ATH_MSG_INFO ("Output MET Container: " << m_outMETCont);

    ATH_MSG_INFO ("Configured to rebuild following MET Terms:");
    if( m_eleColl!=""   ) {
      m_doEle = true;
      ATH_MSG_INFO("  Electrons: " << m_eleColl
		   << " > " << m_eleTerm );
    }
    if( m_gammaColl!="" ) {
      m_doGamma = true;
      ATH_MSG_INFO("  Photons:   " << m_gammaColl
		   << " > " << m_gammaTerm );
    }
    if( m_tauColl!=""   ) {
      m_doTau = true;
      ATH_MSG_INFO("  Taus:      " << m_tauColl
		   << " > " << m_tauTerm );
    }
    if( m_jetColl!=""   ) {
      m_doJet = true;
      ATH_MSG_INFO("  Jets:      " << m_jetColl
		   << " > " << m_jetTerm );
    }
    if( m_muonColl!=""  ) {
      m_doMuon = true;
      ATH_MSG_INFO("  Muons:     " << m_muonColl
		   << " > " << m_muonTerm );
    }
    ATH_MSG_INFO ("  Soft:      " << m_softTerm);

    return StatusCode::SUCCESS;
  }

  StatusCode METRebuilder::finalize()
  {
    ATH_MSG_INFO ("Finalizing " << name() << "...");

    return StatusCode::SUCCESS;
  }

  StatusCode METRebuilder::execute()
  {
    ATH_MSG_DEBUG ( name() << " in execute...");

    const MissingETComponentMap* metMap = 0;
    if( evtStore()->retrieve(metMap, m_inputMap).isFailure() ) {
      ATH_MSG_WARNING("Unable to retrieve MissingETComponentMap: " << m_inputMap);
      return StatusCode::SUCCESS;
    }

    // Create a MissingETContainer with its aux store
    MissingETContainer* outCont = new MissingETContainer();
    if( evtStore()->record(outCont, m_outMETCont).isFailure() ) {
      ATH_MSG_WARNING("Unable to record MissingETContainer: " << m_outMETCont);
      return StatusCode::SUCCESS;
    }
    MissingETAuxContainer* metAuxCont = new MissingETAuxContainer();
    if( evtStore()->record(metAuxCont, m_outMETCont+"Aux.").isFailure() ) {
      ATH_MSG_WARNING("Unable to record MissingETAuxContainer: " << m_outMETCont+"Aux.");
      return StatusCode::SUCCESS;
    }
    outCont->setStore(metAuxCont);

    MissingET* metFinal = new MissingET(0.,0.,0.);
    metFinal->setName("Final");

    const MissingET* oldSoft = xAOD::MissingETComposition::getMissingET(metMap,m_softTerm);
    ATH_MSG_DEBUG("Original MET Soft --"
		  << " mpx: " << oldSoft->mpx()
		  << " mpy: " << oldSoft->mpy());
    // need -1 because this assumes adding a particle, not a MET
    // copy constructor needs correcting.
    MissingET* metSoft = new MissingET(-1*oldSoft->mpx(),
				       -1*oldSoft->mpy(),
				       -1*oldSoft->sumet(),
				       oldSoft->name(),oldSoft->source());

    bool doTracks = MissingETBase::Source::isTrackTerm(metSoft->source());
    
    m_usedTracks.clear();
   
    const MissingETComponent* metComp = 0;
    if(m_doEle) {
      MissingET* metEle = new MissingET(0.,0.,0.);
      ATH_MSG_DEBUG("Rebuild electron MET term");
      outCont->push_back(metEle);
      //
      const xAOD::ElectronContainer* elec = 0;
      if( evtStore()->retrieve(elec, m_eleColl).isFailure() ) {
	ATH_MSG_WARNING("Unable to retrieve ElectronContainer: " << m_eleColl);
	return StatusCode::SUCCESS;
      }
      metComp = xAOD::MissingETComposition::getComponent(metMap,m_eleTerm);
      if(!metComp) {
	ATH_MSG_WARNING("Could not find current METComponent in MET Map!");
	return StatusCode::SUCCESS;
      }
      metEle->setName(m_eleTerm);
      ATH_CHECK( rebuildMET(metEle, elec, metComp, doTracks) );
      (*metFinal) += *metEle;
    }

    if(m_doGamma) {
      MissingET* metGamma = new MissingET(0.,0.,0.);
      ATH_MSG_DEBUG("Rebuild photon MET term");
      outCont->push_back(metGamma);
      //
      const xAOD::PhotonContainer* gamma = 0;
      if( evtStore()->retrieve(gamma, m_gammaColl).isFailure() ) {
	ATH_MSG_WARNING("Unable to retrieve GammaContainer: " << m_gammaColl);
	return StatusCode::SUCCESS;
      }
      metComp = xAOD::MissingETComposition::getComponent(metMap,m_gammaTerm);
      if(!metComp) {
	ATH_MSG_WARNING("Could not find current METComponent in MET Map!");
	return StatusCode::SUCCESS;
      }
      metGamma->setName(m_gammaTerm);
      ATH_CHECK( rebuildMET(metGamma, gamma, metComp, doTracks) );
      (*metFinal) += *metGamma;
    }

    if(m_doTau) {
      MissingET* metTau = new MissingET(0.,0.,0.);
      ATH_MSG_DEBUG("Rebuild tau MET term");
      outCont->push_back(metTau);
      //
      const xAOD::TauJetContainer* taujet = 0;
      if( evtStore()->retrieve(taujet, m_tauColl).isFailure() ) {
	ATH_MSG_WARNING("Unable to retrieve TauJetContainer: " << m_tauColl);
	return StatusCode::SUCCESS;
      }
      metComp = xAOD::MissingETComposition::getComponent(metMap,m_tauTerm);
      if(!metComp) {
	ATH_MSG_WARNING("Could not find current METComponent in MET Map!");
	return StatusCode::SUCCESS;
      }
      metTau->setName(m_tauTerm);
      ATH_CHECK( rebuildMET(metTau, taujet, metComp, doTracks) );
      (*metFinal) += *metTau;
    }

    if(m_doMuon) {
      // May need implementation of Eloss correction
      // Place in separate tool (?)
      MissingET* metMuon = new MissingET(0.,0.,0.);
      ATH_MSG_DEBUG("Rebuild muon MET term");
      outCont->push_back(metMuon);
      //
      const xAOD::MuonContainer* muon = 0;
      if( evtStore()->retrieve(muon, m_muonColl).isFailure() ) {
	ATH_MSG_WARNING("Unable to retrieve MuonContainer: " << m_muonColl);
	return StatusCode::SUCCESS;
      }
      metComp = xAOD::MissingETComposition::getComponent(metMap,m_muonTerm);
      if(!metComp) {
	ATH_MSG_WARNING("Could not find current METComponent in MET Map!");
	return StatusCode::SUCCESS;
      }
      metMuon->setName(m_muonTerm);
      ATH_CHECK( rebuildMET(metMuon, muon, metComp, doTracks) );
      (*metFinal) += *metMuon;
    }

    if(m_doJet) {
      // Needs implementation of the jet/soft term rebuilding too.
      // Place in separate tool (?)
      MissingET* metJet = new MissingET(0.,0.,0.);
      ATH_MSG_DEBUG("Rebuild jet and soft MET terms");
      outCont->push_back(metJet);
      //
      const xAOD::JetContainer* jet = 0;
      if( evtStore()->retrieve(jet, m_jetColl).isFailure() ) {
	ATH_MSG_WARNING("Unable to retrieve JetContainer: " << m_jetColl);
	return StatusCode::SUCCESS;
      }

      metComp = xAOD::MissingETComposition::getComponent(metMap,m_jetTerm);
      if(!metComp) {
	ATH_MSG_WARNING("Could not find current METComponent in MET Map!");
	return StatusCode::SUCCESS;
      }
      metJet->setName(m_jetTerm);
      ATH_CHECK( rebuildJetMET(metJet, metSoft, jet, metComp, doTracks) );
      (*metFinal) += *metJet;
    }

    ATH_MSG_DEBUG("Add MET soft term");
    outCont->push_back(metSoft);
    (*metFinal) += *metSoft;

    ATH_MSG_VERBOSE( "Rebuilt MET soft --"
		     << " mpx: " << metSoft->mpx()
		     << " mpy: " << metSoft->mpy()
		     );

    ATH_MSG_VERBOSE( "Rebuilt MET Final --"
		     << " mpx: " << metFinal->mpx()
		     << " mpy: " << metFinal->mpy()
		     );

    outCont->push_back(metFinal);

    return StatusCode::SUCCESS;
  }

  StatusCode METRebuilder::rebuildMET(xAOD::MissingET* met,
				      const xAOD::IParticleContainer* collection,
				      const xAOD::MissingETComponent* component,
				      bool doTracks) {

    if(component->size()==0) return StatusCode::SUCCESS;

    ATH_MSG_VERBOSE("Rebuilding MET term " << component->metObject()->name());

    const IParticleContainer* testCollection = dynamic_cast<const IParticleContainer*>(component->objects().front()->container());
    bool originalInputs = (testCollection == collection);
    bool matchCollection = true;
    if(collection->size()>0) {
      // Consistency test: check that the collection supplied is the original one
      // used for MET reconstruction, or is derived from this collection
      if(!originalInputs) {
	const IParticle* pObj = collection->front();
	if(!m_objLinkAcc.isAvailable(*pObj)) {
	  ATH_MSG_WARNING("Modified container provided without originalObjectLink -- cannot proceed.");
	  matchCollection = false;
	} else {
	  const IParticleContainer* sourceCollection = dynamic_cast<const IParticleContainer*>((*m_objLinkAcc(*pObj))->container());
	  matchCollection = (sourceCollection == testCollection);
	}
      }
      if(!matchCollection) {
	ATH_MSG_WARNING("Problem with input object container -- skipping this term.");
	return StatusCode::SUCCESS;
      }
    }

    // Method flow:
    // 1. Loop over the objects in the collection
    // 2. Find them or their originals in the METComponent
    // 3. Add to the MET term with appropriate weights
    
    for( IParticleContainer::const_iterator iObj=collection->begin();
	 iObj!=collection->end(); ++iObj ) {

      const IParticle* pObj = *iObj;
      // check if this is a copy - if so, get the original object pointer
      if(!originalInputs) pObj = *m_objLinkAcc(*pObj);

      if(component->findIndex(pObj) != MissingETBase::Numerical::invalidIndex()) {
	MissingETBase::Types::weight_t objWeight = component->weight(*iObj);

	if(doTracks) {
	  associateTracks(*iObj);
	}

	met->add((*iObj)->pt()*cos((*iObj)->phi())*objWeight.wpx(),
		 (*iObj)->pt()*sin((*iObj)->phi())*objWeight.wpy(),
		 (*iObj)->pt()*objWeight.wet());
      } // used object in MET
    }

    ATH_MSG_VERBOSE( "Original MET --"
		     << " mpx: " << component->metObject()->mpx()
		     << " mpy: " << component->metObject()->mpy()
		     );
    ATH_MSG_VERBOSE( "Rebuilt MET  --"
		     << " mpx: " << met->mpx()
		     << " mpy: " << met->mpy()
		     );

    return StatusCode::SUCCESS;
  }

  StatusCode METRebuilder::rebuildJetMET(xAOD::MissingET* metJet,
					 xAOD::MissingET* metSoft,
					 const xAOD::JetContainer* jets,
					 const xAOD::MissingETComponent* component,
					 bool doTracks) {

    if(component->size()==0) return StatusCode::SUCCESS;

    const xAOD::VertexContainer* vtx = 0;
    if(m_trk_doPVsel && doTracks) {
      if( evtStore()->retrieve( vtx, m_vtxColl).isFailure() ) {
        ATH_MSG_WARNING("Unable to retrieve input primary vertex container");
        return StatusCode::SUCCESS;
      }
      if(vtx->size()>0) {
	ATH_MSG_DEBUG("Main primary vertex has z = " << (*vtx)[0]->z());
      } else{
	ATH_MSG_WARNING("Event has no primary vertices");
	return StatusCode::SUCCESS;
      }
    }

    const IParticleContainer* testCollection = dynamic_cast<const IParticleContainer*>(component->objects().front()->container());
    const IParticleContainer* collcast = dynamic_cast<const IParticleContainer*>(jets);
    bool originalInputs = (testCollection == collcast);
    ATH_MSG_INFO("ting " << testCollection << " tong " << collcast);
    bool matchCollection = true;
    if(jets->size()>0) {
      // Consistency test: check that the collection supplied is the original one
      // used for MET reconstruction, or is derived from this collection
      if(!originalInputs) {
	const IParticle* pJet = jets->front();
	if(!m_objLinkAcc.isAvailable(*pJet)) {
	  ATH_MSG_WARNING("Modified container provided without originalObjectLink -- cannot proceed.");
	  matchCollection = false;
	} else {
	  const IParticleContainer* sourceCollection = dynamic_cast<const IParticleContainer*>((*m_objLinkAcc(*pJet))->container());
	  matchCollection = (sourceCollection == testCollection);
	}
      }
    }
    if(!matchCollection) {
      ATH_MSG_WARNING("Problem with input object container -- skipping these terms.");
      return StatusCode::SUCCESS;
    }
    // 1. Loop over the jets in the collection
    // 2. Find them or their originals in the METComponent
    // 3. Add to the MET term with appropriate weights
    for( JetContainer::const_iterator iJet=jets->begin();
	 iJet!=jets->end(); ++iJet ) {

      const xAOD::IParticle* pJet = *iJet;
      if(!originalInputs) pJet = *m_objLinkAcc(*pJet);

      if(component->findIndex(pJet) != MissingETBase::Numerical::invalidIndex()) {

	MissingETBase::Types::weight_t jetWeight = component->weight(pJet);
	if((*iJet)->pt()>m_jetPtCut) {
	  
	  ATH_MSG_VERBOSE("Retain jet with pt " << (*iJet)->pt() << " at full scale.");
	  
	  metJet->add((*iJet)->px()*jetWeight.wpx(),
		      (*iJet)->py()*jetWeight.wpy(),
		      (*iJet)->pt()*jetWeight.wet());
	} // minimum pt cut for jet calibration
	else {
	  if(doTracks) {
	    ATH_MSG_VERBOSE("Add tracks from jet with pt " << (*iJet)->pt());
	    vector<const TrackParticle*> jettracks = (*iJet)->getAssociatedObjects<TrackParticle>(xAOD::JetAttribute::GhostTrack);
	    for(vector<const TrackParticle*>::const_iterator iTrk = jettracks.begin();
		iTrk!=jettracks.end(); ++iTrk) {
	      // duplicate ST track selection -- should be in a tool
	      if(fabs((*iTrk)->pt())>500/*MeV*/ && fabs((*iTrk)->eta())<2.5) {
		uint8_t nPixHits(0), nSctHits(0);
		(*iTrk)->summaryValue(nPixHits,xAOD::numberOfPixelHits);
		(*iTrk)->summaryValue(nSctHits,xAOD::numberOfSCTHits);
		if(nPixHits>=1 && nSctHits>=6) {
		  bool badTrack = false;
		  if( (fabs((*iTrk)->eta())<1.5 && (*iTrk)->pt()>200e3) ||
		      (fabs((*iTrk)->eta())>=1.5 && (*iTrk)->pt()>120e3) ) {
		    // Get relative error on qoverp
		    float Rerr = Amg::error((*iTrk)->definingParametersCovMatrix(),4)/fabs((*iTrk)->qOverP());
		    // Simplified cut -- remove tracks that are more energetic than the jet
		    if(Rerr>0.4 || (*iTrk)->pt()>2*(*iJet)->pt()) badTrack = true;
		  } // additional cuts against high pt mismeasured tracks
		  bool uniqueTrack = true;
		  for(vector<const xAOD::TrackParticle*>::const_iterator jTrk=m_usedTracks.begin();
		      jTrk!=m_usedTracks.end(); ++jTrk) {
		    if(*iTrk==*jTrk) {
		      ATH_MSG_INFO("Track already used.");
		      uniqueTrack = false; break;
		    }
		  } // check for track usage
		  if(!badTrack && uniqueTrack &&
		     (!m_trk_doPVsel||isPVTrack(*iTrk,(*vtx)[0]))) {
		    ATH_MSG_VERBOSE("  + track with pt " << (*iTrk)->pt());
		    metSoft->add((*iTrk)->pt()*cos((*iTrk)->phi()),
				(*iTrk)->pt()*sin((*iTrk)->phi()),
				(*iTrk)->pt());
		  }
		} // PIX & SCT hits
	      } // pt, eta
	    } // track loop
	  } // track-based soft term
	  else {
	    // just add the weighted constituent-scale jet
	    xAOD::JetFourMom_t jetP = (*iJet)->jetP4(xAOD::JetConstitScaleMomentum);
	    ATH_MSG_VERBOSE("Add jet with pt " << (*iJet)->pt()
			    << " at constituent scale (pt = " << jetP.Pt() << ").");
	    metSoft->add(jetP.Px()*jetWeight.wpx(),
			 jetP.Py()*jetWeight.wpy(),
			 jetP.Pt()*jetWeight.wet());	    
	  } // cluster-based soft term
	} // jets below threshold should be added to the soft terms
      } // used jet in MET
    }

    ATH_MSG_VERBOSE( "Original jet MET --"
		     << " mpx: " << component->metObject()->mpx()
		     << " mpy: " << component->metObject()->mpy()
		     );
    ATH_MSG_VERBOSE( "Rebuilt jet MET  --"
		     << " mpx: " << metJet->mpx()
		     << " mpy: " << metJet->mpy()
		     );

    return StatusCode::SUCCESS;
  }


  /////////////////////////////////////////////////////////////////// 
  // Const methods: 
  ///////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////// 
  // Non-const methods: 
  /////////////////////////////////////////////////////////////////// 

  /////////////////////////////////////////////////////////////////// 
  // Protected methods: 
  /////////////////////////////////////////////////////////////////// 

  // Implement for now, but should move to common tools when possible
  bool METRebuilder::isPVTrack(const xAOD::TrackParticle* trk,
			       const xAOD::Vertex* pv) const
  {

    if(trk->d0()>m_trk_d0Max) return false;
    if(fabs(trk->z0() - pv->z()) > m_trk_z0Max) return false;

    return true;
  }

  void METRebuilder::associateTracks(const xAOD::IParticle* obj) {
    
    if(obj->type()==xAOD::Type::Electron) {
      const xAOD::Electron* el = dynamic_cast<const xAOD::Electron*>(obj);
      for(size_t iTrk=0; iTrk<el->nTrackParticles(); ++iTrk) {
	const TrackParticle* eltrk = xAOD::EgammaHelpers::getOriginalTrackParticleFromGSF(el->trackParticle(iTrk));
	m_usedTracks.push_back(eltrk);
      }
    }
    if(obj->type()==xAOD::Type::Photon) {
      const xAOD::Photon* ph = dynamic_cast<const xAOD::Photon*>(obj);
      for(size_t iVtx=0; iVtx<ph->nVertices(); ++iVtx) {
	const xAOD::Vertex* phvx = ph->vertex(iVtx);
	for(size_t iTrk=0; iTrk<phvx->nTrackParticles(); ++iTrk) {
	  const TrackParticle* phtrk = xAOD::EgammaHelpers::getOriginalTrackParticleFromGSF(phvx->trackParticle(iTrk));
	  m_usedTracks.push_back(phtrk);
	}
      }
    }
    if(obj->type()==xAOD::Type::Tau) {
      const xAOD::TauJet* tau = dynamic_cast<const xAOD::TauJet*>(obj);
      // now find associated tracks
      for(size_t iTrk=0; iTrk<tau->nTracks(); ++iTrk) {
	m_usedTracks.push_back(tau->track(iTrk));
      }
    }
    if(obj->type()==xAOD::Type::Muon) {
      const xAOD::Muon* mu = dynamic_cast<const xAOD::Muon*>(obj);
      if(mu->inDetTrackParticleLink().isValid())
	m_usedTracks.push_back(*mu->inDetTrackParticleLink());
    }
  }

  /////////////////////////////////////////////////////////////////// 
  // Const methods: 
  ///////////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////////////////// 
  // Non-const methods: 
  /////////////////////////////////////////////////////////////////// 

} //> end namespace met
