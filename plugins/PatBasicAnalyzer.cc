#include <map>
#include <string>
#include <vector>
#include <memory>
#include <TLorentzVector.h>

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

class TTree;

class PatBasicAnalyzer : public edm::EDAnalyzer {

public:
  /// default constructor
  explicit PatBasicAnalyzer(const edm::ParameterSet&);
  /// default destructor
  ~PatBasicAnalyzer();
  
private:

  /// everything that needs to be done before the event loop
  virtual void beginJob() ;
  /// everything that needs to be done during the event loop
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  /// everything that needs to be done after the event loop
  virtual void endJob() ;
  
  // input tags  

  // ----------member data ---------------------------
  TTree *myTree;
  double xSecNorm_, weight_zero, weight_one, weight_two, weight_three;

  int event,realdata,run,lumi,bxnumber; 
  double EvtInfo_NumVtx,PU_npT,PU_npIT;
  //particles 
  
  edm::InputTag elecSrc_;

  std::vector<double> EvtWeights;

  std::vector<double> GLepDr01Pt;
  std::vector<double> GLepDr01Eta;
  std::vector<double> GLepDr01Phi;
  std::vector<double> GLepDr01E;
  std::vector<double> GLepDr01M;
  std::vector<double> GLepDr01Id;
  std::vector<double> GLepDr01Status;

  std::vector<double> GLepBarePt;
  std::vector<double> GLepBareEta;
  std::vector<double> GLepBarePhi;
  std::vector<double> GLepBareE;
  std::vector<double> GLepBareM;
  std::vector<double> GLepBareId;
  std::vector<double> GLepBareStatus;

  std::vector<double> GMETPt;
  std::vector<double> GMETEta;
  std::vector<double> GMETPhi;
  std::vector<double> GMETE;
  std::vector<double> GMETM;
  std::vector<double> GMETId;
  std::vector<double> GMETStatus;


  std::vector<double> St03Pt;
  std::vector<double> St03Eta;
  std::vector<double> St03Phi;
  std::vector<double> St03E;
  std::vector<double> St03M;
  std::vector<double> St03MotherId;
  std::vector<double> St03Id;
  std::vector<double> St03Status;
  std::vector<double> St03PhotonNumberMom;
  std::vector<double> GLepClosePhotPt;
  std::vector<double> GLepClosePhotEta;
  std::vector<double> GLepClosePhotPhi;
  std::vector<double> GLepClosePhotE;
  std::vector<double> GLepClosePhotM;
  std::vector<double> GLepClosePhotId;
  std::vector<double> GLepClosePhotMomId;
  std::vector<double> GLepClosePhotNumberMom;
  std::vector<double> GLepClosePhotStatus;

  std::vector<double> GJetAk04Pt;
  std::vector<double> GJetAk04Eta;
  std::vector<double> GJetAk04Phi;
  std::vector<double> GJetAk04E;
  std::vector<double> GJetAk04Px;
  std::vector<double> GJetAk04Py;
  std::vector<double> GJetAk04Pz;
  std::vector<double> GJetAk04ChFrac;

  std::vector<double> pdfInfo_;

};
//Electron includes

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "FWCore/Framework/interface/GenericHandle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/JetCorrector.h"






using namespace std;
using namespace reco;

int ccnevent=0;

PatBasicAnalyzer::PatBasicAnalyzer(const edm::ParameterSet& iConfig):
  elecSrc_(iConfig.getUntrackedParameter<edm::InputTag>("electronSrc"))
{
}

PatBasicAnalyzer::~PatBasicAnalyzer()
{
}

void
PatBasicAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

using namespace edm; //ADD  
//cout << "event#" << ccnevent << endl;
 ++ccnevent;

  // PAT trigger event


  edm::Handle<GenParticleCollection> genParticles_h;
  iEvent.getByLabel("genParticles", genParticles_h);
  const GenParticleCollection* genParticles  = genParticles_h.failedToGet () ? 0 : &*genParticles_h;


///Clear vector//////

  EvtWeights.clear();

  GLepDr01Pt.clear();
  GLepDr01Eta.clear();
  GLepDr01Phi.clear();
  GLepDr01E.clear();
  GLepDr01M.clear();
  GLepDr01Id.clear();
  GLepDr01Status.clear();

  GLepBarePt.clear();
  GLepBareEta.clear();
  GLepBarePhi.clear();
  GLepBareE.clear();
  GLepBareM.clear();
  GLepBareId.clear();
  GLepBareStatus.clear();


  GMETPt.clear();
  GMETEta.clear();
  GMETPhi.clear();
  GMETE.clear();
  GMETM.clear();
  GMETId.clear();
  GMETStatus.clear();


  St03Pt.clear();
  St03Eta.clear();
  St03Phi.clear();
  St03E.clear();
  St03M.clear();
  St03MotherId.clear();
  St03Id.clear();
  St03Status.clear();
  St03PhotonNumberMom.clear();

  GLepClosePhotPt.clear();
  GLepClosePhotEta.clear();
  GLepClosePhotPhi.clear();
  GLepClosePhotE.clear();
  GLepClosePhotM.clear();
  GLepClosePhotId.clear();
  GLepClosePhotMomId.clear();
  GLepClosePhotNumberMom.clear();
  GLepClosePhotStatus.clear();

  GJetAk04Pt.clear();
  GJetAk04Eta.clear();
  GJetAk04Phi.clear();
  GJetAk04E.clear();
  GJetAk04Px.clear();
  GJetAk04Py.clear();
  GJetAk04Pz.clear();
  GJetAk04ChFrac.clear();
  pdfInfo_.clear();

  edm::Handle<HepMCProduct> evt;
  iEvent.getByLabel("source","", evt);
  const HepMC::WeightContainer& w = evt->GetEvent()->weights();
//  w.print(std::cout);
  //EvtWeights = w.size() > 0 ? EvtWeights.push_back(w[0]) : EvtWeights.push_back(1);
  w.size() > 0 ? EvtWeights.push_back(w[0]) : EvtWeights.push_back(1);
  if(w.size() > 3){
    //cout<<xSecNorm_<<endl;
    xSecNorm_ += w[3];
    weight_zero = w[0];
    weight_one =  w[1];
    weight_two =  w[2];
    weight_three =  w[3];
  } 	
 
  const HepMC::GenEvent * mc = evt->GetEvent();
  if( mc == 0 )
    throw edm::Exception( edm::errors::InvalidReference )
      << "HepMC has null pointer to GenEvent" << endl;
  //const size_t size = mc->particles_size();
  

  const HepMC::PdfInfo *pdf = mc->pdf_info();

  
  if ( pdf ) {
    pdfInfo_.push_back(pdf->id1());
    pdfInfo_.push_back(pdf->id2());
    pdfInfo_.push_back(pdf->x1());
    pdfInfo_.push_back(pdf->x2());
    pdfInfo_.push_back(pdf->scalePDF());
    
    
    // std::cout << q <<"   "<< x1 <<"  " << x2 <<"   " << size <<std::endl; 
  }
  

  const std::vector<reco::GenParticle> & gen = *genParticles_h;
  for (size_t i=0; i<genParticles->size(); ++i){
    TLorentzVector genR1DressLep1(0,0,0,0);
    //      TLorentzVector genPho(0,0,0,0); 
    int st = gen[i].status();
    int id = gen[i].pdgId();
    //    if(gen[i].numberOfMothers()){
    if(true){
      //if (st!=3 && fabs(id)!=13&& fabs(id)!=11 && fabs(id)!=22 && fabs(id)!=23) continue;
      
      if (st==3 && gen[i].pt() > 0.1 && fabs(gen[i].eta())<3.0){
	TLorentzVector genLep3(0,0,0,0);
	genLep3.SetPtEtaPhiE(gen[i].pt(),gen[i].eta(),gen[i].phi(),gen[i].energy());
	St03Pt.push_back(genLep3.Pt());
	St03Eta.push_back(genLep3.Eta());
	St03Phi.push_back(genLep3.Phi());
	St03E.push_back(genLep3.Energy());
	St03M.push_back(genLep3.M());
	St03MotherId.push_back(gen[i].mother()->pdgId());
	St03PhotonNumberMom.push_back(gen[i].numberOfMothers());
	St03Id.push_back(id);
	St03Status.push_back(st);
      }

      //      if(gen[i].numberOfMothers() ==1 && gen[i].mother()->pdgId() != id) continue;
      
      //neutrinos:
      if(st==1 && (abs(id)==12 || abs(id)==14 || abs(id)==16)){
	TLorentzVector genGMET(0,0,0,0);
	genGMET.SetPtEtaPhiE(gen[i].pt(),gen[i].eta(),gen[i].phi(),gen[i].energy());
	GMETPt.push_back(genGMET.Pt());
	GMETEta.push_back(genGMET.Eta());
	GMETPhi.push_back(genGMET.Phi());
	GMETE.push_back(genGMET.Energy());
	GMETM.push_back(genGMET.M());
	GMETId.push_back(id);
	GMETStatus.push_back(st);
      }
	
      if (st==1 && (abs(id)==13||abs(id)==11) && gen[i].pt() > 0.1 && fabs(gen[i].eta())<3.0){
	TLorentzVector genLep1(0,0,0,0);
	genLep1.SetPtEtaPhiE(gen[i].pt(),gen[i].eta(),gen[i].phi(),gen[i].energy());
	TLorentzVector genR1Pho1(0,0,0,0);
	
	edm::Handle<std::vector<reco::GenParticle> > genpart2;//DONT know why we Need to handle another collection
	iEvent.getByLabel("genParticles", genpart2);
	const std::vector<reco::GenParticle> & gen2 = *genpart2;

	//LOOP over photons//
	for(unsigned int j=0; j<genpart2->size(); ++j){
	  if(gen2[j].numberOfMothers()){
	    if( gen2[j].status()!=1 || gen2[j].pdgId()!=22 || gen2[j].energy()<0.000001 /*|| fabs(MomId2)!=fabs(id)*/) continue;
	    TLorentzVector thisPho1(0,0,0,0);
	    thisPho1.SetPtEtaPhiE(gen2[j].pt(),gen2[j].eta(),gen2[j].phi(),gen2[j].energy());
	    double dR = genLep1.DeltaR(thisPho1);
	    if(dR<0.1){
	      genR1Pho1+=thisPho1;
	    }
	    
	    if(dR<0.2){
	      GLepClosePhotPt.push_back(thisPho1.Pt());
	      GLepClosePhotEta.push_back(thisPho1.Eta());
	      GLepClosePhotPhi.push_back(thisPho1.Phi());
	      GLepClosePhotE.push_back(thisPho1.Energy());
	      GLepClosePhotM.push_back(thisPho1.M());
	      GLepClosePhotId.push_back(gen2[j].pdgId());
	      GLepClosePhotMomId.push_back(fabs(gen2[j].mother()->pdgId()));
	      GLepClosePhotNumberMom.push_back(gen2[j].numberOfMothers());
	      GLepClosePhotStatus.push_back(gen2[j].status());
	    } //dR<0.2
	  }//gen2[j].numberOfMothers()
	}//next j, loop on genpart2 for dressing
	
	genR1DressLep1=genLep1+genR1Pho1;
	GLepDr01Pt.push_back(genR1DressLep1.Pt());
	GLepDr01Eta.push_back(genR1DressLep1.Eta());
	GLepDr01Phi.push_back(genR1DressLep1.Phi());
	GLepDr01E.push_back(genR1DressLep1.Energy());
	GLepDr01M.push_back(genR1DressLep1.M());
	GLepDr01Id.push_back(id);
	GLepDr01Status.push_back(st);
	
	GLepBarePt.push_back(genLep1.Pt());
	GLepBareEta.push_back(genLep1.Eta());
	GLepBarePhi.push_back(genLep1.Phi());
	GLepBareE.push_back(genLep1.Energy());
	GLepBareM.push_back(genLep1.M());
	GLepBareId.push_back(id);
	GLepBareStatus.push_back(st);
	
      }
    }
  }
  edm::Handle<reco::GenJetCollection> genjetColl;
  iEvent.getByLabel("ak5GenJets", genjetColl);
  const reco::GenJetCollection & genjet = *genjetColl;
  for(unsigned int k=0; k<genjetColl->size(); ++k){
    if(genjet[k].pt()<=10 || fabs(genjet[k].eta())>4.7)continue;
    GJetAk04Pt.push_back(genjet[k].pt());
    GJetAk04Eta.push_back(genjet[k].eta());
    GJetAk04Phi.push_back(genjet[k].phi());
    GJetAk04E.push_back(genjet[k].energy());
    GJetAk04Px.push_back(genjet[k].px());
    GJetAk04Py.push_back(genjet[k].py());
    GJetAk04Pz.push_back(genjet[k].pz());
    double isChargedJet=false;
    double chargedFraction = 0.;
    std::vector<const GenParticle*> mcparticles = genjet[k].getGenConstituents();
    for(std::vector <const GenParticle*>::const_iterator thepart =mcparticles.begin();thepart != mcparticles.end(); ++ thepart ) {
          if ( (**thepart).charge()!=0 ){
            isChargedJet=true;
            chargedFraction += (**thepart).pt();
	  }
    }
    if ( chargedFraction == 0 ) cout << " is chargeid: " << isChargedJet << "   " << chargedFraction/genjet[k].pt()<< endl;
    GJetAk04ChFrac.push_back(chargedFraction/genjet[k].pt());	   
  }

  /////////////////////////////////////////////////////////
    myTree->Fill();

}

void 
PatBasicAnalyzer::beginJob()
{
  

  // register to the TFileService
  edm::Service<TFileService> fs;
  
  TFileDirectory TestDir = fs->mkdir("test");
  myTree = new TTree("EventTree","EventTree");

  myTree->Branch("EvtWeights",&EvtWeights);
  myTree->Branch("xSecNorm_",&xSecNorm_);
  myTree->Branch("weight_zero",&weight_zero);
  myTree->Branch("weight_one",&weight_one);
  myTree->Branch("weight_two",&weight_two);
  myTree->Branch("weight_three",&weight_three);

  myTree->Branch("GLepDr01Pt",&GLepDr01Pt);
  myTree->Branch("GLepDr01Eta",&GLepDr01Eta);
  myTree->Branch("GLepDr01Phi",&GLepDr01Phi);
  myTree->Branch("GLepDr01E",&GLepDr01E);
  myTree->Branch("GLepDr01M",&GLepDr01M);
  myTree->Branch("GLepDr01Id",&GLepDr01Id);
  myTree->Branch("GLepDr01Status",&GLepDr01Status);
  myTree->Branch("GLepBarePt",&GLepBarePt);
  myTree->Branch("GLepBareEt",&GLepBareEta);
  myTree->Branch("GLepBarePhi",&GLepBarePhi);
  myTree->Branch("GLepBareE",&GLepBareE);
  myTree->Branch("GLepBareM",&GLepBareM);
  myTree->Branch("GLepBareId",&GLepBareId);
  myTree->Branch("GLepBareStatus",&GLepBareStatus);

  myTree->Branch("GMETPt",&GMETPt);
  myTree->Branch("GMETEt",&GMETEta);
  myTree->Branch("GMETPhi",&GMETPhi);
  myTree->Branch("GMETE",&GMETE);
  myTree->Branch("GMETM",&GMETM);
  myTree->Branch("GMETId",&GMETId);
  myTree->Branch("GMETStatus",&GMETStatus);

  myTree->Branch("St03Pt",&St03Pt);
  myTree->Branch("St03Eta",&St03Eta);
  myTree->Branch("St03Phi",&St03Phi);
  myTree->Branch("St03E",&St03E);
  myTree->Branch("St03M",&St03M);
  myTree->Branch("St03Id",&St03Id);
  myTree->Branch("St03Status",&St03Status);
  myTree->Branch("St03MotherId",&St03MotherId);
  myTree->Branch("St03PhotonNumberMom",&St03PhotonNumberMom);


  myTree->Branch("GLepClosePhotPt",&GLepClosePhotPt);
  myTree->Branch("GLepClosePhotEta",&GLepClosePhotEta);
  myTree->Branch("GLepClosePhotPhi",&GLepClosePhotPhi);
  myTree->Branch("GLepClosePhotE",&GLepClosePhotE);
  myTree->Branch("GLepClosePhotM",&GLepClosePhotM);
  myTree->Branch("GLepClosePhotId",&GLepClosePhotId);
  myTree->Branch("GLepClosePhotMomId",&GLepClosePhotMomId);
  myTree->Branch("GLepClosePhotNumberMom",&GLepClosePhotNumberMom);
  myTree->Branch("GLepClosePhotStatus",&GLepClosePhotStatus);


  myTree->Branch("GJetAk04Pt",&GJetAk04Pt);
  myTree->Branch("GJetAk04Eta",&GJetAk04Eta);
  myTree->Branch("GJetAk04Phi",&GJetAk04Phi);
  myTree->Branch("GJetAk04E",&GJetAk04E);
  myTree->Branch("GJetAk04Px",&GJetAk04Px);
  myTree->Branch("GJetAk04Py",&GJetAk04Py);
  myTree->Branch("GJetAk04Pz",&GJetAk04Pz);
  myTree->Branch("GJetAk04ChFrac",&GJetAk04ChFrac);
  myTree->Branch("pdfInfo_",&pdfInfo_);


}

void 
PatBasicAnalyzer::endJob() 
{

 myTree->Print();
}

#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(PatBasicAnalyzer);
