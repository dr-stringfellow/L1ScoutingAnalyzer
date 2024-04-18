#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/global/EDAnalyzer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Span.h"

// root include files
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TLorentzVector.h"

// L1 scouting 
#include "DataFormats/L1Scouting/interface/L1ScoutingMuon.h"
#include "DataFormats/L1Scouting/interface/L1ScoutingCalo.h"
#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "L1TriggerScouting/Utilities/interface/conversion.h"
#include "L1TriggerScouting/Utilities/interface/convertToL1TFormat.h"

// l1trigger
#include "DataFormats/L1Trigger/interface/Muon.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EGamma.h"
#include "DataFormats/L1Trigger/interface/Tau.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include <memory>
#include <utility>
#include <vector>
#include <iostream>
#include <map>
#include <string>

using namespace l1ScoutingRun3;

class DemoAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DemoAnalyzer(const edm::ParameterSet&);
  ~DemoAnalyzer() {}
  static void fillDescriptions(edm::ConfigurationDescriptions&);

private:
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void beginJob() override;
  void endJob() override;

  void processDataBx(
    unsigned orbitNum,
    unsigned bx,
    const edm::Handle<MuonOrbitCollection>& muonsCollection,
    const edm::Handle<JetOrbitCollection>& jetsCollection,
    const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
    const edm::Handle<TauOrbitCollection>& tausCollection,
    const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection
  );

  // tokens for scouting data
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Muon>> muonsTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Jet>> jetsTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::EGamma>> eGammasTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::Tau>> tausTokenData_;
  edm::EDGetTokenT<OrbitCollection<l1ScoutingRun3::BxSums>> bxSumsTokenData_;


  // l1t standard data format
  std::vector<l1t::Jet> l1jets_;
  std::vector<l1t::EGamma> l1egs_;
  std::vector<l1t::EtSum> l1sums_;
  std::vector<l1t::Tau> l1taus_;
  std::vector<l1t::Muon> l1muons_;

  // map containing TH1D histograms
  std::map<std::string, TH1D*> m_1dhist_;

  TTree* tree;

  Int_t nJet;
  Int_t orbit;
  Int_t bxid;
  vector<Float16_t> Jet_pt;
  vector<Float16_t> Jet_eta;
  vector<Float16_t> Jet_phi;
  vector<Float16_t> Jet_e;  
  Float16_t Dijet_deta;
  Float16_t Dijet_pt;
  Float16_t Dijet_eta;
  Float16_t Dijet_m;
  vector<Int_t> Jet_qual;
  
 };

DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iPSet)
  : muonsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("muonsTag"))),
    jetsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("jetsTag"))),
    eGammasTokenData_(consumes(iPSet.getParameter<edm::InputTag>("eGammasTag"))),
    tausTokenData_(consumes(iPSet.getParameter<edm::InputTag>("tausTag"))),
    bxSumsTokenData_(consumes(iPSet.getParameter<edm::InputTag>("bxSumsTag")))
  {

  // the root file service to handle the output file
  edm::Service<TFileService> fs;

  // Create the TTree
  tree = fs->make<TTree>("events" , "events");

  tree->Branch("orbit",      &orbit,         "orbit/i");
  tree->Branch("bx",         &bxid,             "bx/i");
  tree->Branch("nJet",       &nJet,           "nJet/i");
  tree->Branch("Jet_pt",     &Jet_pt                  );
  tree->Branch("Jet_eta",    &Jet_eta                 );
  tree->Branch("Jet_phi",    &Jet_phi                 );
  tree->Branch("Jet_e",      &Jet_e                   );
  tree->Branch("Jet_qual",   &Jet_qual                );
  tree->Branch("Dijet_deta", &Dijet_deta              );
  tree->Branch("Dijet_pt",   &Dijet_pt                );
  tree->Branch("Dijet_eta",  &Dijet_eta               );
  tree->Branch("Dijet_m",    &Dijet_m                 );

  // init internal containers for l1 objects
  l1muons_.reserve(8);
  l1jets_.reserve(12);
  l1egs_.reserve(12);
  l1taus_.reserve(12);
  l1sums_.reserve(12);
  
  // create histogram subdir
  TFileDirectory histoSubDir = fs->mkdir("histograms");

  // Init histograms
  m_1dhist_["MuonBxOcc"] = histoSubDir.make<TH1D>("MuonBxOcc", "BX in orbit with at least one muon", 3566, -0.5, 3565.5);
  m_1dhist_["Muon_JetVeto_BxOcc"] = histoSubDir.make<TH1D>("Muon_JetVeto_BxOcc", "BX in orbit with at least one muon and no jet in BX-1", 3566, -0.5, 3565.5);
  m_1dhist_["MuonPt"] = histoSubDir.make<TH1D>("MuonPt", "MuonPt", 512, 0, 256);
  m_1dhist_["numJetsBx_wMuon"] = histoSubDir.make<TH1D>("numJetsBx_wMuon", "Jet multiplicity in BX with >1Muon", 13, -0.5, 12.5);

}

void DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  edm::Handle<OrbitCollection<l1ScoutingRun3::Muon>> muonsCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::Jet>> jetsCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::EGamma>> eGammasCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::Tau>> tausCollection; 
  edm::Handle<OrbitCollection<l1ScoutingRun3::BxSums>> bxSumsCollection; 

  iEvent.getByToken(muonsTokenData_, muonsCollection); 
  iEvent.getByToken(jetsTokenData_, jetsCollection); 
  iEvent.getByToken(eGammasTokenData_, eGammasCollection); 
  iEvent.getByToken(tausTokenData_, tausCollection); 
  iEvent.getByToken(bxSumsTokenData_, bxSumsCollection); 

  // get orbit number orbit
  unsigned orbitNum = iEvent.id().event();

  // process all BX in orbit containing at least a Muon
  // getFilledBxs() returns the list of filled BX in the muon orbit collection
  for (const unsigned& bx : muonsCollection->getFilledBxs()) {
    processDataBx(
	orbitNum,
        bx,
        muonsCollection,
        jetsCollection,
        eGammasCollection,
        tausCollection,
        bxSumsCollection
    ); 
  }
  
 }

void DemoAnalyzer::processDataBx(
    unsigned orbitNum,
    unsigned bx,
    const edm::Handle<MuonOrbitCollection>& muonsCollection,
    const edm::Handle<JetOrbitCollection>& jetsCollection,
    const edm::Handle<EGammaOrbitCollection>& eGammasCollection,
    const edm::Handle<TauOrbitCollection>& tausCollection,
    const edm::Handle<BxSumsOrbitCollection>& bxSumsCollection
  ) {


    orbit = orbitNum;
    bxid = bx;
  
    // get iterator for the current BX
    const auto& jets = jetsCollection->bxIterator(bx);
    const auto& eGammas = eGammasCollection->bxIterator(bx);
    const auto& taus = tausCollection->bxIterator(bx);
    const auto& bxSums = bxSumsCollection->bxIterator(bx);
    const auto& muons = muonsCollection->bxIterator(bx);

    // convert scouting objects to l1t::objects for semplicity
    // Note: Scouting objects are stored in hw quantities, if only a subset
    // of the features is needed, the function in L1TriggerScouting/Utilities/interface/conversion.h
    // can be used to get physical quantities.
    // for example, the momentum of a muon can be obtained with ugmt::fPt(hwPt);
    l1jets_.clear();
    l1egs_.clear();
    l1taus_.clear();
    l1sums_.clear();
    l1muons_.clear();

    Jet_pt.clear();
    Jet_eta.clear();
    Jet_phi.clear();
    Jet_e.clear();
    Jet_qual.clear();
    
    for (const auto& muon : muons) {
      l1muons_.emplace_back(getL1TMuon(muon));
    }

    int jet_counter = 0;
    nJet = jets.size();
      
    for (const auto& jet : jets) {
      l1jets_.emplace_back(getL1TJet(jet));
      //std::cout << getL1TJet(jet) << std::endl;
      Jet_pt.push_back(getL1TJet(jet).pt());
      Jet_eta.push_back(getL1TJet(jet).eta());
      Jet_phi.push_back(getL1TJet(jet).phi());
      Jet_e.push_back(getL1TJet(jet).energy());
      Jet_qual.push_back(jet.hwQual());
      jet_counter++;
      if (jet_counter > 1) break;
    }
    
    for (const auto& egamma : eGammas) {
      l1egs_.emplace_back(getL1TEGamma(egamma));
    }
    for (const auto& tau : taus) {
      l1taus_.emplace_back(getL1TTau(tau));
    }    

    // store some of the sums
    if (bxSums.size()>0){
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTotalEt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTotalHt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kMissingEt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kMissingHt));
      l1sums_.emplace_back(getL1TEtSum(bxSums[0], l1t::EtSum::EtSumType::kTowerCount));
    }

    // fill histograms
    m_1dhist_["MuonBxOcc"]->Fill(bx);

    for (const auto& muon: l1muons_){
      m_1dhist_["MuonPt"]->Fill(muon.pt());
    }

    if (nJet>1){
      //TLorentzVector j1 = TLorentzVector();
      math::PtEtaPhiELorentzVectorF j1 = math::PtEtaPhiELorentzVectorF(Jet_pt[0],Jet_eta[0],Jet_phi[0],Jet_e[0]); 
      math::PtEtaPhiELorentzVectorF j2 = math::PtEtaPhiELorentzVectorF(Jet_pt[1],Jet_eta[1],Jet_phi[1],Jet_e[1]); 
      Dijet_pt = (j1+j2).Pt();
      Dijet_eta = (j1+j2).Eta();
      Dijet_m = (j1+j2).M();
      Dijet_deta = Jet_eta[0]-Jet_eta[1];
      std::cout << Dijet_pt << std::endl;
      /*
      TLorentzVector j1; j1.SetPtEtaPhiE(Jet_pt[0],Jet_eta[0],Jet_phi[0],Jet_e[0]);
      TLorentzVector j2; j2.SetPtEtaPhiE(Jet_pt[1],Jet_eta[1],Jet_phi[1],Jet_e[1]);
      Dijet_pt = (j1+j2).Pt();
      Dijet_eta = (j1+j2).Eta();
      Dijet_m = (j1+j2).M();
      */
      tree->Fill();
    }
    
    // collections are sorted based on the object Pt. For exampel, the leading muon Pt
    // can be obtained with l1muons_[0].pt()
    
    // number of jets in bx
    m_1dhist_["numJetsBx_wMuon"]->Fill(jets.size());

    // require no jets in the previous BX
    // check that we are not in the first BX. There could be halo muons from non colliding
    // bunches / cosmics.
    if (bx>1){
      // we can access the jets (if any) with jetsCollection->bxIterator(bx-1);
      // for this example, we will just check how many of them there are in bx-1
      if (jetsCollection->getBxSize(bx-1)>0){
        m_1dhist_["Muon_JetVeto_BxOcc"]->Fill(bx); 

        // fill more hists, do something else....
      }
    }

}

void DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
 
// ------------ method called once each job just before starting event loop  ------------
void DemoAnalyzer::beginJob() {
}

// ------------ method called once each job just after ending the event loop  ------------
void DemoAnalyzer::endJob() {
}

DEFINE_FWK_MODULE(DemoAnalyzer);
