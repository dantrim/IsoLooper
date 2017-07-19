#include "IsoLooper/IsoLooper.h"

// SusyNtuple
#include "SusyNtuple/KinematicTools.h"
#include "SusyNtuple/SusyDefs.h"
using namespace Susy; // everything in SusyNtuple is in this namespace

//ROOT

// std/stl
#include <iomanip> // setw
#include <iostream>
#include <string>
#include <sstream> // stringstream, ostringstream
using namespace std;

//////////////////////////////////////////////////////////////////////////////
IsoLooper::IsoLooper() :
    m_dbg(0),
    m_input_chain(nullptr),
    m_mc_weight(1.0),
    m_lep_sf(1.0),
    m_weight(1.0),
    m_res_mass(0),
    m_lep_type(DiLepEvtType::ET_Unknown)
{
}
//////////////////////////////////////////////////////////////////////////////
void IsoLooper::Begin(TTree* /*tree*/)
{
    // call base class' Begin method
    SusyNtAna::Begin(0);
    if(dbg()) cout << "IsoLooper::Begin" << endl;

    m_outfile = new TFile("iso_tester.root", "RECREATE");
    initialize_histos();

    return;
}
//////////////////////////////////////////////////////////////////////////////
void IsoLooper::initialize_histos()
{
    // electron ptvarcone20 for EE channel
    h_ele_ptvarcone20_ee = new TH1F("h_ele_ptvarcone20_ee", "Electron ptvarcone20 (EE);ptvarcone20;Entries",20,0,5);
    h_ele_ptvarcone20_ee->SetLineColor(kBlue);
    histos.push_back(h_ele_ptvarcone20_ee);

    // electron ptvarcone20 for EM channel, leading E
    h_ele_ptvarcone20_em = new TH1F("h_ele_ptvarcone20_em", "Electron ptvarcone20 (EM);ptvarcone20;Entries",100,0,10);
    h_ele_ptvarcone20_em->SetLineColor(kGreen);
    histos.push_back(h_ele_ptvarcone20_em);

    // electron ptvarcone20 for EM channel, leading M
    h_ele_ptvarcone20_me = new TH1F("h_ele_ptvarcone20_me", "Electron ptvarcone20 (ME);ptvarcone20;Entries",100,0,10);
    h_ele_ptvarcone20_me->SetLineColor(kMagenta);
    histos.push_back(h_ele_ptvarcone20_me);

    // muon ptvarcone20 for MM channel
    h_muo_ptvarcone20_mm = new TH1F("h_muo_ptvarcone20_mm", "Muon ptvarcone20 (MM);ptvarcone20;Entries",20, 0, 5);
    h_muo_ptvarcone20_mm->SetLineColor(kRed);
    histos.push_back(h_muo_ptvarcone20_mm);

    // muon ptvarcone20 for EM channel, leading E
    h_muo_ptvarcone20_em = new TH1F("h_muo_ptvarcone20_em", "Muon ptvarcone20 (EM);ptvarcone20;Entries",100,0,10);
    h_muo_ptvarcone20_em->SetLineColor(kRed);
    histos.push_back(h_muo_ptvarcone20_em);

    // muon ptvarcone20 for EM channel, leading M
    h_muo_ptvarcone20_me = new TH1F("h_muo_ptvarcone20_me", "Muon ptvarcone20 (ME);ptvarcone20;Entries",100,0,10);
    h_muo_ptvarcone20_me->SetLineColor(kMagenta);
    histos.push_back(h_muo_ptvarcone20_me);


    // RELATIVE
    h_ele_ptvarcone20_pt_ee = new TH1F("h_ele_ptvarcone20_pt_ee", "Electron ptvarcone20/pt (EE);ptvarcone20/pt;Entries",40, 0, 0.125);
    h_ele_ptvarcone20_pt_ee->SetLineColor(kBlue);
    histos.push_back(h_ele_ptvarcone20_pt_ee);

    h_muo_ptvarcone20_pt_mm = new TH1F("h_muo_ptvarcone20_pt_mm", "Muon ptvarcone20/pt (MM);ptvarcone20/pt;Entries", 40, 0, 0.125);
    h_muo_ptvarcone20_pt_mm->SetLineColor(kRed);
    histos.push_back(h_muo_ptvarcone20_pt_mm);

    // dphi
    h_dphi = new TH1F("h_dphi", "#Delta#phi_{ll};#Delta#phi_{ll};Entries",32, 0, 3.2);
    h_dphi->SetLineColor(kBlack);
    histos.push_back(h_dphi);

    // 2D
    h2_dphi_iso_ee = new TH2F("h2_dphi_iso_ee", "#Delta#phi_{ll} vs. Iso (EE);ptvarcone20/pt;#Delta#phi",100,0,0.1, 32,0,3.2);

    h2_dphi_iso_mm = new TH2F("h2_dphi_iso_mm", "#Delta#phi_{ll} vs. Iso (MM);ptvarcone20/pt;#Delta#phi",100,0,0.1, 32,0,3.2);

    // delta R
    h_dr_ll = new TH1F("h_dr_ll", "#DeltaR_{ll};#DeltaR_{ll};Entries", 100, 0, 5);

    // detla R for ee channel
    h_dr_ll_ee = new TH1F("h_dr_ll_ee", "#DeltaR_{ll} (EE);#DeltaR_{ll};Entries", 100, 0, 5);

    h_dr_ll_mm = new TH1F("h_dr_ll_mm", "#DeltaR_{ll} (MM);#DeltaR_{ll};Entries", 100, 0, 5);

    h_dr_ll_em = new TH1F("h_dr_ll_em", "#DeltaR_{ll} (EM);#DeltaR_{ll};Entries", 100, 0, 5);

    // bjet dr
    h_dr_l0_b0 = new TH1F("h_dr_l0_b0", "#DeltaR_{l0,b0};#DeltaR_{l0,b0};Entries", 100, 0, 5);
    h_dr_l1_b0 = new TH1F("h_dr_l1_b0", "#DeltaR_{l1,b0};#DeltaR_{l1,b0};Entries", 100, 0, 5);
    h_dr_e0_b0 = new TH1F("h_dr_e0_b0", "#DeltaR_{e0,b0};#DeltaR_{e0,b0};Entries", 100, 0, 5);
    h_dr_e1_b0 = new TH1F("h_dr_e1_b0", "#DeltaR_{e1,b0};#DeltaR_{e1,b0};Entries", 100, 0, 5);
    h_dr_m0_b0 = new TH1F("h_dr_m0_b0", "#DeltaR_{m0,b0};#DeltaR_{m0,b0};Entries", 100, 0, 5);
    h_dr_m1_b0 = new TH1F("h_dr_m1_b0", "#DeltaR_{m1,b0};#DeltaR_{m1,b0};Entries", 100, 0, 5);

    h_dr_l0_b1 = new TH1F("h_dr_l0_b1", "#DeltaR_{l0,b1};#DeltaR_{l0,b1};Entries", 100, 0, 5);
    h_dr_l1_b1 = new TH1F("h_dr_l1_b1", "#DeltaR_{l1,b1};#DeltaR_{l1,b1};Entries", 100, 0, 5);
    h_dr_e0_b1 = new TH1F("h_dr_e0_b1", "#DeltaR_{e0,b1};#DeltaR_{e0,b1};Entries", 100, 0, 5);
    h_dr_e1_b1 = new TH1F("h_dr_e1_b1", "#DeltaR_{e1,b1};#DeltaR_{e1,b1};Entries", 100, 0, 5);
    h_dr_m0_b1 = new TH1F("h_dr_m0_b1", "#DeltaR_{m0,b1};#DeltaR_{m0,b1};Entries", 100, 0, 5);
    h_dr_m1_b1 = new TH1F("h_dr_m1_b1", "#DeltaR_{m1,b1};#DeltaR_{m1,b1};Entries", 100, 0, 5);
}
//////////////////////////////////////////////////////////////////////////////
Bool_t IsoLooper::Process(Long64_t entry)
{

    // calling "GetEntry" loads into memory the susyNt class objects for this event
    GetEntry(entry);
    SusyNtAna::clearObjects(); // clear the previous event's objects

    // increment the chain entry (c.f. SusyNtuple/SusyNtAna.h)
    m_chainEntry++;

    // evt() provides pointer to the SusyNt::Event class object for this event
    int run_number = nt.evt()->run;
    int event_number = nt.evt()->eventNumber;

    if(dbg() || m_chainEntry%1000==0) {
        cout << "IsoLooper::Process    **** Processing entry " << setw(6) << m_chainEntry
                << "  run " << run_number << "  event " << event_number << " **** " << endl;
    }

    // SusyNtAna::selectObject fills the baseline and signal objects
    // for the given AnalysisType
    // m_preX    = objects before any selection (as they are in susyNt)
    // m_baseX   = objects with the Analysis' baseline selection AND overlap removal applied
    // m_signalX = objects with the Analysis' signal selection applied (and baseline AND overlap removal)
    SusyNtAna::selectObjects();

    // get the MC weight using the inherited MCWeighter object
    // (c.f. SusyNtuple/MCWeighter.h)
    if(nt.evt()->isMC) {
        float lumi = 100000; // normalize the MC to 100 fb-1
        m_mc_weight = SusyNtAna::mcWeighter().getMCWeight(nt.evt(), lumi, NtSys::NOM);
    }
    else {
        m_mc_weight = 1.; // don't re-weight data
    }

    // check that the event passes the standard ATLAS event cleaning cuts
    if(!passEventCleaning(m_preMuons, m_baseMuons, m_baseJets)) return false;

    LeptonVector leptons;
    if(!passSignalSelection(m_baseLeptons, m_signalJets, leptons)) return false;
    std::sort(leptons.begin(), leptons.end(), comparePt);
    m_lep_type = getDiLepEvtType(leptons);

    m_lep_sf = nttools().leptonEffSF(leptons);

    JetVector bjets;
    get_bjets(m_signalJets, bjets); 
    std::sort(bjets.begin(), bjets.end(), comparePt);

    fill_histos(leptons, bjets);

    return kTRUE;
}
//////////////////////////////////////////////////////////////////////////////
bool IsoLooper::passEventCleaning(const MuonVector& preMuons, const MuonVector& baseMuons,
            const JetVector& baseJets)
{
    int flags = nt.evt()->cutFlags[NtSys::NOM];

    if(!nttools().passGRL(flags))           return false;

    if(!nttools().passLarErr(flags))        return false;

    if(!nttools().passTileErr(flags))       return false;

    if(!nttools().passTTC(flags))           return false;

    if(!nttools().passSCTErr(flags))        return false;

    if(!nttools().passGoodVtx(flags))       return false;


    ///////////////////////////////////////////////////////
    // for bad muon, cosmic moun, and jet cleaning the
    // cuts depend on the baseline object defintion
    // (and in thec ase of the cosmic muon cut, it also
    // depends on the analysis' overlap removal
    // procedure) -- so we do not use the cutFlags but
    // rather use the objects that have passed the various
    // analysis selections to do the checks
    ///////////////////////////////////////////////////////
    if(!nttools().passBadMuon(preMuons))    return false;

    if(!nttools().passCosmicMuon(baseMuons)) return false;

    if(!nttools().passJetCleaning(baseJets)) return false;

    return true;
}
//////////////////////////////////////////////////////////////////////////////
bool IsoLooper::passSignalSelection(const LeptonVector& leptons, const JetVector& jets, LeptonVector& outleptons)
{
    int n_electrons = 0;
    int n_muons = 0;

    for(auto & lep : leptons) {
        if(lep->isEle()) {
            Susy::Electron* el = dynamic_cast<Susy::Electron*>(lep);
            if(is_signal_ele(el)) { n_electrons++; outleptons.push_back(el); }
        }
        else {
            Susy::Muon* mu = dynamic_cast<Susy::Muon*>(lep);
            if(is_signal_mu(mu)) { n_muons++; outleptons.push_back(mu); }
        }
    }

    if( (n_electrons + n_muons)==2) return true;
    else { return false; }
}
//////////////////////////////////////////////////////////////////////////////
bool IsoLooper::is_signal_ele(Electron* el)
{
    bool pass_id = (el->mediumLLH);
    bool passIP = ( fabs(el->d0sigBSCorr) < 5.0 && fabs(el->z0SinTheta()) < 0.5 ); 
    return (pass_id && passIP);
}
//////////////////////////////////////////////////////////////////////////////
bool IsoLooper::is_signal_mu(Muon* mu)
{
    bool pass_id = (mu->medium);
    bool passIP = ( fabs(mu->d0sigBSCorr) < 3.0 && fabs(mu->z0SinTheta()) < 0.5 );
    return (pass_id && passIP);
}

//////////////////////////////////////////////////////////////////////////////
void IsoLooper::get_bjets(const JetVector& injets, JetVector& bjets)
{
    for(auto & j : injets) {
        if(is_bjet(j)) bjets.push_back(j);
    }
}
bool IsoLooper::is_bjet(const Jet* j)
{
    bool pass_pt = (j->Pt() > 20.0);
    bool pass_eta = (fabs(j->Eta()) < 2.5);
    bool pass_mv2 = (j->mv2c10 > nttools().jetSelector().mv2c10_85efficiency());
    bool pass_jvt = (nttools().jetSelector().passJvt(j));

    return (pass_pt && pass_eta && pass_mv2 && pass_jvt);
}
//////////////////////////////////////////////////////////////////////////////
void IsoLooper::fill_histos(const LeptonVector& leptons, const JetVector& bjets)
{
    float iso0 = leptons.at(0)->ptvarcone20;
    float iso1 = leptons.at(1)->ptvarcone20;

    float pt0 = leptons.at(0)->Pt();
    float pt1 = leptons.at(1)->Pt();

    float dphi = fabs(leptons.at(0)->DeltaPhi(*leptons.at(1)));

    float dr = leptons.at(0)->DeltaR(*leptons.at(1));

    h_dr_ll->Fill(dr, w());
    if(ee()) {
        h_dr_ll_ee->Fill(dr, w());
    }
    else if(mm()) {
        h_dr_ll_mm->Fill(dr, w());
    }
    else if(em() || me()) {
        h_dr_ll_em->Fill(dr, w());
    }


    if(leptons.at(0)->isEle()) {
        float iso = iso0/pt0;
        if(iso!=0)
        h2_dphi_iso_ee->Fill(iso, dphi);
    }
    else {
        float iso = iso0/pt0;
        if(iso!=0)
        h2_dphi_iso_mm->Fill(iso, dphi);
    }


    if(ee()) {
        if(iso0!=0) {
            h_ele_ptvarcone20_ee->Fill(iso0);
            h_ele_ptvarcone20_pt_ee->Fill(iso0/pt0);
        }
        if(iso1!=0) {
            h_ele_ptvarcone20_ee->Fill(iso1);
            h_ele_ptvarcone20_pt_ee->Fill(iso1/pt1);
        }
    }
    else if(em()) {
        if(iso0!=0) {
            h_ele_ptvarcone20_ee->Fill(iso0);
            h_ele_ptvarcone20_pt_ee->Fill(iso0/pt0);
            //h_ele_ptvarcone20_em->Fill(iso0);
        }
        if(iso1!=0) {
            //h_muo_ptvarcone20_em->Fill(iso1);
            h_muo_ptvarcone20_mm->Fill(iso1);
            h_muo_ptvarcone20_pt_mm->Fill(iso1/pt1);
        }
    }
    else if(me()) {
        if(iso0!=0) {
            h_muo_ptvarcone20_mm->Fill(iso0);
            h_muo_ptvarcone20_pt_mm->Fill(iso0/pt0);
            //h_ele_ptvarcone20_me->Fill(iso0);
        }
        if(iso1!=0) {
            h_ele_ptvarcone20_ee->Fill(iso1);
            h_ele_ptvarcone20_pt_ee->Fill(iso1/pt1);
            //h_muo_ptvarcone20_me->Fill(iso1);
        }
    }
    else if(mm()) {
        if(iso0!=0) {
            h_muo_ptvarcone20_mm->Fill(iso0);
            h_muo_ptvarcone20_pt_mm->Fill(iso0/pt0);
            //h_muo_ptvarcone20_mm->Fill(iso0);
        }
    }

    if(bjets.size()>0) {
        float dr0 = leptons.at(0)->DeltaR(*bjets.at(0));
        float dr1 = leptons.at(1)->DeltaR(*bjets.at(0));
        h_dr_l0_b0->Fill(dr0, w());
        h_dr_l1_b0->Fill(dr1, w());

        if(leptons.at(0)->isEle()) {
            h_dr_e0_b0->Fill(dr0, w());
        }
        else {
            h_dr_m0_b0->Fill(dr0, w());
        }

        if(leptons.at(1)->isEle()) {
            h_dr_e1_b0->Fill(dr1, w());
        }
        else {
            h_dr_m1_b0->Fill(dr1, w());
        }
    }
    if(bjets.size()>1) {
        float dr0 = leptons.at(0)->DeltaR(*bjets.at(1));
        float dr1 = leptons.at(1)->DeltaR(*bjets.at(1));
        h_dr_l0_b1->Fill(dr0, w());
        h_dr_l1_b1->Fill(dr1, w());

        if(leptons.at(0)->isEle()) {
            h_dr_e0_b1->Fill(dr0, w());
        }
        else {
            h_dr_m0_b1->Fill(dr0, w());
        }

        if(leptons.at(1)->isEle()) {
            h_dr_e1_b1->Fill(dr1, w());
        }
        else {
            h_dr_m1_b1->Fill(dr1, w());
        }
    }


}
//////////////////////////////////////////////////////////////////////////////
void IsoLooper::save_histos()
{
    m_outfile->cd();

    for(auto & h : histos) h->Write();


    cout << "------------------------------------" << endl;
    cout << "Mean dR_ll      : " << h_dr_ll->GetMean() << endl;
    cout << "Mean dR_ll (EE) : " << h_dr_ll_ee->GetMean() << endl;
    cout << "Mean dR_ll (MM) : " << h_dr_ll_mm->GetMean() << endl;
    cout << "Mean dR_ll (EM) : " << h_dr_ll_em->GetMean() << endl;
    cout << "------------------------------------" << endl;

    stringstream sx;

    cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
    int n_bins_ll = h_dr_ll->GetNbinsX();
    h_dr_ll->Scale(1/h_dr_ll->Integral());
    for(int i = 0; i < n_bins_ll; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " LL " << i << "  " << h_dr_ll->GetBinContent(bin) << endl;
    }
    int n_bins_ee = h_dr_ll_ee->GetNbinsX();
    h_dr_ll_ee->Scale(1/h_dr_ll_ee->Integral());
    for(int i = 0; i < n_bins_ee; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " EE " << i << "  " << h_dr_ll_ee->GetBinContent(bin) << endl;
    }
    int n_bins_mm = h_dr_ll_mm->GetNbinsX();
    h_dr_ll_mm->Scale(1/h_dr_ll_mm->Integral());
    for(int i = 0; i < n_bins_mm; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " MM " << i << "  " << h_dr_ll_mm->GetBinContent(bin) << endl;
    }
    int n_bins_em = h_dr_ll_em->GetNbinsX();
    h_dr_ll_em->Scale(1/h_dr_ll_em->Integral());
    for(int i = 0; i < n_bins_em; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " EM " << i << "  " << h_dr_ll_em->GetBinContent(bin) << endl;
    }

    /*
    cout << "++++++++++++++++++++++++++++++++++++++++" << endl;
    int n_bins = h_dr_l0_b0->GetNbinsX();
    h_dr_l0_b0->Scale(1/h_dr_l0_b0->Integral());
    for(int i = 0; i < n_bins; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " L0 " << i << "  " << h_dr_l0_b0->GetBinContent(bin) << endl;
    }
    n_bins = h_dr_l1_b0->GetNbinsX();
    h_dr_l1_b0->Scale(1/h_dr_l1_b0->Integral());
    for(int i = 0; i < n_bins; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " L1 " << i << "  " << h_dr_l1_b0->GetBinContent(bin) << endl;
    }

    n_bins = h_dr_e0_b0->GetNbinsX();
    h_dr_e0_b0->Scale(1/h_dr_e0_b0->Integral());
    for(int i = 0; i < n_bins; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " E0 " << i << "  " << h_dr_e0_b0->GetBinContent(bin) << endl;
    }

    n_bins = h_dr_e1_b0->GetNbinsX();
    h_dr_e1_b0->Scale(1/h_dr_e1_b0->Integral());
    for(int i = 0; i < n_bins; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " E1 " << i << "  " << h_dr_e1_b0->GetBinContent(bin) << endl;
    }

    n_bins = h_dr_m0_b0->GetNbinsX();
    h_dr_m0_b0->Scale(1/h_dr_m0_b0->Integral());
    for(int i = 0; i < n_bins; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " M0 " << i << "  " << h_dr_m0_b0->GetBinContent(bin) << endl;
    }

    n_bins = h_dr_m1_b0->GetNbinsX();
    h_dr_m1_b0->Scale(1/h_dr_m1_b0->Integral());
    for(int i = 0; i < n_bins; i++) {
        int bin = i+1;
        cout << "DRLL DATA X" << xmass() << " M1 " << i << "  " << h_dr_m1_b0->GetBinContent(bin) << endl;
    }
    */

    m_outfile->Write();
    cout << "IsoLooper::save_histos    Histo file saved to: " << m_outfile->GetName() << endl;

}
//////////////////////////////////////////////////////////////////////////////
void IsoLooper::Terminate()
{
    // close SusyNtAna and print timers
    SusyNtAna::Terminate();
    save_histos();

    return;
}
//////////////////////////////////////////////////////////////////////////////
