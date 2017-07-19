#ifndef SusyNtuple_IsoLooper_h
#define SusyNtuple_IsoLooper_h

//ROOT
#include "TTree.h"
#include "TChain.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"

//SusyNtuple
#include "SusyNtuple/SusyNtAna.h"
#include "SusyNtuple/SusyNtTools.h"

//std/stl
#include <fstream>
#include <vector>

/////////////////////////////////////////////////////////////
//
// IsoLooper
// Class auto-generated with SusyNtuple/make_susy_skeleton on 2017-07-05 11:48
//
/////////////////////////////////////////////////////////////

// for TSelector analysis loopers processing susyNt you MUST inherit from SusyNtAna
// in order to pick up the susyNt class objects
class IsoLooper : public SusyNtAna
{

    public :
        IsoLooper();
        virtual ~IsoLooper() {};

        void set_debug(int dbg) { m_dbg = dbg; }
        int dbg() { return m_dbg; }

        void set_chain(TChain* chain) { m_input_chain = chain; }
        TChain* chain() { return m_input_chain; }

        void set_xmass(int xmass) { m_res_mass = xmass; }
        int xmass() { return m_res_mass; }

        ////////////////////////////////////////////
        // analysis methods
        ////////////////////////////////////////////

        void initialize_histos();

        float w() { return m_mc_weight * m_lep_sf; }

        bool ee() { return (m_lep_type == ET_ee); }
        bool mm() { return (m_lep_type == ET_mm); }
        bool em() { return (m_lep_type == ET_em); }
        bool me() { return (m_lep_type == ET_me); }

        // standard ATLAS event cleaning
        bool passEventCleaning(const MuonVector& preMuons, const MuonVector& baseMuons,
                const JetVector& baseJets);

        bool passSignalSelection(const LeptonVector& leptons, const JetVector& jets, LeptonVector& outleptons);
        bool is_signal_ele(Electron* el);
        bool is_signal_mu(Muon* mu);

        void get_bjets(const JetVector& injets, JetVector& outjets);
        bool is_bjet(const Jet* j);

        void fill_histos(const LeptonVector& leptons, const JetVector& bjets);

        void save_histos();

        ////////////////////////////////////////////
        // TSelector methods override
        ////////////////////////////////////////////
        virtual void Begin(TTree* tree); // Begin is called before looping on entries
        virtual Bool_t Process(Long64_t entry); // Main event loop function called on each event
        virtual void Terminate(); // Terminate is called after looping has finished


    private :
        int m_dbg;
        TChain* m_input_chain; // the TChain object we are processing
        float m_mc_weight;
        float m_lep_sf;
        float m_weight;
        float m_res_mass;
        DiLepEvtType m_lep_type;

        // file
        TFile* m_outfile;

        // histos
        TH1F* h_ele_ptvarcone20_ee;
        TH1F* h_ele_ptvarcone20_em;
        TH1F* h_ele_ptvarcone20_me;
        TH1F* h_muo_ptvarcone20_mm;
        TH1F* h_muo_ptvarcone20_em;
        TH1F* h_muo_ptvarcone20_me;
        TH1F* h_ele_ptvarcone20_pt_ee;
        TH1F* h_ele_ptvarcone20_pt_em;
        TH1F* h_ele_ptvarcone20_pt_me;
        TH1F* h_muo_ptvarcone20_pt_mm;
        TH1F* h_muo_ptvarcone20_pt_em;
        TH1F* h_muo_ptvarcone20_pt_me;
        TH1F* h_dphi;
        TH2F* h2_dphi_iso_ee;
        TH2F* h2_dphi_iso_mm;
        TH1F* h_dr_ll;
        TH1F* h_dr_ll_ee;
        TH1F* h_dr_ll_mm;
        TH1F* h_dr_ll_em;

        TH1F* h_dr_l0_b0;
        TH1F* h_dr_l1_b0;
        TH1F* h_dr_e0_b0;
        TH1F* h_dr_e1_b0;
        TH1F* h_dr_m0_b0;
        TH1F* h_dr_m1_b0;

        TH1F* h_dr_l0_b1;
        TH1F* h_dr_l1_b1;
        TH1F* h_dr_e0_b1;
        TH1F* h_dr_e1_b1;
        TH1F* h_dr_m0_b1;
        TH1F* h_dr_m1_b1;

        std::vector<TH1F*> histos;

}; //class


#endif
