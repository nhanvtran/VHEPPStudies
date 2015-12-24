#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include "TROOT.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include <cmath>
#include "TTree.h"

#include <fastjet/ClusterSequence.hh>
#include <fastjet/GhostedAreaSpec.hh>
#include <fastjet/ClusterSequenceArea.hh>
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"
#include "fastjet/tools/MassDropTagger.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "fastjet/PseudoJet.hh"
#include "fastjet/JetDefinition.hh"

#include "fastjet/contrib/Nsubjettiness.hh"
#include "fastjet/contrib/Njettiness.hh"
#include "fastjet/contrib/NjettinessPlugin.hh"
#include "fastjet/contrib/EnergyCorrelator.hh"
// #include "fastjet/contrib/ModifiedMassDropTagger.hh"
#include "fastjet/contrib/SoftDrop.hh"

//#ifdef __MAKECINT__
//#pragma link C++ class vector<float>+;
//#endif

using namespace std;
using namespace fastjet;
using namespace fastjet::contrib;

/*
 
 TO COMPILE:
 
 export ROOTSYS=~/Desktop/root
 export PATH=$ROOTSYS/bin:$PATH
 export LD_LIBRARY_PATH=$ROOTSYS/lib:$LD_LIBRARY_PATH
 export DYLD_LIBRARY_PATH=$ROOTSYS/lib:$DYLD_LIBRARY_PATH 
 
 c++ -o anaSubstructure `root-config --glibs --cflags` `/Users/ntran/Research/CMS/PhysicsAnalysis/boost2013/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins` -lvectorDict -lEnergyCorrelator anaSubstructure.cpp
 
 TO RUN:     
 
 ./analysis01 ../madevent_3/Events/pietro_test14_unweighted_events.lhe ../Z2j/run_01_unweighted_events.lhe 
 
 */

//! ========================================================================================

////////////////////-----------------------------------------------
// Global variables
ifstream fin;
ifstream finMC;
ifstream finCalo;

int evtCtr;

int njets;
std::vector<float> jpt;
std::vector<float> jeta;
std::vector<float> jphi;
std::vector<float> jmass;
std::vector<float> jmultiplicity;
std::vector<float> jisleptag;
std::vector<float> jmass_sd;

std::vector<float> jpt_calo;

////////////////////-----------------------------------------------
void readEvent( std::vector< fastjet::PseudoJet > &allParticles );
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles );
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles );

void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles, std::vector < fastjet::PseudoJet > caloclusters);

////////////////////////////////////////////////////////////////////////////////////////
int main (int argc, char **argv) {
    
    // std::cout << "hello world" << std::endl;
    std::string type = argv[1];   // type "gg" or "qq"

    std::vector < fastjet::PseudoJet > allParticles;
    std::vector < fastjet::PseudoJet > allParticlesCalo;
    std::vector < fastjet::PseudoJet > allParticlesMC;

    char fname[150];
    sprintf( fname, "dat/%s.dat", type.c_str() );
    fin.open(fname);

    char fnameMC[150];
    sprintf( fnameMC, "dat/of_status3.dat" );
    finMC.open(fnameMC);

    char fnameCalo[150];
    sprintf( fnameCalo, "dat/of_CaloHits.dat" );
    finCalo.open(fnameCalo);

    char outName[192];
    sprintf( outName, "dat/%s.root", type.c_str() );
    TFile *f = TFile::Open(outName,"RECREATE");
    TTree *t = new TTree("t","Tree with vectors");
    t->Branch("njets"          , &njets      );
    t->Branch("jpt"            , &jpt        );
    t->Branch("jeta"           , &jeta       );
    t->Branch("jphi"           , &jphi       );
    t->Branch("jmass"          , &jmass      );    
    t->Branch("jmass_sd"       , &jmass_sd      );    
    t->Branch("jmultiplicity"  , &jmultiplicity      );    
    t->Branch("jisleptag"      , &jisleptag      );    

    t->Branch("jpt_calo"            , &jpt_calo        );

    int ctr = 0;
    while(true){
        
        //  INIT
        jpt.clear();
        jeta.clear();
        jphi.clear();        
        jmass.clear();
        jmass_sd.clear();        
        jmultiplicity.clear();
        jisleptag.clear();
        jpt_calo.clear();

        readEvent( allParticles );
        readEventMC( allParticlesMC );
        readEventCalo( allParticlesCalo );
        std::cout << "size of collection = " << allParticles.size() << "," << allParticlesMC.size() << ", " << allParticlesCalo.size() << std::endl;


        if (ctr > 0){
            analyzeEvent( allParticles, 0.4, allParticlesMC, allParticlesCalo );
            t->Fill();
        }

        allParticles.clear();
        allParticlesMC.clear();
        allParticlesCalo.clear();
        
        ctr++;
        std::cout << "ctr = " << ctr << std::endl;

        if(fin.eof()) break;
    }

    f->cd();
    t->Write();
    f->Close();

}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEvent( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 5.;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        fin  >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
    
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(fin.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventMC( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 5.;
    
    float npart, px, py, pz, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finMC  >> pdgid >> px >> py >> pz >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        if (pdgid == -99){
            return;
        }        
    
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( px, py, pz, e );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finMC.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void readEventCalo( std::vector< fastjet::PseudoJet > &allParticles ){
    
    // hard-coded! 
    float etaMax = 5.;
    
    float npart, x, y, z, e, pdgid, isCh, isPU = 0;
    
    int ctr = 0;
    while(true){
        
        finCalo  >> pdgid >> x >> y >> z >> e;         
        // std::cout << "pdgid = " << pdgid << ", " << px << ", " << py << std::endl;
    
        TVector3 caloposition(x,y,z);
        double curPhi = caloposition.Phi();
        double curEta = caloposition.Eta();
        double curPt = sin(caloposition.Theta())*e;
        TLorentzVector cp4;
        cp4.SetPtEtaPhiE(curPt,curEta,curPhi,e);

        if (pdgid == -99){
            return;
        }        
    
        // fill vector of pseudojets
        fastjet::PseudoJet curPseudoJet( cp4.Px(), cp4.Py(), cp4.Pz(), cp4.E() );
        curPseudoJet.set_user_index(pdgid);
        if (fabs(curPseudoJet.eta()) < etaMax){
            allParticles.push_back( curPseudoJet );
        }
        
        if(finCalo.eof()) break;
        ctr++;
        // std::cout << "ctr2 = " << ctr << std::endl;
    }
    
}
////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////
void analyzeEvent(std::vector < fastjet::PseudoJet > particles, float rVal, std::vector < fastjet::PseudoJet > MCparticles, std::vector< fastjet::PseudoJet > caloclusters){

    std::cout << "analyzing event..." << particles.size() << std::endl;
    double rParam = rVal;
    fastjet::JetDefinition jetDef(fastjet::antikt_algorithm, rParam);    

    int activeAreaRepeats = 1;
    double ghostArea = 0.01;
    double ghostEtaMax = 7.0;
    
    fastjet::GhostedAreaSpec fjActiveArea(ghostEtaMax,activeAreaRepeats,ghostArea);
    fastjet::AreaDefinition fjAreaDefinition( fastjet::active_area, fjActiveArea );
    
    fastjet::ClusterSequenceArea* thisClustering = new fastjet::ClusterSequenceArea(particles, jetDef, fjAreaDefinition);
    fastjet::ClusterSequenceArea* caloClustering = new fastjet::ClusterSequenceArea(caloclusters, jetDef, fjAreaDefinition);
    std::vector<fastjet::PseudoJet> out_jets = sorted_by_pt(thisClustering->inclusive_jets(25.0));
    std::vector<fastjet::PseudoJet> calo_jets = sorted_by_pt(caloClustering->inclusive_jets(25.0));

    double beta_sd = 1.0;
    double zcut_sd = 0.1;
    double mu_sd   = 1.0;
    fastjet::contrib::SoftDrop soft_drop_mmdt(0.0, zcut_sd, mu_sd);

    njets = out_jets.size();
    std::cout << "number of high pT jets! = " << njets << std::endl;
    for (unsigned int i = 0; i < out_jets.size(); i++){

        double isleptag = 0;
        double minDR = 9999.;
        for (unsigned int j = 0; j < MCparticles.size(); j++){
            // std::cout << "MCparticles = " << MCparticles[j].pt() << "," << MCparticles[j].user_index() << std::endl;
            double pdgid =  fabs(MCparticles[j].user_index());
            double dr = sqrt( pow(MCparticles[j].phi()-out_jets[i].phi(),2) + pow(MCparticles[j].eta()-out_jets[i].eta(),2) );
            // std::cout << "dr = " << dr << "," << pdgid << std::endl;
            if (minDR > dr && (pdgid >= 11) ){ minDR = dr; }
        }
        // std::cout<< "minDR = " << minDR << std::endl;
        if (minDR < 0.8){ isleptag = 1.; }

        jpt.push_back( out_jets[i].pt() );
        jeta.push_back( out_jets[i].eta() );
        jphi.push_back( out_jets[i].phi() );
        jmass.push_back( out_jets[i].m() );  
        jmass_sd.push_back( soft_drop_mmdt( out_jets.at(i) ).m() );        
        jmultiplicity.push_back( (float) out_jets.at(i).constituents().size() );
        jisleptag.push_back( isleptag );
      
    }

    // calo jets...
    int njets_calo = calo_jets.size();
    std::cout << "number of high pT calo jets! = " << njets_calo << std::endl;
    for (unsigned int i = 0; i < calo_jets.size(); i++){

        jpt_calo.push_back( calo_jets[i].pt() );
      
    }    

}

