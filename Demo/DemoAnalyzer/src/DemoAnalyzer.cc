// -*- C++ -*-
//
// Package:    DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/src/DemoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  
//         Created:  Mon May  4 15:24:13 CEST 2015
// $Id$
// ..
//

// *** UNDER CONSTRUCTION ! ***
// *** not (yet) for public use ! *** 

// ***************************************************************************
// version of DEMO setup provided by CMS open data access team               *
// expanded/upgraded to contain validation examples for                      *
// - general reference distributions (also technical)                        *
// - minumum bias track multiplicities and pt/eta distributions (QCD-10-006) *
// - J/psi mass peak (BPH-10-002)                                            *
// - dimuon mass spectrum (MUO-10-004)                                       *
// - Z mass peak (EWK-10-)                                                   *
// - inelastic proton-proton cross section (FWD-11-001)                      *
//                                                                           *
//                         I.Dutta and A.Geiser, 30.6.15                     *
// ***************************************************************************

// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


//------ ADD EXTRA HEADER FILES--------------------//
#include"math.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Ref.h"
// for histogramming
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TStyle.h"
// for tracking information
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
// for vertex information 
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
// for muon information
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
//for beamspot information
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
//for electron informaton
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

//for MET informaton

#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETFwd.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METFwd.h"

// for jet information
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/JetReco/interface/TrackJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/JetExtendedAssociation.h"
#include "DataFormats/JetReco/interface/JetID.h"

#include <iostream>
using namespace std;
// class declaration
//

class DemoAnalyzer : public edm::EDAnalyzer 
{
   public:
      explicit DemoAnalyzer(const edm::ParameterSet&);
      ~DemoAnalyzer();

  // declare variables
  double deltaz;
  double firstz;

   private:
      virtual void beginJob() ;
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

  // ----------member data ---------------------------
  // declare histograms
  TH1D *histset[1000];
  TH2D *hxhy[1000];
};

//
// constants, enums and typedefs
//
double s1,s2,s,scorr,s3,s4,pz,rap,w,M, pt_jpsi;
double sqm1 = pow(0.105658,2);

//(0.105658)*(0.105658);
double mumass=0.105658;
double jpsim=3.097;
double Kmass=0.49367;
double Pimass=0.13957;
double MD0Actual=1.86484;  // [PDG]
double MDstarActual=2.01022;

double binsy[43]={1};
double binsz[101];



//
// static data member definitions
//

//
// constructors and destructor


DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig)

{
  //now do what ever initialization is needed
  edm::Service<TFileService> fs;
//Loop for
    for (int i=0; i!=43; i++){
        
        if (i>=0 && i<=19) binsy[i]= 0.1*i;
        else if (i>19 && i<=30) binsy[i]=2+(i-20)*0.5;
        else if (i>=31 && i<=38) binsy[i]=8+(i-31);
        else if (i==39) binsy[i]=17;
        else if (i==40) binsy[i]=20;
        else if (i==41) binsy[i]=25;
        else if (i==42) binsy[i]=30;
        //cout<< binsy[i] << endl;
    }

  ////////------------------------------------ Book Histograms ------------------------------------////////
    {
        
    //-------------------------Minimum Bias---------------------------------//
    
  histset[6]  = fs->make<TH1D>("tracks" , "Tracks" , 300 , 0 , 300 ); 					// Track Multiplicity

  histset[34] = fs->make<TH1D>("momentum", "Momentum",200, 0,20.0); 					// Track momentum
  histset[35] = fs->make<TH1D>("pt", "PT",200, 0,4.0); 							// Track pt
  histset[36] = fs->make<TH1D>("posx", "Position X",100, 0.0,.2); 					// Track position x
  histset[37] = fs->make<TH1D>("posy", "Position Y",100, -.1,.1);					// track position y
  histset[7]  = fs->make<TH1D>("posz", "Position Z",100, -25.0,25.0); 					// Track position z
  histset[38] = fs->make<TH1D>("eta", "Eta",100, -3.5,3.5); 						// Track eta
  histset[81] = fs->make<TH1D>("track_phi", "Track_Phi",314, -3.15,3.15); 				// Track phi
  histset[8]  = fs->make<TH1D>("propvertex" , "ProperVertices" , 20 , 0 , 20 ); 			// Proper vertices, i.e.after checking the number of tracks for the vertex
  histset[9]  = fs->make<TH1D>("0.2pt" , "0.2pt" , 100 , 0 , 4.0 );					// Tracks with absolute value of eta less than 0.2
  histset[10] = fs->make<TH1D>("0.4pt" , "0.4pt" , 100 , 0 , 4.0 );					// Tracks with absolute value mof eta between 0.2 and 0.4
  histset[11] = fs->make<TH1D>("0.6pt" , "0.6pt" , 100 , 0 , 4.0 );					// and so on...
  histset[12] = fs->make<TH1D>("0.8pt" , "0.8pt" , 100 , 0 , 4.0 );
  histset[13] = fs->make<TH1D>("1.0pt" , "1.0pt" , 100 , 0 , 4.0 );
  histset[14] = fs->make<TH1D>("1.2pt" , "1.2pt" , 100 , 0 , 4.0 );
  histset[15] = fs->make<TH1D>("1.4pt" , "1.4pt" , 100 , 0 , 4.0 );
  histset[16] = fs->make<TH1D>("1.6pt" , "1.6pt" , 100 , 0 , 4.0 );
  histset[17] = fs->make<TH1D>("1.8pt" , "1.8pt" , 100 , 0 , 4.0 );
  histset[18] = fs->make<TH1D>("2.0pt" , "2.0pt" , 100 , 0 , 4.0 );
  histset[19] = fs->make<TH1D>("2.2pt" , "2.2pt" , 100 , 0 , 4.0 );
  histset[20] = fs->make<TH1D>("2.4pt" , "2.4pt" , 100 , 0 , 4.0 );
  histset[21] = fs->make<TH1D>("2.6pt" , "2.6pt" , 100 , 0 , 4.0 );
  histset[22] = fs->make<TH1D>("deltaz" , "DeltaZ" , 300 , 0 , 30.0 );	// the absolute value of the distance in the z coordinate of all vertices from Primvertex->begin()
//histset[23] = fs->make<TH1D>("track_vertex_x" , "TrackVertexX" , 100 , -30.0 , 30.0 );
//histset[24] = fs->make<TH1D>("track_vertex_y" , "TrackVertexY" , 100 , -30.0 , 30.0 );
//histset[25] = fs->make<TH1D>("track_vertex_z" , "TrackVertexZ" , 100 , -30.0 , 30.0 );
  histset[26] = fs->make<TH1D>("vertex" , "Vertices" , 20 , 0 , 20 ); 					// Just vertex multiplicity, without any corrections

  histset[27] = fs->make<TH1D>("VertexTrack_pt", "Vertex.Track_Pt",200,0,4.0); 				// Pt of tracks associated with a particular vertex
  histset[28] = fs->make<TH1D>("VertexTrack_momentum", "VertexTrack_Momentum",200, 0,20.0); 		// Track momentum
  histset[82] = fs->make<TH1D>("Vertex.Track_phi", "Vertex.Track_Phi",100, -3.5,3.5); 			//Track phi
  
// set axis labels
    {
  histset[34]->GetXaxis()->SetTitle("Momentum (in GeV/c)");
  histset[34]->GetYaxis()->SetTitle("Number of Events");

  histset[35]->GetXaxis()->SetTitle("Transverse Momentum (in GeV/c)");
  histset[35]->GetYaxis()->SetTitle("Number of Events");

  histset[36]->GetXaxis()->SetTitle("Position X of Vertices (in cm)");
  histset[36]->GetYaxis()->SetTitle("Number of Events");

  histset[37]->GetXaxis()->SetTitle("Position Y of Vertices (in cm)");
  histset[37]->GetYaxis()->SetTitle("Number of Events");

  histset[7]->GetXaxis()->SetTitle("Position z of Vertices (in cm)");
  histset[7]->GetYaxis()->SetTitle("Number of Events");

  histset[38]->GetXaxis()->SetTitle("Eta (in radians)");
  histset[38]->GetYaxis()->SetTitle("Number of Events");

  histset[8]->GetXaxis()->SetTitle("Number of vertices with a non zero number of tracks associated to it");
  histset[8]->GetYaxis()->SetTitle("Number of Events");

  histset[9]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<0.2(in GeV/c)");
  histset[9]->GetYaxis()->SetTitle("Number of Events");
  histset[9]->SetLineColor(kRed);

  histset[10]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<0.4 and |eta|>=0.2(in GeV/c)");
  histset[10]->GetYaxis()->SetTitle("Number of Events");

  histset[11]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<0.6 and |eta|>=0.4(in GeV/c)");
  histset[11]->GetYaxis()->SetTitle("Number of Events");

  histset[12]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<0.8 and |eta|>=0.6(in GeV/c)");
  histset[12]->GetYaxis()->SetTitle("Number of Events");

  histset[13]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<1.0 and |eta|>=0.8(in GeV/c)");
  histset[13]->GetYaxis()->SetTitle("Number of Events");

  histset[14]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<1.2 and |eta|>=1.0(in Gev/c)");
  histset[14]->GetYaxis()->SetTitle("Number of Events");

  histset[15]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<1.4 and |eta|>=1.2(in GeV/c)");
  histset[15]->GetYaxis()->SetTitle("Number of Events");

  histset[16]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<1.6 and |eta|>=1.4(in GeV/c)");
  histset[16]->GetYaxis()->SetTitle("Number of Events");

  histset[17]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<1.8 and |eta|>=1.6(in Gev/c)");
  histset[17]->GetYaxis()->SetTitle("Number of Events");

  histset[18]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<2.0 and |eta|>=1.8(in GeV/c)");
  histset[18]->GetYaxis()->SetTitle("Number of Events");

  histset[19]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<2.2 and |eta|>=2.0(in Gev/c)");
  histset[19]->GetYaxis()->SetTitle("Number of Events");

  histset[20]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<2.4 and |eta|>=2.2(in GeV/c)");
  histset[20]->GetYaxis()->SetTitle("Number of Events");

  histset[21]->GetXaxis()->SetTitle("Transverse Momentum for |eta|<2.6 and |eta|>=2.4(in GeV/c)");
  histset[21]->GetYaxis()->SetTitle("Number of Events");

  histset[22]->GetXaxis()->SetTitle("Absolute value of the z distance of Vertices from the first primary vertex (in cm)");
  histset[22]->GetYaxis()->SetTitle("Number of Events");

  histset[26]->GetXaxis()->SetTitle("Vertex multiplicity");
  histset[26]->GetYaxis()->SetTitle("Number of Events");

  histset[28]->GetXaxis()->SetTitle("Momentum of tracks associated with a vertex (in GeV/c)");
  histset[28]->GetYaxis()->SetTitle("Number of Events");

  histset[27]->GetXaxis()->SetTitle("Transverse Momentum of tracks associated with a vertex (in GeV/c)");
  histset[27]->GetYaxis()->SetTitle("Number of Events");
}

 //------------------------------Muons and Onia--------------------------------------//
 // use either global muons from inclusive Mu sample 
 // or tracker muons from Onia sample (preselected 2.5<m<4.1 GeV) 


  histset[1] = fs->make<TH1D>("GMmomentum" , "GM_Momentum" , 240 , 0. , 120. );				//TrackCollection_GMuon momentum
  histset[2] = fs->make<TH1D>("GM_Transverse_momentum" , "TransverseMomentum" , 240 , 0. , 120. );	//TrackCollection_GMuon Transverse_momentum
  histset[3] = fs->make<TH1D>("GM_eta" , "GM_Eta" , 140 ,-3.5 , 3.5 );					//TrackColletion_GMuon eta
  histset[83]= fs->make<TH1D>("GM_phi", "GM_phi",314, -3.15,3.15);
  histset[4] = fs->make<TH1D>("GMmultiplicty" , "GMmultiplicity" , 8 , 0 , 8 );				//TrackCollection_GMuon multiplicity
  histset[5] = fs->make<TH1D>("GMmass" , "GMmass" ,120 , 0. , 12. );					//TrackCollection_GMuon mas
  histset[44]= fs->make<TH1D>("GMmass_extended" , "GMmass" ,120 , 0. , 120. );				//TrackCollection_GMuon mass
  histset[29]= fs->make<TH1D>("Muon_momentum" , "Muon_momentum" ,100 , 0. , 20. );			//MuonCollection_ momentum
  histset[30]= fs->make<TH1D>("Muon_pt" , "Muon_PT" ,100 , 0. , 20. );					//MuonCollection_pt
  histset[31]= fs->make<TH1D>("Muon_eta" , "Muon_eta" ,140 , -3.5 , 3.5 );				//MuonCollection_ eta
  histset[84]= fs->make<TH1D>("Muon_phi", "Muon_phi",314, -3.15,3.15);
  histset[32]= fs->make<TH1D>("TM_mass" , "TM_mass" ,120 , 0. , 12. );					//MuonCollection_TMuon mass
  histset[33]= fs->make<TH1D>("NMuons" , "Muon_Size" ,10 , 0. ,10  );					//MuonCollection_Muon size

  histset[39]= fs->make<TH1D>("TM_mass1.2_j_psi" , "TM_mass" ,90 ,2.6  , 3.5 );				//MuonCollection_TMuon mass_j/psi
  histset[40]= fs->make<TH1D>("TM_mass1.6_j_psi" , "TM_mass" ,90 , 2.6 , 3.5 );				//MuonCollection_TMuon mass_j/psi
  histset[41]= fs->make<TH1D>("TM_mass2.4_j_psi" , "TM_mass" ,90 ,2.6 , 3.5 );				//MuonCollection_TMuon mass_j/psi

  histset[42]= fs->make<TH1D>("TM_mass1_upsilon" , "TM_mass" ,80 ,8.  , 12. );				//MuonCollection_TMuon mass_upsilon

  histset[43]= fs->make<TH1D>("TM_mass2.4_upsilon" , "TM_mass" ,80 ,8. ,12. );				//MuonCollection_TMuon mass_ upsilon
  histset[45]= fs->make<TH1D>("TM_mass_likeCharges" , "TM_mass_Like Charges" ,120 , 0. , 12. );		//MuonCollection_TMuon mass_ FOR LIKE CHARGES
  histset[46]= fs->make<TH1D>("GM_mass_likeCharges" , "GM_mass_Like Charges" ,120 , 0. , 12. );		//TrackCollection_GMuon mass_ FOR LIKE CHARGES


  histset[47]= fs->make<TH1D>("GM_mass1.2_j_psi" , "GM_mass" ,90 ,2.6  , 3.5 );				//TrackCollection_GMuon mass_j/psi
  histset[48]= fs->make<TH1D>("GM_mass1.6_j_psi" , "GM_mass" ,90 , 2.6 , 3.5 );				//TrackCollection_GMuon mass_j/psi
  histset[49]= fs->make<TH1D>("GM_mass2.4_j_psi" , "GM_mass" ,90 ,2.6 , 3.5 );				//TrackCollection_gMuon mass_j/psi

  histset[50]= fs->make<TH1D>("TM_chi2" , "TM_Chi2" ,500 ,0 , 100 );					//MuonCollection
  histset[51]= fs->make<TH1D>("TM_Ndof" , "TM_Ndof" ,60 ,0 , 60 );					//MuonCollection
  histset[52]= fs->make<TH1D>("TM_NormalizedChi2" , "TM_NormalizedChi2" ,200 ,0 , 20 );			//MuonCollection

  histset[53]= fs->make<TH1D>("GM_chi2" , "GM_Chi2" ,300 ,0 , 150 );					//TrackCollection
  histset[54]= fs->make<TH1D>("GM_ndof" , "GM_ndof" ,100 ,0 , 100 );					//TrackCollection
  histset[55]= fs->make<TH1D>("GM_normalizedchi2" , "GM_normalizedChi2" ,200 ,0 , 20 );			//TrackCollection

  histset[56]= fs->make<TH1D>("Muon_chi2" , "Muon_Chi2" ,500 ,0 , 100 );				//MuonCollection
  histset[57]= fs->make<TH1D>("Muon_Ndof" , "Muon_Ndof" ,60 ,0 , 60 );					//MuonCollection
  histset[58]= fs->make<TH1D>("Muon_NormalizedChi2" , "Muon_NormalizedChi2" ,200 ,0 , 20 );		//MuonCollection

  histset[59]= fs->make<TH1D>("GM_HitsOK" , "GM_Tracks with good Hits" ,10 ,0 , 10 );			//TrackCollection

  histset[60]= fs->make<TH1D>("GM_validhits" , "GM_ValidHits" ,100, 0. ,100 );				//TrackCollection
  histset[61]= fs->make<TH1D>("GM_pixelhits" , "GM_pixelhits" ,14 , 0. ,14 );				//TrackCollection

  histset[62]= fs->make<TH1D>("TM_HitsOK" , "TM_Tracks with good Hits" ,10 ,0 , 10 );			//MuonCollection

  histset[63]= fs->make<TH1D>("TM_validhits" , "TM_ValidHits" ,40, 0. ,40 );				//MuonCollection
  histset[64]= fs->make<TH1D>("TM_pixelhits" , "TM_pixelhits" ,14 , 0. ,14 );				//MuonCollection


  histset[65]= fs->make<TH1D>("TM_mass_C" , "TM_mass" ,120 , 0. , 12. );				//MuonCollection_TMuon mass_corrected
  histset[66]= fs->make<TH1D>("GMmass_C" , "GMmass" ,120 , 0. , 12. );					//TrackCollection_GMuon mass_corrected
  histset[67]= fs->make<TH1D>("GMmass_extended_C" , "GMmass" ,120 , 0. , 120. );			//TrackCollection_GMuon mass_corrected

  histset[68]= fs->make<TH1D>("GM_mass1.2_j_psi_C" , "GM_mass" ,90 ,2.6  , 3.5 );			//TrackCollection_GMuon mass_j/psi
  histset[69]= fs->make<TH1D>("GM_mass1.6_j_psi_C" , "GM_mass" ,90 , 2.6 , 3.5 );			//TrackCollection_GMuon mass_j/psi
  histset[70]= fs->make<TH1D>("GM_mass2.4_j_psi_C" , "GM_mass" ,90 ,2.6 , 3.5 );			//TrackCollection_gMuon mass_j/psi

  histset[71]= fs->make<TH1D>("TM_mass1_upsilon_C" , "TM_mass" ,80 ,8.  , 12. );			//MuonCollection_TMuon mass_upsilon

  histset[72]= fs->make<TH1D>("TM_mass2.4_upsilon_C" , "TM_mass" ,80 ,8. ,12. );			//MuonCollection_TMuon mass_ upsilon


  histset[73]= fs->make<TH1D>("TM_mass1.2_j_psi_C" , "TM_mass" ,90 ,2.6  , 3.5 );			//MuonCollection_TMuon mass_j/psi
  histset[74]= fs->make<TH1D>("TM_mass1.6_j_psi_C" , "TM_mass" ,90 , 2.6 , 3.5 );			//MuonCollection_TMuon mass_j/psi
  histset[75]= fs->make<TH1D>("TM_mass2.4_j_psi_C" , "TM_mass" ,90 ,2.6 , 3.5 );			//MuonCollection_TMuon mass_j/psi

  histset[76]= fs->make<TH1D>("GM_mass1_upsilon" , "GM_mass" ,80 ,8.  , 12. );				//TrackCollection_TMuon mass_upsilon

  histset[77]= fs->make<TH1D>("GM_mass2.4_upsilon" , "GM_mass" ,80 ,8. ,12. );				//TrackCollection_TMuon mass_ upsilon

  histset[78]= fs->make<TH1D>("GM_mass1_upsilon_C" , "GM_mass" ,80 ,8.  , 12. );			//TrackCollection_TMuon mass_upsilon

  histset[79]= fs->make<TH1D>("GM_mass2.4_upsilon_C" , "GM_mass" ,80 ,8. ,12. );			//TrackCollection_TMuon mass_ upsilon

  histset[80]=fs->make<TH1D>("GM_Zmass", "GM_Zmass",30, 60.0,120.0); 
  histset[94]=fs->make<TH1D>("GM_Zmass_cuts", "GM_Zmass_cuts",30, 60.0,120.0); 

  histset[95]=fs->make<TH1D>("GM_Zmass_cuts_like_charges", "GM_Zmass_cuts_like_charges",30, 60.0,120.0); 

  histset[85]= fs->make<TH1D>("eta_no_mu", "Eta_with no Muon",100, -3.5,3.5); 				//Track eta if Muon not found

  histset[86]= fs->make<TH1D>("tracks_per_vertex_using_counter" , "Tracks per Vertex using counter" , 200 , 0 , 200 );// Track Multiplicity
  histset[87]= fs->make<TH1D>("tracks_per_vertex" , "Tracks per Vertex" , 200 , 0 , 200 );		// Track Multiplicity

  histset[88]=fs->make<TH1D>("dx_muon" , "dx_muon" , 100 , 0 , 1 );
  histset[89]=fs->make<TH1D>("dy_muon" , "dy_muon" , 100 , 0 , 1 );
  histset[90]=fs->make<TH1D>("dz_muon" , "dz_muon" , 200 , 0 , 2 );
  histset[91]=fs->make<TH1D>("dxy_Gmuon" , "dxy_gmuon" , 100 , 0 , 1 );
  histset[92]=fs->make<TH1D>("goodMuonChamberHit_Gmuon" , "goodMuonChamberHit_gmuon" , 40,0 , 40 );
  histset[93]=fs->make<TH1D>("MuonStations_Tmuon" , "MuonStations_Tmuon" , 2,0 , 2 );

  histset[96]= fs->make<TH1D>("Electron_momentum" , "Electron_momentum" ,400 , 0. , 200. );		//ElectronCollection_ momentum
  histset[97]= fs->make<TH1D>("Electron_pt" , "Electron_PT" ,400 , 0. , 200. );				//ElectronCollection_pt
  histset[98]= fs->make<TH1D>("Electron_eta" , "Electron_eta" ,140 , -3.5 , 3.5 );			//ElectronCollection_ eta
  histset[99]= fs->make<TH1D>("Electron_phi", "Electron_phi",314, -3.17,3.17);
  histset[104]= fs->make<TH1D>("Z_electron_mass" , "electron_mass" ,30 , 60. , 120. );			//ElectronCollection mass
  histset[105]= fs->make<TH1D>("Nelectrons" , "Electron_Size" ,10 , 0. ,10  );				//ElectronCollection_Muon size
  histset[106]= fs->make<TH1D>("Super Cluster_eta" , "Super Cluster_eta" ,140 , 0 , 3.5 );		//Electron super cluster_ eta
  histset[107]= fs->make<TH1D>("Super Cluster_rawenergy" , "Super Cluster_rawenergy" ,200 , 0 , 200 );	//Electron super cluster_ energy
  histset[118]= fs->make<TH1D>("Electron Et" , "Electron Et" ,200 , 0 , 200 );				//Electron et

  histset[119]= fs->make<TH1D>("Z_electron_mass_cuts" , "electron_mass_cuts" ,30 , 60. , 120. );	//ElectronCollection mass
  histset[108]= fs->make<TH1D>("Barrel_Electron_track isolation" , "Barrel_Electron Track Isolation" ,500 , 0. , 0.5 );
  histset[109]= fs->make<TH1D>("Barrel_Electron_Ecal isolation" , "Barrel_Electron ECAL Isolation" ,500 , 0. ,0.5 );
  histset[110]= fs->make<TH1D>("Barrel_Electron_hcal isolation" , "Barrel_Electron HCAL Isolation" ,500 , 0. ,0.5 );
  histset[111]= fs->make<TH1D>("Electron_track missing hit" , "Electron Track missing hits" ,5 , 0. ,5);
  histset[112]= fs->make<TH1D>("Electron_Dcot" , "Electron Dcot" ,1000 , -0.05 ,0.05);
  histset[113]= fs->make<TH1D>("Electron_Dist" , "Electron Dist" ,1000 , -0.05 ,0.05);
  histset[114]= fs->make<TH1D>("Barrel_Electron_sigmaIetaIeta" , "Barrel_Electron sigmaIetaIeta" ,300 , 0. ,0.03);
  histset[115]= fs->make<TH1D>("Barrel_Electron_deltaphi" , "Barrel_Electron deltaphi" ,10000 , -0.5 ,0.5);
  histset[116]= fs->make<TH1D>("Barrel_Electron_deltaeta" , "Barrel_Electron deltaeta" ,10000 , -0.05 ,0.05);
  histset[117]= fs->make<TH1D>("Barrel_Electron_H/E" , "Barrel_Electron H/E" ,1200 , 0. ,0.2);


  histset[120]= fs->make<TH1D>("Endcap_Electron_track isolation" , "Endcap_Electron Track Isolation" ,500 , 0. , 0.5 );
  histset[121]= fs->make<TH1D>("Endcap_Electron_Ecal isolation" , "Endcap_Electron ECAL Isolation" ,500 , 0. ,0.5 );
  histset[122]= fs->make<TH1D>("Endcap_Electron_hcal isolation" , "Endcap_Electron HCAL Isolation" ,500 , 0. ,0.5 );
  histset[123]= fs->make<TH1D>("Endcap_Electron_sigmaIetaIeta" , "Endcap_Electron sigmaIetaIeta" ,1000 , 0. ,0.1);
  histset[124]= fs->make<TH1D>("Endcap_Electron_deltaphi" , "Endcap_Electron deltaphi" ,1000 , -0.5,0.5);
  histset[125]= fs->make<TH1D>("Endcap_Electron_deltaeta" , "Endcap_Electron deltaeta" ,400 , -0.01 ,0.01);
  histset[126]= fs->make<TH1D>("Endcap_Electron_H/E" , "Endcap_Electron H/E" ,1000 , 0. ,0.1);
  histset[127]= fs->make<TH1D>("Z_electron_mass_likeCharges" , "electron_mass_likeCharges" ,30 , 60. , 120. );		//ElectronCollection mass
  histset[128]= fs->make<TH1D>("Z_electron_mass_cuts_likeCharges" , "electron_mass_cuts_likeCharges" ,30 , 60. , 120. );//ElectronCollection mass
  histset[129]= fs->make<TH1D>("W_muon_transverseMass" , "muon_transverseMass" ,30 , 0. , 120. );			//MuonCollection mass
  histset[130]= fs->make<TH1D>("W_muon_transverseMass_cuts" , "muon_transverseMass_cuts" ,30 , 0. , 120. );		//MuonCollection mass
  histset[131]= fs->make<TH1D>("W_electron_transverseMass" , "electron_transverseMass" ,30 , 0. , 120. );		//ElectronCollection mass
  histset[132]= fs->make<TH1D>("W_electron_transverseMass_cuts" , "electron_transverseMass_cuts" ,30 , 0. , 120. );	//ElectronCollection mass
  histset[133]= fs->make<TH1D>("Event_PFMET" , "Event_PFMET" ,28 , 0. , 70. );						//ElectronCollection MET for W's
  histset[134]= fs->make<TH1D>("Event_CaloMET" , "Event_CaloMET" ,28 , 0. , 70. );					//ElectronCollection MET for W's
  histset[135]= fs->make<TH1D>("muonCorrCaloMET" , "muonCorrCaloMET" ,28 , 0. , 70. );					//ElectronCollection MET for W's
  histset[136]= fs->make<TH1D>("Event_PFMET_electroncuts" , "Event_PFMET" ,28 , 0. , 70. );				//ElectronCollection MET for W's
  histset[137]= fs->make<TH1D>("Event_PFMET_muoncuts" , "Event_PFMET" ,28 , 0. , 70. );					//ElectronCollection MET for W's
  histset[138]= fs->make<TH1D>("Electron_Z_Dcot" , "Electron Dcot" ,1000 ,-0.5 ,0.5);
  histset[139]= fs->make<TH1D>("Electron_Z_Dist" , "Electron Dist" ,1000 ,-0.5 ,0.5);





// change binning to correspond to log(0.3) - log(500), 200 bins/log10 unit
  histset[100]=fs->make<TH1D>("GM_mass_log", "GM mass log", 644, -.52, 2.7 );					//TrackCollection_GMuon mass
  histset[101]=fs->make<TH1D>("TM_mass_log", "TM mass log", 644, -.52, 2.7 );					//MuonCollection_TMuon mass
  histset[102]=fs->make<TH1D>("GM_mass_log_likecharges", "GM mass log like", 644, -.52, 2.7 );			//TrackCollection_GMuon mass
  histset[103]=fs->make<TH1D>("TM_mass_log_likecharges", "TM mass log like", 644, -.52, 2.7 );			//MuonCollection_TMuon mass


// set axis labels
    {
  histset[1]->GetXaxis()->SetTitle("Global Muon Momentum from track collection(in GeV/c)");
  histset[1]->GetYaxis()->SetTitle("Number of Events");

  histset[2]->GetXaxis()->SetTitle("Transverse Momentum of global muons from track collection(in GeV/c)");
  histset[2]->GetYaxis()->SetTitle("Number of Events");

  histset[3]->GetXaxis()->SetTitle("Eta of global muons from track collection (in radians)");
  histset[3]->GetYaxis()->SetTitle("Number of Events");
  histset[4]->GetXaxis()->SetTitle("Number of Global Muons");
  histset[4]->GetYaxis()->SetTitle("Number of Events");

  histset[5]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[5]->GetYaxis()->SetTitle("Number of Events");

  histset[6]->GetXaxis()->SetTitle("Number of Tracks");
  histset[6]->GetYaxis()->SetTitle("Number of Events");

  histset[29]->GetXaxis()->SetTitle("Momentum (in GeV/c)");
  histset[29]->GetYaxis()->SetTitle("Number of Events");

  histset[30]->GetXaxis()->SetTitle("Transverse Momentum (in GeV/c)");
  histset[30]->GetYaxis()->SetTitle("Number of Events");

  histset[31]->GetXaxis()->SetTitle("Eta (in radians)");
  histset[31]->GetYaxis()->SetTitle("Number of Events"); 

  histset[33]->GetXaxis()->SetTitle("Number of Muons");
  histset[33]->GetYaxis()->SetTitle("Number of Events");

  histset[39]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 with|eta|<1.2(in GeV/c^2)");
  histset[39]->GetYaxis()->SetTitle("Number of Events");

  histset[40]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 with |eta|>1.2 and |eta|<1.6 (in GeV/c^2)");
  histset[40]->GetYaxis()->SetTitle("Number of Events");

  histset[41]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 with |eta|>1,6 and |eta|<2.4 (in GeV/c^2)");
  histset[41]->GetYaxis()->SetTitle("Number of Events");

  histset[32]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[32]->GetYaxis()->SetTitle("Number of Events");
  histset[44]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[44]->GetYaxis()->SetTitle("Number of Events");
  histset[42]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[42]->GetYaxis()->SetTitle("Number of Events");
  histset[43]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[43]->GetYaxis()->SetTitle("Number of Events");

  histset[45]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[45]->GetYaxis()->SetTitle("Number of Events");
  histset[46]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[46]->GetYaxis()->SetTitle("Number of Events");

  histset[47]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[47]->GetYaxis()->SetTitle("Number of Events");

  histset[48]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[48]->GetYaxis()->SetTitle("Number of Events");
  histset[49]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[49]->GetYaxis()->SetTitle("Number of Events");

  histset[50]->GetXaxis()->SetTitle("Chi2 values");
  histset[50]->GetYaxis()->SetTitle("Number of Events");
  histset[51]->GetXaxis()->SetTitle("Ndof values");
  histset[51]->GetYaxis()->SetTitle("Number of Events");
  histset[52]->GetXaxis()->SetTitle("NormalizedChi2 values");
  histset[52]->GetYaxis()->SetTitle("Number of Events");

  histset[53]->GetXaxis()->SetTitle("Chi2 values");
  histset[53]->GetYaxis()->SetTitle("Number of Events");
  histset[54]->GetXaxis()->SetTitle("Ndof values");
  histset[54]->GetYaxis()->SetTitle("Number of Events");
  histset[55]->GetXaxis()->SetTitle("NormalizedChi2 values");
  histset[55]->GetYaxis()->SetTitle("Number of Events");

  histset[56]->GetXaxis()->SetTitle("Chi2 values");
  histset[56]->GetYaxis()->SetTitle("Number of Events");
  histset[57]->GetXaxis()->SetTitle("Ndof values");
  histset[57]->GetYaxis()->SetTitle("Number of Events");
  histset[58]->GetXaxis()->SetTitle("NormalizedChi2 values");
  histset[58]->GetYaxis()->SetTitle("Number of Events");

  histset[59]->GetXaxis()->SetTitle("Good Track multiplicity");
  histset[59]->GetYaxis()->SetTitle("Number of Events");

  histset[60]->GetXaxis()->SetTitle("Number of valid hits");
  histset[60]->GetYaxis()->SetTitle("Number of Events");

  histset[61]->GetXaxis()->SetTitle("Munber of pixel hits");
  histset[61]->GetYaxis()->SetTitle("Number of Events");
  histset[62]->GetXaxis()->SetTitle("Good Track multiplicity");
  histset[62]->GetYaxis()->SetTitle("Number of Events");

  histset[63]->GetXaxis()->SetTitle("Number of valid hits");
  histset[63]->GetYaxis()->SetTitle("Number of Events");

  histset[64]->GetXaxis()->SetTitle("Munber of pixel hits");
  histset[64]->GetYaxis()->SetTitle("Number of Events");

  histset[65]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[65]->GetYaxis()->SetTitle("Number of Events");

  histset[66]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[66]->GetYaxis()->SetTitle("Number of Events");

  histset[67]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[67]->GetYaxis()->SetTitle("Number of Events");
  
  histset[68]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[68]->GetYaxis()->SetTitle("Number of Events");

  histset[69]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[69]->GetYaxis()->SetTitle("Number of Events");
  histset[70]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[70]->GetYaxis()->SetTitle("Number of Events");
  histset[72]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[72]->GetYaxis()->SetTitle("Number of Events");
  histset[71]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[71]->GetYaxis()->SetTitle("Number of Events");

  histset[73]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[73]->GetYaxis()->SetTitle("Number of Events");

  histset[74]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[74]->GetYaxis()->SetTitle("Number of Events");
  histset[75]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[75]->GetYaxis()->SetTitle("Number of Events");

  histset[76]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[76]->GetYaxis()->SetTitle("Number of Events");
  histset[77]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[77]->GetYaxis()->SetTitle("Number of Events");
  histset[78]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[78]->GetYaxis()->SetTitle("Number of Events");
  histset[79]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[79]->GetYaxis()->SetTitle("Number of Events");
  histset[80]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[80]->GetYaxis()->SetTitle("Number of Events");

  histset[81]->GetXaxis()->SetTitle("Phi");
  histset[81]->GetYaxis()->SetTitle("Number of Events");
  histset[82]->GetXaxis()->SetTitle("Phi");
  histset[82]->GetYaxis()->SetTitle("Number of Events");
  histset[83]->GetXaxis()->SetTitle("Phi");
  histset[83]->GetYaxis()->SetTitle("Number of Events");
  histset[84]->GetXaxis()->SetTitle("Phi");
  histset[84]->GetYaxis()->SetTitle("Number of Events");

  histset[85]->GetXaxis()->SetTitle("Eta (in radians)");
  histset[85]->GetYaxis()->SetTitle("Number of Events"); 

  histset[86]->GetXaxis()->SetTitle("Number of Tracks");
  histset[86]->GetYaxis()->SetTitle("Number of Vertices");
  histset[87]->GetXaxis()->SetTitle("Number of Tracks");
  histset[87]->GetYaxis()->SetTitle("Number of Vertices");

  histset[88]->GetXaxis()->SetTitle("difference in x");
  histset[88]->GetYaxis()->SetTitle("Number of Events");
  histset[89]->GetXaxis()->SetTitle("difference in y");
  histset[89]->GetYaxis()->SetTitle("Number of Events");
  histset[90]->GetXaxis()->SetTitle("difference in z");
  histset[90]->GetYaxis()->SetTitle("Number of Events");
  histset[91]->GetXaxis()->SetTitle("Transverse impact paramemter w.r.t. BS");
  histset[91]->GetYaxis()->SetTitle("Number of Events");
  histset[92]->GetXaxis()->SetTitle("Number of valid hits");
  histset[92]->GetYaxis()->SetTitle("Number of Events");
  histset[93]->GetXaxis()->SetTitle("Muons matching with at least two muon stations");
  histset[93]->GetYaxis()->SetTitle("Number of Events");

  histset[94]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[94]->GetYaxis()->SetTitle("Number of Events");
  histset[95]->GetXaxis()->SetTitle("Invariant Mass for Nmuon>=2 (in GeV/c^2)");
  histset[95]->GetYaxis()->SetTitle("Number of Events");

  histset[96]->GetXaxis()->SetTitle("Momentum (in GeV/c)");
  histset[96]->GetYaxis()->SetTitle("Number of Events");
  histset[97]->GetXaxis()->SetTitle("Transverse Momentum (in GeV/c)");
  histset[97]->GetYaxis()->SetTitle("Number of Events");
  histset[98]->GetXaxis()->SetTitle("Eta (in radians)");
  histset[98]->GetYaxis()->SetTitle("Number of Events");
  histset[99]->GetXaxis()->SetTitle("Eta (in radians)");
  histset[99]->GetYaxis()->SetTitle("Number of Events"); 
  histset[104]->GetXaxis()->SetTitle("Invariant Mass for Nelectron>=2 (in GeV/c^2)");
  histset[104]->GetYaxis()->SetTitle("Number of Events");
  histset[105]->GetXaxis()->SetTitle("Number of Electrons");
  histset[105]->GetYaxis()->SetTitle("Number of Events");

  histset[106]->GetXaxis()->SetTitle("Super Cluster Eta(in radians)");
  histset[106]->GetYaxis()->SetTitle("Number of Events");
  histset[107]->GetXaxis()->SetTitle("Super Cluster Energy");
  histset[107]->GetYaxis()->SetTitle("Number of Events");
  histset[108]->GetXaxis()->SetTitle("Track Isolation");
  histset[108]->GetYaxis()->SetTitle("Number of Events");
  histset[109]->GetXaxis()->SetTitle("Ecal Isolation");
  histset[109]->GetYaxis()->SetTitle("Number of Events");
  histset[110]->GetXaxis()->SetTitle("Hcal Isolation");
  histset[110]->GetYaxis()->SetTitle("Number of Events");
  histset[111]->GetXaxis()->SetTitle("gsfTrack Hit type");
  histset[111]->GetYaxis()->SetTitle("Number of Events");
  histset[112]->GetXaxis()->SetTitle("Dcot");
  histset[112]->GetYaxis()->SetTitle("Number of Events");
  histset[113]->GetXaxis()->SetTitle("Dist");
  histset[113]->GetYaxis()->SetTitle("Number of Events");
  histset[114]->GetXaxis()->SetTitle("sigmaIetaIeta");
  histset[114]->GetYaxis()->SetTitle("Number of Events");
  histset[115]->GetXaxis()->SetTitle("DeltaPhi");
  histset[115]->GetYaxis()->SetTitle("Number of Events");
  histset[116]->GetXaxis()->SetTitle("DeltaEta");
  histset[116]->GetYaxis()->SetTitle("Number of Events");
  histset[117]->GetXaxis()->SetTitle("H/E");
  histset[117]->GetYaxis()->SetTitle("Number of Events");
  histset[118]->GetXaxis()->SetTitle("Electron Et");
  histset[118]->GetYaxis()->SetTitle("Number of Events");
  histset[119]->GetXaxis()->SetTitle("Invariant Mass for Nelectron>=2 (in GeV/c^2)");
  histset[119]->GetYaxis()->SetTitle("Number of Events");
  histset[120]->GetXaxis()->SetTitle("Track Isolation");
  histset[120]->GetYaxis()->SetTitle("Number of Events");
  histset[121]->GetXaxis()->SetTitle("Ecal Isolation");
  histset[121]->GetYaxis()->SetTitle("Number of Events");
  histset[122]->GetXaxis()->SetTitle("Hcal Isolation");
  histset[122]->GetYaxis()->SetTitle("Number of Events");

  histset[123]->GetXaxis()->SetTitle("sigmaIetaIeta");
  histset[123]->GetYaxis()->SetTitle("Number of Events");
  histset[124]->GetXaxis()->SetTitle("DeltaPhi");
  histset[124]->GetYaxis()->SetTitle("Number of Events");
  histset[125]->GetXaxis()->SetTitle("DeltaEta");
  histset[125]->GetYaxis()->SetTitle("Number of Events");
  histset[126]->GetXaxis()->SetTitle("H/E");
  histset[126]->GetYaxis()->SetTitle("Number of Events");
  histset[127]->GetXaxis()->SetTitle("Invariant Mass for Nelectron>=2 (in GeV/c^2)");
  histset[127]->GetYaxis()->SetTitle("Number of Events");
  histset[128]->GetXaxis()->SetTitle("Invariant Mass for Nelectron>=2 (in GeV/c^2)");
  histset[128]->GetYaxis()->SetTitle("Number of Events");

  histset[129]->GetXaxis()->SetTitle("Transverse Mass for muons+MET (in GeV/c^2)");
  histset[129]->GetYaxis()->SetTitle("Number of Events");
  histset[130]->GetXaxis()->SetTitle("Transverse Mass for muons+MET (in GeV/c^2)");
  histset[130]->GetYaxis()->SetTitle("Number of Events");
  histset[131]->GetXaxis()->SetTitle("Transverse Mass for electrons+MET (in GeV/c^2)");
  histset[131]->GetYaxis()->SetTitle("Number of Events");
  histset[132]->GetXaxis()->SetTitle("Transverse Mass for electrons+MET (in GeV/c^2)");
  histset[132]->GetYaxis()->SetTitle("Number of Events");
  histset[133]->GetXaxis()->SetTitle("MET (in GeV)");
  histset[133]->GetYaxis()->SetTitle("Number of Events");
  histset[134]->GetXaxis()->SetTitle("MET (in GeV)");
  histset[134]->GetYaxis()->SetTitle("Number of Events");
  histset[135]->GetXaxis()->SetTitle("MET (in GeV)");
  histset[135]->GetYaxis()->SetTitle("Number of Events");
  histset[136]->GetXaxis()->SetTitle("MET (in GeV)");
  histset[136]->GetYaxis()->SetTitle("Number of Events");
  histset[137]->GetXaxis()->SetTitle("MET (in GeV)");
  histset[137]->GetYaxis()->SetTitle("Number of Events");





  histset[100]->GetXaxis()->SetTitle("Invariant Log10(Mass) for Nmuon>=2 (in GeV/c^2)");
  histset[100]->GetYaxis()->SetTitle("Number of Events/GeV");
  histset[101]->GetXaxis()->SetTitle("Invariant Log10(Mass) for Nmuon>=2 (in GeV/c^2)");
  histset[101]->GetYaxis()->SetTitle("Number of Events/GeV");
  histset[102]->GetXaxis()->SetTitle("Invariant Log10(Mass) for Nmuon>=2 (in GeV/c^2)");
  histset[102]->GetYaxis()->SetTitle("Number of Events/GeV");
  histset[103]->GetXaxis()->SetTitle("Invariant Log10(Mass) for Nmuon>=2 (in GeV/c^2)");
  histset[103]->GetYaxis()->SetTitle("Number of Events/GeV");

    }


//                event property histograms                //

// histograms to check and simulate JSON
  histset[200]=fs->make<TH1D>("Run number", "Run number", 3100, 146400, 149500);
  histset[201]=fs->make<TH1D>("Run number, no JSON", "Run number, no JSON, should be empty!", 3100, 146400, 149500);
  histset[202]=fs->make<TH1D>("Event number", "Event number", 2000, 0, 2000000000);
  histset[203]=fs->make<TH1D>("Lumi section", "Lumi section", 300, 0, 3000);
  histset[204]=fs->make<TH1D>("Run number, JSON", "Run number, JSON, should be filled", 3100, 146400, 149500);
} //Fold for I's Histograms
    
                //// ---------- Check InvMass Function ---------- ////
    
                    histset[207]  =   fs->make<TH1D>("h_Jspimass" , "" ,90 ,2.6 , 3.5 );
    
                //// ---------- J/Psi Acceptance ---------- ////
    
                    hxhy[1]       =   fs->make<TH2D>("h_Paircount", "", 100,0, 2.4, 100, 0,30);
                    histset[205]  =   fs->make<TH1D>("h_JpsiPT", "", 100, 0, 30);
                    hxhy[2]       =   fs->make<TH2D>("h_Paircount2", "", 24, 0, 2.4, 42, binsy );
                //// ------------- D* Anaylsis ------------ ////
    
                    histset[206]  =   fs->make<TH1D>("h_D0mass", "", 200, 1.5, 2); //D0 mass histogram
                    histset[208]  =   fs->make<TH1D>("h_deltaM", "", 100, 0.138, 0.17); //Mass difference histogram for 'right charge' paring
                    histset[209]  =   fs->make<TH1D>("h_deltaMwrongcharge", "", 100, 0.14, 0.17); //Mass difference histogram for 'wrong charge' pairing
                    histset[210]  =   fs->make<TH1D>("h_z","",100,0,0.5);
                    histset[211]  =   fs->make<TH1D>("h_n","",100, 0 ,200);
    
                //// ------------- Y cut implementation ------------ ////
    
                histset[250]  =  fs->make<TH1D>("h_Upsilon",80, 8, 12);
                histset[251]  =  fs->make<TH1D>("h_Upsilon_eta",80, 8, 12);
    
    
}


DemoAnalyzer::~DemoAnalyzer()
{
 
   // do anything here that needs to be done at destruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

inline double invMass(double p1, double p1x, double p1y, double p1z, double m1 , double p2, double p2x, double p2y, double p2z, double m2);

// ------------ method called for each event  ------------
void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;
   using namespace reco;
   using namespace std;
   //using reco::TrackCollection;
   //using reco::VertexCollection;
   //using reco::MuonCollection;

Handle<reco::BeamSpot> beamSpotHandle;
Handle<reco::TrackCollection> tracks;
Handle<reco::TrackCollection> gmuons;
Handle<reco::VertexCollection> Primvertex;
Handle<reco::MuonCollection> muons;
Handle<reco::GsfElectronCollection> electrons;
 edm:: Handle<edm::View<reco::PFMET> > pfmets;// Handle<reco::PFMET> mets;  didn't work for some reason
edm:: Handle<edm::View<reco::CaloMET> > calomets;
edm:: Handle<edm::View<reco::CaloMET> > muCorrmets;

  #ifdef THIS_IS_AN_EVENT_EXAMPLE 
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
   #endif

   // should we use this?   
   #ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
   #endif
 

  
   //                   fill basic event information                  //
   //    *** still under very preliminary construction   ****
    // JSON Checkers
   //LogInfo("Demo")<<"Event number"<<iEvent.id()<<" Run number"<<iEvent.run();
    {
 histset[200]->Fill(iEvent.run());
 histset[202]->Fill((iEvent.id()).event());
 histset[203]->Fill(iEvent.luminosityBlock());
 // check JSON "by hand"; according to documentation, JSON provides a
 // "list of good lumi sections"
 // if badJSON nonzero event should not be used for analysis
 int badJSON = 0;
 if (iEvent.run() == 146428) {
  if ((iEvent.luminosityBlock() > 2 && iEvent.luminosityBlock() < 11) ||
      (iEvent.luminosityBlock() > 52 && iEvent.luminosityBlock() < 55) ||
      (iEvent.luminosityBlock() > 92))
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146430) {
  if ((iEvent.luminosityBlock() > 12 && iEvent.luminosityBlock() < 18) ||
      (iEvent.luminosityBlock() > 47 && iEvent.luminosityBlock() < 50) ||
      (iEvent.luminosityBlock() > 62 && iEvent.luminosityBlock() < 65) ||
      (iEvent.luminosityBlock() > 90))
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146431) {
  if (iEvent.luminosityBlock() > 23)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146436) {
  if (iEvent.luminosityBlock() > 532)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146437) {
  if (iEvent.luminosityBlock() > 798)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146510) {
  if (iEvent.luminosityBlock() > 510)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146511) {
  if (iEvent.luminosityBlock() == 64 || iEvent.luminosityBlock() > 778)
    badJSON = iEvent.run(); 
 }
else if (iEvent.run() == 146513) {
  if (iEvent.luminosityBlock() == 2 || iEvent.luminosityBlock() > 15)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146514) {
  if (iEvent.luminosityBlock() == 546 || iEvent.luminosityBlock() == 547 || iEvent.luminosityBlock() > 871)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146589) {
  if (iEvent.luminosityBlock() < 34 || iEvent.luminosityBlock() > 248)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146644) {
   if (iEvent.luminosityBlock() < 89 || iEvent.luminosityBlock() == 118 || iEvent.luminosityBlock() == 566 || iEvent.luminosityBlock() == 868 || iEvent.luminosityBlock() == 1033 || (iEvent.luminosityBlock() > 2171 && iEvent.luminosityBlock() < 2369) || iEvent.luminosityBlock() > 2465)
   badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146698) {
   if (iEvent.luminosityBlock() < 155 || (iEvent.luminosityBlock() > 180 && iEvent.luminosityBlock() < 186) || iEvent.luminosityBlock() > 189)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146710) {
   if (iEvent.luminosityBlock() == 116 || (iEvent.luminosityBlock() > 49 && iEvent.luminosityBlock() < 114) || iEvent.luminosityBlock() > 214)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146712) {
   if (iEvent.luminosityBlock() > 69)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146713) {
   if (iEvent.luminosityBlock() == 49 || iEvent.luminosityBlock() > 256)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146715) {
   if (iEvent.luminosityBlock() > 125)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146721) {
   if (iEvent.luminosityBlock() > 6)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146724) {
   if (iEvent.luminosityBlock() == 107 || iEvent.luminosityBlock() == 108 || (iEvent.luminosityBlock() > 109 && iEvent.luminosityBlock() < 112) || iEvent.luminosityBlock() == 151 || iEvent.luminosityBlock() > 159)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146804) {
   if (iEvent.luminosityBlock() < 111 || iEvent.luminosityBlock() == 149 || iEvent.luminosityBlock() == 521 || iEvent.luminosityBlock() == 790 || iEvent.luminosityBlock() == 823 || iEvent.luminosityBlock() > 905)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146807) {
   if (iEvent.luminosityBlock() < 132 || iEvent.luminosityBlock() == 363 || iEvent.luminosityBlock() == 364 || iEvent.luminosityBlock() == 421 || iEvent.luminosityBlock() > 469)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 146944) {
   if (iEvent.luminosityBlock() < 177 || iEvent.luminosityBlock() > 669)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147043) {
   if (iEvent.luminosityBlock() < 161 || iEvent.luminosityBlock() > 500)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147048) {
   if (iEvent.luminosityBlock() < 94 || iEvent.luminosityBlock() > 484)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147114) {
   if (iEvent.luminosityBlock() < 180 || (iEvent.luminosityBlock() > 187 && iEvent.luminosityBlock() < 227) || (iEvent.luminosityBlock() > 240 && iEvent.luminosityBlock() < 247) || iEvent.luminosityBlock() > 667)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147115) {
   if (iEvent.luminosityBlock() > 546)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147116) {
   if (iEvent.luminosityBlock() > 54)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147196) {
   if (iEvent.luminosityBlock() > 90)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147214) {
   if (iEvent.luminosityBlock() > 79)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147216) {
   if (iEvent.luminosityBlock() > 63)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147217) {
   if (iEvent.luminosityBlock() > 193)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147218) {
   if (iEvent.luminosityBlock() > 45)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147219) {
   if ((iEvent.luminosityBlock() > 293 && iEvent.luminosityBlock() < 309) || iEvent.luminosityBlock() > 320)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147222) {
   if (iEvent.luminosityBlock() > 444)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147284) {
   if (iEvent.luminosityBlock() < 32 || iEvent.luminosityBlock() > 306)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147390) {
   if (iEvent.luminosityBlock() == 479 || iEvent.luminosityBlock() > 837)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147450) {
   if (iEvent.luminosityBlock() < 80 || iEvent.luminosityBlock() > 166)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147451) {
   if (iEvent.luminosityBlock() == 117 || iEvent.luminosityBlock() > 129)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147452) {
   if (iEvent.luminosityBlock() > 44)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147453) {
   if (iEvent.luminosityBlock() > 146)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147454) {
   if (iEvent.luminosityBlock() > 97)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147754) {
   if (iEvent.luminosityBlock() == 168 || iEvent.luminosityBlock() == 169 || iEvent.luminosityBlock() > 377)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147755) {
   if (iEvent.luminosityBlock() < 81 || iEvent.luminosityBlock() > 231)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147757) {
   if (iEvent.luminosityBlock() > 359)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147926) {
   if (iEvent.luminosityBlock() < 77 || iEvent.luminosityBlock() > 548)
   badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147927) {
   if (iEvent.luminosityBlock() > 152)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 147929) {
   if ((iEvent.luminosityBlock() > 266 && iEvent.luminosityBlock() < 272) || iEvent.luminosityBlock() == 619 || iEvent.luminosityBlock() > 643)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148002) {
   if (iEvent.luminosityBlock() < 92 || iEvent.luminosityBlock() > 203)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148029) {
   if (iEvent.luminosityBlock() < 50 || iEvent.luminosityBlock() == 484 || iEvent.luminosityBlock() == 570 || iEvent.luminosityBlock() > 571)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148031) {
   if ((iEvent.luminosityBlock() > 341 && iEvent.luminosityBlock() < 472) || iEvent.luminosityBlock() == 758 || iEvent.luminosityBlock() > 855)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148032) {
  if (iEvent.luminosityBlock() > 199)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148058) {
   if (iEvent.luminosityBlock() > 97)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148822) {
   if (iEvent.luminosityBlock() > 446)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148829) {
   if ((iEvent.luminosityBlock() > 240 && iEvent.luminosityBlock() < 244) || iEvent.luminosityBlock() == 74 || iEvent.luminosityBlock() > 303)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148860) {
   if (iEvent.luminosityBlock() > 39)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148862) {
   if (iEvent.luminosityBlock() == 19 || iEvent.luminosityBlock() == 109 || iEvent.luminosityBlock() == 150 || (iEvent.luminosityBlock() > 165 && iEvent.luminosityBlock() < 224) || (iEvent.luminosityBlock() > 258 && iEvent.luminosityBlock() < 262) || iEvent.luminosityBlock() == 298 || iEvent.luminosityBlock() == 367 || (iEvent.luminosityBlock() > 504 && iEvent.luminosityBlock() < 512) || iEvent.luminosityBlock() > 679)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148864) {
   if (iEvent.luminosityBlock() == 32 || (iEvent.luminosityBlock() > 141 && iEvent.luminosityBlock() < 224) || iEvent.luminosityBlock() == 237 || iEvent.luminosityBlock() == 477 || iEvent.luminosityBlock() > 680)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148952) {
   if (iEvent.luminosityBlock() < 70 || iEvent.luminosityBlock() > 257)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 148953) {
   if (iEvent.luminosityBlock() > 100)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149003) {
   if (iEvent.luminosityBlock() < 84 || iEvent.luminosityBlock() > 238)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149011) {
   if (iEvent.luminosityBlock() == 342 || iEvent.luminosityBlock() > 706)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149058) {
   if (iEvent.luminosityBlock() > 65)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149063) {
   if (iEvent.luminosityBlock() > 102)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149181) {
   if (iEvent.luminosityBlock() < 229 || (iEvent.luminosityBlock() > 1840 && iEvent.luminosityBlock() < 1844) || iEvent.luminosityBlock() > 1920)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149291) {
   if (iEvent.luminosityBlock() < 79 || (iEvent.luminosityBlock() > 79 && iEvent.luminosityBlock() < 82) || iEvent.luminosityBlock() == 787 || iEvent.luminosityBlock() == 789 || (iEvent.luminosityBlock() > 790 && iEvent.luminosityBlock() < 794) || iEvent.luminosityBlock() > 794)
    badJSON = iEvent.run(); 
 }
 else if (iEvent.run() == 149294) {
  if (iEvent.luminosityBlock() > 171)
    badJSON = iEvent.run(); 
 }
 else {
   badJSON = iEvent.run(); 
 }

// signal presence of bad JSON event
if (badJSON > 0) {
  LogInfo("Demo")<<"bad JSON Event id"<<iEvent.id();
  histset[201]->Fill(iEvent.run()); 
  // find out how to exit here ...
 }
else {
  histset[204]->Fill(iEvent.run());
 }
}
   // LogInfo("Demo")<<"Event number"<<(iEvent.id()).event()<<" Run number"<<iEvent.run()<<"  Lumi:"<<iEvent.luminosityBlock();


 //                 event is to be kept                      //
 //                 load relevant event information          //

 iEvent.getByLabel("muons", muons);
 iEvent.getByLabel("globalMuons", gmuons);
 iEvent.getByLabel("generalTracks", tracks); 
 iEvent.getByLabel("offlineBeamSpot", beamSpotHandle);
 iEvent.getByLabel("gsfElectrons",electrons);
 iEvent.getByLabel("pfMet",pfmets);
 iEvent.getByLabel("met",calomets);
 iEvent.getByLabel("corMetGlobalMuons",muCorrmets);

  histset[133]-> Fill((pfmets->front()).et());
 histset[134]-> Fill((calomets->front()).et());
 histset[135]-> Fill((muCorrmets->front()).et());
 reco::BeamSpot beamSpot;
 if ( beamSpotHandle.isValid() )
{
    beamSpot = *beamSpotHandle;

} else
{
    edm::LogInfo("Demo")
      << "No beam spot available from EventSetup \n";
}

 // choose primary verices with/without beam spot
 
 iEvent.getByLabel("offlinePrimaryVerticesWithBS",Primvertex);
 //iEvent.getByLabel("offlinePrimaryVertices",Primvertex);
 
 // fill basic track information (all tracks in event)

 histset[6]->Fill(tracks->size());
 histset[26]->Fill(Primvertex->size());

 // loop over tracks

  for( reco::TrackCollection::const_iterator it = tracks->begin(); it !=
 tracks->end(); it++) {
    // LogInfo("Demo")<<"track p"<<it->p()<< "  track reference position"<<it->referencePoint()<< "   track vertex position"<<it->vertex();
    histset[34]->Fill(it->p());
    histset[35]->Fill(it->pt());
    histset[38]->Fill(it->eta());
    histset[81]->Fill(it->phi());

      

  //  if (it->eta() < 0) {hxhy[0]->Fill((-1)*(it->eta()), it->pt());}
//	else {hxhy[0]->Fill(it->eta(), it->pt());}

  
   // histset[23]->Fill(it->vx());
   //histset[24]->Fill(it->vy());
   // histset[25]->Fill(it->vz());

   // fill histograms for QCD-10-006

   if(fabs(it->eta())<=0.200){
     
     histset[9]->Fill(it->pt());}

if(fabs(it->eta())>0.200 && fabs(it->eta())<0.400){
   
     histset[10]->Fill(it->pt());}
 if(fabs(it->eta())>0.400 && fabs(it->eta())<=0.600){
     
     histset[11]->Fill(it->pt());}
if(fabs(it->eta())<0.600 && fabs(it->eta())<=0.800){
 
     histset[12]->Fill(it->pt());}
if(fabs(it->eta())>0.800 && fabs(it->eta())<=1.00){
  
     histset[13]->Fill(it->pt());}
if(fabs(it->eta())>1.00 && fabs(it->eta())<=1.20){
   
     histset[14]->Fill(it->pt());}
if(fabs(it->eta())>1.200 && fabs(it->eta())<=1.40){
                                      
     histset[15]->Fill(it->pt());}
if(fabs(it->eta())>1.400 && fabs(it->eta())<=1.60){
    
     histset[16]->Fill(it->pt());}
if(fabs(it->eta())>1.600 && fabs(it->eta())<=1.80){
  
     histset[17]->Fill(it->pt());}
if(fabs(it->eta())>1.800 && fabs(it->eta())<=2.00){
  
     histset[18]->Fill(it->pt());}
if(fabs(it->eta())>2.000 && fabs(it->eta())<=2.20){

     histset[19]->Fill(it->pt());}
if(fabs(it->eta())>2.200 && fabs(it->eta())<=2.40){
    
     histset[20]->Fill(it->pt());}
if(fabs(it->eta())>2.400 && fabs(it->eta())<=2.60){

     histset[21]->Fill(it->pt());}

  }//end of track collection loop


  // fill basic vertex information

  int size; // variable to store the size of vertices in an event
  size=Primvertex->size();

  // loop over primary vertices

  for(reco::VertexCollection::const_iterator ite = Primvertex->begin(); ite !=Primvertex->end(); ite++) {
    
  int Muflag=0; //  flag indicating whether a muon candidate is associated to this vertex

  // store z position of first vertex 
  if(ite== Primvertex->begin())firstz=ite->z();
  
  // loop over muon collection to check for muon tracks associated to this vertex, i.e. track refernce position close to the vertex
  for (MuonCollection::const_iterator itMuon = muons->begin(); itMuon != muons->end(); ++itMuon) {
    if((itMuon->track()).isNonnull()){
if(fabs(ite->x()-(itMuon->track())->vx())<=0.1 && fabs(ite->y()-(itMuon->track())->vy())<=0.1 && fabs(ite->z()-(itMuon->track())->vz())<=1.0)
 Muflag=1; // The vertex has a muon associated to it

    }
  }  // end of loop over muon collection

    if(ite->tracksSize()==0){ //This is done so that there are no vertices with zero tracks stored in the histogram.
      size=size-1; }

    /// store distance of vertex to first vertex
 if(ite->z()!=firstz)
    {
      deltaz=fabs(firstz-(ite->z()));
      histset[22]->Fill(deltaz);
      }

    if(ite->tracksSize()!=0 && Muflag==0)

// only store the coordinates of those vertices which are valid and have tracks && have no muons (= "next-to-miminum bias")
// very first attempt in context of FWD-11-001

      { histset[36]->Fill(ite->x());
        histset[37]->Fill(ite->y());
        histset[7]->Fill(ite->z());
        histset[87]->Fill(ite->tracksSize());

              //-----to get track details for valid vertices-------------//

    //  LogInfo("Demo")<<"Vertex position"<<ite->position()<<endl;

    //  loop over all tracks from this vertex

    int counter=0;   
 for(reco::Vertex::trackRef_iterator iTrack  =ite->tracks_begin(); iTrack != ite->tracks_end();++iTrack) { 
// for dealing with tracks associated with a particular vertex

   // why *iTrack ??
      histset[27]->Fill((*iTrack)->pt()); // vertex.track pt
      histset[28]->Fill((*iTrack)->p()); // vertex.track p
      histset[85]->Fill((*iTrack)->eta());
      histset[82]->Fill((*iTrack)->phi());

      // LogInfo("Demo")<<"Vertex track p"<<(*iTrack)->p()<<"  Vertex trackposition"<<(*iTrack)->referencePoint();

      counter++;
 }
 // fill number of tracks associated to this vertex (is it Ok to fill int?)
 histset[86]->Fill(counter);
    

      }// end of if(ite->tracksSize()!=0 && Muflag==0)
  }//end of vertex collection loop

  
  histset[8]->Fill(size); //only valid vertices get stored.
 

   //------------------analysing global muons with the track collection-------------------------// 


   histset[4]->Fill(gmuons->size());
   
  
  for( reco::TrackCollection::const_iterator it= gmuons->begin(); it !=gmuons->end(); it++) {
   
   histset[1]->Fill(it->p());
   histset[2]->Fill(it->pt());
   histset[3]->Fill(it->eta());
   histset[53]->Fill(it->chi2());
   histset[54]->Fill(it->ndof());
   histset[55]->Fill(it->normalizedChi2());
   histset[83]->Fill(it->phi());
 // LogInfo("Demo")  <<"global  muon track pointer "<<it;
   // LogInfo("Demo")<<"global muon track p"<<it->p()<<"  global muon track pos"<<it->referencePoint()<<" global muon track vertex"<<it->vertex();

       //-----------------now prepare quality cuts---------------------//
   
   int size1=0;
   int ValidHits=0, PixelHits=0;
       const reco::HitPattern& p = it->hitPattern();
      
       // loop over the hits of the track
       for (int i=0; i<p.numberOfHits(); i++) {
	 uint32_t hit = p.getHitPattern(i);
	 
	 // if the hit is valid and in pixel 
	 if (p.validHitFilter(hit) && p.pixelHitFilter(hit)) PixelHits++;
	 
	 if (p.validHitFilter(hit))ValidHits++;
	 
       }
       if(ValidHits>=10 && PixelHits>=1) size1++; // For Mu sample
	 
	   
       if(size1)histset[59]->Fill(size1);// Fill hitsOK histo
       histset[60]->Fill(ValidHits);
       histset[61]->Fill(PixelHits);

      
       // loop over muons satisfying quality cuts //
       // as applied for plots used in presentation at CMS physics meetring //

       if(gmuons->size()>=2 && ValidHits>=12 && PixelHits>=2 && it->normalizedChi2()<4.0)
	 //       if(gmuons->size()>=2)
     {

       reco::TrackCollection::const_iterator i=it;
       i++;
       for(;i!=gmuons->end();i++){
	 
	 int ValidHits1=0, PixelHits1=0;
	 const reco::HitPattern& p1 = i->hitPattern();
	
       // loop over the hits of the track
	 for (int n=0; n<p1.numberOfHits(); n++) {
	 uint32_t hit = p1.getHitPattern(n);

	 // if the hit is valid and in pixel 
	 if (p1.validHitFilter(hit) && p1.pixelHitFilter(hit)) PixelHits1++;
	 
	 if (p1.validHitFilter(hit))ValidHits1++;
	 
	 }
	 

	 if(it->charge()==-(i->charge())) // unlike charges
	 {
	   //----------to calculate invariant mass-----------------//

	   s1=sqrt(((it->p())*(it->p())+sqm1)*((i->p())*(i->p())+sqm1));
	   s2=it->px()*i->px()+it->py()*i->py()+it->pz()*i->pz();
	   s=sqrt(2.0*(sqm1+(s1-s2)));
	   histset[5]->Fill(s);
	   histset[44]->Fill(s);
	   histset[80]->Fill(s);
	   // apply weight 200/(ln10*m/GeV) to convert to events/GeV (see sheet)
	   w=200/log(10)/s;
	   histset[100]->Fill(log10(s),w);
	    
           // scorr needs to be fixed ... use average pseudorapidity for approximate correction, should use rapidity

	   scorr=(1.0+0.00038+0.0003*(i->eta()+it->eta())*(i->eta()+it->eta())/4)*s;
		  histset[66]->Fill(scorr);
		  histset[67]->Fill(scorr);

	   //----------to calculate rapidity---------------------//

	     s3=(it->p())*(it->p())+(i->p())*(i->p());
		  s4= sqrt(s3+2*s2+s*s);
		  pz=it->pz()+i->pz();
		  rap=0.5*log((s4+pz)/(s4-pz));
		  if(fabs(rap)<1.2) {histset[47]->Fill(s);histset[68]->Fill(scorr);}
		  else if(fabs(rap)<1.6){ histset[48]->Fill(s);histset[69]->Fill(scorr);}
		  else if(fabs(rap)<2.4){ histset[49]->Fill(s);histset[70]->Fill(scorr);}

		   //--------------------------upsilon ranges ----------------------------------------------------------//
		  
		  if(fabs(it->eta())<2.4 && fabs(i->eta())<2.4){ histset[77]->Fill(s);histset[79]->Fill(scorr);}
		  if(fabs(it->eta())<1. && fabs(i->eta())<1.) {  histset[76]->Fill(s);histset[78]->Fill(scorr);}
	 }

	 if(it->charge()==(i->charge()))   // like charges invariant mass
	 {
	   s1=sqrt(((it->p())*(it->p())+sqm1)*((i->p())*(i->p())+sqm1));
	   s2=it->px()*i->px()+it->py()*i->py()+it->pz()*i->pz();
	   s=sqrt(2.0*(sqm1+(s1-s2)));
	   histset[46]->Fill(s);
	   w=200/log(10)/s;
	   histset[102]->Fill(log10(s),w);
	  
	 }


       }//end of for(;i!=gmuons....)

     }//end of if(gmuons>=2.....)
  }//end of reco ::TrackCollection
  


  //--------------------------analysing muons with MuonCollection----------------------------------------//
  //                          for J/psi (BPH-10-002) and Z (EWK-10-002)                                  //
 
 
  histset[33]->Fill(muons->size()); // only the total number of muons
 
for (MuonCollection::const_iterator itMuon = muons->begin(); itMuon != muons->end(); ++itMuon) {

  double Mt;
 int  ValidHits=0; int PixelHits=0;//TM_J/Psi
 int  GM_ValidHits=0; int GM_PixelHits=0;//GM_Z
 int size2=0;
 int jpsi_flag1=0;
 
 int Zflag1=0;
 
 int Wflag1=0;
 int Wflag2=1;
 int goodhit=0;//GM_Z Muon Chamber hits
 double relIso=1.;// relative Isolation for muons(Z paper);
 
    
 int Yflag=0;

 math::XYZPoint point(beamSpot.position());

  // LogInfo("Demo")  <<" muon p "<<itMuon->momentum()<<" p "<<itMuon->p()<<" pt "<<itMuon->pt();
  // LogInfo("Demo") <<" global "<<itMuon->isGlobalMuon()<<" tracker "<<itMuon->isTrackerMuon();

  histset[29]->Fill(itMuon->p());histset[30]->Fill(itMuon->pt());histset[31]->Fill(itMuon->eta());  histset[84]->Fill(itMuon->phi());
  
 
  if((itMuon->track()).isNonnull()){// some muons might not have valid track references

      histset[56]->Fill(((itMuon->track()))->chi2());
      histset[57]->Fill(((itMuon->track()))->ndof());
      histset[58]->Fill(((itMuon->track()))->normalizedChi2());
    
// LogInfo("Demo")  <<" muon track p "<<(itMuon->track())->p()<<"  muon track position"<<(itMuon->track())->referencePoint()<<"  Muon Vertex Position"<<itMuon->vertex();


	  //-----------------------------------------------------------------//
       //-----------------now putting the cuts for j/psi---------------------//
	  //--------------------------------------------------------------//

	  const reco::HitPattern& p = (itMuon->track())->hitPattern();
	  

       // loop over the hits of the track
       for (int i=0; i<p.numberOfHits(); i++) {
	 uint32_t hit = p.getHitPattern(i);
	  
	
	 // if the hit is valid and in pixel
	 if (p.validHitFilter(hit) && p.pixelHitFilter(hit)){
	   PixelHits++;
           }
	 if (p.validHitFilter(hit)){
	   ValidHits++;
	 }
       }
       if(ValidHits>=12 && PixelHits>=2){
	 size2++; }
       histset[63]->Fill(ValidHits);
       histset[64]->Fill(PixelHits);
       if(size2)histset[62]->Fill(size2);// Fill hitsOK histo


	if (ValidHits>=12 && PixelHits>=2 && (itMuon->track())->normalizedChi2()<4.0 && muon::isGoodMuon(*itMuon,muon::TMLastStationAngTight) && muon::isGoodMuon(*itMuon,muon::TrackerMuonArbitrated) && itMuon->isTrackerMuon())
   {
     if(fabs(itMuon->eta())<1.3 && itMuon->pt()>3.3)jpsi_flag1=1;
     else if(fabs(itMuon->eta())>1.3 && fabs(itMuon->eta())<2.2 && itMuon->p()>2.9)jpsi_flag1=1;
     else if(fabs(itMuon->eta())>2.2  && fabs(itMuon->eta())<2.4 && itMuon->pt()>0.8)jpsi_flag1=1;
   }
 

		//--------------------------------------------------------------------------------//
		//----------------------------Z & W cuts----------------------------------------------//
		//--------------------------------------------------------------------------------//

       if(itMuon->isGlobalMuon() && (itMuon->globalTrack()).isNonnull())
	 {
	   //relative isolation
	  relIso=((itMuon->isolationR03()).sumPt+(itMuon->isolationR03()).emEt+(itMuon->isolationR03()).hadEt)/itMuon->pt();



	  // checking hit pattern info
	  const reco::HitPattern& GM_p = (itMuon->globalTrack())->hitPattern();
	  goodhit=GM_p.numberOfValidMuonHits();
	  for (int i=0; i<GM_p.numberOfHits(); i++) {
	  uint32_t hit = GM_p.getHitPattern(i); 
	 // if the hit is valid and in pixel
	 if (GM_p.validHitFilter(hit) && GM_p.pixelHitFilter(hit)){
	   GM_PixelHits++;
           }
	 if (GM_p.validHitFilter(hit)){
	   GM_ValidHits++; 
	  }
	  }// for (int i=0; i<GM_p.numberOfHits(); i++) ends
           histset[92]->Fill(goodhit);


	   bool result= muon::isGoodMuon(*itMuon, muon::TMLastStationTight);
           histset[93]->Fill(result);
	   histset[91]->Fill((itMuon->globalTrack())->dxy(point)); // Transverse impact parameter w.r.t. BeamSpot


	   if(GM_ValidHits>=10 && GM_PixelHits>=1 && goodhit>=1 && (itMuon->globalTrack())->normalizedChi2()<10.0 && fabs((itMuon->globalTrack())->dxy(point))<0.2 && result && relIso<0.15)
	     {
	       if(itMuon->pt()>20 && fabs(itMuon->eta())<2.1)
		 {
		   Zflag1=1;
		   Wflag1=1;
		 }
	      
	     }
	 }//if(itMuon->isGlobalMuon()..........) ends

      }//isNonnull() ends
 
     

 
  //--------------------------look for dimuons----------------------------------------------------------------------------//
  //------------------------------------------------------------------------------------------------------//
  if(muons->size()>=2)
    {
      
      if(itMuon->isTrackerMuon())
	{
	  //---------------only for tracker muons--------------//
          //   fill for first muon of multimuon candidate      //
          //   *** logic to be reconsidered ***  //        

  histset[50]->Fill(((itMuon->track()))->chi2());
  histset[51]->Fill(((itMuon->track()))->ndof());
  histset[52]->Fill(((itMuon->track()))->normalizedChi2());
	}

  //----------------------------------2nd tracker muon loop ---------------------------------------------//

//-------------------------------------------------------------------//
	  MuonCollection::const_iterator itM=itMuon;
	  itM++;

	  for(;itM!=muons->end();itM++)
	    {
	     int jpsi_flag2=0;
	     int Zflag2=0;
	     Wflag2=1;// to reject events with 2 muons
	     int goodhit1=0;//GM_Z Muon Chamber hits
	     double relIso1=1.0;// relative Isolation for muons(Z paper);
	      if((itM->track()).isNonnull() && (itMuon->track()).isNonnull()){ 
		
 
		histset[88]->Fill(fabs((itMuon->track())->vx()-(itM->track())->vx()));// for j_psi
		histset[89]->Fill(fabs((itMuon->track())->vy()-(itM->track())->vy()));//for j_psi
		histset[90]->Fill(fabs((itMuon->track())->vz()-(itM->track())->vz()));//for_j/psi
		int ValidHits1=0, PixelHits1=0;
		int GM_ValidHits1=0, GM_PixelHits1=0;
		double dx,dy,dz;
		dx=(itMuon->track())->vx()-(itM->track())->vx();// for j_psi
		dy=(itMuon->track())->vy()-(itM->track())->vy();// for j_psi
		dz=(itMuon->track())->vz()-(itM->track())->vz();// for j_psi
    
        //Upisilon Varriables
            double dref[3]={abs((itMuon->track())->vx()-(itM->track())->vx()),abs((itMuon->track())->vy()-(itM->track())->vy()),abs((itMuon->track())->vz()-(itM->track())->vz())};
        
            
            
            
            
		const reco::HitPattern& p1 = (itM->track())->hitPattern();
	 
		// loop over the hits of the track
		for (int n=0; n<p1.numberOfHits(); n++) {
		  uint32_t hit = p1.getHitPattern(n);

		  // if the hit is valid and in pixel 
		  if (p1.validHitFilter(hit) && p1.pixelHitFilter(hit)) PixelHits1++;
	 
		  if (p1.validHitFilter(hit))ValidHits1++;
		  
		}


		// for W 
        
		if(itM->isGlobalMuon() &&(itM->globalTrack()).isNonnull()){
		  if(itM->pt()>10.0) Wflag2=0;// to reject events for W mass


		  //for  Z
		  relIso1=((itM->isolationR03()).sumPt+(itM->isolationR03()).emEt+(itM->isolationR03()).hadEt)/itM->pt();
		  bool result1= muon::isGoodMuon(*itMuon, muon::TMLastStationTight);
		  const reco::HitPattern& GM_p1 = (itM->globalTrack())->hitPattern();
		  goodhit1=GM_p1.numberOfValidMuonHits();
		  for (int n=0; n<GM_p1.numberOfHits(); n++) {
		    uint32_t hit = GM_p1.getHitPattern(n);
		     
		    // if the hit is valid and in pixel 
		    if (GM_p1.validHitFilter(hit) && GM_p1.pixelHitFilter(hit)) GM_PixelHits1++;
	 
		    if (GM_p1.validHitFilter(hit))GM_ValidHits1++;
	 
		  }


		  if( GM_ValidHits1>=10 && GM_PixelHits1>=1 && goodhit1>=1 && (itM->globalTrack())->normalizedChi2()<10.0 && fabs((itM->globalTrack())->dxy(point))<0.2 && result1 && relIso1<0.15)
		    {
		      if(itM->pt()>20 && fabs(itM->eta())<2.1) Zflag2=1;
		      
		    }


		}// if(itM->isGlobalMuon().....) ends





     // for j_psi
       if (ValidHits1>=12 && PixelHits1>=2 && (itM->track())->normalizedChi2()<4.0 && dx<=0.1 && dy<=0.1 && dz<=0.3 && muon::isGoodMuon(*itM,muon::TMLastStationAngTight) && muon::isGoodMuon(*itM,muon::TrackerMuonArbitrated) && itM->isTrackerMuon())
	 {
	    if(fabs(itM->eta())<1.3 && itM->pt()>3.3)jpsi_flag2=1;
	    else if(fabs(itM->eta())>1.3 && fabs(itM->eta())<2.2 && itM->p()>2.9)jpsi_flag2=1;
	    else if(fabs(itM->eta())>2.2 && fabs(itM->eta())<2.4 && itM->pt()>0.8)jpsi_flag2=1;
	 }
       
    //Cuts for Upsilon [arXiv 1012.5545 §3.3]
    if (ValidHits1>=12 && PixelHits1>=1 && (itM->track())->normalizedChi2()<5 && sqrt(dref[0]*dref[0] + dref[1]*dref[1])<=0.2 && dref[2]<=5 && muon::isGoodMuon(*itM,muon::TMLastStationAngTight) && muon::isGoodMuon(*itM,muon::TrackerMuonArbitrated) && itM->isTrackerMuon())
              {
                  if (abs (itM->eta()) < 1.6 && abs(itMuon->eta()) <1.6 && itM->pt()>3.5 && itMuon->pt()>3.5 || abs(abs(itM->eta())-2) < 0.4  && abs(abs(itMuon->eta())-2) < 0.4  && itM->pt()>2.5 && itMuon->pt()>2.5 ) Yflag=1;
                  
              }

              
              
              
              
              
	      }//isNonnull() ends
       if(itM->charge()==-itMuon->charge() )// unlike charges
		{
		  
		  //------------invariant mass------------------------------//

		  s1=sqrt(((itMuon->p())*(itMuon->p())+sqm1)*((itM->p())*(itM->p())+sqm1));
		  s2=itMuon->px()*itM->px()+itMuon->py()*itM->py()+itMuon->pz()*itM->pz();
		  s=sqrt(2.0*(sqm1+(s1-s2)));
		 
		  scorr=(1.0+0.00038+0.0003*(itM->eta()+itMuon->eta())*(itM->eta()+itMuon->eta())/4)*s;
		  w=200/log(10)/s;
		  histset[101]->Fill(log10(s),w);
		  if(jpsi_flag2==1 && jpsi_flag1==1){//fill the jpsi hists only if these criteria are fulfilled
		  histset[32]->Fill(s);
		  histset[65]->Fill(scorr);
          M=invMass(itM->p(),itM->px(),itM->py(),itM->pz(),mumass,itMuon->p(),itMuon->px(),itMuon->py(),itMuon->pz(),mumass);
          histset[207]->Fill(M);
		  }//end of if(jpsi_flag2)
		  if(Zflag1==1 && Zflag2==1)// fill only if  z flags are fine
		    {
		      histset[94]->Fill(s);
		      
		    }


            
		  //------------------rapidity------------------------------//

		  s3=(itMuon->p())*(itMuon->p())+(itM->p())*(itM->p());
		  s4= sqrt(s3+2*s2+s*s);
		  pz=itM->pz()+itMuon->pz();
		  rap=0.5*log((s4+pz)/(s4-pz));

		   if(jpsi_flag2 && jpsi_flag1){//fill the jpsi hists only if these criteria are fulfilled
		  if(fabs(rap)<1.2){ histset[39]->Fill(s); histset[73]->Fill(scorr);}
		  else if(fabs(rap)<1.6){ histset[40]->Fill(s);histset[74]->Fill(scorr);}
		  else if(fabs(rap)<2.4){ histset[41]->Fill(s);histset[75]->Fill(scorr);}
		   }
            
            
            //-------------------------------Acceptance histogram----------------------------------------------//
            
            if (abs(s-jpsim)<.1 && jpsi_flag2==1 && jpsi_flag1==1 ){
                pt_jpsi=sqrt( pow( itM->px() + itMuon->px(), 2 ) + pow( itM->py() + itMuon->py(), 2 ) );
                histset[205]->Fill(pt_jpsi);
                hxhy[1]->Fill(  rap, pt_jpsi );
                hxhy[2]->Fill(  rap, pt_jpsi );
            }
            //--------------------------upsilon range--------------------------------------------------------//
            
            
            if (Yflag=1 && abs(rap) < 2 && abs(M-11)<3 ) {
                histset[250]->Fill(M);
                if (abs(itM->eta())<1 && abs(itMuon->eta())<1) histset[251]->Fill(M);
            }
            
            
            
		  if(fabs(itM->eta())<2.4 && fabs(itMuon->eta())<2.4){ histset[43]->Fill(s);histset[72]->Fill(scorr);}
		  if(fabs(itM->eta())<1. && fabs(itMuon->eta())<1.) {  histset[42]->Fill(s);histset[71]->Fill(scorr);}
		}

       if(itM->charge()==itMuon->charge() )// like charges invariant mass
		{
		  s1=sqrt(((itMuon->p())*(itMuon->p())+sqm1)*((itM->p())*(itM->p())+sqm1));
		  s2=itMuon->px()*itM->px()+itMuon->py()*itM->py()+itMuon->pz()*itM->pz();
		  s=sqrt(2.0*(sqm1+(s1-s2)));
		  
		  w=200/log(10)/s;
		  histset[103]->Fill(log10(s),w);
		  if(jpsi_flag2 && jpsi_flag1){// only if the jpsi cuts are fulfilled
		    histset[45]->Fill(s);
		    histset[82]->Fill(s);
		  }
		  if(Zflag1==1 && Zflag2==1)//only if Z cuts are satisfied
		    {
		      
		      histset[95]->Fill(s);
		      
		    }
		}
	    }//for(;itM!=muons;...)
	  	 
  
    }//(muons->size()>=2) ends
  	  //-------------------------------------------W Block--------------------------------------------------------------------//
    Mt=sqrt(2.0*(itMuon->pt())*((pfmets->front()).pt())*(1.0-cos((itMuon->phi())-((pfmets->front()).phi()))));
    histset[129]->Fill(Mt);
	  if(Wflag1==1 && Wflag2==1)
	    {
	      histset[130]->Fill(Mt);
	      histset[137]->Fill((pfmets->front()).pt());
	    }

}// Muon Collection for loop ends




/*------------------------------------------------------------------------------------------------------------------------------------------------
  -----------------------------------------------ELECTRON COLLECTION--------------------------------------------------------------------------------
  ---------------------------------------------------------------------------------------------------------------------------------------------------*/






histset[105]->Fill(electrons->size());
 for( reco::GsfElectronCollection::const_iterator it = electrons->begin(); it !=
 electrons->end(); it++) {
   int Z_ElecFlag1=0;
   int W_ElecFlag1=0;
    int Z_ElecFlag2=0;
    int W_ElecFlag2=1;
   
   double Mt;
   int misshits, misshits1;
   histset[96]->Fill(it->p());
    histset[97]->Fill(it->pt());
    histset[98]->Fill(it->eta());
    histset[99]->Fill(it->phi());
    histset[106]->Fill(fabs((it->superCluster())->eta()));
    histset[118]->Fill((it->et()));
    histset[107]->Fill(fabs((it->superCluster())->rawEnergy()));
    if(it->isEB()){
    
      histset[108]->Fill(it->dr03TkSumPt()/(it->pt()));//for electron, pt and et are nearly the same.
      histset[109]->Fill(it->dr03EcalRecHitSumEt()/(it->pt()));
      histset[110]->Fill(it->dr03HcalTowerSumEt()/(it->pt()));
    }
    if(it->isEE()){
   
      histset[120]->Fill(it->dr03TkSumPt()/(it->pt()));
      histset[121]->Fill(it->dr03EcalRecHitSumEt()/(it->pt()));
      histset[122]->Fill(it->dr03HcalTowerSumEt()/(it->pt()));
    }
    misshits=((it->gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
    histset[111]->Fill(misshits); // number of missing hits
    histset[113]->Fill(it->convDist());
    histset[112]->Fill(it->convDcot());

    if(it->isEB()){
    histset[114]->Fill(it->sigmaIetaIeta());
    histset[115]->Fill(it->deltaPhiSuperClusterTrackAtVtx());
    histset[116]->Fill(it->deltaEtaSuperClusterTrackAtVtx());
    histset[117]->Fill(it->hcalOverEcal());
    }
    if(it->isEE()){
    histset[123]->Fill(it->sigmaIetaIeta());
    histset[124]->Fill(it->deltaPhiSuperClusterTrackAtVtx());
    histset[125]->Fill(it->deltaEtaSuperClusterTrackAtVtx());
    histset[126]->Fill(it->hcalOverEcal());
    }
    
    //checking if the electron superCluster is within ECAL acceptance
   
    if(((fabs((it->superCluster())->eta()))<1.4442) ||( (fabs((it->superCluster())->eta()))<2.5 && (fabs((it->superCluster())->eta()))>1.5666))
 {
       
	if(it->et()>20.)

	  {
	    //starting with WP80 cuts!!!

	    if(it->isEB())//for barrel
	      {
		if((it->dr03TkSumPt()/(it->pt()))<0.09 && (it->dr03EcalRecHitSumEt()/(it->pt()))<0.07 && (it->dr03HcalTowerSumEt()/(it->pt()))<0.10 && misshits<1 && it->convDist()<0.02 && it->convDcot()<0.02)
		  {
		    if(it->sigmaIetaIeta()<0.01 && fabs(it->deltaPhiSuperClusterTrackAtVtx())<0.06 && fabs(it->deltaEtaSuperClusterTrackAtVtx())<0.004 && it->hcalOverEcal()<0.04)
		      { Z_ElecFlag1=1;
			W_ElecFlag1=1;
			
		      }
		    
		  }
	      }
	    else if(it->isEE())//for endcap
	      {
		if((it->dr03TkSumPt()/(it->pt()))<0.04 && (it->dr03EcalRecHitSumEt()/(it->pt()))<0.05 && (it->dr03HcalTowerSumEt()/(it->pt()))<0.025 && misshits<1 && it->convDist()<0.02 && it->convDcot()<0.02)
		  {
		    if(it->sigmaIetaIeta()<0.03 && fabs(it->deltaPhiSuperClusterTrackAtVtx())<0.03 && it->hcalOverEcal()<0.025)
		      { Z_ElecFlag1=1;
			W_ElecFlag1=1;
			
		      }
		  }	
	      }
	  }//	if(it->et()>20) ends

	  
      }//if(((fabs(it->superCluster())->eta())<1.4442) ||( (fabs(it->superCluster())->eta())<2.5 && (fabs(it->superCluster())->eta())>1.5666)) ends

    if(electrons->size()>=2){
    GsfElectronCollection::const_iterator ite=it;
     ite++;
    
    
      for(;ite!=electrons->end();ite++)
	{
	   Z_ElecFlag2=0;
	   W_ElecFlag2=1;// to reject events with a second electron satisfying some loose criteria
	  misshits1=((ite->gsfTrack())->trackerExpectedHitsInner()).numberOfHits();
	  if((fabs((ite->superCluster())->eta())<1.4442) ||( fabs((ite->superCluster())->eta())<2.5 && fabs((ite->superCluster())->eta())>1.5666))
	 
      {
	
	if(ite->et()>20)
	  { 
	   
	  //--------------------------------------------------------Z block------------------------------------------------------------------------------//
	  //starting with WP80 cuts!!!

	    if( Z_ElecFlag1==1)
	      {// to check for the 2nd electron only if the 1st flag is true
		
	    if(ite->isEB() )
	      {
		if((ite->dr03TkSumPt()/(ite->pt()))<0.09 && (ite->dr03EcalRecHitSumEt()/(ite->pt()))<0.07 && (ite->dr03HcalTowerSumEt()/(ite->pt()))<0.10 && misshits1<1 && ite->convDist()<0.02 && ite->convDcot()<0.02)
		  {
		    
		    if(ite->sigmaIetaIeta()<0.01 && fabs(ite->deltaPhiSuperClusterTrackAtVtx())<0.06 && fabs(ite->deltaEtaSuperClusterTrackAtVtx())<0.004 && ite->hcalOverEcal()<0.04)
		      {  
			Z_ElecFlag2=1;
			
		      }
		  }
	      }
	    else if(ite->isEE())//for endcap
	      {
		if((ite->dr03TkSumPt()/(ite->pt()))<0.04 && (ite->dr03EcalRecHitSumEt()/(ite->pt()))<0.05 && (ite->dr03HcalTowerSumEt()/(ite->pt()))<0.025 && misshits1<1 && ite->convDist()<0.02 && ite->convDcot()<0.02)
		  {
		   
		    if(ite->sigmaIetaIeta()<0.03 && fabs(ite->deltaPhiSuperClusterTrackAtVtx())<0.03 && ite->hcalOverEcal()<0.025)
		      {
			Z_ElecFlag2=1;
		
		      }
		  }	
	      }
	    } //if( Z_ElecFlag1==1) ends

	    //----------------------------------------------------W block---------------------------------------------------------------------------------//

	    //starting with WP95 cuts
	    
	    //for rejecting events with a second electron, for calculating the W peak
	    if(W_ElecFlag1==1)
	      {// to check for the 2nd electron only if the 1st flag is true

		
	     if(ite->isEB() )
	      {
		if((ite->dr03TkSumPt()/(ite->pt()))<0.15 && (ite->dr03EcalRecHitSumEt()/(ite->pt()))<2.0 && (ite->dr03HcalTowerSumEt()/(ite->pt()))<0.12 && misshits1<2 )
		  {
		    if(ite->sigmaIetaIeta()<0.01  && fabs(ite->deltaEtaSuperClusterTrackAtVtx())<0.007 && ite->hcalOverEcal()<0.15)
		      { W_ElecFlag2=0;
		      }
		  }
	      }
	    else if(ite->isEE())//for endcap
	      {
		if((ite->dr03TkSumPt()/(ite->pt()))<0.08 && (ite->dr03EcalRecHitSumEt()/(ite->pt()))<0.06 && (ite->dr03HcalTowerSumEt()/(ite->pt()))<0.05 && misshits1<2)
		  {
		    if(ite->sigmaIetaIeta()<0.03 && ite->hcalOverEcal()<0.07)
		      { W_ElecFlag2=0;}
		  }	
	      }

	      }//if( W_ElecFlag1==1) ends


	  }//if(ite->et()>20) ends

	  
      }//if(((fabs(ite->superCluster())->eta())<1.4442) ||( (fabs(ite->superCluster())->eta())<2.5 && (fabs(ite->superCluster())->eta())>1.5666))  ends


	  //---------------------------------------------Z  Block----------------------------------------------------------------//
	  if(it->charge()==-ite->charge() )// unlike charges
		{
		  
		  //------------invariant mass------------------------------//

		  s1=sqrt(((ite->p())*(ite->p())+sqm1)*((it->p())*(it->p())+sqm1));
		  s2=it->px()*ite->px()+it->py()*ite->py()+it->pz()*ite->pz();
		  s=sqrt(2.0*(sqm1+(s1-s2)));
		  histset[104]->Fill(s);
		  if(Z_ElecFlag1==1 && Z_ElecFlag2==1){
		    histset[119]->Fill(s);
		    if(s>=85. && s<=95.)
		      {
			histset[138]->Fill(it->convDcot());
			histset[139]->Fill(it->convDist());
		      }
		  }
		  
		}//if((it->charge()==-ite->charge() ) ends


	  if(it->charge()==ite->charge() )// like charges
		{
		  
		  //------------invariant mass------------------------------//

		  s1=sqrt(((ite->p())*(ite->p())+sqm1)*((it->p())*(it->p())+sqm1));
		  s2=it->px()*ite->px()+it->py()*ite->py()+it->pz()*ite->pz();
		  s=sqrt(2.0*(sqm1+(s1-s2)));
		  histset[127]->Fill(s);
		  if(Z_ElecFlag1==1 && Z_ElecFlag2==1)histset[128]->Fill(s);
		}//if((it->charge()==-ite->charge() ) ends
	  

	  
	}//for(;ite!=electrons....)ends
    }//if(electrons->size()>=2) ends


	  //-------------------------------------------W Block--------------------------------------------------------------------//
    Mt=sqrt(2.0*(it->pt())*((pfmets->front()).pt())*(1.0-cos((it->phi())-((pfmets->front()).phi()))));
    histset[131]->Fill(Mt);
	  if(W_ElecFlag1==1 && W_ElecFlag2==1)
	    {
	      histset[132]->Fill(Mt);
	      histset[136]->Fill((pfmets->front()).pt());
	    }

 }//for(reco::GsfElectronCollection........) ends
 
    
//////////////// ------------------------------ D* Meson Analysis ------------------------------ ////////////////
    
    //Loop over all non zero tracks with an iteractor K
    
    if (tracks->size()>=3) { //check for track collection with more than 3 tracks
        
       // cout<< "Analysing " << tracks->size() << "tracks" << endl; // checking no of tracks
 
        
        //Now loop over tracks assuming they are kaons
            for( reco::TrackCollection::const_iterator itK1 = tracks->begin(); itK1 !=tracks->end(); itK1++) {

               
                if ( itK1->pt() > 0.5) {
                     //display Kaon momentum for validation
                    //cout << " 'Kaon' Momenta values are " << itK1->pt() << " " << itK1->p() << " GeV/c " << endl;
                    for (reco::TrackCollection::const_iterator itP2 = tracks->begin(); itP2 !=tracks->end() ; itP2 ++){
                        
                        //check P2 is not also K1 and apply Pt cut
                        if ( itP2 !=itK1 && itP2->pt() >0.5 )  {
                            //display pion momentum for validation
                            //cout <<  " 'Pion' Momenta values are " << itP2->pt()<< " " << itP2->p()  << " GeV/c " <<endl;
                    
                            // check distance of track origins,
                            //Calculate tracks origin seperations in x, y and z
                            double vc[3] ={abs(itK1->vx()-itP2->vx()), abs(itK1->vy()- itP2->vy()), abs(itK1->vz()-itP2->vz())};
                            
                                if ( sqrt (vc[0]*vc[0] + vc[1]*vc[1])<0.1 && vc[2] <0.1){ //Check they are seperated by a radius of 1mm in x and y and a distance of 0.1 in z
                        
                                    //print them out
                                    //cout << "X seperation:" << vc[0] << " | Y seperation: " << vc[1] << " | Z seperation: " << vc[2] << " cm"<< endl;
                                    
                                    //Assign flags to right and wrong charge
                                    int cf=0; // 1 if right charge, -1 if wrong charge
                                    
                                        if(itP2->charge()!= itK1->charge() )      cf = 1; //right charge
                                        else if (itP2->charge() == itK1->charge() ) cf =-1; //wrong charge
                                    //cout << "Kaon charge: " << itK1->charge() << "| Pion Charge:  " << itP2->charge() <<endl;
                                
                                        //Calculate D0 invariant mass,
                                    double MD0 = invMass( itK1->p(),itK1->px(),itK1->py(),itK1->pz(),Kmass, itP2->p(),itP2->px(),itP2->py(),itP2->pz(), Pimass );
                                       
                                        if (cf==1) {
                                            //cout <<"BINGO: " << MD0 <<endl;
                                            histset[206]->Fill(MD0);
                                        }
                                    
                                //----Begin a loop over a third track to find D*
                                    //First check D0 mass is within 50 MeV of the Actual D0 mass
                                    if (abs(MD0- MD0Actual)<.05) { // check it is a reasonable D0 candidate
                                        
                                        //Calculate D px, py, pz, p and pt
                                        double paxisD0[3]={itK1->px()+itP2->px(), itK1->py() +itP2->py(), itK1->pz() + itP2->pz() };
                                        double pD0 = sqrt(paxisD0[0]*paxisD0[0] +paxisD0[1]*paxisD0[1] +paxisD0[2]*paxisD0[2]);
                                        double ptD0= sqrt(paxisD0[0]*paxisD0[0] +paxisD0[1]*paxisD0[1]);
                                         
                                        //locate D0 vertex using average coordinate of tracks
                                        double vcD0[3] ={0.5*(itK1->vx()+itP2->vx()), 0.5*(itK1->vy()- itP2->vy()), 0.5*(itK1->vz()-itP2->vz())};
                                            
                                            
                                            //Now start a third loop for the slow pion
                                            for (reco::TrackCollection::const_iterator itPS3 = tracks->begin(); itPS3 !=tracks->end() ; itPS3 ++){
                                                if (itPS3->charge()!=itK1->charge() && itPS3!=itP2 && itPS3!=itK1 ) { //check charge and that it not K1 or P2
                                            
                                                    double vc2[3] ={abs(vcD0[0]-itPS3->vx()), abs(vcD0[1]- itPS3->vy()), abs(vcD0[2]-itPS3->vz())};
                                                    
                                                        if (sqrt(vc2[0]*vc2[0] + vc2[1]*vc2[1])<0.1 && vc2[2] <0.1){ // Check the slow pion originates from the same region as the D0
                                                            
                                                            //Now Calculate D* mass
                                                            double MDstar= invMass(pD0,paxisD0[0],paxisD0[1],paxisD0[2],MD0, itPS3->p(), itPS3->px(), itPS3->py(), itPS3->pz(),Pimass);
                                                            
                                                            //D* mass diff cut <0.17 to reduce overflow of histogram
                                                            if (MDstar-MD0 <0.17) {
                                                                
                                                                
                                                             //Start of z cut
                                                                
                                                                double sumpt=0;
                                                                int n=0;     //declare varriable
                                                                //Loop over all tracks
                                                                for( reco::TrackCollection::const_iterator itSum = tracks->begin(); itSum!=tracks->end(); itSum++) {
                                                                    
                                                                //Check track origins
                                                                double vc3[3] ={abs(vcD0[0]-itSum->vx()), abs(vcD0[1]- itSum->vy()), abs(vcD0[2]-itSum->vz())};
                                                                    if (sqrt(vc3[0]*vc3[0] + vc3[1]*vc3[1])<0.1 && vc3[2] <0.1){
                                                                        sumpt += abs(itSum->pt()); // sum pt for all tracks
                                                                        n++; //Fill a counter histogram to see how many tracks are present here?!
                                                                    }// end of Sumpt vertex check
                                                                    
                                                                    
                                                                }
                                                                //end of sum loop
                                                                double z=ptD0/sumpt; // calculate z after loop
                                                                histset[210]->Fill(z);
                                                                histset[211]->Fill(n);
                                                                //Fill deltaM histograms depending on Charge flag cf and z values
                                                                if      (cf==1 && z>0.05)  histset[208]->Fill(MDstar-MD0);
                                                                else if (cf==-1 && z>0.05) histset[209]->Fill(MDstar-MD0);
                                                               
                                                            }
                                       
                                                    }//end of Pi Slow vertex check.
                                                }//end of PS3 vs K1 charge and iterators check
                                            }//end of PS3 loop
                                    }//End of charge check
                            }//vertex seperation cut
                        } //Pt and K1 cut
                    }//end of P2 loop
                }//end of K1 Pt cut
            }//end of K loop
    }//end of track size cut

    

}//DemoAnalyzer ends


inline double invMass(double p1,double p1x, double p1y, double p1z, double m1 , double p2, double p2x, double p2y, double p2z, double m2)
{
    s= sqrt(   m1*m1 + m2*m2 + 2.0*( sqrt( ( m1*m1 + p1 * p1 ) * ( m2*m2 + p2 * p2 ) ) - ( p1x * p2x + p1y * p2y + p1z * p2z ) ) ) ;
    return s;
}

// ------------ method called once each job just before starting event loop  ------------
void 
DemoAnalyzer::beginJob()
{

}

// ------------ method called once each job just after ending the event loop  ------------
void 
DemoAnalyzer::endJob() 
{
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);