#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TVector3.h>

#include "Constants.h"

using namespace std;
using namespace Constants;

bool inFVVector(TVector3 vector) {

  if(vector.X() < (FVx - borderx) && (vector.X() > borderx) && (vector.Y() < (FVy/2. - bordery)) && (vector.Y() > (-FVy/2. + bordery)) && 
     (vector.Z() < (FVz - borderz)) && (vector.Z() > borderz)) return true;
  else return false;
}

void NeutrinoSelectionFilter::Loop() {

  //--------------------------------------------------//

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nbytes = 0, nb = 0;

   TH1D::SetDefaultSumw2();
   TH2D::SetDefaultSumw2();

   //--------------------------------------------------//

   // Output file

   TString FileName = "MicroBooNE_Reco.root";
   TFile* OutputFile = new TFile(FileName,"recreate");
   std::cout << std::endl << "File " << FileName << " to be created"<< std::endl << std::endl;

   //--------------------------------------------------//

   // Plot declaration

   // Max: extend the plots to have a 2nd dimension aka [NInte][3]
   // 2nd index = 0/1/2 (all events/true CC1p0pi/true non-CC1p0pi)

   // Max: introduce delta pt, alpha alphat, pmiss mag & cos theta plots
   // Max: beyond the 1D plots, introduce the 2D plots connecting the true vs backtracked values

   TH1D* RecoMuonCosThetaPlot[NInte];

   TH1D* BacktrackedMuonCosThetaPlot[NInte];

   //--------------------------------------------------//

   // Loop over the interaction processes

   for (int inte = 0; inte < NInte; inte++) {

     //--------------------------------------------------//

     // 1D analysis

     RecoMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"RecoMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

     BacktrackedMuonCosThetaPlot[inte] = new TH1D(InteractionLabels[inte]+"BacktrackedMuonCosThetaPlot",LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);

     //--------------------------------------------------//

   } // End of the loop over the interaction processes

   //--------------------------------------------------//

   // Loop over the events

   for (Long64_t jentry=0; jentry<nentries;jentry++) {

     //--------------------------------------------------//

      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;

      //--------------------------------------------------//

      // MC weight to scale events to data pot

      double event_weight = (data_pot / mc_pot) * weightSplineTimesTune;

      //--------------------------------------------------//

      // Max: introduce neutron counter & vector

      int TrueMuonCounter = 0, TrueProtonCounter = 0, TrueChargedPionCounter = 0, TruePi0Counter = 0;
      int TrueNeutronCounter = 0, TrueHeavierMesonCounter = 0;
      
      std::vector<int> VectorTrueMuonIndex; VectorTrueMuonIndex.clear();
      std::vector<int> VectorTrueProtonIndex; VectorTrueProtonIndex.clear();
      std::vector<int> VectorTrueNeutronIndex; VectorTrueNeutronIndex.clear();

      int NMCParticles = mc_pdg->size();

      for (int WhichMCParticle = 0; WhichMCParticle < NMCParticles; WhichMCParticle++) {

	// MC truth information for the final-state primary particles

	// CC numu events

	if ( ccnc == 0 && nu_pdg == 14) {

	  TVector3 MCParticle(mc_px->at(WhichMCParticle),mc_py->at(WhichMCParticle),mc_pz->at(WhichMCParticle));
	  double MCParticleMomentum = MCParticle.Mag();
	  int MCParticlePdg = mc_pdg->at(WhichMCParticle);

	  if ( MCParticlePdg == MuonPdg && MCParticleMomentum >= ArrayNBinsMuonMomentum[0] && MCParticleMomentum <= ArrayNBinsMuonMomentum[NBinsMuonMomentum]) 
	    { TrueMuonCounter++;  VectorTrueMuonIndex.push_back(WhichMCParticle); }

	  if ( MCParticlePdg == ProtonPdg && MCParticleMomentum >= ArrayNBinsProtonMomentum[0] && MCParticleMomentum <= ArrayNBinsProtonMomentum[NBinsProtonMomentum] ) 
	    { TrueProtonCounter++; VectorTrueProtonIndex.push_back(WhichMCParticle); }

	  if ( fabs(MCParticlePdg) == AbsChargedPionPdg && MCParticleMomentum >= ChargedPionMomentumThres ) 
	    { TrueChargedPionCounter++; }

	  if (MCParticlePdg == NeutralPionPdg) { TruePi0Counter++; }


	} // End of the demand stable final state particles and primary interactions

      } // end of the loop over the MCParticles

      //--------------------------------------------------//      

      // CC1p0pi signal events

      bool CC1p0pi = false;

      if (TrueMuonCounter == 1 && TrueProtonCounter == 1 && TrueChargedPionCounter == 0 && TruePi0Counter == 0 && TrueHeavierMesonCounter == 0) {

	CC1p0pi = true;

      }

      // Unlike the truth level study, we just flag events as signal or not but we don't continue
      //if (!CC1p0pi) { continue; }

      //--------------------------------------------------//

      // Max: check for true CC1p0pi events inside fiducial volume of interest
      // if true vertex outside, set CC1p0pi = false 
      // but do NOT "continue"!

      //--------------------------------------------------//

      // Max: Need if-statement for true CC1p0pi events
      // else: fill in the vectors with (-999.,-999.,-999.)
      // and still calculate the quantities below


      // Max: Define muon / proton vectors, deltapt, delta alphaT & pmiss

      //--------------------------------------------------//

      // Classify events based on interaction

      int genie_mode = -1;

      if (interaction == 0) { genie_mode = 1; } // QE
      else if (interaction == 10) { genie_mode = 2; } // MEC
      else if (interaction == 1) { genie_mode = 3; } // RES
      else if (interaction == 2) { genie_mode = 4; } // DIS
      else { genie_mode = 5; } // COH / other

      //--------------------------------------------------//
      //--------------------------------------------------//

      // Reco level

      // Loop over the PFParticles and keep track of the number of tracks and showers 

      int reco_shower_count = 0;
      int reco_track_count = 0;

      std::vector<int> CandidateIndex;

      for ( int p = 0; p < n_pfps; ++p ) {

	// Only check direct neutrino daughters (generation == 2)

	unsigned int generation = pfp_generation_v->at( p );
	if ( generation != 2u ) continue;

	float tscore = trk_score_v->at( p );
	if ( tscore <= TRACK_SCORE_CUT ) { ++reco_shower_count; }
	else { ++reco_track_count; CandidateIndex.push_back(p); }
  
      }

      // Max: add track & shower multiplicity plots to see what
      // these distributions look like on MicroBooNE
      // Use the weight to fill in the plots

      //--------------------------------------------------//

      // Requirement for exactly 2 tracks and 0 showers for a reco CC1p0pi selection 

      if (reco_shower_count != 0) { continue; }
      if (reco_track_count != 2) { continue; }

      // Requirement for exactly one neutrino slice

      if (nslice != 1) { continue; }

      //--------------------------------------------------//

      // Identify candidate muon & proton
      // Muon = the one with the highest LLR PID Score
      // Proton = the one with the lowest LLR PID Score

      if (CandidateIndex.size() != 2) { continue; }
      float first_pid_score = trk_llr_pid_score_v->at( CandidateIndex.at(0) );
      float second_pid_score = trk_llr_pid_score_v->at( CandidateIndex.at(1) );

      int CandidateMuonIndex = -1.;
      int CandidateProtonIndex = -1.;

      if (first_pid_score > second_pid_score) { 
	
	CandidateMuonIndex = CandidateIndex.at(0); 
	CandidateProtonIndex = CandidateIndex.at(1); 
	
      } else { 
	
	CandidateMuonIndex = CandidateIndex.at(1); 
	CandidateProtonIndex = CandidateIndex.at(0); 
	
      }

      // MuonPdg in this case / Pandora = Track-like object

      if (pfpdg->at(CandidateMuonIndex) != MuonPdg) { continue; }
      if (pfpdg->at(CandidateProtonIndex) != MuonPdg) { continue; }

      //--------------------------------------------------//

      // Reconstructed Vertex

      TVector3 reco_vertex(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);

      // Max: demand that the reco vertex is contained in the fiducial volume

      //--------------------------------------------------//

      // Muon kinematics

      double CandidateMuonTrackTheta = trk_theta_v->at(CandidateMuonIndex); // rad
      double CandidateMuonTrackPhi = trk_phi_v->at(CandidateMuonIndex); // rad
      double CandidateMuMom = trk_range_muon_mom_v->at(CandidateMuonIndex); // GeV/c
      
      double MuonTrackStartX = trk_sce_start_x_v->at(CandidateMuonIndex);
      double MuonTrackStartY = trk_sce_start_y_v->at(CandidateMuonIndex);
      double MuonTrackStartZ = trk_sce_start_z_v->at(CandidateMuonIndex);
      
      double MuonTrackEndX = trk_sce_end_x_v->at(CandidateMuonIndex);
      double MuonTrackEndY = trk_sce_end_y_v->at(CandidateMuonIndex);
      double MuonTrackEndZ = trk_sce_end_z_v->at(CandidateMuonIndex);
      
      TVector3 CandidateMuonTrackStart(MuonTrackStartX,MuonTrackStartY,MuonTrackStartZ);
      TVector3 CandidateMuonTrackEnd(MuonTrackEndX,MuonTrackEndY,MuonTrackEndZ);
      bool CandidateMuonTrackStartContainment = inFVVector(CandidateMuonTrackStart);
      bool CandidateMuonTrackEndContainment = inFVVector(CandidateMuonTrackEnd);

      TVector3 TVector3CandidateMuon(-1,-1,-1);
      TVector3CandidateMuon.SetMag(CandidateMuMom); // GeV/c
      TVector3CandidateMuon.SetTheta(CandidateMuonTrackTheta); // rad
      TVector3CandidateMuon.SetPhi(CandidateMuonTrackPhi); // rad

      //--------------------------------------------------//

      // Max: declare the proton kinematics

      //--------------------------------------------------//

      // Momenta ranges for the muon candidate

      if (CandidateMuMom < ArrayNBinsMuonMomentum[0]) { continue; }
      if (CandidateMuMom > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }

      // Max: apply the momentum the proton candidate

      //--------------------------------------------------//

      // Fully contained muons

      if (!CandidateMuonTrackStartContainment) { continue; }
      if (!CandidateMuonTrackEndContainment) { continue; }

      // Max: do the same for proton candidate

      //--------------------------------------------------//

      // Start-To-Start & End-To-End Distance to avoid flipped tracks

      // Max: add the 3 lines below in when you have the proton candidate equivalents above declared

      //double LocalStartToStartDistance = (CandidateMuonTrackStart - CandidateProtonTrackStart).Mag();
      //double LocalEndToEndDistance = (CandidateMuonTrackEnd - CandidateProtonTrackEnd).Mag();
      //if (LocalStartToStartDistance > LocalEndToEndDistance) { continue; }

      //--------------------------------------------------//

      // Vertex - Track & relative differences

      // Max: add the 5 lines below in after you declare the proton start/end points

      //double CandidateMuStartVertexMag = (reco_vertex - CandidateMuonTrackStart).Mag();
      //double CandidatePStartVertexMag = (reco_vertex - CandidateProtonTrackStart).Mag();
      //double CandidateMuEndVertexMag = (reco_vertex - CandidateMuonTrackEnd).Mag();
      //double CandidatePEndVertexMag = (reco_vertex - CandidateProtonTrackEnd).Mag();
      //if (CandidateMuStartVertexMag > CandidateMuEndVertexMag) { continue; }

      // Max: do the same for the proton

      //--------------------------------------------------//

      // Max: declare delta pt & delta alphat for reco CC1p0pi candidates
      // Max: plot the x/y/z components of the reco vertex for the CC1p0pi candidates

      //--------------------------------------------------//
      //--------------------------------------------------//

      // Backtracking: matching the reconstructed objects to true ones
      // Not the same as the truth-level study
      // where things don't need to be reconstructed

      //--------------------------------------------------//

      // Muon candidate backtracking

      double CandidateMuonPx = backtracked_px->at(CandidateMuonIndex);
      double CandidateMuonPy = backtracked_py->at(CandidateMuonIndex);
      double CandidateMuonPz = backtracked_pz->at(CandidateMuonIndex);
      TVector3 TrueCandidateMuonP(CandidateMuonPx,CandidateMuonPy,CandidateMuonPz);
      
      double TrueCandidateMuonTrackPhi = TrueCandidateMuonP.Phi(); // rad
      double TrueCandidateMuonTrackTheta = TrueCandidateMuonP.Theta(); // rad
      double TrueCandidateMuonPMag = TrueCandidateMuonP.Mag();

      TVector3 True_TVector3CandidateMuon(-1,-1,-1);
      True_TVector3CandidateMuon.SetMag(TrueCandidateMuonPMag);
      True_TVector3CandidateMuon.SetTheta(TrueCandidateMuonTrackTheta);
      True_TVector3CandidateMuon.SetPhi(TrueCandidateMuonTrackPhi);

      // Max: do the same for the proton backtracking

      //--------------------------------------------------//

      // If it is a true CC1p0pi event
      // Check that it is a backtracked CC1p0pi event as well

      if (CC1p0pi) {

	if ( !(backtracked_pdg->at(CandidateMuonIndex) == MuonPdg && TrueCandidateMuonPMag > ArrayNBinsMuonMomentum[0] && TrueCandidateMuonPMag < ArrayNBinsMuonMomentum[NBinsMuonMomentum]) ) 
	  { CC1p0pi = false; }

	// Max: do the same for the proton

      }

      //--------------------------------------------------//

      // Max: define delta pt, delta alphat, and pmiss for backtracked vectors

      //--------------------------------------------------//

      // Fill in the plots

      // All events, 2nd index = 0


      // CC1p0pi events, 2nd index = 1


      // CC1p0pi events, 2nd index = 2


      //--------------------------------------------------//

      // Loop over the blips

      //--------------------------------------------------//

      // Reco CC1p0pi candidate event selected, print out details
      if (CC1p0pi) { cout << "run = " << run << ", sub = " << sub << ", evt = " << evt << endl; }

      //--------------------------------------------------//

   } // End of the loop over the events

   //--------------------------------------------------//

   // Divide by bin width with Reweight function

   //--------------------------------------------------//

} // End of the program
