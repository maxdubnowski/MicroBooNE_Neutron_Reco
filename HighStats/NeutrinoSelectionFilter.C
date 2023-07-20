#define NeutrinoSelectionFilter_cxx
#include "NeutrinoSelectionFilter.h"
#include <TH1D.h>
#include <TH2D.h>
#include <TString.h>
#include <TVector3.h>

#include "../Constants.h"

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

   // Max: introduce delta pt, delta alphat, pmiss mag & cos theta plots, neutron multiplicity*at backtracked level*
   // Max: beyond the 1D plots, introduce the 2D plots connecting the true vs backtracked values

   TH1D* RecoMuonCosThetaPlot[NInte][3];
   TH1D* RecoDeltaPtPlot[NInte][3];
   TH1D* RecoDeltaAlphaTPlot[NInte][3];
   TH1D* RecoPMissMagPlot[NInte][3];
   TH1D* RecoPMissDirPlot[NInte][3];


   TH1D* BacktrackedMuonCosThetaPlot[NInte][3];
   TH1D* BacktrackedNeutronMultiplicityPlot[NInte][3];
   TH1D* BacktrackedDeltaPtPlot[NInte][3];
   TH1D* BacktrackedDeltaAlphaTPlot[NInte][3];
   TH1D* BacktrackedPMissMagPlot[NInte][3];
   TH1D* BacktrackedPMissDirPlot[NInte][3];

   TH2D* BacktrackedvsRecoDeltaPtPlot[NInte][3];
   TH2D* BacktrackedvsRecoDeltaAlphaTPlot[NInte][3];
   TH2D* BacktrackedvsRecoPMissMagPlot[NInte][3];
   TH2D* BacktrackedvsRecoPMissDirPlot[NInte][3];
   

   TH1D* TrueVertexXPlot[NInte][3];
   TH1D* TrueVertexYPlot[NInte][3];
   TH1D* TrueVertexZPlot[NInte][3];
   TH1D* TruePMissMagPlot[NInte][3];
   TH1D* TruePMissDirPlot[NInte][3];
   TH1D* TrueMuonEnergyPlot[NInte][3];
   TH1D* TrueNeutronMultiplicityPlot[NInte][3];


   TH1D* RecoVertexXPlot[NInte][3];
   TH1D* RecoVertexYPlot[NInte][3];
   TH1D* RecoVertexZPlot[NInte][3];


   TH1D* TrackMultiplicity[3];
   TH1D* ShowerMultiplicity[3];
   TH2D* checkMomentumSame[3]; //For the three directions (0-X, 1-Y, 2-Z)
   
   TH1D* BlipPDG[NInte][3]; //Interaction mechanisms, 3 interaction types (All, CC1p0pi, nonCC1p0pi)
   TH1D* BlipMultiplicity[NInte][3]; //3 interaction types
   TH1D* AssociatedBlipMultiplicity[NInte][3];
   TH1D* BlipLocation[NInte][3][3]; //Interaction mechanisms, 3 directions (X,Y,Z) and 3 interaction types
   TH1D* BlipVertexDist[NInte][3];
   TH1D* ProtonBlipDist[NInte][3];

   TH1D* BlipVertexDistNeut[NInte][3][NNeut+1];
   TH1D* BlipProxTrkDistNeut[NInte][3][NNeut+1];
   
   TH2D* AssociatedBlipMultVsVertexDist[NInte][3][NNeut+1];
   TH2D* AssociatedBlipMultVsProxTrkDist[NInte][3][NNeut+1];

   TH2D* AssociatedBlipVertexDistVsProxTrkDist[NInte][3];

   //TH1D* BlipCharge[NInte][3];
   TH1D* BlipProxTrkDist[NInte][3];
   TH1D* BlipTime[NInte][3];
   TH1D* BlipLength[NInte][3];
   TH1D* BlipMaxWireSpan[NInte][3];
   TH1D* BlipEnergy[NInte][3];
   //TH1D* ProtonBlipCharge[NInte][3];
   TH1D* ProtonBlipProxTrkDist[NInte][3];
   TH1D* ProtonBlipTime[NInte][3];
   TH1D* ProtonBlipLength[NInte][3];
   TH1D* ProtonBlipMaxWireSpan[NInte][3];
   TH1D* ProtonBlipEnergy[NInte][3];

   TH2D* BlipProxDistVsPMissMag[NInte][3]; 
   TH2D* BlipVertexDistVsPMissMag[NInte][3]; 
   TH2D* BlipProxDistVsPMissDir[NInte][3];
   TH2D* BlipVertexDistVsPMissDir[NInte][3];
   TH2D* BlipMultiplicityVsVertexDist[NInte][3];
   TH2D* BlipMultiplicityVsProxDist[NInte][3];
   TH2D* ProxDistVsNeutronMultiplicity[NInte][3];
   TH2D* VertexDistVsNeutronMultiplicity[NInte][3];

   TH2D* BlipVertexDistVsEnergy[3];
   TH2D* ProtonBlipVertexDistVsEnergy[3];
   TH2D* BlipProxDistVsEnergy[3];
   TH2D* ProtonBlipProxDistVsEnergy[3];
   
   TH1D* PostBlipCutRecoMuonCosThetaPlot[NInte][3][NCuts];
   TH1D* PostBlipCutRecoDeltaPtPlot[NInte][3][NCuts];
   TH1D* PostBlipCutRecoDeltaAlphaTPlot[NInte][3][NCuts];
   TH1D* PostBlipCutTrueNeutronMultiplicityPlot[NInte][3][NCuts];


    


   // Initialize Histograms
   for(int type=0; type<3; type++){
     TrackMultiplicity[type] = new TH1D("TrackMultiplicity"+CC1p0piLabels[type], ";Track Multiplicity; Weighted Events", 10, -0.5, 9.5);
     ShowerMultiplicity[type] = new TH1D("ShowerMultiplicity"+CC1p0piLabels[type], ";Shower Multiplicity; Weighted Events", 10, -0.5, 9.5);
     BlipVertexDistVsEnergy[type] = new TH2D("BlipVertexDistVsEnergy"+ CC1p0piLabels[type], ";Vertex Distance; Energy Deposited" , 20, 0, 400, 20, 0, 25);
     ProtonBlipVertexDistVsEnergy[type] =new  TH2D("ProtonBlipVertexDistVsEnergy"+ CC1p0piLabels[type], ";Vertex Distance; Energy Deposited" , 20, 0, 400, 20, 0, 25);
     BlipProxDistVsEnergy[type] = new TH2D("BlipProxDistVsEnergy"+ CC1p0piLabels[type], ";Prox Trk Distance; Energy Deposited" , 20, 0, 150, 20, 0, 25);
     ProtonBlipProxDistVsEnergy[type] =new  TH2D("ProtonBlipProxDistVsEnergy"+ CC1p0piLabels[type], ";Prox Trk Distance; Energy Deposited" , 20, 0, 150, 20, 0, 25);
   }

   checkMomentumSame[0]= new TH2D("CheckMomentumXDirection", ";True ;Truer", 15, 0,1, 15, 0,1 );
   checkMomentumSame[1]= new TH2D("CheckMomentumYDirection", ";True ;Truer", 15, 0,1, 15, 0,1 );
   checkMomentumSame[2]= new TH2D("CheckMomentumZDirection", ";True ;Truer", 15, 0,1, 15, 0,1 );


   // Loop over the interaction processes
   for (int inte = 0; inte < NInte; inte++) {
     // 1D analysis 
     for (int type=0; type < NType; type++){
       RecoMuonCosThetaPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoMuonCosThetaPlot"+ CC1p0piLabels[type],LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
       RecoDeltaPtPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoDeltaPtPlot"+ CC1p0piLabels[type],LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
       RecoDeltaAlphaTPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoDeltaAlphaTPlot"+ CC1p0piLabels[type],LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
       RecoPMissMagPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoPMissMagPlot"+ CC1p0piLabels[type],";P_{miss} ;Weighted Events" , 10 ,0 ,1 );
       RecoPMissDirPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoPMissDirPlot"+ CC1p0piLabels[type],";cos(#theta_{miss}) ;Weighted Events" , 10 ,-1 ,1 );


       for (int cut=0; cut <NCuts; cut++){
	 PostBlipCutRecoMuonCosThetaPlot[inte][type][cut] = new TH1D(InteractionLabels[inte]+"PostBlipCutRecoMuonCosThetaPlot"+ BlipCuts[cut]+ CC1p0piLabels[type],LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
	
	 PostBlipCutRecoDeltaPtPlot[inte][type][cut] = new TH1D(InteractionLabels[inte]+"PostBlipCutRecoDeltaPtPlot"+ BlipCuts[cut]+ CC1p0piLabels[type],LabelXAxisDeltaPT,NBinsDeltaPT,ArrayNBinsDeltaPT);
	 
	 PostBlipCutRecoDeltaAlphaTPlot[inte][type][cut] = new TH1D(InteractionLabels[inte]+"PostBlipCutRecoDeltaAlphaTPlot"+ BlipCuts[cut]+ CC1p0piLabels[type],LabelXAxisDeltaAlphaT,NBinsDeltaAlphaT,ArrayNBinsDeltaAlphaT);
	 
	 PostBlipCutTrueNeutronMultiplicityPlot[inte][type][cut] = new TH1D(InteractionLabels[inte] +"PostBlipCutTrueNeutronMultiplicityPlot" + BlipCuts[cut]+ CC1p0piLabels[type], ";Neutron Multiplicity ; Weighted Events" , 6, -0.5, 5.5);
       } //End of selection cut naming



       BacktrackedMuonCosThetaPlot[inte][type] = new TH1D(InteractionLabels[inte]+"BacktrackedMuonCosThetaPlot"+ CC1p0piLabels[type],LabelXAxisMuonCosTheta,NBinsMuonCosTheta,ArrayNBinsMuonCosTheta);
       BacktrackedDeltaPtPlot[inte][type] = new TH1D(InteractionLabels[inte]+"BacktrackedDeltaPtPlot"+ CC1p0piLabels[type], LabelXAxisDeltaPT, NBinsDeltaPT, ArrayNBinsDeltaPT);
       BacktrackedDeltaAlphaTPlot[inte][type] = new TH1D(InteractionLabels[inte]+"BacktrackedDeltaAlphaTPlot"+ CC1p0piLabels[type], LabelXAxisDeltaAlphaT, NBinsDeltaAlphaT, ArrayNBinsDeltaAlphaT);
       BacktrackedPMissMagPlot[inte][type] = new TH1D(InteractionLabels[inte]+"BacktrackedPMissMagPlot"+ CC1p0piLabels[type], ";P_{miss} ;Weighted Events" , 10 ,0 ,1 );
       BacktrackedPMissDirPlot[inte][type] = new TH1D(InteractionLabels[inte]+"BacktrackedPMissDirPlot"+ CC1p0piLabels[type], ";cos(#theta_{miss}) ;Weighted Events" , 10 ,-1 ,1 );
       
       BacktrackedvsRecoDeltaPtPlot[inte][type] = new TH2D(InteractionLabels[inte]+"BacktrackedvsRecoDeltaPtPlot"+ CC1p0piLabels[type], ";Backtracked #deltap_{T} [GeV/c] ; Reco #deltap_{T} [GeV/c]",NBinsDeltaPT, ArrayNBinsDeltaPT, NBinsDeltaPT, ArrayNBinsDeltaPT); 
       BacktrackedvsRecoDeltaAlphaTPlot[inte][type] = new TH2D(InteractionLabels[inte]+"BacktrackedvsRecoDeltaAlphaTPlot"+ CC1p0piLabels[type], ";Backtracked #delta#alpha_{T} [deg] ; Reco #delta#alpha_{T} [deg]",NBinsDeltaAlphaT, ArrayNBinsDeltaAlphaT, NBinsDeltaAlphaT, ArrayNBinsDeltaAlphaT); 
       BacktrackedvsRecoPMissMagPlot[inte][type] = new TH2D(InteractionLabels[inte]+"BacktrackedvsRecoPMissMagPlot"+ CC1p0piLabels[type], ";Backtracked P_{miss} [GeV/c] ; Reco P_{miss} [GeV/c]",20,0,1.5,20,0,1.5); 
       BacktrackedvsRecoPMissDirPlot[inte][type] = new TH2D(InteractionLabels[inte]+"BacktrackedvsRecoPMissDirPlot"+ CC1p0piLabels[type], ";Backtracked cos#theta_{miss}  ; Reco cos#theta_{miss} ",20,-1,1,20,-1,1); 

       TrueVertexXPlot[inte][type] = new TH1D(InteractionLabels[inte]+"TrueVertexXPlot"+ CC1p0piLabels[type],";Vertex x [cm] ;Weighted Events" , NBinsVertexX , MinVertexX , MaxVertexX);
       TrueVertexYPlot[inte][type] = new TH1D(InteractionLabels[inte]+"TrueVertexYPlot"+ CC1p0piLabels[type],";Vertex y [cm] ;Weighted Events" , NBinsVertexY , MinVertexY , MaxVertexY);
       TrueVertexZPlot[inte][type] = new TH1D(InteractionLabels[inte]+"TrueVertexZPlot"+ CC1p0piLabels[type],";Vertex z [cm] ;Weighted Events" , NBinsVertexZ , MinVertexZ , MaxVertexZ);
       TruePMissMagPlot[inte][type] = new TH1D(InteractionLabels[inte]+"TruePMissMagPlot"+ CC1p0piLabels[type],";P_{miss} Mag [GeV/c] ;Weighted Events" , 10 , 0 ,1);
       TruePMissDirPlot[inte][type] = new TH1D(InteractionLabels[inte]+"TruePMissDirPlot"+ CC1p0piLabels[type],";cos(#theta_{miss}) ;Weighted Events" , 10 , -1 , 1);
       TrueMuonEnergyPlot[inte][type] = new TH1D(InteractionLabels[inte]+"TrueMuonEnergyPlot"+ CC1p0piLabels[type],";Energy [GeV] ;Weighted Events" , 15 , 0 , 1.5);
       TrueNeutronMultiplicityPlot[inte][type] = new TH1D(InteractionLabels[inte] +"TrueNeutronMultiplicityPlot" + CC1p0piLabels[type], ";Neutron Multiplicity ; Weighted Events" , 6, -0.5, 5.5);

       RecoVertexXPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoVertexXPlot"+ CC1p0piLabels[type],";Vertex x [cm] ;Weighted Events" , NBinsVertexX , MinVertexX , MaxVertexX);
       RecoVertexYPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoVertexYPlot"+ CC1p0piLabels[type],";Vertex y [cm] ;Weighted Events" , NBinsVertexY , MinVertexY , MaxVertexY);
       RecoVertexZPlot[inte][type] = new TH1D(InteractionLabels[inte]+"RecoVertexZPlot"+ CC1p0piLabels[type],";Vertex z [cm] ;Weighted Events" , NBinsVertexZ , MinVertexZ , MaxVertexZ);
            
       BlipPDG[inte][type] = new TH1D(InteractionLabels[inte]+"BlipPDG" + CC1p0piLabels[type], ";Blip PDG ;Weighted Events", 6000,-3000.5,2999.5);
       BlipMultiplicity[inte][type] = new TH1D(InteractionLabels[inte]+"BlipMultiplicity" + CC1p0piLabels[type], ";Blip Multiplicity ;Weighted Events", 50,29.5,179.5);
       AssociatedBlipMultiplicity[inte][type] = new TH1D(InteractionLabels[inte]+"AssociatedBlipMultiplicity" + CC1p0piLabels[type], ";Blip Multiplicity ;Weighted Events", 20,-0.5,19.5);
       BlipLocation[inte][0][type] = new TH1D(InteractionLabels[inte]+"XAxisBlipLocation" + CC1p0piLabels[type], ";X - Axis [cm]; Weighted Events",  NBinsVertexX , MinVertexX , MaxVertexX);
       BlipLocation[inte][1][type] = new TH1D(InteractionLabels[inte]+"YAxisBlipLocation" + CC1p0piLabels[type], ";Y - Axis [cm]; Weighted Events",  NBinsVertexY , MinVertexY , MaxVertexY);
       BlipLocation[inte][2][type] = new TH1D(InteractionLabels[inte]+"ZAxisBlipLocation" + CC1p0piLabels[type], ";Z - Axis [cm]; Weighted Events",  NBinsVertexZ , MinVertexZ , MaxVertexZ);

       BlipVertexDist[inte][type] = new TH1D(InteractionLabels[inte]+"BlipVertexDist" + CC1p0piLabels[type], ";Distance [cm] ; Weighted Events" , 20,0, 800);
       ProtonBlipDist[inte][type] = new TH1D(InteractionLabels[inte]+"ProtonBlipDist" + CC1p0piLabels[type], ";Distance [cm] ; Weighted Events" , 20,0, 800);

       //BlipCharge[inte][type] =  new TH1D(InteractionLabels[inte]+"BlipCharge" + CC1p0piLabels[type], ";Charge ; Weighted Events" , 20,0, 1000000);
       BlipProxTrkDist[inte][type] =  new TH1D(InteractionLabels[inte]+"BlipProxTrkDist" + CC1p0piLabels[type], ";Nearest Track Distance [cm]; Weighted Events" , 20,0, 300);
       BlipTime[inte][type] =  new TH1D(InteractionLabels[inte]+"BlipTime" + CC1p0piLabels[type], ";Time [ticks]; Weighted Events" , 40,0, 5000);
       //BlipLength[inte][type] =  new TH1D(InteractionLabels[inte]+"BlipLength" + CC1p0piLabels[type], ";Length [cm] ; Weighted Events" , 20,60, 80);
       BlipMaxWireSpan[inte][type] =  new TH1D(InteractionLabels[inte]+"BlipMaxWireSpan" + CC1p0piLabels[type], ";Max Wire Span ; Weighted Events" , 15,-0.5, 14.5);
       BlipEnergy[inte][type] =  new TH1D(InteractionLabels[inte]+"BlipEnergy" + CC1p0piLabels[type], ";Energy ; Weighted Events" , 20,0, 50);

       //ProtonBlipCharge[inte][type] =  new TH1D(InteractionLabels[inte]+"ProtonBlipCharge" + CC1p0piLabels[type], ";Charge ; Weighted Events" , 20,0, 1000000);
       ProtonBlipProxTrkDist[inte][type] =  new TH1D(InteractionLabels[inte]+"ProtonBlipProxTrkDist" + CC1p0piLabels[type], ";Nearest Track Distance [cm]; Weighted Events" , 20,0, 300);
       ProtonBlipTime[inte][type] =  new TH1D(InteractionLabels[inte]+"ProtonBlipTime" + CC1p0piLabels[type], ";Time [ticks]; Weighted Events" , 20,0, 5000);
       //ProtonBlipLength[inte][type] =  new TH1D(InteractionLabels[inte]+"ProtonBlipLength" + CC1p0piLabels[type], ";Length [cm] ; Weighted Events" , 20,0, 400);
       ProtonBlipMaxWireSpan[inte][type] =  new TH1D(InteractionLabels[inte]+"ProtonBlipMaxWireSpan" + CC1p0piLabels[type], ";Max Wire Span ; Weighted Events" , 15,-0.5, 14.5);
       ProtonBlipEnergy[inte][type] =  new TH1D(InteractionLabels[inte]+"ProtonBlipEnergy" + CC1p0piLabels[type], ";Energy ; Weighted Events" , 20,0, 50);



       BlipProxDistVsPMissMag[inte][type] = new TH2D(InteractionLabels[inte]+"BlipProxDistVsPMissMag" + CC1p0piLabels[type], ";Nearest Track Distance ;P_{miss} Magnitude [GeV/c] " , 20,0, 150, 20,0,1.5);
       BlipVertexDistVsPMissMag[inte][type] = new TH2D(InteractionLabels[inte]+"BlipVertexDistVsPMissMag" + CC1p0piLabels[type], ";Vertex Distance ;P_{miss} Magnitude [GeV/c] " , 20,0, 400, 20,0,1.5);
       BlipProxDistVsPMissDir[inte][type] = new TH2D(InteractionLabels[inte]+"BlipProxDistVsPMissDir" + CC1p0piLabels[type], ";Nearest Track Distance ;cos(#theta_{miss})" , 20,0, 150, 20,-1,1);
       BlipVertexDistVsPMissDir[inte][type] = new TH2D(InteractionLabels[inte]+"BlipVertexDistVsPMissDir" + CC1p0piLabels[type], ";Vertex Track Distance ;cos(#theta_{miss})" , 20,0, 400, 20,-1,1);
       BlipMultiplicityVsProxDist[inte][type] = new TH2D(InteractionLabels[inte]+"BlipMultiplicityVsProxDist" + CC1p0piLabels[type], ";Nearest Track Distance ;Blip Multiplicity" , 20,0,300, 20,30,150);
       BlipMultiplicityVsVertexDist[inte][type] = new TH2D(InteractionLabels[inte]+"BlipMultiplicityVsVertexDist" + CC1p0piLabels[type], ";Vertex Distance ;Blip Multiplicity" , 20,0, 800, 20,30,150);
       ProxDistVsNeutronMultiplicity[inte][type] = new TH2D(InteractionLabels[inte]+"ProxDistVsNeutronMultiplicity" + CC1p0piLabels[type], ";Nearest Track Distance ;Neutron Multiplicity" , 20,0, 300, 21,-0.5,20.5);
       VertexDistVsNeutronMultiplicity[inte][type] = new TH2D(InteractionLabels[inte]+"VertexDistVsNeutronMultiplicity" + CC1p0piLabels[type], ";Vertex Distance ;Neutron Multiplicity" , 20,0, 800, 21,-0.5,20.5);
       
       

       AssociatedBlipVertexDistVsProxTrkDist[inte][type] = new TH2D(InteractionLabels[inte]+"AssociatedBlipVertexDistVsProxTrkDist" + CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Vertex Distance [cm]" , 20,0, 300, 20,0 , 600);


       BlipVertexDistNeut[inte][type][NNeut-1] =  new TH1D(InteractionLabels[inte]+"BlipVertexDist"+Form("OverNeutron%d", NNeut-1) + CC1p0piLabels[type], ";Distance [cm] ; Weighted Events" , 20,0, 600);
       BlipProxTrkDistNeut[inte][type][NNeut-1] =  new TH1D(InteractionLabels[inte]+"BlipProxTrkDist" +Form("OverNeutron%d", NNeut-1)+ CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Weighted Events" , 20,0, 300);
       AssociatedBlipMultVsVertexDist[inte][type][NNeut-1] = new TH2D(InteractionLabels[inte]+"AssociatedBlipMultVsVertexDist" +Form("OverNeutron%d", NNeut-1)+ CC1p0piLabels[type], ";Distance [cm] ; Associated Blip Multiplicity" , 20,0, 600, 20,-0.5,19.5);
       
       AssociatedBlipMultVsProxTrkDist[inte][type][NNeut-1] = new TH2D(InteractionLabels[inte]+"AssociatedBlipMultVsProxTrkDist" +Form("OverNeutron%d", NNeut-1)+ CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Associated Blip Multiplicity" , 20,0, 300, 20,-0.5,19.5);


       BlipVertexDistNeut[inte][type][NNeut] =  new TH1D(InteractionLabels[inte]+"BlipVertexDist"+"AllNeutron" + CC1p0piLabels[type], ";Distance [cm] ; Weighted Events" , 20,0, 600);
       BlipProxTrkDistNeut[inte][type][NNeut] =  new TH1D(InteractionLabels[inte]+"BlipProxTrkDist" +"AllNeutron"+ CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Weighted Events" , 20,0, 300);
       AssociatedBlipMultVsVertexDist[inte][type][NNeut] = new TH2D(InteractionLabels[inte]+"AssociatedBlipMultVsVertexDist" +"AllNeutron"+ CC1p0piLabels[type], ";Distance [cm] ; Associated Blip Multiplicity" , 20,0, 600, 20,-0.5,19.5);
       
       AssociatedBlipMultVsProxTrkDist[inte][type][NNeut] = new TH2D(InteractionLabels[inte]+"AssociatedBlipMultVsProxTrkDist" +"AllNeutron"+ CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Associated Blip Multiplicity" , 20,0, 300, 20,-0.5,19.5);

       
       
       for (int neut=0; neut<NNeut-1; neut++){
	 BlipVertexDistNeut[inte][type][neut] =  new TH1D(InteractionLabels[inte]+"BlipVertexDist" +Form("Neutron%d", neut)+ CC1p0piLabels[type], ";Distance [cm] ; Weighted Events" , 20,0, 600);
	 BlipProxTrkDistNeut[inte][type][neut] =  new TH1D(InteractionLabels[inte]+"BlipProxTrkDist"+Form("Neutron%d", neut) + CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Weighted Events" , 20,0, 300);
	 
	 AssociatedBlipMultVsVertexDist[inte][type][neut] = new TH2D(InteractionLabels[inte]+"AssociatedBlipMultVsVertexDist" +Form("Neutron%d", neut)+ CC1p0piLabels[type], ";Distance [cm] ; Associated Blip Multiplicity" , 20,0, 600, 20,-0.5,19.5);
	 AssociatedBlipMultVsProxTrkDist[inte][type][neut] = new TH2D(InteractionLabels[inte]+"AssociatedBlipMultVsProxTrkDist" +Form("Neutron%d", neut)+ CC1p0piLabels[type], ";Nearest Track Distance [cm] ; Associated Blip Multiplicity" , 20,0, 300, 20,-0.5,19.5);
	 
       } // End of loop over the neutron multiplicity


     } // End of For loop over the type of signal
   } // End of the loop over the interaction processes
   











   int counter =0,  trueCounter=0, noBlipCounter=0;
   // Loop over the events
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
     //for (Long64_t jentry=0; jentry<20000;jentry++) {


      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if (jentry%1000 == 0) std::cout << jentry/1000 << " k " << std::setprecision(2) << double(jentry)/nentries*100. << " %"<< std::endl;



      // MC weight to scale events to data pot
      //if (fabs(weightSplineTimesTune) != weightSplineTimesTune){ continue;}
      if (weightSplineTimesTune <= 0 || weightSplineTimesTune > 30) { continue; } // bug fix weight
      double event_weight = (data_pot / mc_pot_highStat) * weightSplineTimesTune;



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

	  if ( MCParticlePdg == NeutronPdg ) 
	    { TrueNeutronCounter++; VectorTrueNeutronIndex.push_back(WhichMCParticle); }

	  if ( fabs(MCParticlePdg) == AbsChargedPionPdg && MCParticleMomentum >= ChargedPionMomentumThres ) 
	    { TrueChargedPionCounter++; }

	  if (MCParticlePdg == NeutralPionPdg) { TruePi0Counter++; }

	} // End of the demand stable final state particles and primary interactions
      } // end of the loop over the MCParticles




      // CC1p0pi signal events

      bool CC1p0pi = false;
      if (TrueMuonCounter == 1 && TrueProtonCounter == 1 && TrueChargedPionCounter == 0 && TruePi0Counter == 0 && TrueHeavierMesonCounter == 0) {
	//trueCounter++;
	CC1p0pi = true; 
      }
      // Unlike the truth level study, we just flag events as signal or not but we don't continue. Do Not: if (!CC1p0pi) { continue; }
      
      
      // Max: check for true CC1p0pi events inside fiducial volume of interest
      // if true vertex outside, set CC1p0pi = false 
      // but do NOT "continue"!
      
      TVector3 trueVertex;
      if(CC1p0pi == true) {trueVertex.SetXYZ(mc_vx->at(VectorTrueMuonIndex[0]), mc_vy->at(VectorTrueMuonIndex[0]), mc_vz->at(VectorTrueMuonIndex[0]) );}
      else {trueVertex.SetXYZ(-999., -999., -999.); }

      
      if ( inFVVector(trueVertex) == false) { CC1p0pi = false; }
      


      //else {std::cout << "trueVertexInFV" << std::endl;}
      // else {trueCounter++;}



      
      // Max: Need if-statement for true CC1p0pi events
      // else: fill in the vectors with (-999.,-999.,-999.)
      // and still calculate the quantities below
      
      TVector3 TrueMuonMomentumV;
      TVector3 TrueProtonMomentumV;
      double TrueMuonEnergy, TrueProtonEnergy;
      if (CC1p0pi == false){
	TrueMuonMomentumV.SetXYZ(-999.,-999.,-999.);
	TrueProtonMomentumV.SetXYZ(-999.,-999.,-999.);
	TrueMuonEnergy = -999.0;
	TrueProtonEnergy =-999.0;
      }
      
      else{
	TrueMuonMomentumV.SetXYZ(mc_px->at(VectorTrueMuonIndex[0]), mc_py->at(VectorTrueMuonIndex[0]), mc_pz->at(VectorTrueMuonIndex[0]));
	TrueProtonMomentumV.SetXYZ(mc_px->at(VectorTrueProtonIndex[0]), mc_py->at(VectorTrueProtonIndex[0]), mc_pz->at(VectorTrueProtonIndex[0]));	
	TrueMuonEnergy = mc_E->at(VectorTrueMuonIndex[0]);
	TrueProtonEnergy = mc_E->at(VectorTrueProtonIndex[0]);
      }
      
      // Max: Define muon / proton vectors, deltapt, delta alphaT & pmiss
      TVector3 TrueMuonTransverseMom(TrueMuonMomentumV.X(), TrueMuonMomentumV.Y(), 0);
      TVector3 TrueProtonTransverseMom(TrueProtonMomentumV.X(), TrueProtonMomentumV.Y(), 0);
      TVector3 TrueDeltaPtV = TrueMuonTransverseMom +TrueProtonTransverseMom;
      double TrueDeltaPtMag = TrueDeltaPtV.Mag();
      
      double TrueMuonCosTheta = TrueMuonMomentumV.CosTheta();
      double TrueProtonCosTheta = TrueProtonMomentumV.CosTheta();
      
      double TruecosAlphaT = (-1.0*(TrueMuonTransverseMom.Dot(TrueDeltaPtV))) / (TrueMuonTransverseMom.Mag() * TrueDeltaPtMag);
      double TrueDeltaAlphaT = TMath::ACos(TruecosAlphaT)*( 180.0 / 3.1415926 );
      
      double TrueprotonKE = TrueProtonEnergy - ProtonMass_GeV;
      double TrueCalEne = TrueMuonEnergy + TrueprotonKE + 0.04;
      TVector3 TrueNuMomentumVector(0,0,TrueCalEne);
      
      TVector3 TruePMissing = TrueNuMomentumVector - TrueProtonMomentumV - TrueMuonMomentumV;
      // cout << "Muon Kinematics:  X = " << TrueMuonMomentumV.X() << "  Y = " <<  TrueMuonMomentumV.Y() << "  Z = " << TrueMuonMomentumV.Z() << endl;
      // cout << "Proton Kinematics:  X = " << TrueProtonMomentumV.X() << "  Y = " <<  TrueProtonMomentumV.Y() << "  Z = " << TrueProtonMomentumV.Z() << endl;
      // cout << "PMiss Kinematics:  X = " << TruePMissing.X() << "  Y = " <<  TruePMissing.Y() << "  Z = " << TruePMissing.Z() << endl;
      double TruePMissingMag = TruePMissing.Mag();
      double TruePMissingDir = TruePMissing.CosTheta();
      
      // Classify events based on interaction

      int genie_mode = -1;

      if (interaction == 0) { genie_mode = 1; } // QE
      else if (interaction == 10) { genie_mode = 2; } // MEC
      else if (interaction == 1) { genie_mode = 3; } // RES
      else if (interaction == 2) { genie_mode = 4; } // DIS
      else { genie_mode = 5; } // COH / other


      
      TrueVertexXPlot[0][0]->Fill(trueVertex.X(),event_weight);
      TrueVertexYPlot[0][0]->Fill(trueVertex.Y(), event_weight);
      TrueVertexZPlot[0][0]->Fill(trueVertex.Z(), event_weight);
      TruePMissMagPlot[0][0]->Fill(TruePMissingMag, event_weight);
      TruePMissDirPlot[0][0]->Fill(TruePMissingDir, event_weight);
      TrueMuonEnergyPlot[0][0]->Fill(TrueMuonEnergy, event_weight);


      if (CC1p0pi == true){ 
      	TrueVertexXPlot[0][1]->Fill(trueVertex.X(), event_weight);
      	TrueVertexYPlot[0][1]->Fill(trueVertex.Y(), event_weight);
	TruePMissMagPlot[0][1]->Fill(TruePMissingMag, event_weight);
	TruePMissDirPlot[0][1]->Fill(TruePMissingDir, event_weight);
      	TrueVertexZPlot[0][1]->Fill(trueVertex.Z(), event_weight);	  
	TrueMuonEnergyPlot[0][1]->Fill(TrueMuonEnergy, event_weight);

      }
	
      else {   
      	TrueVertexXPlot[0][2]->Fill(trueVertex.X(), event_weight);
      	TrueVertexYPlot[0][2]->Fill(trueVertex.Y(), event_weight);
      	TrueVertexZPlot[0][2]->Fill(trueVertex.Z(), event_weight);
	TruePMissMagPlot[0][2]->Fill(TruePMissingMag, event_weight);
	TruePMissDirPlot[0][2]->Fill(TruePMissingDir, event_weight);
	TrueMuonEnergyPlot[0][2]->Fill(TrueMuonEnergy, event_weight);

      }


      
      TrueVertexXPlot[genie_mode][0]->Fill(trueVertex.X(), event_weight);
      TrueVertexYPlot[genie_mode][0]->Fill(trueVertex.Y(), event_weight);
      TrueVertexZPlot[genie_mode][0]->Fill(trueVertex.Z(), event_weight);
      TruePMissMagPlot[genie_mode][0]->Fill(TruePMissingMag, event_weight);
      TruePMissDirPlot[genie_mode][0]->Fill(TruePMissingDir, event_weight);
      TrueMuonEnergyPlot[genie_mode][0]->Fill(TrueMuonEnergy, event_weight);

      if (CC1p0pi == true){ 
      	TrueVertexXPlot[genie_mode][1]->Fill(trueVertex.X(), event_weight);
      	TrueVertexYPlot[genie_mode][1]->Fill(trueVertex.Y(), event_weight);
      	TrueVertexZPlot[genie_mode][1]->Fill(trueVertex.Z(), event_weight);	  
	TruePMissMagPlot[genie_mode][1]->Fill(TruePMissingMag, event_weight);
	TruePMissDirPlot[genie_mode][1]->Fill(TruePMissingDir, event_weight);
	TrueMuonEnergyPlot[genie_mode][1]->Fill(TrueMuonEnergy, event_weight);

      }
	
      else {   
      	TrueVertexXPlot[genie_mode][2]->Fill(trueVertex.X(), event_weight);
      	TrueVertexYPlot[genie_mode][2]->Fill(trueVertex.Y(), event_weight);
      	TrueVertexZPlot[genie_mode][2]->Fill(trueVertex.Z(), event_weight);
	TruePMissMagPlot[genie_mode][2]->Fill(TruePMissingMag, event_weight);
	TruePMissDirPlot[genie_mode][2]->Fill(TruePMissingDir, event_weight);
	TrueMuonEnergyPlot[genie_mode][2]->Fill(TrueMuonEnergy, event_weight);

      }
      




      // ------------------------  Reco level  ------------------------------- //
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
      TrackMultiplicity[0]->Fill(reco_track_count, event_weight);
      ShowerMultiplicity[0]->Fill(reco_shower_count, event_weight);
      


      // Requirement for exactly 2 tracks and 0 showers for a reco CC1p0pi selection 
      if (reco_shower_count != 0) { continue; }
      if (reco_track_count != 2) { continue; }



      // Requirement for exactly one neutrino slice
      if (nslice != 1) { continue; }





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
      }
      
      else { 
	CandidateMuonIndex = CandidateIndex.at(1); 
	CandidateProtonIndex = CandidateIndex.at(0); 
      }
      
      //cout <<"Muon Index: " <<CandidateMuonIndex << "  Proton Index: " << CandidateProtonIndex << endl;

      // MuonPdg in this case / Pandora = Track-like object
      if (pfpdg->at(CandidateMuonIndex) != MuonPdg) { continue; }
      if (pfpdg->at(CandidateProtonIndex) != MuonPdg) { continue; }



      // Reconstructed Vertex
      TVector3 reco_vertex(reco_nu_vtx_sce_x,reco_nu_vtx_sce_y,reco_nu_vtx_sce_z);

      // Max: demand that the reco vertex is contained in the fiducial volume
      if (inFVVector(reco_vertex)==false) { continue; }
      // else {std::cout << "recoVertexInFV" << std::endl;}

      // Muon kinematics

      double CandidateMuonTrackTheta = trk_theta_v->at(CandidateMuonIndex); // rad
      double CandidateMuonTrackPhi = trk_phi_v->at(CandidateMuonIndex); // rad
      double CandidateMuMom = trk_range_muon_mom_v->at(CandidateMuonIndex); // GeV/c
      
      double CandidateMuE_GeV = TMath::Sqrt( TMath::Power(CandidateMuMom ,2) + TMath::Power(MuonMass_GeV ,2) );
      double CandidateMuKE_GeV =  CandidateMuE_GeV - MuonMass_GeV;

      
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
      
      TVector3 CandidateMuon(-1,-1,-1);
      CandidateMuon.SetMag(CandidateMuMom); // GeV/c
      CandidateMuon.SetTheta(CandidateMuonTrackTheta); // rad
      CandidateMuon.SetPhi(CandidateMuonTrackPhi); // rad



      // Max: declare the proton kinematics
      double CandidateProtonTrackTheta = trk_theta_v->at(CandidateProtonIndex); // rad
      double CandidateProtonTrackPhi = trk_phi_v->at(CandidateProtonIndex); // rad
      
      double CandidatePKE_GeV = trk_energy_proton_v->at(CandidateProtonIndex); // GeV // Watch out, kinetic energy not energy
      double CandidatePE_GeV = CandidatePKE_GeV + ProtonMass_GeV; // GeV 
      double CandidateProtonMom = TMath::Sqrt( TMath::Power(CandidatePE_GeV,2.) - TMath::Power(ProtonMass_GeV,2.)); // GeV/c

      double ProtonTrackStartX = trk_sce_start_x_v->at(CandidateProtonIndex);
      double ProtonTrackStartY = trk_sce_start_y_v->at(CandidateProtonIndex);
      double ProtonTrackStartZ = trk_sce_start_z_v->at(CandidateProtonIndex);
      
      double ProtonTrackEndX = trk_sce_end_x_v->at(CandidateProtonIndex);
      double ProtonTrackEndY = trk_sce_end_y_v->at(CandidateProtonIndex);
      double ProtonTrackEndZ = trk_sce_end_z_v->at(CandidateProtonIndex);
      
      TVector3 CandidateProtonTrackStart(ProtonTrackStartX,ProtonTrackStartY,ProtonTrackStartZ);
      TVector3 CandidateProtonTrackEnd(ProtonTrackEndX,ProtonTrackEndY,ProtonTrackEndZ);
      bool CandidateProtonTrackStartContainment = inFVVector(CandidateProtonTrackStart);
      bool CandidateProtonTrackEndContainment = inFVVector(CandidateProtonTrackEnd);

      TVector3 CandidateProton(-1,-1,-1);
      CandidateProton.SetMag(CandidateProtonMom); // GeV/c
      CandidateProton.SetTheta(CandidateProtonTrackTheta); // rad
      CandidateProton.SetPhi(CandidateProtonTrackPhi); // rad


      if (trk_llr_pid_score_v->at( CandidateProtonIndex ) > ProtonLLRPIDScore ) {continue;}




      // Momenta ranges for the muon candidate

      if (CandidateMuMom < ArrayNBinsMuonMomentum[0]) { continue; }
      if (CandidateMuMom > ArrayNBinsMuonMomentum[NBinsMuonMomentum]) { continue; }

      // Max: apply the momentum the proton candidate
      if (CandidateProtonMom < ArrayNBinsProtonMomentum[0]) { continue; }
      if (CandidateProtonMom > ArrayNBinsProtonMomentum[NBinsProtonMomentum]) { continue; }
      //else {std::cout << "Momentum In Range" << std::endl;}


      // Fully contained muons
      if (!CandidateMuonTrackStartContainment) { continue; }
      if (!CandidateMuonTrackEndContainment) { continue; }

      // Max: do the same for proton candidate
      if (!CandidateProtonTrackStartContainment) { continue; }
      if (!CandidateProtonTrackEndContainment) { continue; }
      //else {std::cout << "Tracks Fully Contained" << std::endl;}



      // Start-To-Start & End-To-End Distance to avoid flipped tracks
      double LocalStartToStartDistance = (CandidateMuonTrackStart - CandidateProtonTrackStart).Mag();
      double LocalEndToEndDistance = (CandidateMuonTrackEnd - CandidateProtonTrackEnd).Mag();
      if (LocalStartToStartDistance >= LocalEndToEndDistance) { continue; }
      


      // Vertex - Track & relative differences
      double CandidateMuStartVertexMag = (reco_vertex - CandidateMuonTrackStart).Mag();
      double CandidatePStartVertexMag = (reco_vertex - CandidateProtonTrackStart).Mag();
      double CandidateMuEndVertexMag = (reco_vertex - CandidateMuonTrackEnd).Mag();
      double CandidatePEndVertexMag = (reco_vertex - CandidateProtonTrackEnd).Mag();
      if (CandidateMuStartVertexMag > CandidateMuEndVertexMag) { continue; }

      // Max: do the same for the proton
      if (CandidatePStartVertexMag > CandidatePEndVertexMag) { continue; }



      // Max: declare delta pt & delta alphat for reco CC1p0pi candidates
      TVector3 recoMuonTransvMom(CandidateMuon.X(), CandidateMuon.Y(), 0);
      TVector3 recoProtonTransvMom(CandidateProton.X(), CandidateProton.Y(), 0);
      TVector3 recoNuMom(0,0,  CandidateMuE_GeV + CandidatePKE_GeV + 0.04);
      TVector3 recoPMissing = recoNuMom - CandidateProton -CandidateMuon;


      double recoDeltaPT = (recoMuonTransvMom + recoProtonTransvMom).Mag();
      double recoDeltaAlphaT = TMath::ACos( (-1.0 * (recoMuonTransvMom.Dot(recoMuonTransvMom+recoProtonTransvMom)) ) / (recoMuonTransvMom.Mag() * recoDeltaPT) ) * (180.0 / 3.1415926);
      double recoPMissingMag = recoPMissing.Mag();
      double recoPMissingDir = recoPMissing.CosTheta();
      

      //PMissing Magnitude Cuts
      if (PMissingCut == true){
	if (recoPMissingMag > 0.3) { continue; }
      } //End of application of cut


      // Max: plot the x/y/z component of the reco vertex for the CC1p0pi candidates
      RecoVertexXPlot[0][0]->Fill(reco_vertex.X()) ;
      RecoVertexYPlot[0][0]->Fill(reco_vertex.Y()) ;
      RecoVertexZPlot[0][0]->Fill(reco_vertex.Z()) ;      

      RecoVertexXPlot[genie_mode][0]->Fill(reco_vertex.X()) ;
      RecoVertexYPlot[genie_mode][0]->Fill(reco_vertex.Y()) ;
      RecoVertexZPlot[genie_mode][0]->Fill(reco_vertex.Z()) ;      
      
      RecoVertexXPlot[0][1]->Fill(reco_vertex.X()) ;
      RecoVertexYPlot[0][1]->Fill(reco_vertex.Y()) ;
      RecoVertexZPlot[0][1]->Fill(reco_vertex.Z()) ;      

      RecoVertexXPlot[genie_mode][1]->Fill(reco_vertex.X()) ;
      RecoVertexYPlot[genie_mode][1]->Fill(reco_vertex.Y()) ;
      RecoVertexZPlot[genie_mode][1]->Fill(reco_vertex.Z()) ;      
      
      RecoVertexXPlot[0][2]->Fill(reco_vertex.X()) ;
      RecoVertexYPlot[0][2]->Fill(reco_vertex.Y()) ;
      RecoVertexZPlot[0][2]->Fill(reco_vertex.Z()) ;      

      RecoVertexXPlot[genie_mode][2]->Fill(reco_vertex.X()) ;
      RecoVertexYPlot[genie_mode][2]->Fill(reco_vertex.Y()) ;
      RecoVertexZPlot[genie_mode][2]->Fill(reco_vertex.Z()) ;      
      


      








      // ------------------------- Backtracking --------------------------- //

      // Backtracking: matching the reconstructed objects to true ones
      // Not the same as the truth-level study
      // where things don't need to be reconstructed   

      // Muon candidate backtracking
      double CandidateMuonPx = backtracked_px->at(CandidateMuonIndex);
      double CandidateMuonPy = backtracked_py->at(CandidateMuonIndex);
      double CandidateMuonPz = backtracked_pz->at(CandidateMuonIndex);
      TVector3 TrueCandidateMuonP(CandidateMuonPx,CandidateMuonPy,CandidateMuonPz);
      
      double TrueCandidateMuonTrackPhi = TrueCandidateMuonP.Phi(); // rad
      double TrueCandidateMuonTrackTheta = TrueCandidateMuonP.Theta(); // rad
      double TrueCandidateMuonPMag = TrueCandidateMuonP.Mag();
      
      TVector3 True_CandidateMuon(-1,-1,-1);
      True_CandidateMuon.SetMag(TrueCandidateMuonPMag);
      True_CandidateMuon.SetTheta(TrueCandidateMuonTrackTheta);
      True_CandidateMuon.SetPhi(TrueCandidateMuonTrackPhi);
      
      checkMomentumSame[0]->Fill(TrueCandidateMuonP.X(),True_CandidateMuon.X() , event_weight);
      checkMomentumSame[1]->Fill(TrueCandidateMuonP.Y(),True_CandidateMuon.Y() , event_weight);
      checkMomentumSame[2]->Fill(TrueCandidateMuonP.Z(),True_CandidateMuon.Z() , event_weight);

   
      
      // Max: do the same for the proton backtracking

      double CandidateProtonPx = backtracked_px->at(CandidateProtonIndex);
      double CandidateProtonPy = backtracked_py->at(CandidateProtonIndex);
      double CandidateProtonPz = backtracked_pz->at(CandidateProtonIndex);
      TVector3 TrueCandidateProtonP(CandidateProtonPx,CandidateProtonPy,CandidateProtonPz);
      
      double TrueCandidateProtonTrackPhi = TrueCandidateProtonP.Phi(); // rad
      double TrueCandidateProtonTrackTheta = TrueCandidateProtonP.Theta(); // rad
      double TrueCandidateProtonPMag = TrueCandidateProtonP.Mag();
      
      TVector3 True_CandidateProton(-1,-1,-1);
      True_CandidateProton.SetMag(TrueCandidateProtonPMag);
      True_CandidateProton.SetTheta(TrueCandidateProtonTrackTheta);
      True_CandidateProton.SetPhi(TrueCandidateProtonTrackPhi);
      

      

      // If it is a true CC1p0pi event
      // Check that it is a backtracked CC1p0pi event as well
      if (CC1p0pi==true) {	
	//put back in
	if ( !(backtracked_pdg->at(CandidateMuonIndex) == MuonPdg && TrueCandidateMuonPMag > ArrayNBinsMuonMomentum[0] && TrueCandidateMuonPMag < ArrayNBinsMuonMomentum[NBinsMuonMomentum]) ) {
	  CC1p0pi = false;
	}

	
	if ( !(backtracked_pdg->at(CandidateProtonIndex) == ProtonPdg && TrueCandidateProtonPMag > ArrayNBinsProtonMomentum[0] && TrueCandidateProtonPMag < ArrayNBinsProtonMomentum[NBinsProtonMomentum]) ) {
	  CC1p0pi = false;
	}
      }



      // Max: define delta pt, delta alphat, and pmiss for backtracked vectors
      TVector3 backtrackedMuonTransVect(TrueCandidateMuonP.X(), TrueCandidateMuonP.Y(), 0);
      TVector3 backtrackedProtonTransVect(TrueCandidateProtonP.X(), TrueCandidateProtonP.Y(), 0);
      double backtrackedDeltaPt = (backtrackedMuonTransVect + backtrackedProtonTransVect).Mag();

      double backtrackedDeltaAlphaT = TMath::ACos(-1.0 *(backtrackedMuonTransVect.Dot(backtrackedMuonTransVect + backtrackedProtonTransVect)) / (backtrackedMuonTransVect.Mag() *(backtrackedMuonTransVect + backtrackedProtonTransVect).Mag() ) ) * (180.0 / 3.1415926) ;
      

      double backtrackedProtonKE = TMath::Sqrt( TMath::Power(TrueCandidateProtonP.Mag(),2) + TMath::Power(ProtonMass_GeV, 2)) - ProtonMass_GeV; // backtracked_e->at(CandidateProtonIndex);
      double backtrackedMuonEnergy = TMath::Sqrt( TMath::Power(TrueCandidateMuonP.Mag(),2) + TMath::Power(MuonMass_GeV, 2));
      TVector3 backtrackedNumuMomVect(0,0, backtrackedMuonEnergy + backtrackedProtonKE + 0.04);
      TVector3 backtrackedPMiss = backtrackedNumuMomVect - TrueCandidateProtonP - TrueCandidateMuonP;
      double backtrackedPMissMag = backtrackedPMiss.Mag();
      double backtrackedPMissDir = backtrackedPMiss.CosTheta();
      

      // Fill in the plots
      // All events, 2nd index = 0 
      RecoMuonCosThetaPlot[0][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
      RecoDeltaPtPlot[0][0]->Fill(recoDeltaPT,event_weight);
      RecoDeltaAlphaTPlot[0][0]->Fill(recoDeltaAlphaT,event_weight);
      RecoPMissMagPlot[0][0]->Fill(recoPMissingMag,event_weight);
      RecoPMissDirPlot[0][0]->Fill(recoPMissingDir,event_weight);
      BacktrackedvsRecoDeltaPtPlot[0][0]->Fill(backtrackedDeltaPt ,recoDeltaPT , event_weight);
      BacktrackedvsRecoDeltaAlphaTPlot[0][0]->Fill(backtrackedDeltaAlphaT ,recoDeltaAlphaT , event_weight);
      BacktrackedDeltaPtPlot[0][0]->Fill(backtrackedDeltaPt , event_weight);
      BacktrackedDeltaAlphaTPlot[0][0]->Fill(backtrackedDeltaAlphaT , event_weight);
      BacktrackedPMissMagPlot[0][0]->Fill(backtrackedPMissMag, event_weight);
      BacktrackedPMissDirPlot[0][0]->Fill(backtrackedPMissDir, event_weight);
      BacktrackedvsRecoPMissMagPlot[0][0]->Fill(backtrackedPMissMag, recoPMissingMag, event_weight);
      BacktrackedvsRecoPMissDirPlot[0][0]->Fill(backtrackedPMissDir, recoPMissingDir ,event_weight);
      
      BacktrackedvsRecoPMissMagPlot[genie_mode][0]->Fill(backtrackedPMissMag, recoPMissingMag, event_weight);
      BacktrackedvsRecoPMissDirPlot[genie_mode][0]->Fill(backtrackedPMissDir, recoPMissingDir ,event_weight);	
      RecoMuonCosThetaPlot[genie_mode][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
      RecoDeltaPtPlot[genie_mode][0]->Fill(recoDeltaPT,event_weight);
      RecoDeltaAlphaTPlot[genie_mode][0]->Fill(recoDeltaAlphaT,event_weight);
      RecoPMissMagPlot[genie_mode][0]->Fill(recoPMissingMag,event_weight);
      RecoPMissDirPlot[genie_mode][0]->Fill(recoPMissingDir,event_weight);
      BacktrackedvsRecoDeltaPtPlot[genie_mode][0]->Fill(backtrackedDeltaPt ,recoDeltaPT , event_weight);
      BacktrackedvsRecoDeltaAlphaTPlot[genie_mode][0]->Fill(backtrackedDeltaAlphaT ,recoDeltaAlphaT , event_weight);
      BacktrackedDeltaPtPlot[genie_mode][0]->Fill(backtrackedDeltaPt , event_weight);
      BacktrackedDeltaAlphaTPlot[genie_mode][0]->Fill(backtrackedDeltaAlphaT , event_weight);
      BacktrackedPMissMagPlot[genie_mode][0]->Fill(backtrackedPMissMag, event_weight);
      BacktrackedPMissDirPlot[genie_mode][0]->Fill(backtrackedPMissDir, event_weight);

      TrueNeutronMultiplicityPlot[0][0]->Fill(TrueNeutronCounter, event_weight);
      TrueNeutronMultiplicityPlot[genie_mode][0]->Fill(TrueNeutronCounter, event_weight);
      

      // CC1p0pi events, 2nd index = 1
      if (CC1p0pi){
	RecoMuonCosThetaPlot[0][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	RecoDeltaPtPlot[0][1]->Fill(recoDeltaPT,event_weight);
	RecoDeltaAlphaTPlot[0][1]->Fill(recoDeltaAlphaT,event_weight);
	RecoPMissMagPlot[0][1]->Fill(recoPMissingMag,event_weight);
	RecoPMissDirPlot[0][1]->Fill(recoPMissingDir,event_weight);
	BacktrackedvsRecoDeltaPtPlot[0][1]->Fill(backtrackedDeltaPt ,recoDeltaPT , event_weight);
	BacktrackedvsRecoDeltaAlphaTPlot[0][1]->Fill(backtrackedDeltaAlphaT ,recoDeltaAlphaT , event_weight);
	BacktrackedDeltaPtPlot[0][1]->Fill(backtrackedDeltaPt , event_weight);
	BacktrackedDeltaAlphaTPlot[0][1]->Fill(backtrackedDeltaAlphaT , event_weight);
	BacktrackedPMissMagPlot[0][1]->Fill(backtrackedPMissMag, event_weight);
	BacktrackedPMissDirPlot[0][1]->Fill(backtrackedPMissDir, event_weight);
	BacktrackedvsRecoPMissMagPlot[0][1]->Fill(backtrackedPMissMag, recoPMissingMag, event_weight);
	BacktrackedvsRecoPMissDirPlot[0][1]->Fill(backtrackedPMissDir, recoPMissingDir ,event_weight);
	
	BacktrackedvsRecoPMissMagPlot[genie_mode][1]->Fill(backtrackedPMissMag, recoPMissingMag, event_weight);
	BacktrackedvsRecoPMissDirPlot[genie_mode][1]->Fill(backtrackedPMissDir, recoPMissingDir ,event_weight);	
	RecoMuonCosThetaPlot[genie_mode][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	RecoDeltaPtPlot[genie_mode][1]->Fill(recoDeltaPT,event_weight);
	RecoDeltaAlphaTPlot[genie_mode][1]->Fill(recoDeltaAlphaT,event_weight);
	RecoPMissMagPlot[genie_mode][1]->Fill(recoPMissingMag,event_weight);
	RecoPMissDirPlot[genie_mode][1]->Fill(recoPMissingDir,event_weight);
	BacktrackedvsRecoDeltaPtPlot[genie_mode][1]->Fill(backtrackedDeltaPt ,recoDeltaPT , event_weight);
	BacktrackedvsRecoDeltaAlphaTPlot[genie_mode][1]->Fill(backtrackedDeltaAlphaT ,recoDeltaAlphaT , event_weight);
	BacktrackedDeltaPtPlot[genie_mode][1]->Fill(backtrackedDeltaPt , event_weight);
	BacktrackedDeltaAlphaTPlot[genie_mode][1]->Fill(backtrackedDeltaAlphaT , event_weight);
	BacktrackedPMissMagPlot[genie_mode][1]->Fill(backtrackedPMissMag, event_weight);
	BacktrackedPMissDirPlot[genie_mode][1]->Fill(backtrackedPMissDir, event_weight);

	TrueNeutronMultiplicityPlot[0][1]->Fill(TrueNeutronCounter, event_weight);
	TrueNeutronMultiplicityPlot[genie_mode][1]->Fill(TrueNeutronCounter, event_weight);
	

      } //End of if statement for plotting CC1p0pi events
      
      //non-CC1p0pi events, 2nd index = 2
      else {
	RecoMuonCosThetaPlot[0][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	RecoDeltaPtPlot[0][2]->Fill(recoDeltaPT,event_weight);
	RecoDeltaAlphaTPlot[0][2]->Fill(recoDeltaAlphaT,event_weight);
	RecoPMissMagPlot[0][2]->Fill(recoPMissingMag,event_weight);
	RecoPMissDirPlot[0][2]->Fill(recoPMissingDir,event_weight);
	BacktrackedvsRecoDeltaPtPlot[0][2]->Fill(backtrackedDeltaPt ,recoDeltaPT , event_weight);
	BacktrackedvsRecoDeltaAlphaTPlot[0][2]->Fill(backtrackedDeltaAlphaT ,recoDeltaAlphaT , event_weight);
	BacktrackedvsRecoPMissMagPlot[0][2]->Fill(backtrackedPMissMag, recoPMissingMag, event_weight);
	BacktrackedvsRecoPMissDirPlot[0][2]->Fill(backtrackedPMissDir, recoPMissingDir ,event_weight);
	BacktrackedDeltaPtPlot[0][2]->Fill(backtrackedDeltaPt , event_weight);
	BacktrackedDeltaAlphaTPlot[0][2]->Fill(backtrackedDeltaAlphaT , event_weight);
	BacktrackedPMissMagPlot[0][2]->Fill(backtrackedPMissMag, event_weight);
	BacktrackedPMissDirPlot[0][2]->Fill(backtrackedPMissDir, event_weight);

	RecoMuonCosThetaPlot[genie_mode][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	RecoDeltaPtPlot[genie_mode][2]->Fill(recoDeltaPT,event_weight);
	RecoDeltaAlphaTPlot[genie_mode][2]->Fill(recoDeltaAlphaT,event_weight);
	RecoPMissMagPlot[genie_mode][2]->Fill(recoPMissingMag,event_weight);
	RecoPMissDirPlot[genie_mode][2]->Fill(recoPMissingDir,event_weight);
	BacktrackedvsRecoDeltaPtPlot[genie_mode][2]->Fill(backtrackedDeltaPt ,recoDeltaPT , event_weight);
	BacktrackedvsRecoPMissMagPlot[genie_mode][2]->Fill(backtrackedPMissMag, recoPMissingMag, event_weight);
	BacktrackedvsRecoPMissDirPlot[genie_mode][2]->Fill(backtrackedPMissDir, recoPMissingDir ,event_weight);
	BacktrackedvsRecoDeltaAlphaTPlot[genie_mode][2]->Fill(backtrackedDeltaAlphaT ,recoDeltaAlphaT , event_weight);
	BacktrackedDeltaPtPlot[genie_mode][2]->Fill(backtrackedDeltaPt , event_weight);
	BacktrackedDeltaAlphaTPlot[genie_mode][2]->Fill(backtrackedDeltaAlphaT , event_weight);      
	BacktrackedPMissMagPlot[genie_mode][2]->Fill(backtrackedPMissMag, event_weight);
	BacktrackedPMissDirPlot[genie_mode][2]->Fill(backtrackedPMissDir, event_weight);

	TrueNeutronMultiplicityPlot[0][2]->Fill(TrueNeutronCounter, event_weight);
	TrueNeutronMultiplicityPlot[genie_mode][2]->Fill(TrueNeutronCounter, event_weight);
	
      } //End of else statement for plotting nonCC1p0pi events

      




      // Loop over the blips
      // if (blip_ID->size() ==0) { cout << "No Blip" << endl; noBlipCounter++; }
      int numBlips = blip_ID->size();
      
      BlipMultiplicity[0][0]->Fill(numBlips, event_weight);
      BlipMultiplicity[genie_mode][0]->Fill(numBlips, event_weight);
      bool protonBlip =false;
      bool regularBlip = false;
      if (CC1p0pi == true){
	BlipMultiplicity[0][1]->Fill(numBlips, event_weight);
	BlipMultiplicity[genie_mode][1]->Fill(numBlips, event_weight);
      }
      
      else {
	BlipMultiplicity[0][2]->Fill(numBlips, event_weight); 
	BlipMultiplicity[genie_mode][2]->Fill(numBlips, event_weight);	      
      }


      int blipCounter =0;
      vector<int> AssociatedBlipID;
      vector<double> blipVertex, blipProxTrk;
      bool blipVertexTooSmall = false;
      if (blip_ID->size() >0 ){
		
	for (int blip =0; blip < TMath::Abs(numBlips); blip++){
	  if (blip_pdg->at(blip) !=-9){

	    TVector3 blip_vertex(blip_X->at(blip),blip_Y->at(blip),blip_Z->at(blip));
	    double blipVertexDistance = (blip_vertex - reco_vertex).Mag();
	    regularBlip =true;
	    //cout << " Index: " << blip_ProxTrkID->at(blip) << endl;
	    
	    if (blip_ProxTrkID->at(blip) == CandidateMuonIndex || blip_ProxTrkID->at(blip) == CandidateProtonIndex){
	      //cout << "Blip in Proximity" << endl;
	      //cout << "ID:  "  << endl;
	      blipVertexTooSmall = false;
	      blipCounter++;
	      AssociatedBlipID.push_back(blipCounter-1);
	      blipVertex.push_back(blipVertexDistance);
	      blipProxTrk.push_back(blip_ProxTrkDist->at(blip));
	      
	      // cout << " Blip Vertex Distance: " << blipVertexDistance << "  Blip Prox Track Distance: " << blip_ProxTrkDist->at(blip) << "  ID: " << blip  << "  Associated Blip ID: " << AssociatedBlipID.at(blipCounter-1)<< endl;
	      // cout << " Blip Vertex Distance: " << blipVertexDistance << "  Associated Blip Vertex Distance: " << blipVertex.at(blipCounter-1) <<  endl;


	      BlipProxDistVsPMissMag[0][0]->Fill(blip_ProxTrkDist->at(blip),recoPMissingMag ,event_weight);
	      BlipVertexDistVsPMissMag[0][0]->Fill(blipVertexDistance,recoPMissingMag ,event_weight);
	      BlipProxDistVsPMissDir[0][0]->Fill(blip_ProxTrkDist->at(blip),recoPMissingDir ,event_weight);
	      BlipVertexDistVsPMissDir[0][0]->Fill(blipVertexDistance,recoPMissingDir ,event_weight);
	      BlipMultiplicityVsProxDist[0][0]->Fill(blip_ProxTrkDist->at(blip),TMath::Abs(numBlips) ,event_weight);
	      BlipMultiplicityVsVertexDist[0][0]->Fill(blipVertexDistance, TMath::Abs(numBlips) ,event_weight);
	      ProxDistVsNeutronMultiplicity[0][0]->Fill(blip_ProxTrkDist->at(blip) ,TrueNeutronCounter , event_weight);
	      VertexDistVsNeutronMultiplicity[0][0]->Fill(blipVertexDistance ,TrueNeutronCounter , event_weight);
	      AssociatedBlipVertexDistVsProxTrkDist[0][0]->Fill(blip_ProxTrkDist->at(blip),blipVertexDistance , event_weight);
	      
	      AssociatedBlipVertexDistVsProxTrkDist[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip),blipVertexDistance , event_weight);
	      BlipProxDistVsPMissMag[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip),recoPMissingMag ,event_weight);
	      BlipVertexDistVsPMissMag[genie_mode][0]->Fill(blipVertexDistance,recoPMissingMag ,event_weight);
	      BlipProxDistVsPMissDir[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip),recoPMissingDir ,event_weight);
	      BlipVertexDistVsPMissDir[genie_mode][0]->Fill(blipVertexDistance,recoPMissingDir ,event_weight);
	      BlipMultiplicityVsProxDist[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip),TMath::Abs(numBlips) ,event_weight);
	      BlipMultiplicityVsVertexDist[genie_mode][0]->Fill(blipVertexDistance, TMath::Abs(numBlips) ,event_weight);
	      ProxDistVsNeutronMultiplicity[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip) ,TrueNeutronCounter , event_weight);
	      VertexDistVsNeutronMultiplicity[genie_mode][0]->Fill(blipVertexDistance ,TrueNeutronCounter , event_weight);
	     
	      
	      
	      BlipLocation[0][0][0]->Fill(blip_X->at(blip), event_weight);
	      BlipLocation[0][1][0]->Fill(blip_Y->at(blip), event_weight);
	      BlipLocation[0][2][0]->Fill(blip_Z->at(blip), event_weight);
	      BlipLocation[genie_mode][0][0]->Fill(blip_X->at(blip), event_weight);
	      BlipLocation[genie_mode][1][0]->Fill(blip_Y->at(blip), event_weight);
	      BlipLocation[genie_mode][2][0]->Fill(blip_Z->at(blip), event_weight);
	      
	      
	      BlipVertexDist[0][0]->Fill(blipVertexDistance,event_weight);
	      BlipVertexDist[genie_mode][0]->Fill(blipVertexDistance,event_weight);
	      
	      //BlipCharge[0][0]->Fill(blip_Charge->at(blip), event_weight);
	      BlipProxTrkDist[0][0]->Fill(blip_ProxTrkDist->at(blip), event_weight);
	      BlipTime[0][0]->Fill(blip_Time->at(blip), event_weight);
	    //BlipLength[0][0]->Fill(blip_Length->at(blip), event_weight);
	      BlipMaxWireSpan[0][0]->Fill(blip_MaxWireSpan->at(blip), event_weight);
	      BlipEnergy[0][0]->Fill(blip_Energy->at(blip), event_weight);
	      
	      //BlipCharge[genie_mode][0]->Fill(blip_Charge->at(blip), event_weight);
	      BlipProxTrkDist[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip), event_weight);
	      BlipTime[genie_mode][0]->Fill(blip_Time->at(blip), event_weight);
	      //BlipLength[genie_mode][0]->Fill(blip_Length->at(blip), event_weight);
	      BlipMaxWireSpan[genie_mode][0]->Fill(blip_MaxWireSpan->at(blip), event_weight);
	      BlipEnergy[genie_mode][0]->Fill(blip_Energy->at(blip), event_weight);

	      
	      
	      BlipVertexDistVsEnergy[0]->Fill(blipVertexDistance, blip_Energy->at(blip), event_weight);
	      BlipProxDistVsEnergy[0]->Fill(blip_ProxTrkDist->at(blip), blip_Energy->at(blip), event_weight);
	      

	      BlipVertexDistNeut[0][0][NNeut]->Fill(blipVertexDistance, event_weight);
	      BlipProxTrkDistNeut[0][0][NNeut]->Fill(blip_ProxTrkDist->at(blip), event_weight);

	      BlipVertexDistNeut[genie_mode][0][NNeut]->Fill(blipVertexDistance, event_weight);
	      BlipProxTrkDistNeut[genie_mode][0][NNeut]->Fill(blip_ProxTrkDist->at(blip), event_weight);

	      if (TrueNeutronCounter <NNeut-1){
		BlipVertexDistNeut[0][0][TrueNeutronCounter]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[0][0][TrueNeutronCounter]->Fill(blip_ProxTrkDist->at(blip), event_weight);

		BlipVertexDistNeut[genie_mode][0][TrueNeutronCounter]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[genie_mode][0][TrueNeutronCounter]->Fill(blip_ProxTrkDist->at(blip), event_weight);
	      } //End of true neutron Counter check
	      
	      else{
		BlipVertexDistNeut[0][0][NNeut-1]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[0][0][NNeut-1]->Fill(blip_ProxTrkDist->at(blip), event_weight);

		BlipVertexDistNeut[genie_mode][0][NNeut-1]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[genie_mode][0][NNeut-1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
	      } //End of overflow for neutron statement



	      if (blip_pdg->at(blip) == 2212){
		ProtonBlipDist[0][0]->Fill(blipVertexDistance, event_weight);
		ProtonBlipDist[genie_mode][0]->Fill(blipVertexDistance, event_weight);	      
		
		//ProtonBlipCharge[0][0]->Fill(blip_Charge->at(blip), event_weight);
		ProtonBlipProxTrkDist[0][0]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		ProtonBlipTime[0][0]->Fill(blip_Time->at(blip), event_weight);
		//ProtonBlipLength[0][0]->Fill(blip_Length->at(blip), event_weight);
		ProtonBlipMaxWireSpan[0][0]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		ProtonBlipEnergy[0][0]->Fill(blip_Energy->at(blip), event_weight);
		
		//ProtonBlipCharge[genie_mode][0]->Fill(blip_Charge->at(blip), event_weight);
		ProtonBlipProxTrkDist[genie_mode][0]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		ProtonBlipTime[genie_mode][0]->Fill(blip_Time->at(blip), event_weight);
		//ProtonBlipLength[genie_mode][0]->Fill(blip_Length->at(blip), event_weight);
		ProtonBlipMaxWireSpan[genie_mode][0]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		ProtonBlipEnergy[genie_mode][0]->Fill(blip_Energy->at(blip), event_weight);
		ProtonBlipVertexDistVsEnergy[0]->Fill(blipVertexDistance, blip_Energy->at(blip), event_weight);
		ProtonBlipProxDistVsEnergy[0]->Fill(blip_ProxTrkDist->at(blip), blip_Energy->at(blip), event_weight);
		protonBlip = true;
	      } // End of if statement for proton blips


	      BlipPDG[0][0]->Fill(blip_pdg->at(blip) , event_weight);
	      BlipPDG[genie_mode][0]->Fill(blip_pdg->at(blip) , event_weight);
	    
	      //if (blip_pdg->at(blip) < 0 ) cout << "Blip PDG:  " << blip_pdg->at(blip) << endl;
	    

	      if (CC1p0pi == true) { 
		BlipPDG[0][1]->Fill(blip_pdg->at(blip), event_weight); 
		BlipPDG[genie_mode][1]->Fill(blip_pdg->at(blip), event_weight); 
	      
		BlipLocation[0][0][1]->Fill(blip_X->at(blip), event_weight);
		BlipLocation[0][1][1]->Fill(blip_Y->at(blip), event_weight);
		BlipLocation[0][2][1]->Fill(blip_Z->at(blip), event_weight);

		BlipLocation[genie_mode][0][1]->Fill(blip_X->at(blip), event_weight);
		BlipLocation[genie_mode][1][1]->Fill(blip_Y->at(blip), event_weight);
		BlipLocation[genie_mode][2][1]->Fill(blip_Z->at(blip), event_weight);
	      
		BlipVertexDist[0][1]->Fill(blipVertexDistance,event_weight);
		BlipVertexDist[genie_mode][1]->Fill(blipVertexDistance,event_weight);
		//BlipCharge[0][1]->Fill(blip_Charge->at(blip), event_weight);
		BlipProxTrkDist[0][1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		BlipTime[0][1]->Fill(blip_Time->at(blip), event_weight);
		//BlipLength[0][1]->Fill(blip_Length->at(blip), event_weight);
		BlipMaxWireSpan[0][1]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		BlipEnergy[0][1]->Fill(blip_Energy->at(blip), event_weight);
	      
		//BlipCharge[genie_mode][1]->Fill(blip_Charge->at(blip), event_weight);
		BlipProxTrkDist[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		BlipTime[genie_mode][1]->Fill(blip_Time->at(blip), event_weight);
		// BlipLength[genie_mode][1]->Fill(blip_Length->at(blip), event_weight);
		BlipMaxWireSpan[genie_mode][1]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		BlipEnergy[genie_mode][1]->Fill(blip_Energy->at(blip), event_weight);


		BlipProxDistVsPMissMag[0][1]->Fill(blip_ProxTrkDist->at(blip),recoPMissingMag ,event_weight);
		BlipVertexDistVsPMissMag[0][1]->Fill(blipVertexDistance,recoPMissingMag ,event_weight);
		BlipProxDistVsPMissDir[0][1]->Fill(blip_ProxTrkDist->at(blip),recoPMissingDir ,event_weight);
		BlipVertexDistVsPMissDir[0][1]->Fill(blipVertexDistance,recoPMissingDir ,event_weight);
		BlipMultiplicityVsProxDist[0][1]->Fill(blip_ProxTrkDist->at(blip),TMath::Abs(numBlips) ,event_weight);
		BlipMultiplicityVsVertexDist[0][1]->Fill(blipVertexDistance, TMath::Abs(numBlips) ,event_weight);
		ProxDistVsNeutronMultiplicity[0][1]->Fill(blip_ProxTrkDist->at(blip) ,TrueNeutronCounter , event_weight);
		VertexDistVsNeutronMultiplicity[0][1]->Fill(blipVertexDistance ,TrueNeutronCounter , event_weight);
		AssociatedBlipVertexDistVsProxTrkDist[0][1]->Fill(blip_ProxTrkDist->at(blip),blipVertexDistance , event_weight);
	      
		AssociatedBlipVertexDistVsProxTrkDist[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip),blipVertexDistance , event_weight);
		ProxDistVsNeutronMultiplicity[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip) ,TrueNeutronCounter , event_weight);
		VertexDistVsNeutronMultiplicity[genie_mode][1]->Fill(blipVertexDistance ,TrueNeutronCounter , event_weight);
		BlipProxDistVsPMissMag[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip),recoPMissingMag ,event_weight);
		BlipVertexDistVsPMissMag[genie_mode][1]->Fill(blipVertexDistance,recoPMissingMag ,event_weight);
		BlipProxDistVsPMissDir[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip),recoPMissingDir ,event_weight);
		BlipVertexDistVsPMissDir[genie_mode][1]->Fill(blipVertexDistance,recoPMissingDir ,event_weight);
		BlipMultiplicityVsProxDist[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip),TMath::Abs(numBlips) ,event_weight);
		BlipMultiplicityVsVertexDist[genie_mode][1]->Fill(blipVertexDistance, TMath::Abs(numBlips) ,event_weight);

		BlipVertexDistVsEnergy[1]->Fill(blipVertexDistance, blip_Energy->at(blip), event_weight);
		BlipProxDistVsEnergy[1]->Fill(blip_ProxTrkDist->at(blip), blip_Energy->at(blip), event_weight);
	      

		BlipVertexDistNeut[0][1][NNeut]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[0][1][NNeut]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		
		BlipVertexDistNeut[genie_mode][1][NNeut]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[genie_mode][1][NNeut]->Fill(blip_ProxTrkDist->at(blip), event_weight);
	      		
		if (TrueNeutronCounter <NNeut-1){
		  BlipVertexDistNeut[0][1][TrueNeutronCounter]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[0][1][TrueNeutronCounter]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  
		  BlipVertexDistNeut[genie_mode][1][TrueNeutronCounter]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[genie_mode][1][TrueNeutronCounter]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		} //End of true neutron Counter check
		
		else{
		  BlipVertexDistNeut[0][1][NNeut-1]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[0][1][NNeut-1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  
		  BlipVertexDistNeut[genie_mode][1][NNeut-1]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[genie_mode][1][NNeut-1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		} //End of overflow for neutron statement

		

		if (blip_pdg->at(blip) == 2212){
		  ProtonBlipDist[0][1]->Fill(blipVertexDistance, event_weight);
		  ProtonBlipDist[genie_mode][1]->Fill(blipVertexDistance, event_weight);	      
		
		  //ProtonBlipCharge[0][1]->Fill(blip_Charge->at(blip), event_weight);
		  ProtonBlipProxTrkDist[0][1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  ProtonBlipTime[0][1]->Fill(blip_Time->at(blip), event_weight);
		  //ProtonBlipLength[0][1]->Fill(blip_Length->at(blip), event_weight);
		  ProtonBlipMaxWireSpan[0][1]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		  ProtonBlipEnergy[0][1]->Fill(blip_Energy->at(blip), event_weight);
		
		  //ProtonBlipCharge[genie_mode][1]->Fill(blip_Charge->at(blip), event_weight);
		  ProtonBlipProxTrkDist[genie_mode][1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  ProtonBlipTime[genie_mode][1]->Fill(blip_Time->at(blip), event_weight);
		  //ProtonBlipLength[genie_mode][1]->Fill(blip_Length->at(blip), event_weight);
		  ProtonBlipMaxWireSpan[genie_mode][1]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		  ProtonBlipEnergy[genie_mode][1]->Fill(blip_Energy->at(blip), event_weight);
	      
		  ProtonBlipVertexDistVsEnergy[1]->Fill(blipVertexDistance, blip_Energy->at(blip), event_weight);
		  ProtonBlipProxDistVsEnergy[1]->Fill(blip_ProxTrkDist->at(blip), blip_Energy->at(blip), event_weight);

		} //End of if statement for proton blips
	      } //End of if else statement for plotting CC1p0pi events
	    
	   
	    

	      else { 

	   

		BlipPDG[0][2]->Fill(blip_pdg->at(blip), event_weight);
		BlipLocation[0][0][2]->Fill(blip_X->at(blip), event_weight);
		BlipLocation[0][1][2]->Fill(blip_Y->at(blip), event_weight);
		BlipLocation[0][2][2]->Fill(blip_Z->at(blip), event_weight);

		BlipPDG[genie_mode][2]->Fill(blip_pdg->at(blip), event_weight);
		BlipLocation[genie_mode][0][2]->Fill(blip_X->at(blip), event_weight);
		BlipLocation[genie_mode][1][2]->Fill(blip_Y->at(blip), event_weight);
		BlipLocation[genie_mode][2][2]->Fill(blip_Z->at(blip), event_weight);
	      
		BlipVertexDist[0][2]->Fill(blipVertexDistance,event_weight);
		BlipVertexDist[genie_mode][2]->Fill(blipVertexDistance,event_weight);
	      
		//BlipCharge[0][2]->Fill(blip_Charge->at(blip), event_weight);
		BlipProxTrkDist[0][2]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		BlipTime[0][2]->Fill(blip_Time->at(blip), event_weight);
		//BlipLength[0][2]->Fill(blip_Length->at(blip), event_weight);
		BlipMaxWireSpan[0][2]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		BlipEnergy[0][2]->Fill(blip_Energy->at(blip), event_weight);
	      
		//BlipCharge[genie_mode][2]->Fill(blip_Charge->at(blip), event_weight);
		BlipProxTrkDist[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		BlipTime[genie_mode][2]->Fill(blip_Time->at(blip), event_weight);
		//BlipLength[genie_mode][2]->Fill(blip_Length->at(blip), event_weight);
		BlipMaxWireSpan[genie_mode][2]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		BlipEnergy[genie_mode][2]->Fill(blip_Energy->at(blip), event_weight);



		BlipProxDistVsPMissMag[0][2]->Fill(blip_ProxTrkDist->at(blip),recoPMissingMag ,event_weight);
		BlipVertexDistVsPMissMag[0][2]->Fill(blipVertexDistance,recoPMissingMag ,event_weight);
		BlipProxDistVsPMissDir[0][2]->Fill(blip_ProxTrkDist->at(blip),recoPMissingDir ,event_weight);
		BlipVertexDistVsPMissDir[0][2]->Fill(blipVertexDistance,recoPMissingDir ,event_weight);
		BlipMultiplicityVsProxDist[0][2]->Fill(blip_ProxTrkDist->at(blip),TMath::Abs(numBlips) ,event_weight);
		BlipMultiplicityVsVertexDist[0][2]->Fill(blipVertexDistance, TMath::Abs(numBlips) ,event_weight);
		ProxDistVsNeutronMultiplicity[0][2]->Fill(blip_ProxTrkDist->at(blip) ,TrueNeutronCounter , event_weight);
		VertexDistVsNeutronMultiplicity[0][2]->Fill(blipVertexDistance ,TrueNeutronCounter , event_weight);
		AssociatedBlipVertexDistVsProxTrkDist[0][2]->Fill(blip_ProxTrkDist->at(blip),blipVertexDistance , event_weight);
	      
		AssociatedBlipVertexDistVsProxTrkDist[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip),blipVertexDistance , event_weight);
		ProxDistVsNeutronMultiplicity[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip) ,TrueNeutronCounter , event_weight);
		VertexDistVsNeutronMultiplicity[genie_mode][2]->Fill(blipVertexDistance ,TrueNeutronCounter , event_weight);
		BlipProxDistVsPMissMag[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip),recoPMissingMag ,event_weight);
		BlipVertexDistVsPMissMag[genie_mode][2]->Fill(blipVertexDistance,recoPMissingMag ,event_weight);
		BlipProxDistVsPMissDir[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip),recoPMissingDir ,event_weight);
		BlipVertexDistVsPMissDir[genie_mode][2]->Fill(blipVertexDistance,recoPMissingDir ,event_weight);
		BlipMultiplicityVsProxDist[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip),TMath::Abs(numBlips) ,event_weight);
		BlipMultiplicityVsVertexDist[genie_mode][2]->Fill(blipVertexDistance, TMath::Abs(numBlips) ,event_weight);

		BlipVertexDistVsEnergy[2]->Fill(blipVertexDistance, blip_Energy->at(blip), event_weight);
		BlipProxDistVsEnergy[2]->Fill(blip_ProxTrkDist->at(blip), blip_Energy->at(blip), event_weight);

		BlipVertexDistNeut[0][2][NNeut]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[0][2][NNeut]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		
		BlipVertexDistNeut[genie_mode][2][NNeut]->Fill(blipVertexDistance, event_weight);
		BlipProxTrkDistNeut[genie_mode][2][NNeut]->Fill(blip_ProxTrkDist->at(blip), event_weight);


		if (TrueNeutronCounter <NNeut-1){
		  BlipVertexDistNeut[0][2][TrueNeutronCounter]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[0][2][TrueNeutronCounter]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  
		  BlipVertexDistNeut[genie_mode][2][TrueNeutronCounter]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[genie_mode][2][TrueNeutronCounter]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		} //End of true neutron Counter check
		
		else{
		  BlipVertexDistNeut[0][2][NNeut-1]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[0][2][NNeut-1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  
		  BlipVertexDistNeut[genie_mode][2][NNeut-1]->Fill(blipVertexDistance, event_weight);
		  BlipProxTrkDistNeut[genie_mode][2][NNeut-1]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		} //End of overflow for neutron statement
		


		if (blip_pdg->at(blip) == 2212){
		  ProtonBlipDist[0][2]->Fill(blipVertexDistance, event_weight);
		  ProtonBlipDist[genie_mode][2]->Fill(blipVertexDistance, event_weight);	      
		
		  //ProtonBlipCharge[0][2]->Fill(blip_Charge->at(blip), event_weight);
		  ProtonBlipProxTrkDist[0][2]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  ProtonBlipTime[0][2]->Fill(blip_Time->at(blip), event_weight);
		  //ProtonBlipLength[0][2]->Fill(blip_Length->at(blip), event_weight);
		  ProtonBlipMaxWireSpan[0][2]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		  ProtonBlipEnergy[0][2]->Fill(blip_Energy->at(blip), event_weight);
		
		  //ProtonBlipCharge[genie_mode][2]->Fill(blip_Charge->at(blip), event_weight);
		  ProtonBlipProxTrkDist[genie_mode][2]->Fill(blip_ProxTrkDist->at(blip), event_weight);
		  ProtonBlipTime[genie_mode][2]->Fill(blip_Time->at(blip), event_weight);
		  //ProtonBlipLength[genie_mode][2]->Fill(blip_Length->at(blip), event_weight);
		  ProtonBlipMaxWireSpan[genie_mode][2]->Fill(blip_MaxWireSpan->at(blip), event_weight);
		  ProtonBlipEnergy[genie_mode][2]->Fill(blip_Energy->at(blip), event_weight);

		  ProtonBlipVertexDistVsEnergy[2]->Fill(blipVertexDistance, blip_Energy->at(blip), event_weight);
		  ProtonBlipProxDistVsEnergy[2]->Fill(blip_ProxTrkDist->at(blip), blip_Energy->at(blip), event_weight);
		  
		} //End of if statement for proton blips
	      } //End of else statement for plotting nonCC1p0pi events
	      
	    } //End of if statement for Trk ID == MuonID/ProtonID 
	    
	  } //End of if statement for blip_pdg not default value (-9)
        } //End of for loop over the blips
	

	//cout << blipCounter << endl;
	AssociatedBlipMultiplicity[0][0]->Fill(blipCounter, event_weight);
	AssociatedBlipMultiplicity[genie_mode][0]->Fill(blipCounter, event_weight);

	if (CC1p0pi == true){
	  AssociatedBlipMultiplicity[0][1]->Fill(blipCounter, event_weight);
	  AssociatedBlipMultiplicity[genie_mode][1]->Fill(blipCounter, event_weight);
	}

	else{
	  AssociatedBlipMultiplicity[0][2]->Fill(blipCounter, event_weight);
	  AssociatedBlipMultiplicity[genie_mode][2]->Fill(blipCounter, event_weight);
	}


	//Using the Associated blip Multiplicity, rerun through the blips to make 2D plots of each
	for (int ablip =0; ablip< blipCounter; ablip++){
	  blipVertexTooSmall = false;
	  if (blipVertex.at(ablip) < 150){ blipVertexTooSmall = true; } //Changing flag for a blip too close
	
	  AssociatedBlipMultVsVertexDist[0][0][NNeut]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	  AssociatedBlipMultVsProxTrkDist[0][0][NNeut]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	  
	  AssociatedBlipMultVsVertexDist[genie_mode][0][NNeut]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	  AssociatedBlipMultVsProxTrkDist[genie_mode][0][NNeut]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;


	  if(CC1p0pi ==true){
	    AssociatedBlipMultVsVertexDist[0][1][NNeut]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[0][1][NNeut]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	    
	    AssociatedBlipMultVsVertexDist[genie_mode][1][NNeut]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[genie_mode][1][NNeut]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	  } //End of if statement for true signal events

	  
	  else {
	    AssociatedBlipMultVsVertexDist[0][2][NNeut]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[0][2][NNeut]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	    
	    AssociatedBlipMultVsVertexDist[genie_mode][2][NNeut]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[genie_mode][2][NNeut]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	  }


	  
	  if (TrueNeutronCounter < NNeut-1){
	    AssociatedBlipMultVsVertexDist[0][0][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[0][0][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;

	    //cout << "Associated Blip Vertex Distance: " << blipVertex.at(ablip) << "  Associated Blip Prox Track Distance: " <<blipProxTrk.at(ablip)<< endl;

	    AssociatedBlipMultVsVertexDist[genie_mode][0][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[genie_mode][0][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	 
	    if (CC1p0pi ==true){
	      AssociatedBlipMultVsVertexDist[0][1][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[0][1][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	      
	      AssociatedBlipMultVsVertexDist[genie_mode][1][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[genie_mode][1][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	    } //End of if statement for true signal events

	    else if (CC1p0pi == false){
	      AssociatedBlipMultVsVertexDist[0][2][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[0][2][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	      
	      AssociatedBlipMultVsVertexDist[genie_mode][2][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[genie_mode][2][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;

	    } //End of else if statement for background events
	  } //End of if statement for neutron multiplicity requirements
	  

	  
	  else {
	    

	    AssociatedBlipMultVsVertexDist[0][0][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[0][0][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;

	    //cout << "Associated Blip Vertex Distance: " << blipVertex.at(ablip) << "  Associated Blip Prox Track Distance: " <<blipProxTrk.at(ablip)<< endl;

	    AssociatedBlipMultVsVertexDist[genie_mode][0][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	    AssociatedBlipMultVsProxTrkDist[genie_mode][0][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	 
	    if (CC1p0pi ==true){
	      AssociatedBlipMultVsVertexDist[0][1][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[0][1][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	      
	      AssociatedBlipMultVsVertexDist[genie_mode][1][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[genie_mode][1][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	    } //End of if statement for true signal events

	    else if (CC1p0pi == false){
	      AssociatedBlipMultVsVertexDist[0][2][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[0][2][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;
	      
	      AssociatedBlipMultVsVertexDist[genie_mode][2][TrueNeutronCounter]->Fill(blipVertex.at(ablip) ,blipCounter , event_weight) ;
	      AssociatedBlipMultVsProxTrkDist[genie_mode][2][TrueNeutronCounter]->Fill(blipProxTrk.at(ablip) ,blipCounter , event_weight) ;

	    } //End of else if statement for background events
	  } //End of overflow bin for neutron multiplicity

	 
	}//End of for loop over the associated blips
      } //End of if statement for blip_ID->size() >0 
      



      //Perform the second round of selection cuts using the associated blips and vertex distance

      //Using the Associated blips selection cut
      if (!(blip_ID->size() >0 && blipCounter >5)){
	
	if (blipCounter >5 ) cout << "--------------------------------------- Error 1 Multiplicity ---------------------------------" << endl;

	PostBlipCutRecoMuonCosThetaPlot[0][0][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	PostBlipCutRecoDeltaPtPlot[0][0][0]->Fill(recoDeltaPT,event_weight);
	PostBlipCutRecoDeltaAlphaTPlot[0][0][0]->Fill(recoDeltaAlphaT,event_weight);
	PostBlipCutTrueNeutronMultiplicityPlot[0][0][0]->Fill(TrueNeutronCounter, event_weight);
      
	PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][0][0]->Fill(TrueNeutronCounter, event_weight);      
	PostBlipCutRecoMuonCosThetaPlot[genie_mode][0][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	PostBlipCutRecoDeltaPtPlot[genie_mode][0][0]->Fill(recoDeltaPT,event_weight);
	PostBlipCutRecoDeltaAlphaTPlot[genie_mode][0][0]->Fill(recoDeltaAlphaT,event_weight);

	if (CC1p0pi ==true){
	  PostBlipCutRecoMuonCosThetaPlot[0][1][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[0][1][0]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[0][1][0]->Fill(recoDeltaAlphaT,event_weight);
	  PostBlipCutTrueNeutronMultiplicityPlot[0][1][0]->Fill(TrueNeutronCounter, event_weight);
	
	  PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][1][0]->Fill(TrueNeutronCounter, event_weight);      	
	  PostBlipCutRecoMuonCosThetaPlot[genie_mode][1][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[genie_mode][1][0]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[genie_mode][1][0]->Fill(recoDeltaAlphaT,event_weight);
	} //End of CC1p0pi plotting section

	else {
	  PostBlipCutRecoMuonCosThetaPlot[0][2][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[0][2][0]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[0][2][0]->Fill(recoDeltaAlphaT,event_weight);
	  PostBlipCutTrueNeutronMultiplicityPlot[0][2][0]->Fill(TrueNeutronCounter, event_weight);
      
	  PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][2][0]->Fill(TrueNeutronCounter, event_weight);      
	  PostBlipCutRecoMuonCosThetaPlot[genie_mode][2][0]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[genie_mode][2][0]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[genie_mode][2][0]->Fill(recoDeltaAlphaT,event_weight);
	} //End of nonCC1p0pi plotting section
      } //End of associated blip cuts



      //Applying the vertex distance cut
      if (!(blip_ID->size() >0 && blipVertexTooSmall ==true)){ 

	if (blipVertexTooSmall == true ) cout <<" ------------------------------------- Error 2 Vertex ------------------------------------------" << endl;


	PostBlipCutRecoMuonCosThetaPlot[0][0][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	PostBlipCutRecoDeltaPtPlot[0][0][1]->Fill(recoDeltaPT,event_weight);
	PostBlipCutRecoDeltaAlphaTPlot[0][0][1]->Fill(recoDeltaAlphaT,event_weight);
	PostBlipCutTrueNeutronMultiplicityPlot[0][0][1]->Fill(TrueNeutronCounter, event_weight);
      
	PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][0][1]->Fill(TrueNeutronCounter, event_weight);      
	PostBlipCutRecoMuonCosThetaPlot[genie_mode][0][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	PostBlipCutRecoDeltaPtPlot[genie_mode][0][1]->Fill(recoDeltaPT,event_weight);
	PostBlipCutRecoDeltaAlphaTPlot[genie_mode][0][1]->Fill(recoDeltaAlphaT,event_weight);

	if (CC1p0pi ==true){
	  PostBlipCutRecoMuonCosThetaPlot[0][1][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[0][1][1]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[0][1][1]->Fill(recoDeltaAlphaT,event_weight);
	  PostBlipCutTrueNeutronMultiplicityPlot[0][1][1]->Fill(TrueNeutronCounter, event_weight);
	
	  PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][1][1]->Fill(TrueNeutronCounter, event_weight);      	
	  PostBlipCutRecoMuonCosThetaPlot[genie_mode][1][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[genie_mode][1][1]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[genie_mode][1][1]->Fill(recoDeltaAlphaT,event_weight);
	} //End of CC1p0pi plotting section

	else {
	  PostBlipCutRecoMuonCosThetaPlot[0][2][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[0][2][1]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[0][2][1]->Fill(recoDeltaAlphaT,event_weight);
	  PostBlipCutTrueNeutronMultiplicityPlot[0][2][1]->Fill(TrueNeutronCounter, event_weight);
      
	  PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][2][1]->Fill(TrueNeutronCounter, event_weight);      
	  PostBlipCutRecoMuonCosThetaPlot[genie_mode][2][1]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[genie_mode][2][1]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[genie_mode][2][1]->Fill(recoDeltaAlphaT,event_weight);
	} //End of nonCC1p0pi plotting section
      } //End of vertex distance cut




      //Perform cuts on both vertex distance and associated blip multiplicity      
      if (!(blip_ID->size() >0 && (blipVertexTooSmall ==true || blipCounter >5))){
	
	if (blipVertexTooSmall == true ) cout <<" ----------------------------------------- Error 3 Vertex ------------------------------------------" << endl;
	if (blipCounter >5 ) cout << "--------------------------------------- Error 3 Multiplicity ---------------------------------" << endl;

	PostBlipCutRecoMuonCosThetaPlot[0][0][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	PostBlipCutRecoDeltaPtPlot[0][0][2]->Fill(recoDeltaPT,event_weight);
	PostBlipCutRecoDeltaAlphaTPlot[0][0][2]->Fill(recoDeltaAlphaT,event_weight);
	PostBlipCutTrueNeutronMultiplicityPlot[0][0][2]->Fill(TrueNeutronCounter, event_weight);
      
	PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][0][2]->Fill(TrueNeutronCounter, event_weight);      
	PostBlipCutRecoMuonCosThetaPlot[genie_mode][0][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	PostBlipCutRecoDeltaPtPlot[genie_mode][0][2]->Fill(recoDeltaPT,event_weight);
	PostBlipCutRecoDeltaAlphaTPlot[genie_mode][0][2]->Fill(recoDeltaAlphaT,event_weight);

	if (CC1p0pi ==true){
	  PostBlipCutRecoMuonCosThetaPlot[0][1][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[0][1][2]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[0][1][2]->Fill(recoDeltaAlphaT,event_weight);
	  PostBlipCutTrueNeutronMultiplicityPlot[0][1][2]->Fill(TrueNeutronCounter, event_weight);
	
	  PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][1][2]->Fill(TrueNeutronCounter, event_weight);      	
	  PostBlipCutRecoMuonCosThetaPlot[genie_mode][1][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[genie_mode][1][2]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[genie_mode][1][2]->Fill(recoDeltaAlphaT,event_weight);
	} //End of CC1p0pi plotting section

	else {
	  PostBlipCutRecoMuonCosThetaPlot[0][2][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[0][2][2]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[0][2][2]->Fill(recoDeltaAlphaT,event_weight);
	  PostBlipCutTrueNeutronMultiplicityPlot[0][2][2]->Fill(TrueNeutronCounter, event_weight);
      
	  PostBlipCutTrueNeutronMultiplicityPlot[genie_mode][2][2]->Fill(TrueNeutronCounter, event_weight);      
	  PostBlipCutRecoMuonCosThetaPlot[genie_mode][2][2]->Fill(TMath::Cos(CandidateMuonTrackTheta),event_weight);
	  PostBlipCutRecoDeltaPtPlot[genie_mode][2][2]->Fill(recoDeltaPT,event_weight);
	  PostBlipCutRecoDeltaAlphaTPlot[genie_mode][2][2]->Fill(recoDeltaAlphaT,event_weight);
	} //End of nonCC1p0pi plotting section
      } // End of combination of associated blip cut and vertex distance cut





      // Reco CC1p0pi candidate event selected, print out details
      if (CC1p0pi) {
	//if (protonBlip ==true ){
	cout << "run = " << run << ", sub = " << sub << ", evt = " << evt << ", X-Vertex: " << trueVertex.X() <<", Y-Vertex: " << trueVertex.Y() << ", Z-Vertex: " << trueVertex.Z()  << endl; 
	counter++;
	  //} //End of if statement for proton blips
      } //End of if statement for CC1p0pi events
     

   } // End of the loop over the events
   cout << "Counter: " << counter  << "  True Counter: " << trueCounter << "   No Blip Counter: " << noBlipCounter<< endl;


   // Divide by bin width with Reweight function





   OutputFile->cd();
   OutputFile->Write();
   OutputFile->Close();


} // End of the program
