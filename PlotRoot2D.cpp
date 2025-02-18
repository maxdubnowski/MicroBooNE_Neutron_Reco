#include <TFile.h>
#include <TTree.h>
#include <TString.h>
#include <TPad.h>
#include "Constants.h"
using namespace std;
using namespace Constants;

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <fstream>
#include <stdlib.h>

void PlotRoot2D() {

	//------------------------------//

	TH1D::SetDefaultSumw2();
	TH2D::SetDefaultSumw2();
	gStyle->SetOptStat(0);

	int FontStyle = 132;
	double TextSize = 0.06;			
	double LegendTextSize = 0.03;

	TString OutFilePath = "/uboone/app/users/maxd/MicroBooNE_Neutron_Reco/";


	
	// Event generators

	std::vector<TString> Names; std::vector<TString> Labels; std::vector<int> Colors;

	Names.push_back(OutFilePath+"MicroBooNE_Reco.root"); 
	Labels.push_back(" ");
	Colors.push_back(kBlack);	





	const int NSamples = Names.size();
	std::vector<TFile*> Files; Files.resize(NSamples);
	


	// Plots to overlay

	std::vector<TString> PlotNames;

	PlotNames.push_back("BacktrackedvsRecoDeltaPtPlot_AllEvents");
	PlotNames.push_back("BacktrackedvsRecoDeltaPtPlot_CC1p0piEvent");
	PlotNames.push_back("BacktrackedvsRecoDeltaPtPlot_nonCC1p0piEvent");

	PlotNames.push_back("BacktrackedvsRecoDeltaAlphaTPlot_AllEvents");
	PlotNames.push_back("BacktrackedvsRecoDeltaAlphaTPlot_CC1p0piEvent");
	PlotNames.push_back("BacktrackedvsRecoDeltaAlphaTPlot_nonCC1p0piEvent");
	
	PlotNames.push_back("BacktrackedvsRecoPMissDirPlot_AllEvents");
	PlotNames.push_back("BacktrackedvsRecoPMissDirPlot_CC1p0piEvent");
	PlotNames.push_back("BacktrackedvsRecoPMissDirPlot_nonCC1p0piEvent");

	PlotNames.push_back("BacktrackedvsRecoPMissMagPlot_AllEvents");
	PlotNames.push_back("BacktrackedvsRecoPMissMagPlot_CC1p0piEvent");
	PlotNames.push_back("BacktrackedvsRecoPMissMagPlot_nonCC1p0piEvent");


	PlotNames.push_back("BacktrackedvsRecoPMissMagPlot_nonCC1p0piEvent");
	PlotNames.push_back("BacktrackedvsRecoPMissMagPlot_nonCC1p0piEvent");
	PlotNames.push_back("BacktrackedvsRecoPMissMagPlot_nonCC1p0piEvent");
	
	// PlotNames.push_back("CheckMomentumXDirection");
	// PlotNames.push_back("CheckMomentumYDirection");
	// PlotNames.push_back("CheckMomentumZDirection");


	// PlotNames.push_back("BlipProxDistVsPMissMag_AllEvents");
	// PlotNames.push_back("BlipProxDistVsPMissMag_CC1p0piEvent");
	// PlotNames.push_back("BlipProxDistVsPMissMag_nonCC1p0piEvent");

	// PlotNames.push_back("BlipVertexDistVsPMissMag_AllEvents");
	// PlotNames.push_back("BlipVertexDistVsPMissMag_CC1p0piEvent");
	// PlotNames.push_back("BlipVertexDistVsPMissMag_nonCC1p0piEvent");

	// PlotNames.push_back("BlipProxDistVsPMissDir_AllEvents");
	// PlotNames.push_back("BlipProxDistVsPMissDir_CC1p0piEvent");
	// PlotNames.push_back("BlipProxDistVsPMissDir_nonCC1p0piEvent");

	// PlotNames.push_back("BlipVertexDistVsPMissDir_AllEvents");
	// PlotNames.push_back("BlipVertexDistVsPMissDir_CC1p0piEvent");
	// PlotNames.push_back("BlipVertexDistVsPMissDir_nonCC1p0piEvent");

	// PlotNames.push_back("BlipMultiplicityVsProxDist_AllEvents");
	// PlotNames.push_back("BlipMultiplicityVsProxDist_CC1p0piEvent");
	// PlotNames.push_back("BlipMultiplicityVsProxDist_nonCC1p0piEvent");

	// PlotNames.push_back("BlipMultiplicityVsVertexDist_AllEvents");
	// PlotNames.push_back("BlipMultiplicityVsVertexDist_CC1p0piEvent");
	// PlotNames.push_back("BlipMultiplicityVsVertexDist_nonCC1p0piEvent");
 
	// PlotNames.push_back("BlipVertexDistVsEnergy_AllEvents");
	// PlotNames.push_back("BlipVertexDistVsEnergy_CC1p0piEvent");
	// PlotNames.push_back("BlipVertexDistVsEnergy_nonCC1p0piEvent");

	// PlotNames.push_back("ProtonBlipVertexDistVsEnergy_AllEvents");
	// PlotNames.push_back("ProtonBlipVertexDistVsEnergy_CC1p0piEvent");
	// PlotNames.push_back("ProtonBlipVertexDistVsEnergy_nonCC1p0piEvent");

	// PlotNames.push_back("BlipProxDistVsEnergy_AllEvents");
	// PlotNames.push_back("BlipProxDistVsEnergy_CC1p0piEvent");
	// PlotNames.push_back("BlipProxDistVsEnergy_nonCC1p0piEvent");


	// PlotNames.push_back("ProtonBlipProxDistVsEnergy_AllEvents");
	// PlotNames.push_back("ProtonBlipProxDistVsEnergy_CC1p0piEvent");
	// PlotNames.push_back("ProtonBlipProxDistVsEnergy_nonCC1p0piEvent");




	const int NPlots = PlotNames.size();




	// Loop over the samples to open the files and the TTree

	for (int iSample = 0; iSample < NSamples; iSample++) {

		Files[iSample] = new TFile(Names[iSample],"readonly");

	} // End of the loop over the samples





	// Loop over the plots to be compared


	
	for (int iPlot = 0; iPlot < NPlots; iPlot++) {
	  // Loop over the samples to open the files and to get the corresponding plot
	  for (int iSample = 0; iSample < NSamples; iSample++) {	
	
	    std::vector<TH2D*> Histos; Histos.resize(NSamples);
	    TString CanvasName = "Canvas_" + PlotNames[iPlot];
	    TCanvas* PlotCanvas = new TCanvas(CanvasName,CanvasName,205,34,1024,768);
	    PlotCanvas->cd();
	    PlotCanvas->SetTopMargin(0.12);
	    PlotCanvas->SetLeftMargin(0.18);
	    PlotCanvas->SetBottomMargin(0.15);		
	    PlotCanvas->Draw();	
	    
	    
	    Histos[iSample] = (TH2D*)(Files[iSample]->Get(PlotNames[iPlot]));
	    
	    // Histos[iSample]->SetLineWidth(4);
	    // Histos[iSample]->SetLineColor( Colors.at(iSample) );	
	    Histos[iSample]->SetTitle(Labels[iSample]);
	    Histos[iSample]->GetXaxis()->SetTitleFont(FontStyle);
	    //Histos[iSample]->GetXaxis()->SetTitle("cos(#theta_{miss})");
	    Histos[iSample]->GetXaxis()->SetLabelFont(FontStyle);
	    
	    // Histos[iSample]->GetXaxis()->SetNdivisions(8);
	    Histos[iSample]->GetXaxis()->SetLabelSize(TextSize);
	    Histos[iSample]->GetXaxis()->SetTitleSize(TextSize);	
	    Histos[iSample]->GetXaxis()->SetTitleOffset(1.1);					
	    Histos[iSample]->GetXaxis()->CenterTitle();						
	    
	    Histos[iSample]->GetYaxis()->SetTitleFont(FontStyle);
	    //Histos[iSample]->GetYaxis()->SetTitle("Magnitude p_{missing}");
	    Histos[iSample]->GetYaxis()->SetLabelFont(FontStyle);
	    //Histos[iSample]->GetYaxis()->SetNdivisions(6);
	    
	    Histos[iSample]->GetYaxis()->SetLabelSize(TextSize);
	    Histos[iSample]->GetYaxis()->SetTitleSize(TextSize);
	    Histos[iSample]->GetYaxis()->SetTitleOffset(1.3);
	    Histos[iSample]->GetYaxis()->SetTickSize(0);
	    Histos[iSample]->GetYaxis()->CenterTitle();	
	    
	    PlotCanvas->cd()->SetLogz();
	    
	    
	    PlotCanvas->cd();
	    Histos[iSample]->Draw("colz");
	    Histos[0]->Draw("colz");	
	    
	    // leg->AddEntry(Histos[iSample],Labels[iSample],"l");
	    
	    
	    PlotCanvas->cd();
	    if (PMissingCut == false) PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"_MicroBooNE_Reco.pdf"); 
	    else if (PMissingCut == true) PlotCanvas->SaveAs("myPlots/"+PlotNames[iPlot]+"NeutronCut_MicroBooNE_Reco.pdf"); 
	    //PlotCanvas->SaveAs("myPlots/" PlotNames[iPlot]+"1MillionEvt_MicroBooNE_Reco.pdf");
	 
	  } // End of the loop over the samples grabing the plots	
	} // End of the loop over the plots
	
	
	
	

	

	
	

	

	
	
	
} // End of the program
