/*  ////---------------------------- DRAW ----------------------------////////
 
Description: This is a root macro to draw histograms produced by the CMS Open Data analyzer.
 
/----------------------------  Implementation  ----------------------------/
TEST2
 
It may be executed, via the terminal, with the following commands:
    root -l
    .x draw.cxx 
 
 
This histograms should then appear on the screen. Calculations may also be computed by this macro


To produce a plot one should carry out the following
 
1. Change the file opening string to the output file.
    -alternatively one may add a new line to open additional files
 
            TFile *filehandle = new TFile("outputfilename.root");
 
2. If the root file has sub directories get  a pointer to the directory
 
            TDirectory *directoryname=(TDirectory*) filehandle->Get("directory");
 
 This may need to be repeated
 
3. Add a pointer to the required histogram in the output file
                
            TH2D* newhistoname=(TH[1/2]D*)directoryname->Get("histoname");
 
4. Copy the block below replacing names and options as required:
 
            TCanvas *cN = new TCanvas("cN","",600,400); //Change N to a number
            newhistoname->GetXaxis()->SetTitle("Axis title (Units)");
            newhistoname->GetYaxis()->SetTitle("Axis title (Units)");
            newhistoname->Draw("Some Option");
            // Save the canvas
            c1->SaveAs("imagename.extension");
                    
 
 Further options and customisation details may be found at root.cern.ch
 (eg https://root.cern.ch/root/html/TH1.html )
 ------------------------------------------------------------------------------------------------------------

 */

//Header files
#include "TH1F.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TRandom3.h"
#include "TTree.h"
#include <iostream>
bool

double CrystalBall(double* x, double* par);
double Background(double* x, double* par);
double allFit(double* x, double *par );
void draw(){

  //Open the output file as it is called in line of "demoanalyzer_cfg.py"

            TFile *f = new TFile("OutputD03.root");
            //TFile *f = new TFile("OutputMuoniaUpsilon2.root");
            //f->cd("demo"); //Change to directory containing histograms
    
    
    
    //Get Pointer to main directory.
    TDirectory *dirf=(TDirectory*) f->Get("demo");
    
    //Get pointers to histogram directories
    TDirectory *dirfsubTracks=(TDirectory*) dirf->Get("Tracks");
    TDirectory *dirfsubMuons=(TDirectory*) dirf->Get("Muons");
    TDirectory *dirfsubTMuons=(TDirectory*) dirf->Get("TMuons");
    TDirectory *dirfsubGMuons=(TDirectory*) dirf->Get("GMuons");
    TDirectory *dirfsubPFJets=(TDirectory*) dirf->Get("PFJets");
    TDirectory *dirfsubPrVertex=(TDirectory*) dirf->Get("PrVertex");
    TDirectory *dirfsubUpsilon=(TDirectory*) dirf->Get("Upsilons");
    TDirectory *dirfsubDMesons=(TDirectory*) dirf->Get("DMesons");
    TDirectory *dirfsubAcceptance=(TDirectory*) dirf->Get("Acceptance");
    
    //global style preferences
    gStyle->SetPalette(1); //Set palette for 2D plots
    gStyle->SetOptStat("eou"); //Print name and no of entries
 
    //Define fitting functions


    
    TCanvas *c4 = new TCanvas("c4","",600,400);
    TH1D *hdeltaM=(TH1D*)dirfsubDMesons->Get("h_deltaM2ndcuts");
    TH1D *hdeltaMwrongcharge=(TH1D*)dirfsubDMesons->Get("h_deltaMwrongcharge2ndcuts");
    hdeltaM->SetStats(1);
    hdeltaM->SetLineColor(kBlue);
    hdeltaM->GetXaxis()->SetTitle("M(K#pi#pi)- M(K#pi) (GeV/c^{2})");
    hdeltaM->GetYaxis()->SetTitleOffset(0.05);
    hdeltaM->GetYaxis()->SetTitle("Number of Entries");
    hdeltaM->Draw("E");
    hdeltaMwrongcharge->SetLineColor(kRed);
    hdeltaMwrongcharge->Draw("E same");
    TLegend* legc4 = new TLegend(0.7, 0.1, .9, .3);
    legc4->AddEntry(hdeltaM, "Right Charge", "l");
    legc4->AddEntry(hdeltaMwrongcharge, "Wrong Charge", "l");
    legc4->Draw();
    legc4->SetBorderSize(0);
    TLegend* legc4i= new TLegend(0.6, 0.6, .89, .89);
    legc4i->SetBorderSize(0);
    legc4i->AddEntry((TObject*)0, "CMS,  #sqrt{s} = 7 TeV", "");
    legc4i->AddEntry((TObject*)0, "L = n pb^{-1}", "");
    //legc4i->AddEntry((TObject*)0, "#||{#eta^{#mu}} < 1 ", "");
    legc4i->Draw();

    c4->SaveAs("./Plots/DeltaD0Mass2ndcuts.png");
    
    
    
   // c8->SaveAs("Upsilon.png");
    
    //TCanvas *c9=new TCanvas("c9","", 600,500);
    //
    //cryst2->Draw("same");

  
 return;
}
double CrystalBall(double *x,double *par) {
    
    Double_t t = (x[0]-par[2])/par[3];
    if (par[0] < 0) t = -t;
    
    Double_t absAlpha = fabs((Double_t)par[0]);
    
    if (t >= -absAlpha) {
        return par[4]*exp(-0.5*t*t);
    }
    else {
        Double_t a =  TMath::Power(par[1]/absAlpha,par[1])*exp(-0.5*absAlpha*absAlpha);
        Double_t b= par[1]/absAlpha - absAlpha;
        
        return par[4]*(a/TMath::Power(b - t, par[1]));
    }
}


double Background(double* x, double* par) {
    if ( x[0] < 9.2 || x[0] > 10.5) {
        return exp(-(par[0]*x[0]-par[1]));
    }
    else{
        return 0;
    }
}
double allFit(double* x, double* par){
    return ( CrystalBall(x,par)+ CrystalBall(x,&par[5])+CrystalBall(x,&par[10]) ) ;

}



  /*  ---------------------------- ARCHIVE SECTION ----------------------------
   Use this section to store old/no longer needed code for future use.



   ---------------------------- Muon Pair Count pt/y ----------------------------
   // Plot the number of muons agains rapidity and pt

   
   TCanvas *c1 = new TCanvas("c1","",600,400);
   TH2D* hPaircount=(TH2D*)dirfsubAcceptance->Get("h_Paircount");
   hPaircount->SetStats(1);                           // Add/Remove the statistics box 1||0
   hPaircount->GetXaxis()->SetTitle("y^{J/#psi}");    // add axis titles
   hPaircount->GetYaxis()->SetTitle("p_{T}^{J/#psi} [GeV/c]");
   hPaircount->Draw("colz");// draw with a colour scale
   c1->SaveAs("./Plots/Paircount.png");

   TCanvas *c5 = new TCanvas("c5","",600,400);
   TH2D* hPaircount2=(TH2D*)dirfsubAcceptance->Get("h_Paircount2");
   hPaircount2->SetStats(1);                           // Add/Remove the statistics box 1||0
   hPaircount2->GetXaxis()->SetTitle("y^{J/#psi}");    // add axis titles
   hPaircount2->GetYaxis()->SetTitle("p_{T}^{J/#psi} [GeV/c]");
   hPaircount2->Draw("colz");// draw with a colour scale
   c2->SaveAs("./Plots/Paircount2.png");
   
   TCanvas *c2 = new TCanvas("c2","",600,400);
   h_JpsiPT->SetStats(1);
   h_JpsiPT->GetXaxis()->SetTitle("p_{T}^{J/#psi} [GeV/c]");
   h_JpsiPT->GetYaxis()->SetTitle("Counts");
   h_JpsiPT->Draw();
   c2->SaveAs("./Plots/JpsiPT.png");
   
   ---------------------------- D analysis ----------------------------
   
   TCanvas *c3 = new TCanvas("c3","",600,400);
   h_D0mass->SetStats(1);
   h_D0mass->GetXaxis()->SetTitle("D^{0} Mass [GeV/c^{2}]");
   h_D0mass->GetYaxis()->SetTitle("Number of Entries");
   h_D0mass->Draw();
   c3->SaveAs("./Plots/D0Mass.png");
   

   
   TCanvas *c6=new TCanvas("c6","" ,600,400);
   TH1D *hz=(TH1D*)dirfsubDMesons->Get("h_z");
   hz->SetStats(0);
   hz->GetXaxis()->SetTitle("z");
   hz->GetYaxis()->SetTitleOffset(0.05);
   hz->GetYaxis()->SetTitle("Number of Entries");
   hz->Draw();
   c6->SaveAs("./Plots/D0z.png");
   
   TCanvas *c7=new TCanvas("c7","", 600,400);
   TH1D *hn=(TH1D*)dirfsubDMesons->Get("h_n");
   hn->SetStats(0);
   hn->GetYaxis()->SetTitleOffset(0.05);
   hn->GetXaxis()->SetTitle("Number of Tracks from D0 region");
   hn->GetYaxis()->SetTitle("Number of Entries");
   hn->Draw();
   //cn->SaveAs("./Plots/D0z.png");

   
---------------------------- Upsilon Analysis ----------------------------
   TCanvas *c9=new TCanvas("c9","", 600,500);
   TH1D *hUpsiloneta=(TH1D*)dirfsubUpsilon->Get("h_Upsilon_eta");
   hUpsiloneta->SetStats(1);
   hUpsiloneta->GetYaxis()->SetTitleOffset(1.15);
   hUpsiloneta->GetXaxis()->SetTitle("#mu^{+} #mu^{-} mass (Gev/c^{2})");
   hUpsiloneta->GetYaxis()->SetTitle("Events");
   hUpsiloneta->Draw("E");
   TLegend* legc9 = new TLegend(0.6, 0.6, .89, .89);
   legc9->SetBorderSize(0);
   legc9->AddEntry((TObject*)0, "CMS,  #sqrt{s} = 7 TeV", "");
   legc9->AddEntry((TObject*)0, "L = n pb^{-1}", "");
   legc9->AddEntry((TObject*)0, "#||{#eta^{#mu}} < 1 ", "");
   legc9->Draw();
   c9->SaveAs("Upsiloneta.png");
   
   
   TF1* Fit1=new TF1("Fit1", allFit,8,12,15);
   //Set parameters
   Fit1->SetParameter(0,1.5);
   Fit1->SetParameter(1,0.3);
   Fit1->SetParameter(2,9.46);
   Fit1->SetParameter(3,0.1);
   Fit1->SetParameter(4,1400);
   
   Fit1->SetParameter(5,1);
   Fit1->SetParameter(6,0.22);
   Fit1->SetParameter(7,10.02326);
   Fit1->SetParameter(8,0.1);
   Fit1->SetParameter(9,600);
   
   Fit1->SetParameter(10,1);
   Fit1->SetParameter(11,0.3);
   Fit1->SetParameter(12,10.3552);
   Fit1->SetParameter(13,0.15);
   Fit1->SetParameter(14,550);
   
   
   TCanvas *c8=new TCanvas("c8","", 600,500);
   TH1D *hUpsilon=(TH1D*)dirfsubUpsilon->Get("h_Upsilon");
   hUpsilon->SetStats(0);
   hUpsilon->GetYaxis()->SetTitleOffset(1.15);
   hUpsilon->GetXaxis()->SetTitle("#mu^{+} #mu^{-} mass (Gev/c^{2})");
   hUpsilon->GetYaxis()->SetTitle("Events");
   hUpsilon->Draw("E");
   Fit1->Draw("same");
   //hUpsilon->Fit("Fit1");
   TLegend* legc8 = new TLegend(0.6, 0.6, .89, .89);
   legc8->SetBorderSize(0);
   legc8->AddEntry((TObject*)0, "CMS,  #sqrt{s} = 7 TeV", "");
   legc8->AddEntry((TObject*)0, "L = n pb^{-1}", "");
   legc8->AddEntry((TObject*)0, "#||{#eta^{#mu}} < 2.4 ", "");
   legc8->Draw();
   
   
   
   ----------------------------The DREGS ----------------------------

   
   //Fitting functions
   TF1* cryst1 = new TF1("cryst1",CrystalBall,8,12,5);
   TF1* cryst2 = new TF1("cryst2",CrystalBall,8,12,5);
   TF1 *custgaus = new TF1("custgaus",   "[0]*exp(-0.5*((x-[1])/[2])^2)", 8,12);
   TF1* bkgrnd= new TF1("bkgrnd",Background,0,1000,2);
   TF1* cryst3 = new TF1("cryst3",&CrystalBall,8,12,7);
  */
