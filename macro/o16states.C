#include <string>
#include <iostream>
#include <stdlib.h>

#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>

Int_t o16states(Double_t minEx=19.0, Double_t maxEx=35.0, Double_t step=0.1){
  // gROOT->SetBatch();
  gStyle->SetOptStat(0);

  // file for the experimental spectrum
  TFile *expfile = new TFile("./root/to_adachi_0924.root", "READ");
  // file for the cascade spectrum
  TFile *casfile;

  // 20Ne excitation spectrum
  TH1D* hexp0 = (TH1D*)expfile->Get("h50003")->Clone();
  // 20Ne excited states decay -> 16O(0p6) + alpha
  TH1D* hexp1 = (TH1D*)expfile->Get("h43605")->Clone();

  hexp0->SetName("hexp0");
  hexp1->SetName("hexp1");

  // bin information for CASCADE spectrum
  Double_t binw0 = step;
  Double_t xmin0 = minEx - binw0/2.0;
  Double_t xmax0 = maxEx + binw0/2.0;
  Int_t nbinx0 = (xmax0 - xmin0) / binw0;

  // bin information for experimental spectrum
  Int_t nbinx1 = hexp1->GetNbinsX();
  Double_t binw1 = hexp1->GetBinWidth(1);
  Double_t xmin1 = hexp1->GetBinCenter(1) - binw1/2.0;
  Double_t xmax1 = hexp1->GetBinCenter(nbinx1) + binw1/2.0;

  // Cascade original distribution
  TH1D *hcas0 = new TH1D("hcas0","CASCADE 0p6", nbinx0, xmin0, xmax0);
  // Cascade smeared distribution
  TH1D *hcas1 = new TH1D("hcas1","CASCADE 0p6, Ex-comp.", nbinx1, xmin1, xmax1);
  // Cascade smeared & Ex-compensated distribution
  TH1D *hcas2 = new TH1D("hcas2","CASCADE 0p6, Ex-comp. & smeared", nbinx1, xmin1, xmax1);

  // Cascade evaporation alpha spectrum
  TH1D *halpha;
  Int_t totalevent;
  TTree *tree;

  // cascade root file for evaporation alpha spectrum
  Char_t casname[256];
  Double_t exenergy;
  const Double_t intw = 0.32;
  Int_t leftb, rightb;
  Double_t prob;
  Double_t nevent;
  for(Int_t i=0;i<nbinx0;i++){
    exenergy = minEx + ((Double_t)i)*step;
    sprintf(casname, "./root/out%6.3f.root", exenergy);

    casfile = new TFile( casname, "READ" );
    // file open check
    if( casfile->IsZombie() ) {
      std::cout << " Error to open the cascade root file!: "
		<< casname << std::endl;
      return -1;
    }

    halpha = (TH1D*)casfile->Get("h2")->Clone();
    halpha->SetName("halpha");
    tree = (TTree*)casfile->Get("tree")->Clone();
    tree->SetName("tree");
    tree->SetBranchAddress("totalevent", &totalevent);
    tree->GetEntry(1);
    nevent = totalevent;
    
    leftb = halpha->FindBin( 15.10 - intw );
    rightb = halpha->FindBin( 15.10 + intw );
    prob = halpha->Integral( leftb, rightb ) / nevent;

    hcas0->SetBinContent( hcas0->FindBin(exenergy), prob );

    casfile->Close();
  }

  // smearing by the experimental energy resolution
  const Double_t pi = TMath::Pi();
  const Double_t sigma = 0.12; // smeared by this sigma in MeV
  Double_t count = 0.0;
  Double_t pos = 0.0, center = 0.0, sum = 0.0;
  // hcas1
  for(Int_t i=1;i<=nbinx1;i++){
    pos = hcas1->GetBinCenter( i );
    sum = 0.0;
    for(Int_t j=1;j<=nbinx0;j++){
      count = hcas0->GetBinContent( j );
      center = hcas0->GetBinCenter( j );

      sum += count * binw0
	*TMath::Exp(-(pos-center)*(pos-center)/(2.0*sigma*sigma))
	/TMath::Sqrt(2.0*pi)/sigma;
    }
    hcas1->SetBinContent( i, sum );
  }

  // Ex-compensated spectrum hcas1
  for(Int_t i=0;i<nbinx1;i++){
    Double_t bnumber = hexp0->FindBin( hcas1->GetBinCenter(i+1) );
    Double_t content = hexp0->GetBinContent( bnumber );
    Double_t value = hcas1->GetBinContent( i+1 ) * (content);
    hcas2->SetBinContent( i+1, value);
  }

  hcas0->Scale(50);
  hcas0->Draw("");
  hcas0->GetXaxis()->SetTitle("Excitation energy in {}^{20}Ne (MeV)");
  hcas0->GetYaxis()->SetTitle("Cross section (mb/sr/MeV)");
  hcas1->Scale(50);
  hcas1->SetLineColor(kGreen+2);
  hcas1->Draw("same");
  hcas2->SetLineColor(kRed);
  hcas2->Draw("same");
  hexp0->Scale(0.1);
  hexp0->Draw("same e0 hist");
  gPad->Update();
  
  TFile *ofile = new TFile("./root/ex0p6decay.root", "RECREATE");
  hexp0->Write();
  hexp1->Write();
  hcas0->Write();
  hcas1->Write();
  hcas2->Write();
  ofile->Close();

  return 0;
}
