#include <stdio.h>
#include <iostream>
#include <iomanip>

using namespace std;

Int_t smeared(Double_t exNe=23.5){
  gROOT->SetBatch();
  
  const Double_t pi = TMath::Pi();
  const Double_t sigma = 0.22; // smeared by this sigma in MeV
  Int_t i,j;

  Char_t casname[256], outname[256], ffile[256];
  Double_t exenergy;
  exenergy = exNe;
  sprintf(casname, "./root/out%6.3f.root", exenergy);
  sprintf(outname, "./root/smeared%6.3f.root", exenergy);
  sprintf(ffile, "./root/to_adachi_0924.root");

  TFile *fin = new TFile( casname, "read" );
  TFile *fout = new TFile( outname, "recreate" );
  TFile *ff = new TFile( ffile, "read" );

  ofstream outrfile;
  outrfile.open("./dat/summaryRatio.dat", ofstream::out || ofstream::app);

  fin->cd();
  TH1D *h2org = (TH1D*)fin->Get("h2")->Clone();
  h2org->SetName("h2org");
  TTree *tree = (TTree*)fin->Get("tree")->Clone();
  tree->SetName("tree");
  fout->cd();
   
  Int_t nbinx = h2org->GetNbinsX();
  Double_t binw = h2org->GetBinCenter(2) - h2org->GetBinCenter(1);
  Double_t xmin = h2org->GetBinLowEdge(1);
  Double_t xmax = h2org->GetBinLowEdge(nbinx) + binw;

  Int_t totalevent;
  tree->SetBranchAddress("totalevent", &totalevent);
  tree->GetEntry(1);

  TH1D *h2prob = new TH1D("h2prob","Excitation energy of {}^{16}O", nbinx, xmin, xmax);
  for(i=1;i<nbinx;i++){
    Double_t cont;
    Double_t prob;
    cont = h2org->GetBinContent(i);
    prob = cont / totalevent;

    h2prob->SetBinContent( i, prob);
  }
  h2prob->GetXaxis()->SetTitle("Excitation energy in {}^{16}O (MeV)");
  h2prob->GetYaxis()->SetTitle("Probability");

  // Fujikawa file
  ff->cd();
  TH1F* hexp1 = (TH1F*)ff->Get("h43805")->Clone();
  TH1F* hexp2 = (TH1F*)ff->Get("h43905")->Clone();
  TH1F* hexp3 = (TH1F*)ff->Get("h44005")->Clone();
  hexp1->SetName("hexp1");
  hexp2->SetName("hexp2");
  hexp3->SetName("hexp3");
  hexp1->SetMarkerStyle(0);
  hexp2->SetMarkerStyle(0);
  hexp3->SetMarkerStyle(0);

  hexp1->SetMarkerColor(kBlue+2);
  hexp2->SetMarkerColor(kBlue+2);
  hexp3->SetMarkerColor(kBlue+2);

  hexp1->GetXaxis()->SetTitleFont(132);
  hexp1->GetXaxis()->SetLabelFont(132);
  hexp1->GetYaxis()->SetTitleFont(132);
  hexp1->GetYaxis()->SetLabelFont(132);
  hexp2->GetXaxis()->SetTitleFont(132);
  hexp2->GetXaxis()->SetLabelFont(132);
  hexp2->GetYaxis()->SetTitleFont(132);
  hexp2->GetYaxis()->SetLabelFont(132);
  hexp3->GetXaxis()->SetTitleFont(132);
  hexp3->GetXaxis()->SetLabelFont(132);
  hexp3->GetYaxis()->SetTitleFont(132);
  hexp3->GetYaxis()->SetLabelFont(132);

  // same binning with Fujikawa spectrum
  Int_t nbinxnew = hexp2->GetNbinsX();
  Double_t binwnew = hexp2->GetBinWidth(1);
  Double_t xminnew = hexp2->GetBinLowEdge(1);
  Double_t xmaxnew = hexp2->GetBinLowEdge(nbinxnew) + binwnew;
  TH1D *h2new = new TH1D("h2new", "Excitation energy in {}^{16}O, smeared",
			 nbinxnew, xminnew, xmaxnew);

  Double_t cont = 0.0;
  Double_t pos = 0.0, center = 0.0, sum = 0.0;
  for(i=1;i<=nbinxnew;i++){
    pos = h2new->GetBinCenter( i );
    sum = 0.0;
    for(j=0;j<nbinx;j++){
      cont   = h2prob->GetBinContent( j );
      center = h2prob->GetBinCenter( j );

      sum += cont * binw
	*TMath::Exp(-(pos-center)*(pos-center)/(2.0*sigma*sigma))
	/TMath::Sqrt(2.0*pi)/sigma;
    }

    h2new->SetBinContent( i, sum );
  }
  h2new->GetXaxis()->SetTitle("Excitation energy in {}^{16}O (MeV)");
  h2new->GetYaxis()->SetTitle("Probability");

  h2org->SetLineColor(kRed);
  h2new->SetLineColor(kRed);
  h2new->SetLineWidth(2);


  // Ratio calculation
  // 16O excitated states position in MeV
  const Double_t peak[4] = {0.0, 6.05, 6.92, 15.10};
  Double_t integral = 0.0;
  Double_t intw = 0.32;
  Double_t ratio[4] = {0.0, 0.0, 0.0, 0.0};
  for(i=1;i<=nbinxnew;i++){
    integral += h2new->GetBinContent(i) * binwnew;
  }
  for(i=0;i<4;i++){
    Int_t leftBin, rightBin;
    leftBin = h2new->FindBin( peak[i] - intw );
    rightBin = h2new->FindBin( peak[i] + intw );
    ratio[i] = h2new->Integral( leftBin, rightBin ) * binwnew / integral;
    // cout << setprecision(4) << peak[i] << " MeV : "
    // 	 << setprecision(4) << ratio[i] << endl;
    // printf("%5.2f MeV: %6.4f\n", peak[i], ratio[i]);
  }
  outrfile << setprecision(4) << scientific;
  outrfile << exenergy;
  for(i=0;i<4;i++){
    outrfile << "   " << ratio[i];
  }
  outrfile << endl;

  // file write
  fout->cd();
  h2org->Write();
  h2prob->Write();
  tree->Write();
  h2new->Write();
  hexp1->Write();
  hexp2->Write();
  hexp3->Write();

  fin->Close();
  fout->Close();

  outrfile.close();
  
  return 0;
}
