#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

#include "TROOT.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TMinuit.h"
#include "TF1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLine.h"
#include "TLegend.h"

// Global variables
static TH1D *h1 = nullptr; // For data_1.dat
static TH1D *h2 = nullptr; // For data_2.dat

// We will fit parameters: p0=const, p1=amplitude of Gaussian, p2=mean, p3=sigma.
// f1(x) = p0 + p1 * exp(-(x-p2)^2/(2*p3^2))
// f2(x) = p0

// The combined chi2 to minimize
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    // par[0] = const
    // par[1] = amplitude of Gaussian
    // par[2] = mean
    // par[3] = sigma

    double chi2 = 0.0;
    
    // Chi2 from first histogram (data_1)
    for (int i=1; i<=h1->GetNbinsX(); i++) {
        double x    = h1->GetBinCenter(i);
        double y    = h1->GetBinContent(i);
        double ey   = (y>0)? std::sqrt(y): 1.0; // error = sqrt(N), avoid zero if empty bin
        double fval = par[0] + par[1]*std::exp(-0.5 * (x - par[2])*(x - par[2])/(par[3]*par[3]));
        double diff = (y - fval)/ey;
        chi2 += diff*diff;
    }

    // Chi2 from second histogram (data_2)
    for (int i=1; i<=h2->GetNbinsX(); i++) {
        double x    = h2->GetBinCenter(i);
        double y    = h2->GetBinContent(i);
        double ey   = (y>0)? std::sqrt(y): 1.0;
        // f2(x) = par[0]
        double fval = par[0];
        double diff = (y - fval)/ey;
        chi2 += diff*diff;
    }

    f = chi2;
}

void task10() {
    // Style adjustments
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);

    // Create histograms
    h1 = new TH1D("h1","Data 1",100,500,600);
    h2 = new TH1D("h2","Data 2",100,500,600);

    // Read data_1.dat
    {
        std::ifstream in("data_1.dat");
        if(!in) {
            std::cerr << "Cannot open data_1.dat" << std::endl;
            return;
        }
        double val;
        while (in >> val) {
            if(val>=500 && val<=600) h1->Fill(val);
        }
        in.close();
    }

    // Read data_2.dat
    {
        std::ifstream in("data_2.dat");
        if(!in) {
            std::cerr << "Cannot open data_2.dat" << std::endl;
            return;
        }
        double val;
        while (in >> val) {
            if(val>=500 && val<=600) h2->Fill(val);
        }
        in.close();
    }

    // Create TMinuit object with 4 parameters
    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn);

    // Set initial parameter guesses
    double p0     = 5.0;    // const offset guess
    double p1     = 40.0;   // amplitude guess for Gaussian
    double p2     = 550.0;  // mean guess
    double p3     = 10.0;   // sigma guess

    double step   = 0.1;
    gMinuit->DefineParameter(0,"const", p0, step, 0, 1e6);
    gMinuit->DefineParameter(1,"amplitude", p1, step, 0, 1e6);
    gMinuit->DefineParameter(2,"mean", p2, step, 500, 600);
    gMinuit->DefineParameter(3,"sigma", p3, step, 0.1, 100);

    // Perform minimization
    gMinuit->Command("MIGRAD");
    gMinuit->Command("HESSE");

    // Get fitted parameters
    double par[4], err[4];
    for (int i=0; i<4; i++) {
        gMinuit->GetParameter(i, par[i], err[i]);
    }

    // Compute chi2 with final parameters
    double final_chi2;
    {
        Int_t npar;
        Double_t *gin;
        Int_t iflag=0;
        fcn(npar, gin, final_chi2, par, iflag);
    }

    // Number of degrees of freedom:
    // total bins = 100 + 100 = 200
    // parameters = 4
    // ndof = 200 - 4 = 196
    int ndof = (h1->GetNbinsX() + h2->GetNbinsX()) - 4;
    double chi2_per_ndf = final_chi2/ndof;

    // Integral under Gaussian (number of events under gauss):
    // Integral of Gaussian A*exp(- (x-mu)^2/(2sigma^2)) dx = A * sigma * sqrt(2*pi)
    // Here A = par[1], so integral = par[1] * par[3] * sqrt(2*pi)
    double gauss_integral = par[1]*par[3]*std::sqrt(2*M_PI);

    std::cout << "Fit results:" << std::endl;
    std::cout << "const = " << par[0] << " ± " << err[0] << std::endl;
    std::cout << "amplitude = " << par[1] << " ± " << err[1] << std::endl;
    std::cout << "mean = " << par[2] << " ± " << err[2] << std::endl;
    std::cout << "sigma = " << par[3] << " ± " << err[3] << std::endl;
    std::cout << "Chi2 = " << final_chi2 << " for " << ndof << " d.o.f. => Chi2/ndf = " << chi2_per_ndf << std::endl;
    std::cout << "Events under Gaussian = " << gauss_integral << std::endl;

    // Now draw the histograms and overlay the fit results.
    TCanvas *c = new TCanvas("c","Fits",1200,600);
    c->Divide(2,1);

    c->cd(1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(kBlue);
    h1->SetTitle("Data 1 with fit; x; counts");
    h1->Draw("E");

    // Draw fit function for h1
    // We'll create a TF1 to represent f1
    TF1 *f1 = new TF1("f1","[0] + [1]*exp(-0.5*((x-[2])*(x-[2]))/([3]*[3]))",500,600);
    f1->SetParameters(par[0],par[1],par[2],par[3]);
    f1->SetLineColor(kRed);
    f1->Draw("SAME");

    // Add legend
    {
        TLegend *leg = new TLegend(0.65,0.75,0.9,0.9);
        leg->AddEntry(h1,"Data 1","p");
        leg->AddEntry(f1,"Const+Gauss fit","l");
        leg->Draw();
    }

    c->cd(2);
    h2->SetMarkerStyle(20);
    h2->SetMarkerColor(kBlue);
    h2->SetTitle("Data 2 with fit; x; counts");
    h2->Draw("E");

    // Draw fit function for h2 (just const = par[0])
    TF1 *f2 = new TF1("f2","[0]",500,600);
    f2->SetParameter(0,par[0]);
    f2->SetLineColor(kRed);
    f2->Draw("SAME");

    {
        TLegend *leg2 = new TLegend(0.65,0.75,0.9,0.9);
        leg2->AddEntry(h2,"Data 2","p");
        leg2->AddEntry(f2,"Const fit","l");
        leg2->Draw();
    }

    c->Update();
    // Save canvas to file
    c->SaveAs("fit_results.png");

    // Done
}
