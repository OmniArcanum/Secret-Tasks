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
#include "TLegend.h"
#include "TMath.h"
#include "TGraph.h"

// Глобальные переменные для fcn
static TH1D *g_h1 = nullptr; 
static TH1D *g_h2 = nullptr; 

// FCN для MLE (приближенный подход: mu_i = f(x_center)*binWidth)
void fcn_mle(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double neg2logL = 0.0;

    double cst   = par[0];
    double A     = par[1];
    double mean  = par[2];
    double sigma = par[3];

    if (sigma <= 0) {
        f = 1e20; 
        return;
    }

    double bw1 = (g_h1->GetXaxis()->GetXmax() - g_h1->GetXaxis()->GetXmin()) / g_h1->GetNbinsX();
    for (int i = 1; i <= g_h1->GetNbinsX(); i++) {
        double n = g_h1->GetBinContent(i);
        double x = g_h1->GetBinCenter(i);
        double mu = (cst + A*std::exp(-0.5*(x-mean)*(x-mean)/(sigma*sigma)))*bw1;
        if (mu <= 0) {f=1e20; return;}
        if (n>0) neg2logL += 2*(mu - n*std::log(mu));
        else neg2logL += 2*mu;
    }

    double bw2 = (g_h2->GetXaxis()->GetXmax() - g_h2->GetXaxis()->GetXmin()) / g_h2->GetNbinsX();
    for (int i = 1; i <= g_h2->GetNbinsX(); i++) {
        double n = g_h2->GetBinContent(i);
        double mu = cst*bw2;
        if (mu <=0) {f=1e20; return;}
        if (n>0) neg2logL += 2*(mu - n*std::log(mu));
        else neg2logL += 2*mu;
    }

    f = neg2logL;
}

// Функция для фита и получения N_signal при заданном числе бинов
double FitAndGetNsignal(int nbins, double xlow, double xup) {
    TH1D *h1 = new TH1D(Form("h1_temp_%d",nbins),"Data 1", nbins, xlow, xup);
    TH1D *h2 = new TH1D(Form("h2_temp_%d",nbins),"Data 2", nbins, xlow, xup);

    {
        std::ifstream in1("data_1.dat");
        double val;
        while (in1 >> val) {
            if(val>=xlow && val<=xup) h1->Fill(val);
        }
    }

    {
        std::ifstream in2("data_2.dat");
        double val;
        while (in2 >> val) {
            if(val>=xlow && val<=xup) h2->Fill(val);
        }
    }

    g_h1 = h1;
    g_h2 = h2;

    TMinuit gMinuit(4);
    gMinuit.SetFCN(fcn_mle);

    double p0=5, p1=40, p2=550, p3=10;
    double step=0.1;
    gMinuit.DefineParameter(0,"const", p0, step, 0, 1e6);
    gMinuit.DefineParameter(1,"amplitude", p1, step, 0, 1e6);
    gMinuit.DefineParameter(2,"mean", p2, step, 500, 600);
    gMinuit.DefineParameter(3,"sigma", p3, step, 0.1, 100);

    gMinuit.Command("MIGRAD");
    gMinuit.Command("HESSE");

    double par[4], err[4];
    for (int i=0; i<4; i++) gMinuit.GetParameter(i, par[i], err[i]);

    double N_signal = par[1]*par[3]*std::sqrt(2*M_PI);

    delete h1;
    delete h2;
    return N_signal;
}

void task11() {
    gStyle->SetOptStat(1);
    gStyle->SetOptFit(0);

    double xlow=500, xup=600;
    int nbins=100;

    TH1D *h1 = new TH1D("h1_11","Data 1", nbins, xlow, xup);
    TH1D *h2 = new TH1D("h2_11","Data 2", nbins, xlow, xup);

    {
        std::ifstream in1("data_1.dat");
        double val;
        while (in1 >> val) {
            if(val>=xlow && val<=xup) h1->Fill(val);
        }

        std::ifstream in2("data_2.dat");
        while (in2 >> val) {
            if(val>=xlow && val<=xup) h2->Fill(val);
        }
    }

    g_h1 = h1;
    g_h2 = h2;

    TMinuit gMinuit(4);
    gMinuit.SetFCN(fcn_mle);

    double p0=5, p1=40, p2=550, p3=10;
    double step=0.1;
    gMinuit.DefineParameter(0,"const", p0, step, 0, 1e6);
    gMinuit.DefineParameter(1,"amplitude", p1, step, 0, 1e6);
    gMinuit.DefineParameter(2,"mean", p2, step, 500, 600);
    gMinuit.DefineParameter(3,"sigma", p3, step, 0.1, 100);

    gMinuit.Command("MIGRAD");
    gMinuit.Command("HESSE");

    double par[4], err[4];
    for (int i=0; i<4; i++) gMinuit.GetParameter(i, par[i], err[i]);

    double N_signal = par[1]*par[3]*std::sqrt(2*M_PI);
    std::cout << "N_signal (при nbins=100) = " << N_signal << std::endl;

    // Рисуем контуры ошибок:
    // Сначала контур для большего интервала (например, ~2σ)
    TCanvas *c3 = new TCanvas("c3", "Errors Contour", 800, 600);
    gMinuit.SetErrorDef(2.25); // примерно 2σ уровень
    TGraph *gr2 = (TGraph*)gMinuit.Contour(40, 2, 3); // параметры mean(2), sigma(3)
    gr2->SetTitle("Errors Contour");
    gr2->SetFillColor(42);
    gr2->Draw("A lf");  // оси, линия, закраска

    // Теперь узкий контур (~1σ)
    gMinuit.SetErrorDef(1.0);
    TGraph *gr1 = (TGraph*)gMinuit.Contour(40, 2, 3);
    gr1->Draw("C"); // сглаженная линия поверх

    c3->Update();
    c3->SaveAs("errors_contour.png");

    // Зависимость N_signal от числа бинов
    std::vector<int> nbins_list={50,100,150,200};
    std::vector<double> nsignal_vals;
    for (auto nb : nbins_list) {
        double nsig = FitAndGetNsignal(nb,xlow,xup);
        nsignal_vals.push_back(nsig);
    }

    TCanvas *c_dep = new TCanvas("c_dep","N_signal vs nbins",600,600);
    TGraph *gr = new TGraph((int)nbins_list.size());
    for (size_t i=0;i<nbins_list.size();i++){
        gr->SetPoint((int)i, nbins_list[i], nsignal_vals[i]);
    }
    gr->SetTitle("N_signal vs nbins;nbins;N_signal");
    gr->SetMarkerStyle(20);
    gr->SetMarkerColor(kBlue);
    gr->Draw("AP");
    c_dep->SaveAs("Nsignal_vs_nbins.png");

    // Отрисовка фита
    TCanvas *c = new TCanvas("c","Fit results",1200,600);
    c->Divide(2,1);

    double binWidth = (xup - xlow)/nbins;

    c->cd(1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(kBlue);
    h1->SetTitle("Data 1; x; counts");
    h1->Draw("E");

    TF1 *f1 = new TF1("f1",
      "[0]*[4] + [1]*[4]*exp(-0.5*((x-[2])*(x-[2]))/([3]*[3]))",xlow,xup);
    f1->SetParameters(par[0],par[1],par[2],par[3],binWidth);
    f1->SetLineColor(kRed);
    f1->Draw("SAME");

    {
        TLegend *leg = new TLegend(0.65,0.75,0.9,0.9);
        leg->AddEntry(h1,"Data 1","p");
        leg->AddEntry(f1,"const+Gauss (MLE fit)","l");
        leg->Draw();
    }

    c->cd(2);
    h2->SetMarkerStyle(20);
    h2->SetMarkerColor(kBlue);
    h2->SetTitle("Data 2; x; counts");
    h2->Draw("E");

    TF1 *f2 = new TF1("f2","[0]*[1]", xlow, xup);
    f2->SetParameters(par[0], binWidth);
    f2->SetLineColor(kRed);
    f2->Draw("SAME");

    {
        TLegend *leg2 = new TLegend(0.65,0.75,0.9,0.9);
        leg2->AddEntry(h2,"Data 2","p");
        leg2->AddEntry(f2,"const (MLE fit)","l");
        leg2->Draw();
    }

    c->Update();
    c->SaveAs("fit_results_mle_11_contour.png");

    // Удаляем гистограммы
    delete h1;
    delete h2;
}
