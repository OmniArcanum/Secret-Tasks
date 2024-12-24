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
#include "TMath.h"
#include "TGraph.h"

static TH1D *h1_11 = nullptr; 
static TH1D *h2_11 = nullptr; 
static std::vector<double> data1;
static std::vector<double> data2;

// Функция вычисления -2 ln(L)
void fcn_11(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double neg2lnL = 0.0;

    // Модель для первого спектра: par[0] (фон/бин) + гаусс (par[1] - общее число сигналов)
    // h1_11: M = const + amplitude*(Gauss-плотность)
    for (int i = 1; i <= h1_11->GetNbinsX(); i++) {
        double x = h1_11->GetBinCenter(i);
        int N    = (int)h1_11->GetBinContent(i);

        // Гауссов компонент
        double gauss_val = (par[1]/(par[3]*std::sqrt(2*M_PI))) *
                           std::exp(-0.5*(x - par[2])*(x - par[2])/(par[3]*par[3]));
        // Суммарная модель
        double M = par[0] + gauss_val;
        if (M <= 0) M = 1e-9;

        double p = TMath::Poisson(N, M);
        neg2lnL += (p>0) ? -2.0 * std::log(p) : 1e10;
    }

    // Модель для второго спектра: M = const (нет сигнала)
    for (int i = 1; i <= h2_11->GetNbinsX(); i++) {
        int N = (int)h2_11->GetBinContent(i);
        double M = par[0];
        if (M <= 0) M = 1e-9;

        double p = TMath::Poisson(N, M);
        neg2lnL += (p>0) ? -2.0 * std::log(p) : 1e10;
    }

    f = neg2lnL;
}

// Функция для фитирования с разным количеством бинов
double FitWithBinning_11(int nbins) {
    // Создаем локальные гистограммы с уникальными именами
    TH1D *h1local = new TH1D(Form("h1_11_%d",nbins),"Data 1 local",nbins,500,600);
    TH1D *h2local = new TH1D(Form("h2_11_%d",nbins),"Data 2 local",nbins,500,600);

    for (auto val : data1) {
        if(val>=500 && val<=600) h1local->Fill(val);
    }
    for (auto val : data2) {
        if(val>=500 && val<=600) h2local->Fill(val);
    }

    // Сохраняем глобальные указатели
    TH1D *oldh1 = h1_11;
    TH1D *oldh2 = h2_11;

    // Переключаем на локальные для фита
    h1_11 = h1local;
    h2_11 = h2local;

    TMinuit *gMin = new TMinuit(4);
    gMin->SetFCN(fcn_11);

    double p0=5.0;   // фон/бин
    double p1=40.0;  // общее число сигналов
    double p2=550.0; // среднее гаусса
    double p3=10.0;  // sigma
    double step=0.1;
    gMin->DefineParameter(0,"const", p0, step, 0, 1e6);
    gMin->DefineParameter(1,"amplitude", p1, step, 0, 1e6);
    gMin->DefineParameter(2,"mean", p2, step, 500, 600);
    gMin->DefineParameter(3,"sigma", p3, step, 0.1, 100);

    gMin->Command("MIGRAD");
    gMin->Command("HESSE");

    double par[4], err[4];
    for (int i=0;i<4;i++){
        gMin->GetParameter(i,par[i],err[i]);
    }

    double N_signal = par[1];

    // Возвращаем глобальные гистограммы
    h1_11 = oldh1;
    h2_11 = oldh2;

    // Удаляем временные
    delete h1local;
    delete h2local;
    delete gMin;

    return N_signal;
}


void task11() {
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    // 1) Считываем данные
    {
        std::ifstream in("data_1.dat");
        if(!in) {
            std::cerr << "Не удалось открыть data_1.dat" << std::endl;
            return;
        }
        double val;
        while(in >> val) {
            data1.push_back(val);
        }
        in.close();
    }

    {
        std::ifstream in("data_2.dat");
        if(!in) {
            std::cerr << "Не удалось открыть data_2.dat" << std::endl;
            return;
        }
        double val;
        while(in >> val) {
            data2.push_back(val);
        }
        in.close();
    }

    // 2) Основной фит с 100 бинами
    int nbins_default = 100;
    h1_11 = new TH1D("h1_11","Data 1",nbins_default,500,600);
    h2_11 = new TH1D("h2_11","Data 2",nbins_default,500,600);

    for (auto val : data1) {
        if(val>=500 && val<=600) h1_11->Fill(val);
    }
    for (auto val : data2) {
        if(val>=500 && val<=600) h2_11->Fill(val);
    }

    TMinuit *gMinuit = new TMinuit(4);
    gMinuit->SetFCN(fcn_11);

    double p0=5.0;    // фон/бин
    double p1=40.0;   // общее число сигналов
    double p2=550.0;  // среднее гаусса
    double p3=10.0;   // sigma
    double step=0.1;

    gMinuit->DefineParameter(0,"const", p0, step, 0, 1e6);
    gMinuit->DefineParameter(1,"amplitude", p1, step, 0, 1e6);
    gMinuit->DefineParameter(2,"mean", p2, step, 500, 600);
    gMinuit->DefineParameter(3,"sigma", p3, step, 0.1, 100);

    gMinuit->Command("MIGRAD");
    gMinuit->Command("HESSE");

    double par[4], err[4];
    for (int i=0; i<4; i++){
        gMinuit->GetParameter(i,par[i],err[i]);
    }

    double N_signal     = par[1];
    double N_signal_err = err[1];

    // Пример вычисления общего фона и сигнала (как в программе (1)):
    double bkgTotal = par[0]*h1_11->GetNbinsX(); // фон = const * кол-во бинов
    double sigTotal = N_signal;                  // сигнал = amplitude

    // 3) Вывод результатов в терминал
    std::cout << "================ Fit results ================" << std::endl;
    std::cout << "N_signal (fitted)   = " << N_signal << " ± " << N_signal_err << std::endl;
    std::cout << "Total background    = " << bkgTotal << std::endl;
    std::cout << "=============================================" << std::endl;

    // Можно при желании сохранить в текстовый файл (примерно как в программе (1)):
    /*
    {
       std::ofstream ofile("result_fit.txt");
       ofile << "=== Fit results ===" << std::endl;
       ofile << "N_signal = " << N_signal << " ± " << N_signal_err << std::endl;
       ofile << "Bkg_total= " << bkgTotal  << std::endl;
       ofile.close();
    }
    */

    // 4) Построение гистограмм с подгонкой (как было у вас)
    TCanvas *c = new TCanvas("c","Fits",1200,600);
    c->Divide(2,1);

    // Левая часть: Data 1 + Gauss+Const
    c->cd(1);
    h1_11->SetMarkerStyle(20);
    h1_11->SetMarkerColor(kBlue);
    h1_11->SetTitle("Data 1 with ML fit; x; counts");
    h1_11->Draw("E");

    TF1 *f1 = new TF1("f1","[0] + ([1]/([3]*sqrt(2*pi)))*exp(-0.5*((x-[2])*(x-[2]))/([3]*[3]))",500,600);
    f1->SetParameters(par[0],par[1],par[2],par[3]);
    f1->SetLineColor(kRed);
    f1->Draw("SAME");

    {
        TLegend *leg = new TLegend(0.65,0.75,0.9,0.9);
        leg->AddEntry(h1_11,"Data 1","p");
        leg->AddEntry(f1,"Const+Gauss (MLE)","l");
        leg->Draw();
    }

    // Правая часть: Data 2 + Const
    c->cd(2);
    h2_11->SetMarkerStyle(20);
    h2_11->SetMarkerColor(kBlue);
    h2_11->SetTitle("Data 2 with ML fit; x; counts");
    h2_11->Draw("E");

    TF1 *f2 = new TF1("f2","[0]",500,600);
    f2->SetParameter(0,par[0]);
    f2->SetLineColor(kRed);
    f2->Draw("SAME");

    {
        TLegend *leg2 = new TLegend(0.65,0.75,0.9,0.9);
        leg2->AddEntry(h2_11,"Data 2","p");
        leg2->AddEntry(f2,"Const (MLE)","l");
        leg2->Draw();
    }

    c->Update();
    c->SaveAs("ml_fit_results.png");

    // 5) Строим зависимость N_signal от числа бинов
    std::vector<int> binNumbers = {50,75,100,125,150};
    std::vector<double> Nsignals;
    for (auto nb : binNumbers) {
        double Nsig = FitWithBinning_11(nb);
        Nsignals.push_back(Nsig);
    }

    TCanvas *c2 = new TCanvas("c2","N_signal vs bins",800,600);
    TGraph *gNvsBins = new TGraph((int)binNumbers.size());
    for (size_t i = 0; i < binNumbers.size(); i++) {
        gNvsBins->SetPoint((int)i, binNumbers[i], Nsignals[i]);
    }
    gNvsBins->SetTitle("N_{signal} vs number of bins;number of bins;N_{signal}");
    gNvsBins->SetMarkerStyle(20);
    gNvsBins->Draw("AP");
    c2->SaveAs("Nsignal_vs_bins.png");

    // 6) Построение контуров ошибок
    TCanvas *c3 = new TCanvas("c3","Errors Contour",600,600);

    // Сначала внутренний контур (по умолчанию ErrorDef=1)
    TGraph *gr_inner = (TGraph*)gMinuit->Contour(40,1,2); // (amplitude, mean)

    // Теперь устанавливаем другой ErrorDef для внешнего контура
    gMinuit->SetErrorDef(2.25);
    TGraph *gr_outer = (TGraph*)gMinuit->Contour(40,1,2);

    if (gr_outer) {
        gr_outer->SetTitle("Errors Contour; amplitude; mean"); // Исправили подпись
        gr_outer->SetFillColor(42);
        gr_outer->Draw("ALF");
        if(gr_inner) gr_inner->Draw("C");
        c3->SaveAs("errors_contour.png");
    } else {
        std::cerr << "Не удалось построить внешний контур" << std::endl;
    }
}
