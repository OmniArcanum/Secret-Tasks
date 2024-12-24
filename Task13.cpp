#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>
#include <iostream>
#include <fstream>

// ----------------------------------------------------------------------------
// Две гауссианы, каждая по площади = 1.
// ----------------------------------------------------------------------------
double Resolution(double x, double sigma1, double mu2, double sigma2)
{
   double g1 = TMath::Exp(-0.5 * (x/sigma1)*(x/sigma1))
             / (TMath::Sqrt(2*TMath::Pi()) * sigma1);

   double dx = x - mu2;
   double g2 = TMath::Exp(-0.5 * (dx/sigma2)*(dx/sigma2))
             / (TMath::Sqrt(2*TMath::Pi()) * sigma2);

   return (g1 + g2);
}

// ----------------------------------------------------------------------------
// Функция для фита
// ----------------------------------------------------------------------------
double BWConvolution(double *x, double *par)
{
   double bg = par[0] + par[1]*x[0];
   double A = par[2];
   double sigma1 = par[3];
   double mu2 = par[4];
   double sigma2 = par[5];

   double tMin = 3.0;
   double tMax = 3.2;
   int    Nstep = 200;
   double step = (tMax - tMin) / Nstep;

   double sum = 0.0;
   for (int i = 0; i < Nstep; i++) {
      double t = tMin + (i + 0.5) * step;
      double bwVal = TMath::BreitWigner(t, 3.0969, 0.000093);
      double rVal = Resolution(x[0] - t, sigma1, mu2, sigma2);
      sum += bwVal * rVal;
   }
   sum *= step;
   return bg + A * sum;
}

void task13()
{
    // Настройки стиля
    gStyle->SetOptStat(111); // Включаем только Entries, Mean и Std Dev
    gStyle->SetOptFit(1);    // Включаем только Chi2/NDF и параметры фитинга

    // 1) Открываем файл
    std::ifstream fin("m3piJPSI_cut.dat");
    if(!fin.is_open()){
        std::cerr<<"Нет файла m3piJPSI_cut.dat"<<std::endl;
        return;
    }

    // 2) Гистограмма
    TH1F *hMass = new TH1F("hMass","J/psi -> 3#pi; m (GeV/c^{2}); Events",400,3.0,3.2);
    double massVal;
    while (fin >> massVal) {
        if(!fin.good()) break;
        hMass->Fill(massVal);
    }
    fin.close();

    // 3) Функция для фита
    TF1 *fFit = new TF1("fFit", BWConvolution, 3.0, 3.2, 6);
    fFit->SetNpx(800);

    // Установка параметров
    fFit->SetParName(0,"p0_bg");
    fFit->SetParName(1,"p1_bg");
    fFit->SetParName(2,"A_signal");
    fFit->SetParName(3,"sigma1");
    fFit->SetParName(4,"mu2");
    fFit->SetParName(5,"sigma2");

    fFit->SetParameter(0,  0.0);    
    fFit->SetParameter(1,  0.0);    
    fFit->SetParameter(2,  500.0);  
    fFit->SetParameter(3,  0.003);  
    fFit->SetParameter(4,  0.0);    
    fFit->SetParameter(5,  0.003);  

    fFit->SetParLimits(3, 1e-5, 0.01); 
    fFit->SetParLimits(5, 1e-5, 0.01); 

    hMass->Fit(fFit,"R"); // "R" — ограничить диапазон

    // 5) Рисуем результат
    TCanvas *c1 = new TCanvas("c1","J/psi -> 3pi with BW * Gaus conv",900,700);
    c1->SetLogy(); // Логарифмическая шкала по Y
    hMass->Draw("E"); // Точки с ошибками

    fFit->SetLineColor(kRed);
    fFit->SetLineWidth(2);
    fFit->Draw("same");

    // Отдельно нарисуем фон (чёрная линия): p0 + p1*x
    TF1 *fBg = new TF1("fBg","[0] + [1]*x",3.0,3.2);
    fBg->SetParameters(fFit->GetParameter(0), fFit->GetParameter(1));
    fBg->SetLineColor(kBlack);
    fBg->SetLineStyle(2);
    fBg->Draw("same");

    // Убираем легенду
    delete c1->GetListOfPrimitives()->FindObject("TPave");

    // Сохраняем результат
    c1->SaveAs("task13_corrected_result.png");
}
