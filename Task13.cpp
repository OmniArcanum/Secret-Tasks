////////////////////////////////////////////////////////////
// fitJpsi_BW_conv_noLegend.C
//
// 1) Показываем только статистику гистограммы (Entries, Mean, Std Dev).
// 2) Убираем Legend с подписями линий.
// 3) Численная свёртка узкой BW (из TMath::BreitWigner) и (G1+G2).
////////////////////////////////////////////////////////////

#include <TROOT.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>
#include <iostream>
#include <fstream>

// --------------------------------------------------
// Две гауссианы (сумма), каждая по площади = 1
// --------------------------------------------------
double Resolution(double x, double sigma1, double mu2, double sigma2)
{
   double g1 = TMath::Exp(-0.5 * (x/sigma1)*(x/sigma1))
             / (TMath::Sqrt(2*TMath::Pi())*sigma1);

   double dx = x - mu2;
   double g2 = TMath::Exp(-0.5 * (dx/sigma2)*(dx/sigma2))
             / (TMath::Sqrt(2*TMath::Pi())*sigma2);

   return (g1 + g2);
}

// --------------------------------------------------
// BW свёртка с (G1+G2):
//
//  f(x) = (p0 + p1*x) + A * ∫ BW(t)*Resolution(x-t) dt
// --------------------------------------------------
double BWConvolution(double *x, double *par)
{
   // фон
   double background = par[0] + par[1]*x[0];
   // амплитуда
   double A = par[2];
   // параметры резолюции
   double sigma1 = par[3];
   double mu2    = par[4];
   double sigma2 = par[5];

   // интеграл по t от 3.0 до 3.2
   const double tMin = 3.0;
   const double tMax = 3.2;
   const int    Nstep= 200;
   double step  = (tMax - tMin)/(double)Nstep;

   double sum = 0.0;
   for(int i=0; i<Nstep; i++){
      double t = tMin + (i+0.5)*step;
      double bwVal  = TMath::BreitWigner(t, 3.0969, 0.000093);
      double resVal = Resolution(x[0] - t, sigma1, mu2, sigma2);
      sum += bwVal * resVal;
   }
   sum *= step;

   return background + A * sum;
}

void task13()
{
   // --- Настройка отображения статистики ---
   // Показываем только Entries, Mean, Std Dev (RMS)
   gStyle->SetOptStat(7);   // bits 0..2 => 1+2+4=7
   gStyle->SetOptFit(0);    // не показывать параметры фита

   // Открываем данные
   std::ifstream fin("m3piJPSI_cut.dat");
   if(!fin.is_open()){
      std::cerr<<"Нет файла m3piJPSI_cut.dat"<<std::endl;
      return;
   }

   // Гистограмма
   TH1F *hMass = new TH1F("hMass","J/psi -> 3#pi; m (GeV/c^{2}); Events",
                          400,3.0,3.2);

   double massVal;
   while(fin >> massVal){
      if(!fin.good()) break;
      hMass->Fill(massVal);
   }
   fin.close();

   // Фит-функция
   TF1 *fFit = new TF1("fFit", BWConvolution, 3.0, 3.2, 6);
   fFit->SetNpx(800);

   // Параметры: p0_bg, p1_bg, A, sigma1, mu2, sigma2
   fFit->SetParameter(0,  0.0);   // p0
   fFit->SetParameter(1,  0.0);   // p1
   fFit->SetParameter(2, 500.0);  // A
   fFit->SetParameter(3, 0.003);  // sigma1
   fFit->SetParameter(4, 0.0);    // mu2
   fFit->SetParameter(5, 0.003);  // sigma2

   // Ограничения
   fFit->SetParLimits(3, 1e-5, 0.01); // sigma1 in [0,10 МэВ]
   fFit->SetParLimits(5, 1e-5, 0.01); // sigma2 in [0,10 МэВ]

   // Фит
   hMass->Fit(fFit,"R");

   // Рисуем
   TCanvas *c1 = new TCanvas("c1","J/psi->3pi with BW*Gauss conv",900,700);
   c1->SetLogy();  // логарифмическая шкала

   // Гистограмма (только статистика)
   hMass->Draw("E"); 

   // Общая фит-функция (красная)
   fFit->SetLineColor(kRed);
   fFit->SetLineWidth(2);
   fFit->Draw("same");

   // Доп. кривая: фон (чёрная, пунктир)
   TF1 *fBg = new TF1("fBg","[0] + [1]*x",3.0,3.2);
   fBg->SetParameters(fFit->GetParameter(0),fFit->GetParameter(1));
   fBg->SetLineColor(kBlack);
   fBg->SetLineStyle(2);
   fBg->Draw("same");

   // Всё, легенду не рисуем
}
