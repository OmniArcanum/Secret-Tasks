#include <iostream>
#include <fstream>
#include <TCanvas.h>
#include <TH1D.h>
#include <TF1.h>
#include <TMath.h>
#include <TStyle.h>

// Определение функции фитинга
double fitFunction(double *x, double *par) {
    double poissonPart = par[0] * TMath::Poisson(x[0], par[1]);
    double breitWignerPart = par[4] * TMath::BreitWigner(x[0], par[2], par[3]);
    return poissonPart + breitWignerPart;
}

// Функция для заполнения гистограммы
void FillHistogram(TH1D *hist, const std::string &filename) {
    std::ifstream inFile(filename);
    if (!inFile.is_open()) {
        std::cerr << "Ошибка открытия файла: " << filename << std::endl;
        exit(-1);
    }

    double value;
    while (inFile >> value) {
        hist->Fill(value);
    }
    inFile.close();
}

// Основная функция
void task12() {
    // Настройка стиля
    gStyle->SetOptFit(1111);
    gStyle->SetOptStat(1111);

    // Создание гистограммы
    TH1D *hist = new TH1D("hist", "Distribution data from task10Nov.dat;mm;yields / 0.1 mm", 100, 0, 10);

    // Заполнение гистограммы
    FillHistogram(hist, "task10Nov.dat");

    // Создание функции фитинга
    TF1 *fitFunc = new TF1("fitFunc", fitFunction, 0, 10, 5);

    // Установка параметров
    fitFunc->SetParameter(0, 30);  // Poisson scale
    fitFunc->SetParameter(1, 0.2); // Poisson mean
    fitFunc->SetParameter(2, 5.4); // Breit Wigner mean
    fitFunc->SetParameter(3, 4.0); // Breit Wigner width
    fitFunc->SetParameter(4, 100); // Scale for Breit-Wigner

    // Установка ограничений параметров
    fitFunc->SetParLimits(0, 10, hist->GetSum()); // Poisson scale
    fitFunc->SetParLimits(1, 0.1, 0.3);           // Poisson mean
    fitFunc->SetParLimits(2, 5.0, 5.6);           // Breit Wigner mean
    fitFunc->SetParLimits(3, 3.5, 4.5);           // Breit Wigner width
    fitFunc->SetParLimits(4, 50, 200);            // Scale for Breit-Wigner

    // Фитинг гистограммы
    TFitResultPtr fitRes = hist->Fit(fitFunc, "SME");

    // Создание холста для отображения
    TCanvas *c1 = new TCanvas("c1", "Diagnostic Fit", 1200, 600);
    hist->SetLineColor(kBlue + 2);
    hist->SetLineWidth(2);
    hist->Draw("HIST");

    fitFunc->SetLineColor(kRed);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("SAME");

    // Сохранение графика
    c1->SaveAs("diagnostic_fit.png");

    // Вывод результатов в консоль
    std::cout << "Chi2/NDF = " << fitRes->Chi2() / fitRes->Ndf() << std::endl;
    for (int i = 0; i < 5; ++i) {
        std::cout << "Parameter " << i << ": " << fitFunc->GetParameter(i)
                  << " +/- " << fitFunc->GetParError(i) << std::endl;
    }
}
