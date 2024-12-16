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

// Глобальные переменные для хранения гистограмм
static TH1D *h1 = nullptr; // Гистограмма для data_1.dat
static TH1D *h2 = nullptr; // Гистограмма для data_2.dat

// Прототип функции минимизации
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);

// Функция для минимизации суммы хи-квадрат
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double chi2 = 0.0; // Начальная сумма хи-квадрат

    // Расчёт хи-квадрат для первой гистограммы (data_1)
    for (int i = 1; i <= h1->GetNbinsX(); i++) {
        double x = h1->GetBinCenter(i);  // Центр бина
        double y = h1->GetBinContent(i); // Значение в бине
        double ey = (y > 0) ? std::sqrt(y) : std::sqrt(par[0]); // Ошибка (корень из y или фон)
        // Значение модели: константа + гауссиан
        double fval = par[0] + par[1] * std::exp(-0.5 * (x - par[2]) * (x - par[2]) / (par[3] * par[3]));
        double diff = (y - fval) / ey; // Нормализованное отклонение
        chi2 += diff * diff; // Добавляем вклад в хи-квадрат
    }

    // Расчёт хи-квадрат для второй гистограммы (data_2)
    for (int i = 1; i <= h2->GetNbinsX(); i++) {
        double x = h2->GetBinCenter(i);  // Центр бина
        double y = h2->GetBinContent(i); // Значение в бине
        double ey = (y > 0) ? std::sqrt(y) : std::sqrt(par[0]); // Ошибка (корень из y или фон)
        double fval = par[0]; // Модель: только константа
        double diff = (y - fval) / ey; // Нормализованное отклонение
        chi2 += diff * diff; // Добавляем вклад в хи-квадрат
    }

    f = chi2; // Сохраняем итоговое значение хи-квадрат
}

void task10() {
    // Настройка отображения
    gStyle->SetOptStat(1); // Показ статистики
    gStyle->SetOptFit(0);  // Отключить автоматический вывод фитинга

    // Создание гистограмм
    h1 = new TH1D("h1", "Data 1", 100, 500, 600); // Гистограмма data_1
    h2 = new TH1D("h2", "Data 2", 100, 500, 600); // Гистограмма data_2

    // Загрузка данных из файла data_1.dat
    {
        std::ifstream in("data_1.dat");
        if (!in) {
            std::cerr << "Не удалось открыть data_1.dat" << std::endl;
            return;
        }
        double val;
        while (in >> val) {
            if (val >= 500 && val <= 600) h1->Fill(val); // Заполняем гистограмму
        }
        in.close();
    }

    // Загрузка данных из файла data_2.dat
    {
        std::ifstream in("data_2.dat");
        if (!in) {
            std::cerr << "Не удалось открыть data_2.dat" << std::endl;
            return;
        }
        double val;
        while (in >> val) {
            if (val >= 500 && val <= 600) h2->Fill(val); // Заполняем гистограмму
        }
        in.close();
    }

    // Создание объекта TMinuit для минимизации
    TMinuit *gMinuit = new TMinuit(4); // 4 параметра
    gMinuit->SetFCN(fcn); // Устанавливаем функцию минимизации

    // Задаём начальные параметры
    double p0 = 5.0;    // Константа
    double p1 = 40.0;   // Амплитуда Гаусса
    double p2 = 550.0;  // Среднее значение Гаусса
    double p3 = 10.0;   // Сигма Гаусса

    double step = 0.1; // Шаг для изменения параметров
    gMinuit->DefineParameter(0, "const", p0, step, 0, 1e6);
    gMinuit->DefineParameter(1, "amplitude", p1, step, 0, 1e6);
    gMinuit->DefineParameter(2, "mean", p2, step, 500, 600);
    gMinuit->DefineParameter(3, "sigma", p3, step, 0.1, 100);

    // Выполняем минимизацию
    gMinuit->Command("MIGRAD");
    gMinuit->Command("HESSE");

    // Извлечение параметров после подгонки
    double par[4], err[4];
    for (int i = 0; i < 4; i++) {
        gMinuit->GetParameter(i, par[i], err[i]);
    }

    // Расчёт хи-квадрат с финальными параметрами
    double final_chi2;
    {
        Int_t npar;
        Double_t *gin;
        Int_t iflag = 0;
        fcn(npar, gin, final_chi2, par, iflag);
    }

    // Число степеней свободы
    int ndof = (h1->GetNbinsX() + h2->GetNbinsX()) - 4;
    double chi2_per_ndf = final_chi2 / ndof;

    // Интеграл под Гауссом
    double gauss_integral = par[1] * par[3] * std::sqrt(2 * M_PI);

    // Проверка интеграла с помощью численного метода
    TF1 gauss("gauss", "[0]*exp(-0.5*((x-[1])*(x-[1]))/([2]*[2]))", 500, 600);
    gauss.SetParameters(par[1], par[2], par[3]);
    double numeric_integral = gauss.Integral(500, 600);

    // Вывод результатов
    std::cout << "Результаты подгонки:" << std::endl;
    std::cout << "const = " << par[0] << " ± " << err[0] << std::endl;
    std::cout << "amplitude = " << par[1] << " ± " << err[1] << std::endl;
    std::cout << "mean = " << par[2] << " ± " << err[2] << std::endl;
    std::cout << "sigma = " << par[3] << " ± " << err[3] << std::endl;
    std::cout << "Chi2 = " << final_chi2 << " на " << ndof << " степеней свободы => Chi2/ndf = " << chi2_per_ndf << std::endl;
    std::cout << "Число событий под Гауссом = " << gauss_integral << std::endl;
    std::cout << "Численный интеграл под Гауссом = " << numeric_integral << std::endl;

    // Отрисовка гистограмм и результатов подгонки
    TCanvas *c = new TCanvas("c", "Fits", 1200, 600);
    c->Divide(2, 1);

    c->cd(1);
    h1->SetMarkerStyle(20);
    h1->SetMarkerColor(kBlue);
    h1->SetTitle("Data 1 with fit; x; counts");
    h1->Draw("E");

    TF1 *f1 = new TF1("f1", "[0] + [1]*exp(-0.5*((x-[2])*(x-[2]))/([3]*[3]))", 500, 600);
    f1->SetParameters(par[0], par[1], par[2], par[3]);
    f1->SetLineColor(kRed);
    f1->Draw("SAME");

    TLegend *leg = new TLegend(0.65, 0.75, 0.9, 0.9);
    leg->AddEntry(h1, "Data 1", "p");
    leg->AddEntry(f1, "Const+Gauss fit", "l");
    leg->Draw();

    c->cd(2);
    h2->SetMarkerStyle(20);
    h2->SetMarkerColor(kBlue);
    h2->SetTitle("Data 2 with fit; x; counts");
    h2->Draw("E");

    TF1 *f2 = new TF1("f2", "[0]", 500, 600);
    f2->SetParameter(0, par[0]);
    f2->SetLineColor(kRed);
    f2->Draw("SAME");

    TLegend *leg2 = new TLegend(0.65, 0.75, 0.9, 0.9);
    leg2->AddEntry(h2, "Data 2", "p");
    leg2->AddEntry(f2, "Const fit", "l");
    leg2->Draw();

    c->Update();
    c->SaveAs("fit_results.png");
}
