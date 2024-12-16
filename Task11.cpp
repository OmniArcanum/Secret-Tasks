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

// Глобальные переменные для хранения гистограмм
// Эти переменные используются в функции правдоподобия
static TH1D *g_h1 = nullptr; // Первая гистограмма (Data 1)
static TH1D *g_h2 = nullptr; // Вторая гистограмма (Data 2)

// Функция для минимизации (FCN)
// Используется метод максимального правдоподобия для данных, распределенных по Пуассону
void fcn_poisson(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
    double neg2logL = 0.0; // -2 * логарифм функции правдоподобия

    // Извлечение параметров из массива par
    double cst = par[0];    // Постоянный фон
    double A = par[1];      // Амплитуда нормального распределения
    double mean = par[2];   // Среднее значение нормального распределения
    double sigma = par[3];  // Стандартное отклонение нормального распределения

    // Проверка допустимых значений параметров
    if (sigma <= 0) {
        f = 1e20; // Устанавливаем большое значение функции для исключения
        return;
    }

    // Вычисление ширины бинa для первой гистограммы
    double bw1 = (g_h1->GetXaxis()->GetXmax() - g_h1->GetXaxis()->GetXmin()) / g_h1->GetNbinsX();

    // Проходим по всем бинам первой гистограммы
    for (int i = 1; i <= g_h1->GetNbinsX(); ++i) {
        double n = g_h1->GetBinContent(i); // Количество событий в бине
        double x = g_h1->GetBinCenter(i); // Центр бинa
        double mu = (cst + A * std::exp(-0.5 * (x - mean) * (x - mean) / (sigma * sigma))) * bw1; // Ожидаемое число событий

        // Проверка на корректность mu
        if (mu > 0) {
            neg2logL -= n * std::log(mu) - mu; // Добавляем вклад бинa в логарифм правдоподобия
        } else {
            f = 1e20; // Исключаем некорректное значение
            return;
        }
    }

    // Аналогично вычисляем вклад второй гистограммы
    double bw2 = (g_h2->GetXaxis()->GetXmax() - g_h2->GetXaxis()->GetXmin()) / g_h2->GetNbinsX();
    for (int i = 1; i <= g_h2->GetNbinsX(); ++i) {
        double n = g_h2->GetBinContent(i);
        double mu = cst * bw2;
        if (mu > 0) {
            neg2logL -= n * std::log(mu) - mu;
        } else {
            f = 1e20;
            return;
        }
    }

    // Итоговая функция минимизации
    f = neg2logL * 2.0;
}

// Функция для подгонки и расчета N_signal
// Возвращает число событий под нормальным распределением
double FitAndGetNsignal(int nbins, double xlow, double xup) {
    // Создаем временные гистограммы для данных
    TH1D *h1 = new TH1D(Form("h1_temp_%d", nbins), "Data 1", nbins, xlow, xup);
    TH1D *h2 = new TH1D(Form("h2_temp_%d", nbins), "Data 2", nbins, xlow, xup);

    // Заполняем первую гистограмму из файла data_1.dat
    {
        std::ifstream in1("data_1.dat");
        double val;
        while (in1 >> val) {
            if (val >= xlow && val <= xup) h1->Fill(val);
        }
    }

    // Заполняем вторую гистограмму из файла data_2.dat
    {
        std::ifstream in2("data_2.dat");
        double val;
        while (in2 >> val) {
            if (val >= xlow && val <= xup) h2->Fill(val);
        }
    }

    // Передаем гистограммы в глобальные переменные
    g_h1 = h1;
    g_h2 = h2;

    // Инициализируем объект TMinuit для минимизации
    TMinuit gMinuit(4); // 4 параметра
    gMinuit.SetFCN(fcn_poisson);

    // Устанавливаем начальные значения параметров и шаги
    double p0 = 5, p1 = 40, p2 = 550, p3 = 10; // Начальные значения
    double step = 0.1; // Шаги минимизации
    gMinuit.DefineParameter(0, "const", p0, step, 0, 1e6);
    gMinuit.DefineParameter(1, "amplitude", p1, step, 0, 1e6);
    gMinuit.DefineParameter(2, "mean", p2, step, 500, 600);
    gMinuit.DefineParameter(3, "sigma", p3, step, 0.1, 100);

    // Выполняем минимизацию
    gMinuit.Command("MIGRAD");
    gMinuit.Command("HESSE");

    // Извлекаем параметры подгонки
    double par[4], err[4];
    for (int i = 0; i < 4; ++i) gMinuit.GetParameter(i, par[i], err[i]);

    // Вычисляем N_signal как площадь под кривой гауссовского распределения
    double N_signal = par[1] * par[3] * std::sqrt(2 * M_PI);

    // Очищаем временные гистограммы
    delete h1;
    delete h2;

    return N_signal;
}

// Главная функция
void task11() {
    // Устанавливаем стили графиков
    gStyle->SetOptStat(1); // Отображение статистики
    gStyle->SetOptFit(0);  // Отключение результатов подгонки на графике

    double xlow = 500, xup = 600; // Диапазон гистограмм
    int nbins = 100; // Количество бинов

    // Создаем гистограммы для данных
    TH1D *h1 = new TH1D("h1_11", "Data 1", nbins, xlow, xup);
    TH1D *h2 = new TH1D("h2_11", "Data 2", nbins, xlow, xup);

    // Чтение и заполнение данных для гистограмм
    {
        std::ifstream in1("data_1.dat");
        double val;
        while (in1 >> val) {
            if (val >= xlow && val <= xup) h1->Fill(val);
        }

        std::ifstream in2("data_2.dat");
        while (in2 >> val) {
            if (val >= xlow && val <= xup) h2->Fill(val);
        }
    }

    g_h1 = h1; // Передаем глобальной переменной
    g_h2 = h2;

    // Минимизация и извлечение параметров
    TMinuit gMinuit(4);
    gMinuit.SetFCN(fcn_poisson);

    double p0 = 5, p1 = 40, p2 = 550, p3 = 10;
    double step = 0.1;
    gMinuit.DefineParameter(0, "const", p0, step, 0, 1e6);
    gMinuit.DefineParameter(1, "amplitude", p1, step, 0, 1e6);
    gMinuit.DefineParameter(2, "mean", p2, step, 500, 600);
    gMinuit.DefineParameter(3, "sigma", p3, step, 0.1, 100);

    gMinuit.Command("MIGRAD");
    gMinuit.Command("HESSE");

    double par[4], err[4];
    for (int i = 0; i < 4; ++i) gMinuit.GetParameter(i, par[i], err[i]);

    double N_signal = par[1] * par[3] * std::sqrt(2 * M_PI);
    std::cout << "N_signal (nbins=100) = " << N_signal << std::endl;
}
