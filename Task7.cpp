#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <iostream>

void task7()
{
    // Установка стиля графика
    gStyle->SetOptStat(1111);
    gStyle->SetOptFit(1111);

    // Открываем исходный файл
    TFile *inputFile = new TFile("m3pimc.root", "READ");
    if (!inputFile || inputFile->IsZombie()) {
        std::cout << "Ошибка открытия входного файла!" << std::endl;
        return;
    }

    // Получаем исходное дерево
    TTree *inputTree = (TTree*)inputFile->Get("h10");
    if (!inputTree) {
        std::cout << "Дерево h10 не найдено!" << std::endl;
        return;
    }

    // Создаем новый файл
    TFile *outputFile = new TFile("newroot.root", "RECREATE");

    // Активируем только нужные ветки
    inputTree->SetBranchStatus("*", 0);
    inputTree->SetBranchStatus("nph", 1);
    inputTree->SetBranchStatus("eph", 1);
    inputTree->SetBranchStatus("thetaph", 1);
    inputTree->SetBranchStatus("phiph", 1);
    inputTree->SetBranchStatus("Isrfilter", 1); // Или "lsrfilter", если это опечатка
    inputTree->SetBranchStatus("chi2_3p", 1);

    // Копируем дерево с условием отбора
    TTree *newTree = inputTree->CopyTree("Isrfilter == 1 && chi2_3p < 30");
    if (!newTree) {
        std::cout << "Ошибка при создании нового дерева!" << std::endl;
        return;
    }
    newTree->SetName("MyTree");

    // Деактивируем ветки отбора перед записью
    newTree->SetBranchStatus("Isrfilter", 0);
    newTree->SetBranchStatus("chi2_3p", 0);

    // Сохраняем новое дерево
    outputFile->cd();
    newTree->Write();

    // Получаем размеры файлов
    Long_t id_in, flags_in, modtime_in;
    Long64_t inputSize;
    gSystem->GetPathInfo("m3pimc.root", &id_in, &inputSize, &flags_in, &modtime_in);

    Long_t id_out, flags_out, modtime_out;
    Long64_t outputSize;
    gSystem->GetPathInfo("newroot.root", &id_out, &outputSize, &flags_out, &modtime_out);

    std::cout << "Коэффициент сжатия: " << static_cast<double>(inputSize) / static_cast<double>(outputSize) << std::endl;

    // Объявляем переменные для чтения данных
    Int_t nph;
    Float_t eph[100];  // Увеличиваем размер массива или используем динамическое выделение

    // Устанавливаем адреса веток
    newTree->SetBranchAddress("nph", &nph);
    newTree->SetBranchAddress("eph", eph);

    // Создаем гистограмму
    TH1F *hEnergy = new TH1F("hEnergy", "Energy Distribution of Photons;Energy (GeV);Entries", 100, 0, 10);
    hEnergy->SetLineColor(kBlue);
    hEnergy->SetFillColor(kBlue-10);

    Double_t minEnergy = 1e9;
    Double_t maxEnergy = -1e9;

    // Проверяем количество записей
    Long64_t entries = newTree->GetEntries();
    std::cout << "Количество событий после отбора: " << entries << std::endl;

    // Заполняем гистограмму
    for(Long64_t i = 0; i < entries; i++) {
        newTree->GetEntry(i);
        if (nph > 100) continue; // Пропускаем события с количеством фотонов больше размера массива
        for(Int_t j = 0; j < nph; j++) {
            if (eph[j] > 0) {
                hEnergy->Fill(eph[j]);
                if(eph[j] < minEnergy) minEnergy = eph[j];
                if(eph[j] > maxEnergy) maxEnergy = eph[j];
            }
        }
    }

    std::cout << "Минимальная энергия фотона: " << minEnergy << " GeV" << std::endl;
    std::cout << "Максимальная энергия фотона: " << maxEnergy << " GeV" << std::endl;

    // Создаем канвас
    TCanvas *c1 = new TCanvas("c1", "Energy Distribution", 800, 600);
    c1->SetGrid();
    c1->SetLogy(); // Устанавливаем логарифмическую шкалу по оси Y

    // Рисуем гистограмму
    hEnergy->Draw();

    // Определяем диапазон фитирования на основе данных гистограммы
    Double_t fitMin = 0.1; // Устанавливаем нижний предел фитирования
    Double_t fitMax = 10.0; // Верхний предел фитирования

  // Настраиваем и выполняем подгонку с общим показателем степени
TF1 *fit = new TF1("fit", "[0]*pow(x, [1])*exp(-[2]*x) + [3]*pow(x, [1])*exp(-[4]*x)", fitMin, fitMax);
fit->SetParameters(1e4, -1, 0.5, 1e3, 0.1); // Начальные параметры
fit->SetParNames("A1", "n", "k1", "A2", "k2");
fit->SetLineColor(kRed);

// Ограничиваем параметры, если необходимо
fit->SetParLimits(1, -5, 0); // n между -5 и 0
fit->SetParLimits(2, 0, 10); // k1 положительное
fit->SetParLimits(4, 0, 10); // k2 положительное

hEnergy->Fit(fit, "R"); // Фитирование в заданном диапазоне


    // Обновляем канвас
    c1->Modified();
    c1->Update();

    // Сохраняем канвас
    c1->SaveAs("energy_distribution.png");

    // Закрываем файлы
    //outputFile->Close();
    //inputFile->Close();

    std::cout << "Анализ завершен. Результаты сохранены в newroot.root" << std::endl;
}
