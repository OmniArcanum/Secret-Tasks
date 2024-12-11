#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TMath.h>
#include <iostream>

void task9() {
    // Открываем файл из задачи 7: newroot.root
    TFile *file7 = TFile::Open("newroot.root");
    if (!file7 || file7->IsZombie()) {
        std::cout << "Failed to open file newroot.root" << std::endl;
        return;
    }
    std::cout << "Successfully opened newroot.root from task 7." << std::endl;

    // Открываем файл из задачи 8: results.root
    TFile *file8 = TFile::Open("results.root", "UPDATE");
    if (!file8 || file8->IsZombie()) {
        std::cout << "Failed to open file results.root" << std::endl;
        return;
    }
    std::cout << "Successfully opened results.root from task 8." << std::endl;

    // Предположим, что в results.root у нас есть дерево с уже отобранными кандидатами π0
    // назовем его Pi0CandidatesTree. В нем должны быть ветви theta_pi0 и phi_pi0.
    TTree *candTree = (TTree*)file8->Get("Pi0CandidatesTree");
    if (!candTree) {
        std::cout << "Failed to find Pi0CandidatesTree in results.root" << std::endl;
        return;
    }

    // Подключаем ветви с углами кандидатов π0
    Float_t theta_pi0, phi_pi0;
    candTree->SetBranchAddress("theta_pi0", &theta_pi0);
    candTree->SetBranchAddress("phi_pi0", &phi_pi0);

    // Настраиваем бины по углам:
    Int_t nbins_theta = 18; // 0-π (180°), шаг 10°
    Int_t nbins_phi   = 36; // 0-2π (360°), шаг 10°

    // Создаем гистограммы
    TH1F *h_theta = new TH1F("h_theta", "Number of #pi^{0} Candidates vs Polar Angle; #theta [rad]; Number of Candidates", nbins_theta, 0, TMath::Pi());
    TH1F *h_phi   = new TH1F("h_phi", "Number of #pi^{0} Candidates vs Azimuthal Angle; #phi [rad]; Number of Candidates", nbins_phi, 0, 2*TMath::Pi());

    // Заполняем гистограммы, используя готовых кандидатов из Pi0CandidatesTree
    Long64_t nentries = candTree->GetEntries();
    std::cout << "Total number of pi0 candidates in Pi0CandidatesTree: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        candTree->GetEntry(i);
        h_theta->Fill(theta_pi0);
        h_phi->Fill(phi_pi0);
    }

    // Добавляем легенды, маркеры, рисуем ошибки
    TCanvas *c1 = new TCanvas("c1", "Pi0 Candidates vs Theta", 800, 600);
    h_theta->SetMarkerStyle(21);
    h_theta->SetMarkerColor(kBlue);
    h_theta->SetLineColor(kBlue);
    h_theta->Draw("E");

    TLegend *leg_theta = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg_theta->AddEntry(h_theta, "#pi^{0} candidates vs #theta", "lep");
    leg_theta->Draw();

    TCanvas *c2 = new TCanvas("c2", "Pi0 Candidates vs Phi", 800, 600);
    h_phi->SetMarkerStyle(21);
    h_phi->SetMarkerColor(kRed);
    h_phi->SetLineColor(kRed);
    h_phi->Draw("E");

    TLegend *leg_phi = new TLegend(0.7, 0.8, 0.9, 0.9);
    leg_phi->AddEntry(h_phi, "#pi^{0} candidates vs #phi", "lep");
    leg_phi->Draw();

    // Сохраняем результаты в новый подкаталог Task9 в results.root
    TDirectory *dir9 = file8->mkdir("Task9");
    dir9->cd();

    h_theta->Write();
    h_phi->Write();
    c1->Write("c1_theta");
    c2->Write("c2_phi");

    file8->Close();
    file7->Close();

    std::cout << "Task 9 analysis complete. Results saved in results.root under Task9 directory." << std::endl;
}
