#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TCanvas.h>
#include <iostream>
#include <vector>

void task8() {
    // Открываем файл с данными
    TFile *file = TFile::Open("newroot.root");
    if (!file || file->IsZombie()) {
        std::cout << "Failed to open file newroot.root" << std::endl;
        return;
    }

    // Получаем дерево MyTree
    TTree *tree = (TTree*)file->Get("MyTree");
    if (!tree) {
        std::cout << "Failed to find tree MyTree in the file" << std::endl;
        return;
    }

    // Объявляем переменные для чтения данных
    const Int_t maxnph = 100; // Максимальное количество фотонов в событии
    Int_t nph;
    Float_t eph[maxnph];
    Float_t thetaph[maxnph];
    Float_t phiph[maxnph];

    // Устанавливаем адреса ветвей
    tree->SetBranchAddress("nph", &nph);
    tree->SetBranchAddress("eph", eph);
    tree->SetBranchAddress("thetaph", thetaph);
    tree->SetBranchAddress("phiph", phiph);

    // Создаем гистограммы
    TH1F *h_invmass = new TH1F("h_invmass", "Invariant Mass of #pi^{0} Candidates; M_{#gamma#gamma} [GeV/c^{2}]; Number of Candidates", 100, 0, 0.3);
    TH1F *h_angle = new TH1F("h_angle", "Angle Between Photon Pairs; Angle [rad]; Number of Pairs", 100, 0, TMath::Pi());

    // Переменная для подсчёта общего количества кандидатов (по всем событиям)
    Long64_t total_candidates = 0;

    // Создаём дерево для записи кандидатов π0
    TTree *pi0Tree = new TTree("Pi0CandidatesTree", "Tree with selected pi0 candidates");
    Float_t out_mass, out_theta, out_phi;
    pi0Tree->Branch("mass_pi0", &out_mass, "mass_pi0/F");
    pi0Tree->Branch("theta_pi0", &out_theta, "theta_pi0/F");
    pi0Tree->Branch("phi_pi0", &out_phi, "phi_pi0/F");

    // Цикл по событиям
    Long64_t nentries = tree->GetEntries();
    std::cout << "Total number of events: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Вектор для хранения инвариантных масс кандидатов текущего события
        std::vector<Float_t> inv_masses;
        // Вектор для хранения четырехмерных векторов кандидатов pi0 текущего события
        std::vector<TLorentzVector> pi0_candidates;

        // Цикл по всем парам фотонов (комбинации без повторений)
        for (Int_t j = 0; j < nph - 1; j++) {
            for (Int_t k = j + 1; k < nph; k++) {
                // Энергии фотонов
                Float_t E1 = eph[j];
                Float_t E2 = eph[k];

                // Углы фотонов
                Float_t theta1 = thetaph[j];
                Float_t phi1   = phiph[j];

                Float_t theta2 = thetaph[k];
                Float_t phi2   = phiph[k];

                // Расчет компонент импульса фотонов
                TVector3 p1;
                p1.SetMagThetaPhi(E1, theta1, phi1);
                TVector3 p2;
                p2.SetMagThetaPhi(E2, theta2, phi2);

                // Создание четырехмерных векторов фотонов
                TLorentzVector photon1(p1, E1);
                TLorentzVector photon2(p2, E2);

                // Расчет инвариантной массы пары фотонов
                TLorentzVector pi0_candidate = photon1 + photon2;
                Float_t inv_mass = pi0_candidate.M();

                // Расчет угла между фотонами и заполнение гистограммы углов для всех пар
                Float_t angle = photon1.Angle(photon2.Vect());
                h_angle->Fill(angle);

                // Проверка условия на инвариантную массу
                if (inv_mass >= 0.1 && inv_mass <= 0.2) {
                    // Сохраняем инвариантную массу кандидата
                    inv_masses.push_back(inv_mass);
                    // Сохраняем сам четырехмерный вектор кандидата
                    pi0_candidates.push_back(pi0_candidate);
                }
            }
        }

        // Проверяем, сколько кандидатов мы нашли в этом событии
        if (inv_masses.size() == 2) {
            // Если ровно два кандидата, заполняем гистограмму инвариантной массы
            // и записываем кандидатов в дерево
            for (size_t c = 0; c < inv_masses.size(); c++) {
                Float_t mass = inv_masses[c];
                h_invmass->Fill(mass);
                total_candidates++;

                // Восстанавливаем направление pi0
                TVector3 p_pi0 = pi0_candidates[c].Vect();
                out_mass = mass;
                out_theta = p_pi0.Theta();
                out_phi   = p_pi0.Phi();
                pi0Tree->Fill();
            }
        }
        // Если кандидатов не 2, то мы игнорируем это событие в плане инвариантной массы и не записываем в дерево
    }

    // Вывод общего количества кандидатов
    std::cout << "Total number of pi0 candidates (from all events with exactly 2 candidates): " << total_candidates << std::endl;

    // Сохраняем гистограммы и дерево в файл
    TFile *outfile = new TFile("results.root", "RECREATE");
    h_invmass->Write();
    h_angle->Write();
    pi0Tree->Write();
    outfile->Close();

    std::cout << "Analysis complete. Results saved in file results.root" << std::endl;

    // Отображение гистограмм
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass", 800, 600);
    h_invmass->Draw();

    TCanvas *c2 = new TCanvas("c2", "Angle Between Photons", 800, 600);
    h_angle->Draw();

    // Тест угла, чтобы убедиться, что Angle() работает корректно
    {
        TVector3 p1(1.0, 0.0, 0.0);
        TVector3 p2(0.0, 1.0, 0.0);
        TLorentzVector photon1(p1, 1.0);
        TLorentzVector photon2(p2, 1.0);

        Float_t angle = photon1.Angle(photon2.Vect());
        std::cout << "Angle test: " << angle << " rad, expected ~" << TMath::Pi()/2 << " rad" << std::endl;
    }
}
