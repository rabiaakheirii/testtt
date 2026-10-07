#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TMultiGraph.h>
#include <TLegend.h>
#include <TF1.h>
#include <TStyle.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cmath>

// массы частиц: π+, π–, K+, K–, p, p̄
const double masses[6] = {
    0.13957, 0.13957, 0.49367, 0.49367, 0.93827, 0.93827
};

// Хранение графиков: [эксперимент][частица][центральность]
TGraphErrors* gr[2][6][10]; // допустим, не более 10 центральностей
int   nCent[2][6];          // реальное число центральностей
int   nPoints[2][6];        // число точек по каждой частицы и эксперименту

void ReadSpectra(const char* filename, int iexp) {
    std::ifstream f(filename);
    if (!f.is_open()) {
        std::cerr << "Error opening " << filename << std::endl;
        return;
    }

    for (int part = 0; part < 6; ++part) {
        int N;
        if (!(f >> N)) break;
        nPoints[iexp][part] = N;
        std::string line;
        std::getline(f, line); // Skip the rest of the line

        // Read the first line of the block
        std::getline(f, line);
        std::istringstream iss(line);
        std::vector<double> tokens;
        double x;
        while (iss >> x) tokens.push_back(x);

        // Number of centralities
        int T = tokens.size();
        int M = (T - 1) / 2;
        nCent[iexp][part] = M;

        // Prepare buffers
        std::vector<double> pT(N), mT(N), xerr(N, 0.05);
        std::vector<std::vector<double>> s(M, std::vector<double>(N));
        std::vector<std::vector<double>> serr(M, std::vector<double>(N));

        // First point
        pT[0] = tokens[0];
        mT[0] = sqrt(pT[0] * pT[0] + masses[part] * masses[part]) - masses[part];
        for (int c = 0; c < M; ++c) {
            s[c][0] = tokens[1 + 2 * c];
            serr[c][0] = tokens[1 + 2 * c + 1];
        }

        // Read remaining points
        for (int i = 1; i < N; ++i) {
            std::getline(f, line);
            std::istringstream iss2(line);
            iss2 >> pT[i];
            mT[i] = sqrt(pT[i] * pT[i] + masses[part] * masses[part]) - masses[part];  // Correct mT calculation for all points
            for (int c = 0; c < M; ++c) {
                iss2 >> s[c][i] >> serr[c][i];
            }
        }

        // Create graphs for pT
        for (int c = 0; c < M; ++c) {
            gr[iexp][part][c] = new TGraphErrors(
                N,
                pT.data(), s[c].data(),
                xerr.data(), serr[c].data()
            );
            // Set style
            if (iexp == 0) {
                gr[iexp][part][c]->SetMarkerColor(kRed);
                gr[iexp][part][c]->SetMarkerStyle(28);
                gr[iexp][part][c]->SetLineWidth(0);     
                gr[iexp][part][c]->SetLineColor(kRed);
            } else {
                gr[iexp][part][c]->SetMarkerColor(kBlue);
                gr[iexp][part][c]->SetMarkerStyle(29);
                gr[iexp][part][c]->SetLineWidth(0);
            }
        }
    }
    f.close();
}

void Spectra_PHENIX_STAR() {
    // Read STAR (iexp=1) and PHENIX (iexp=0)
    ReadSpectra("input/PHENIX/AuAu/spectra_PHENIX.txt", 0);
    ReadSpectra("input/STAR/AuAu/spectra_STAR.txt",   1);

    // Create the canvas
    TCanvas* c = new TCanvas("cSpectra", "Spectra PHENIX vs STAR", 1440, 2160);
    c->Divide(2, 3, 0, 0); // 3 columns x 2 rows = 6 pads

    const char* partNames[6] = {"#pi^{+}", "#pi^{-}", "K^{+}", "K^{-}", "p", "#bar{p}"};

    // Fixed limits for all graphs
    const double xmin = 0.0;
    const double xmax = 2; // Extended xmax to include proton/antiproton values
    const double ymin = 1e-3;
    const double ymax = 1e4; // Adjusted to better capture proton/antiproton spectra

    // Define common centralities (indices)
    const std::vector<int> commonCent = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9}; // 0-5%, 5-10%, ..., 70-80%

    // Create the graphs for each particle
    for (int part = 0; part < 6; ++part) {
        c->cd(part+1);
        gPad->SetLogy();

        // Create frame histogram with extended axis limits
        TH1F* frame = gPad->DrawFrame(xmin, ymin, xmax, ymax);
        frame->SetTitle(Form("%s spectra; p_{T} [GeV]; 1/2#pi k_{T} d^{2}N/dp_{T}dy", partNames[part]));
        frame->GetXaxis()->SetTitleSize(0.06);
        frame->GetYaxis()->SetTitleSize(0.06);
        frame->GetXaxis()->SetLabelSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.05);

        TMultiGraph* mg = new TMultiGraph();
        TLegend* leg = new TLegend(0.15, 0.6, 0.55, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.05);
        leg->SetTextFont(42);

        // Iterate through centralities
        int M0 = nCent[0][part]; // PHENIX
        int M1 = nCent[1][part]; // STAR
        for (int cidx : commonCent) {
            if (cidx < M0) {
                mg->Add(gr[0][part][cidx], "P");
            }
            if (cidx < M1) {
                mg->Add(gr[1][part][cidx], "P");
            }
        }

        // Set larger markers for better visibility
        for (int cidx = 0; cidx < M0; ++cidx) {
            gr[0][part][cidx]->SetMarkerSize(1.5);
        }
        for (int cidx = 0; cidx < M1; ++cidx) {
            gr[1][part][cidx]->SetMarkerSize(1.5);
        }

        mg->Draw("P"); // Draw on top of the frame
        leg->Draw();
    }

    c->Update();
    c->SaveAs("output/pics/Spectra_PHENIX_vs_STAR_pT.png");
    c->SaveAs("output/pics/Spectra_PHENIX_vs_STAR_pT.pdf");
}
