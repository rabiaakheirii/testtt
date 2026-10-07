// analyze_f14_centrality_rapidity_overlay_expLike.C
//
// Запуск:
// root -l -q 'analyze_f14_centrality_rapidity_overlay_expLike.C("myrun.f14")'

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <iomanip>

#include "TH1D.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TColor.h"

struct Particle {
    double t, x, y, z;
    double E, px, py, pz, m;
    int ityp;
    int iso3;
    int charge;
    int parent;
    int ncoll;
    int proc;
};

struct Event {
    int eventId = -1;
    int npart = 0;
    double time = 0.0;
    double impact = -1.0;
    int nch = 0;   // multiplicity only in experimental acceptance
    std::vector<Particle> particles;
};

struct CentralityBin {
    std::string label;
    int nchMin;
    int nchMax;
};

struct StoppingObservables
{
    double meanY_proj = 0.0;          // <y> для projectile hemisphere (y > 0)
    double deltaY_proj = 0.0;         // y_beam - <y>_proj

    double midFraction = 0.0;         // N(|y| < yMidCut) / N(total)
    double peakHeight = 0.0;          // высота projectile-side peak
    double midHeight = 0.0;           // высота в midrapidity
    double peakToMid = 0.0;           // peak / mid

    double protonYieldPerEvent = 0.0; // интеграл по (1/Nevt) dN/dy
    int nProtons = -1;                // если нет ненормированной статистики, оставить -1
};


double Rapidity(double E, double pz) {
    const double eps = 1e-12;
    if (E <= std::abs(pz) + eps) return 999.0;
    return 0.5 * std::log((E + pz) / (E - pz));
}

double Pt(double px, double py) {
    return std::sqrt(px * px + py * py);
}

bool IsProton(const Particle& p) {
    return (p.ityp == 1 && p.iso3 == 1);
}

bool IsNeutron(const Particle& p) {
    return (p.ityp == 1 && p.iso3 == -1);
}

bool IsNucleon(const Particle& p) {
    return IsProton(p) || IsNeutron(p);
}

bool IsCharged(const Particle& p) {
    return (p.charge != 0);
}

bool TryParseEventHeader(const std::string& line, int& eventId, double& impact) {
    std::istringstream iss(line);
    int first, ev, Ap, At;
    double b, sqs, sig, elab, plab;
    if (iss >> first >> ev >> Ap >> At >> b >> sqs >> sig >> elab >> plab) {
        if (first == -1) {
            eventId = ev;
            impact = b;
            return true;
        }
    }
    return false;
}

bool TryParseBodyStart(const std::string& line, int& npart, double& time) {
    std::istringstream iss(line);
    if (iss >> npart >> time) return (npart >= 0);
    return false;
}

bool ParseParticleLine(const std::string& line, Particle& p) {
    std::istringstream iss(line);
    if (!(iss >> p.t >> p.x >> p.y >> p.z
              >> p.E >> p.px >> p.py >> p.pz >> p.m
              >> p.ityp >> p.iso3 >> p.charge
              >> p.parent >> p.ncoll >> p.proc)) {
        return false;
    }
    return true;
}

// ------------------------------------------------------------
// ПАРАМЕТРЫ "ЭКСПЕРИМЕНТА"
// ------------------------------------------------------------

// Сдвиг к c.m. системе: y_cm = y_lab - Y_CM_SHIFT
// Для симметричной системы обычно Y_CM_SHIFT = y_beam/2 в lab.
// Для асимметричной системы лучше задать руками по кинематике твоего случая.
// Пока оставляю как настраиваемый параметр.
const double Y_CM_SHIFT = 0; // 1.05;

// Окно, в котором определяется центральность.
// Это "детекторная" область, где ты реально считаешь множественность.
const double CENT_Y_MIN = -0.5;
const double CENT_Y_MAX =  0.5;

// При желании можно добавить порог по pT, чтобы убрать очень мягкие треки
const double CENT_PT_MIN = 0.0;

// Окно, где строятся распределения dN/dy
const double PLOT_Y_MIN = -2.0;
const double PLOT_Y_MAX =  2.0;

// Если хочешь жёстко отрезать остатки пучка,
// можно ввести дополнительное ограничение на нуклоны:
// например, не брать нуклоны с |y_cm| > 1.8
const bool   APPLY_NUCLEON_RAPIDITY_VETO = true;
const double NUCLEON_Y_VETO_ABS = 1.;

// Для observables барионного стоппинга
const double Y_BEAM_CM = 0.99;     // задай сюда beam rapidity в c.m.
const double MID_Y_CUT = 0.5;      // midrapidity window: |y_cm| < 0.5
const double PEAK_SEARCH_MIN = 0.6; // искать пик только вне midrapidity

// ------------------------------------------------------------

double RapidityCM(const Particle& p) {
    double yLab = Rapidity(p.E, p.pz);
    if (std::abs(yLab) > 900) return 999.0;
    return yLab - Y_CM_SHIFT;
}

// Частица попала в "измеряемую" область для центральности
bool IsInCentralityAcceptance(const Particle& p) {
    if (!IsCharged(p)) return false;

    double ycm = RapidityCM(p);
    if (std::abs(ycm) > 900) return false;

    double pt = Pt(p.px, p.py);

    if (ycm < CENT_Y_MIN || ycm > CENT_Y_MAX) return false;
    if (pt < CENT_PT_MIN) return false;

    return true;
}

// Частица попала в область, где строим dN/dy
bool IsInPlotAcceptance(const Particle& p) {
    double ycm = RapidityCM(p);
    if (std::abs(ycm) > 900) return false;
    if (ycm < PLOT_Y_MIN || ycm > PLOT_Y_MAX) return false;

    // Дополнительный "экспериментальный" veto для нуклонов
    // вблизи beam/target rapidity, если понадобится
    if (APPLY_NUCLEON_RAPIDITY_VETO && IsNucleon(p)) {
        if (std::abs(ycm) > NUCLEON_Y_VETO_ABS) return false;
    }

    return true;
}

std::vector<Event> ReadF14(const char* fname) {
    std::ifstream fin(fname);
    if (!fin.is_open()) {
        std::cerr << "Cannot open file: " << fname << std::endl;
        return {};
    }

    std::vector<Event> events;
    std::string line;

    int currentEventId = -1;
    double currentImpact = -1.0;

    while (std::getline(fin, line)) {
        if (line.empty()) continue;

        int evIdTmp;
        double impactTmp;
        if (TryParseEventHeader(line, evIdTmp, impactTmp)) {
            currentEventId = evIdTmp;
            currentImpact = impactTmp;
            continue;
        }

        int npart;
        double time;
        if (!TryParseBodyStart(line, npart, time)) continue;

        std::string countersLine;
        if (!std::getline(fin, countersLine)) break;

        Event ev;
        ev.eventId = currentEventId;
        ev.impact = currentImpact;
        ev.npart = npart;
        ev.time = time;
        ev.particles.reserve(npart);

        for (int i = 0; i < npart; ++i) {
            std::string pline;
            if (!std::getline(fin, pline)) break;

            Particle p;
            if (!ParseParticleLine(pline, p)) {
                std::cerr << "Warning: cannot parse particle line in event "
                          << currentEventId << ", particle " << i << std::endl;
                continue;
            }

            ev.particles.push_back(p);

            // Центральность считаем только по "измеряемым" заряженным частицам
            if (IsInCentralityAcceptance(p)) ev.nch++;
        }

        events.push_back(std::move(ev));
    }

    return events;
}

std::vector<CentralityBin> BuildCentralityBins(const std::vector<Event>& events) {
    std::vector<double> edgesPct = {0, 10, 20, 40, 60, 80};

    std::vector<int> nchs;
    nchs.reserve(events.size());
    for (const auto& ev : events) nchs.push_back(ev.nch);

    std::sort(nchs.begin(), nchs.end(), std::greater<int>());

    auto percentileNch = [&](double pct) -> int {
        if (nchs.empty()) return 0;
        double pos = pct / 100.0 * (nchs.size() - 1);
        size_t idx = std::round(pos);
        if (idx >= nchs.size()) idx = nchs.size() - 1;
        return nchs[idx];
    };

    std::vector<CentralityBin> bins;
    for (size_t i = 0; i + 1 < edgesPct.size(); ++i) {
        double c1 = edgesPct[i];
        double c2 = edgesPct[i + 1];

        int high = percentileNch(c1);
        int low  = percentileNch(c2);

        CentralityBin bin;
        bin.label = Form("%.0f-%.0f%%", c1, c2);
        bin.nchMin = low;
        bin.nchMax = high;

        if (bin.nchMin > bin.nchMax) std::swap(bin.nchMin, bin.nchMax);
        bins.push_back(bin);
    }

    return bins;
}

int FindCentralityBin(int nch, const std::vector<CentralityBin>& bins) {
    for (size_t i = 0; i < bins.size(); ++i) {
        if (nch >= bins[i].nchMin && nch <= bins[i].nchMax) return (int)i;
    }
    return -1;
}

void NormalizePerEventAndWidth(TH1D* h, int nEvents) {
    if (nEvents > 0) h->Scale(1.0 / nEvents, "width");
}

void SetHistStyle(TH1D* h, int color) {
    h->SetLineColor(color);
    h->SetMarkerColor(color);
    h->SetLineWidth(2);
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
}
StoppingObservables ComputeStoppingObservables(TH1D* h,
                                               double yBeam,
                                               double yMidCut = 0.5,
                                               double peakSearchMin = 0.6)
{
    StoppingObservables obs;
    if (!h) return obs;

    double sumAll = 0.0;      // интеграл по всем y
    double sumProj = 0.0;     // интеграл по projectile hemisphere (y > 0)
    double sumYProj = 0.0;    // интеграл y * dN
    double midYield = 0.0;    // интеграл в |y| < yMidCut

    const int nBins = h->GetNbinsX();

    // --- интегральные observables ---
    for (int ib = 1; ib <= nBins; ++ib) {
        const double y     = h->GetBinCenter(ib);
        const double width = h->GetBinWidth(ib);
        const double dNdy  = h->GetBinContent(ib);

        const double yield = dNdy * width;   // средний выход на событие в данном бине
        if (yield <= 0.0) continue;

        sumAll += yield;

        if (std::abs(y) < yMidCut) {
            midYield += yield;
        }

        if (y > 0.0) {
            sumProj += yield;
            sumYProj += y * yield;
        }
    }

    obs.protonYieldPerEvent = sumAll;

    if (sumProj > 0.0) {
        obs.meanY_proj = sumYProj / sumProj;
        obs.deltaY_proj = yBeam - obs.meanY_proj;
    }

    if (sumAll > 0.0) {
        obs.midFraction = midYield / sumAll;
    }

    // nProtons здесь лучше НЕ вычислять из нормированной гистограммы.
    // Если нужно реальное число протонов, его надо считать отдельно до нормировки.
    obs.nProtons = -1;

    // --- высота в midrapidity ---
    const int midBin = h->FindBin(0.0);
    obs.midHeight = h->GetBinContent(midBin);

    // усреднение по 3 центральным бинам для устойчивости
    if (midBin > 1 && midBin < nBins) {
        obs.midHeight =
            (h->GetBinContent(midBin - 1) +
             h->GetBinContent(midBin) +
             h->GetBinContent(midBin + 1)) / 3.0;
    }

    // --- поиск projectile-side fragmentation peak ---
    double peak = 0.0;
    for (int ib = 1; ib <= nBins; ++ib) {
        const double y   = h->GetBinCenter(ib);
        const double val = h->GetBinContent(ib);

        // ищем пик только в projectile hemisphere и вне central region
        if (y < peakSearchMin) continue;
        if (val > peak) peak = val;
    }
    obs.peakHeight = peak;

    obs.peakToMid = (obs.midHeight > 0.0) ? (obs.peakHeight / obs.midHeight) : 0.0;

    return obs;
}


void SaveStoppingObservables(const std::vector<StoppingObservables>& obsVec,
                             const std::vector<CentralityBin>& centBins,
                             const char* outname = "stopping_observables.dat")
{
    std::ofstream fout(outname);
    if (!fout.is_open()) {
        std::cerr << "Cannot open output file: " << outname << std::endl;
        return;
    }

    fout << "# Centrality  <y>_proj  deltaY_proj  midFraction  peakHeight  midHeight  peakToMid  protonYieldPerEvent  nProtons\n";

    std::cout << "\n=== Baryon stopping observables (projectile side, protons) ===\n";
    std::cout << std::setw(10) << "Cent"
              << std::setw(14) << "<y>_proj"
              << std::setw(16) << "deltaY_proj"
              << std::setw(14) << "midFrac"
              << std::setw(14) << "peak"
              << std::setw(14) << "mid"
              << std::setw(14) << "peak/mid"
              << std::setw(18) << "pYield/evt"
              << std::setw(12) << "nProtons"
              << "\n";

    for (size_t i = 0; i < obsVec.size(); ++i) {
        const auto& o = obsVec[i];

        fout << centBins[i].label << "  "
             << o.meanY_proj << "  "
             << o.deltaY_proj << "  "
             << o.midFraction << "  "
             << o.peakHeight << "  "
             << o.midHeight << "  "
             << o.peakToMid << "  "
             << o.protonYieldPerEvent << "  "
             << o.nProtons << "\n";

        std::cout << std::setw(10) << centBins[i].label
                  << std::setw(14) << std::fixed << std::setprecision(4) << o.meanY_proj
                  << std::setw(16) << o.deltaY_proj
                  << std::setw(14) << o.midFraction
                  << std::setw(14) << o.peakHeight
                  << std::setw(14) << o.midHeight
                  << std::setw(14) << o.peakToMid
                  << std::setw(18) << o.protonYieldPerEvent
                  << std::setw(12) << o.nProtons
                  << "\n";
    }

    fout.close();
}

void rapidity_distributions(const char* fname = "/home/dasha/urqmd-3.4/output_XeW/urqmdXeW.f14") {
    gStyle->SetOptStat(0);

    auto events = ReadF14(fname);
    if (events.empty()) {
        std::cerr << "No events read from " << fname << std::endl;
        return;
    }

    std::cout << "Read events: " << events.size() << std::endl;
    std::cout << "Using y_cm = y_lab - " << Y_CM_SHIFT << std::endl;
    std::cout << "Centrality acceptance: "
              << CENT_Y_MIN << " < y_cm < " << CENT_Y_MAX
              << ", pT > " << CENT_PT_MIN << " GeV/c" << std::endl;
    std::cout << "Plot acceptance: "
              << PLOT_Y_MIN << " < y_cm < " << PLOT_Y_MAX << std::endl;

    int nchMaxAll = 0;
    for (const auto& ev : events) {
        if (ev.nch > nchMaxAll) nchMaxAll = ev.nch;
    }

    TH1D* hNch = new TH1D("hNch",
                         "Measured charged-particle multiplicity;N_{ch}^{acc};Events",
                         nchMaxAll + 1, -0.5, nchMaxAll + 0.5);

    for (const auto& ev : events) hNch->Fill(ev.nch);

    auto centBins = BuildCentralityBins(events);

    std::cout << "\nCentrality classes from measured Nch percentiles:\n";
    for (const auto& b : centBins) {
        std::cout << b.label << " : Nch in [" << b.nchMin << ", " << b.nchMax << "]\n";
    }

    const int nYBins = 50;
    const double yMin = PLOT_Y_MIN;
    const double yMax = PLOT_Y_MAX;

    std::vector<TH1D*> hRapP;
    std::vector<TH1D*> hRapCh;
    std::vector<int> nEventsInBin(centBins.size(), 0);

    std::vector<int> colors = {
        kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7,
        kCyan + 1, kViolet + 1, kTeal + 2, kPink + 7
    };

    for (size_t i = 0; i < centBins.size(); ++i) {
        hRapP.push_back(new TH1D(Form("hRapP_%zu", i),
                                 ";y_{cm};(1/N_{evt}) dN/dy_{cm}",
                                 nYBins, yMin, yMax));

        hRapCh.push_back(new TH1D(Form("hRapCh_%zu", i),
                                  ";y_{cm};(1/N_{evt}) dN/dy_{cm}",
                                  nYBins, yMin, yMax));

        SetHistStyle(hRapP[i], colors[i % colors.size()]);
        SetHistStyle(hRapCh[i], colors[i % colors.size()]);
    }

    for (const auto& ev : events) {
        int ic = FindCentralityBin(ev.nch, centBins);
        if (ic < 0) continue;

        nEventsInBin[ic]++;

        for (const auto& p : ev.particles) {
            if (!IsInPlotAcceptance(p)) continue;

            double ycm = RapidityCM(p);
            if (std::abs(ycm) > 900) continue;

            if (IsCharged(p)) hRapCh[ic]->Fill(ycm);
            if (IsProton(p))  hRapP[ic]->Fill(ycm);
        }
    }

    for (size_t i = 0; i < centBins.size(); ++i) {
        NormalizePerEventAndWidth(hRapP[i],  nEventsInBin[i]);
        NormalizePerEventAndWidth(hRapCh[i], nEventsInBin[i]);
    }

      // ===== observables барионного стоппинга =====
    std::vector<StoppingObservables> stoppingObs;
    for (size_t i = 0; i < centBins.size(); ++i) {
        stoppingObs.push_back(
            ComputeStoppingObservables(hRapP[i], Y_BEAM_CM, MID_Y_CUT, PEAK_SEARCH_MIN)
        );
    }

    SaveStoppingObservables(stoppingObs, centBins, "stopping_observables.dat");
    
    TFile* fout = new TFile("of14_centrality_rapidity_overlay_expLike.root", "RECREATE");
    hNch->Write();
    for (auto* h : hRapP)  h->Write();
    for (auto* h : hRapCh) h->Write();
    fout->Close();

    TCanvas* cNch = new TCanvas("cNch", "Nch", 850, 650);
    hNch->SetLineWidth(2);
    hNch->Draw("hist");
    cNch->SaveAs("output/Nch_distribution_expLike.png");

    // ===== Протоны =====
    TCanvas* cP = new TCanvas("cP", "Proton rapidity", 900, 700);

    double maxP = 0.0;
    for (auto* h : hRapP) {
        if (h->GetMaximum() > maxP) maxP = h->GetMaximum();
    }

    hRapP[0]->SetTitle("Protons in experimental acceptance");
    hRapP[0]->SetMaximum(1.15 * maxP);
    hRapP[0]->Draw("hist");
    for (size_t i = 1; i < hRapP.size(); ++i) {
        hRapP[i]->Draw("hist same");
    }

    TLegend* legP = new TLegend(0.70, 0.55, 0.88, 0.88);
    legP->SetBorderSize(0);
    legP->SetFillStyle(0);
    for (size_t i = 0; i < hRapP.size(); ++i) {
        legP->AddEntry(hRapP[i], centBins[i].label.c_str(), "l");
    }
    legP->Draw();

    cP->SaveAs("output/rapidity_protons_overlay_expLike.png");

    // ===== Все заряженные =====
    TCanvas* cCh = new TCanvas("cCh", "Charged rapidity", 900, 700);

    double maxCh = 0.0;
    for (auto* h : hRapCh) {
        if (h->GetMaximum() > maxCh) maxCh = h->GetMaximum();
    }

    hRapCh[0]->SetTitle("Charged particles in experimental acceptance");
    hRapCh[0]->SetMaximum(1.15 * maxCh);
    hRapCh[0]->Draw("hist");
    for (size_t i = 1; i < hRapCh.size(); ++i) {
        hRapCh[i]->Draw("hist same");
    }

    TLegend* legCh = new TLegend(0.70, 0.55, 0.88, 0.88);
    legCh->SetBorderSize(0);
    legCh->SetFillStyle(0);
    for (size_t i = 0; i < hRapCh.size(); ++i) {
        legCh->AddEntry(hRapCh[i], centBins[i].label.c_str(), "l");
    }
    legCh->Draw();

    cCh->SaveAs("output/rapidity_charged_overlay_expLike.png");

    std::cout << "\nEvents per centrality bin:\n";
    for (size_t i = 0; i < centBins.size(); ++i) {
        std::cout << centBins[i].label << " : " << nEventsInBin[i] << " events\n";
    }

    std::cout << "\nSaved:\n"
              << "  f14_centrality_rapidity_overlay_expLike.root\n"
              << "  Nch_distribution_expLike.png\n"
              << "  rapidity_protons_overlay_expLike.png\n"
              << "  rapidity_charged_overlay_expLike.png\n";
}