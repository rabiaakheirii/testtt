#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "io/datapoint.hpp"
#include "io/phenix_reader.hpp"
#include "io/csv_writer.hpp"
#include "particles/table.hpp"

static int failures = 0;

#define CHECK(cond, msg)                                       \
    do {                                                       \
        if (!(cond)) {                                         \
            std::cout << "  FAIL: " << msg << "\n";            \
            ++failures;                                        \
        } else {                                               \
            std::cout << "  ok  : " << msg << "\n";            \
        }                                                      \
    } while (0)

// Write a small synthetic PHENIX file with 6 blocks, 3 centrality bins
static std::string make_test_file() {
    const std::string path = "/tmp/test_phenix_auau.txt";
    std::ofstream f(path);
    // 3 blocks of 2 points each, to keep the test quick.
    // Fill missing blocks with 2-column data.
    for (int blk = 0; blk < 6; ++blk) {
        f << "2\n";
        // 2 data points, 3 centrality columns
        f << "0.5  1.0  0.1  0.8  0.1  0.5  0.05\n";
        f << "1.0  0.5  0.05 0.4  0.05 0.2  0.02\n";
        f << "\n";
    }
    f.close();
    return path;
}

int main() {
    std::cout << "========== Phase 4 unit tests ==========\n";

    const std::string path = make_test_file();
    std::cout << "\n[Test file]\n";
    std::cout << "  " << path << "\n";

    std::cout << "\n[Read all species at centrality 0]\n";
    auto spectra = bwmu::read_phenix(path, 0);
    CHECK(spectra.size() == 6, "6 species read");

    for (const auto& s : spectra) {
        std::cout << "  " << s.species
                  << "  centrality=" << s.centrality
                  << "  n_points=" << s.points.size() << "\n";
    }

    std::cout << "\n[Check pi+ points]\n";
    const auto& pi = spectra[0];
    CHECK(pi.species == "pi+", "first species is pi+");
    CHECK(pi.points.size() == 2, "pi+ has 2 points");

    // pi+ mass ~ 0.13957
    // pT = 0.5 -> mT = sqrt(0.25 + 0.01948) = 0.51912
    // mT - m0 = 0.37955
    const double expected_mT0 = std::sqrt(0.5*0.5 + 0.13957*0.13957) - 0.13957;
    CHECK(std::abs(pi.points[0].mT - expected_mT0) < 1e-6,
          "pi+ point 0: mT-m0 correct");
    CHECK(std::abs(pi.points[0].value - 1.0) < 1e-12,
          "pi+ point 0: value correct");
    CHECK(std::abs(pi.points[0].err   - 0.1) < 1e-12,
          "pi+ point 0: err correct");

    const double expected_mT1 = std::sqrt(1.0*1.0 + 0.13957*0.13957) - 0.13957;
    CHECK(std::abs(pi.points[1].mT - expected_mT1) < 1e-6,
          "pi+ point 1: mT-m0 correct");

    std::cout << "\n[Read different centrality]\n";
    auto spectra_c1 = bwmu::read_phenix(path, 1);
    CHECK(spectra_c1.size() == 6, "6 species at centrality 1");
    // value should be 0.8 for pi+ point 0
    CHECK(std::abs(spectra_c1[0].points[0].value - 0.8) < 1e-12,
          "pi+ centrality 1 value = 0.8");

    std::cout << "\n[Read only pi-]\n";
    auto pi_minus = bwmu::read_phenix_species(path, 2, "pi-");
    CHECK(pi_minus.species == "pi-", "species is pi-");
    CHECK(pi_minus.points.size() == 2, "pi- has 2 points");
    CHECK(std::abs(pi_minus.points[0].value - 0.5) < 1e-12,
          "pi- centrality 2 value = 0.5");

    std::cout << "\n[K- mass check]\n";
    auto k_minus = bwmu::read_phenix_species(path, 0, "K-");
    // Use the exact mass from the particle table, not a rounded constant.
    const double mK = bwmu::get_species("K-").mass;
    const double expected_mT_K =
        std::sqrt(0.5*0.5 + mK*mK) - mK;
    CHECK(std::abs(k_minus.points[0].mT - expected_mT_K) < 1e-10,
          "K- point 0: mT-m0 uses exact K mass");

    std::cout << "\n[CSV writer]\n";
    const std::string csv = "/tmp/test_output.csv";
    bwmu::write_spectrum_csv(csv, pi);
    std::ifstream f(csv);
    CHECK(f.good(), "CSV file was created");
    f.close();

    std::cout << "\n[Missing file handling]\n";
    bool threw = false;
    try {
        bwmu::read_phenix("/tmp/does-not-exist-12345.txt", 0);
    } catch (const std::runtime_error&) {
        threw = true;
    }
    CHECK(threw, "missing file throws");

    std::cout << "\n[Out-of-range centrality]\n";
    auto oob = bwmu::read_phenix(path, 100);  // no such column
    bool all_empty = true;
    for (const auto& s : oob) if (!s.points.empty()) all_empty = false;
    CHECK(all_empty, "centrality=100 gives empty spectra");

    std::cout << "\n========================================\n";
    if (failures == 0) {
        std::cout << " All Phase 4 tests passed.\n";
        return 0;
    } else {
        std::cout << " " << failures << " test(s) failed.\n";
        return 1;
    }
}
