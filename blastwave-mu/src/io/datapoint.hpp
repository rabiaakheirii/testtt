#pragma once
#include <string>
#include <vector>

namespace bwmu {

    // One measured data point.
    //
    //   mT    : transverse mass minus rest mass, in GeV
    //   value : dN/(mT dmT) in GeV^-2
    //   err   : statistical error, same units as value
    //
    struct DataPoint {
        double mT;
        double value;
        double err;
    };

    // All data for one centrality bin of one system.
    //
    //   centrality : integer index (0 = most central)
    //   species    : name of the particle (must match Species::name)
    //   points     : the measurements
    //
    struct Spectrum {
        int    centrality;
        std::string species;
        std::vector<DataPoint> points;
    };

} // namespace bwmu
