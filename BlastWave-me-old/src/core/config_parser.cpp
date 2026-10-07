#include "core/config_parser.hpp"
#include <iostream>

namespace qcd {
    Config load_config(const std::string& filename) {
        Config cfg;
        // TODO: parse YAML in later phases
        std::cout << "[config] using defaults (YAML parsing not yet implemented)\n";
        return cfg;
    }
}