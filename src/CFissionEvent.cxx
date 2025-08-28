#include "CFissionEvent.h"

void CFissionEvent::AddIC(int label, double time_ps, double nrj1, double nrj2) {
    ic_hits.push_back({label, time_ps, nrj1, nrj2});
}

void CFissionEvent::AddGamma(int label, double energy_keV, double time_ps) {
    gamma_hits.push_back({label, energy_keV, time_ps});
}

void CFissionEvent::AddNeutron(double time_after_cathode_us) {
    neutron_times_us.push_back(time_after_cathode_us);
}

double CFissionEvent::GetCathodeTime_ps() const {
    for (const auto& ic : ic_hits) {
        if (ic.label == 1) return ic.time_ps;
    }
    return -1.0;
}

size_t CFissionEvent::GetGammaCount() const {
    return gamma_hits.size();
}

size_t CFissionEvent::GetNeutronCount() const {
    return neutron_times_us.size();
}