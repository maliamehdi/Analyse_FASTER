#ifndef CFISSIONEVENT_H
#define CFISSIONEVENT_H

#include <vector>
#include <map>

class CFissionEvent {
public:
    struct IC_Hit {
        int label;
        double time_ps;
        double energy1;
        double energy2;
    };

    struct GammaHit {
        int label;
        double energy_keV;
        double time_ps;
    };

    std::vector<IC_Hit> ic_hits;
    std::vector<GammaHit> gamma_hits;
    std::vector<double> neutron_times_us;

    void AddIC(int label, double time_ps, double nrj1, double nrj2);
    void AddGamma(int label, double energy_keV, double time_ps);
    void AddNeutron(double time_after_cathode_us);

    double GetCathodeTime_ps() const;
    size_t GetGammaCount() const;
    size_t GetNeutronCount() const;
};

#endif