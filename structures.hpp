// structures
#ifndef STRUCTURES_HPP
#define STRUCTURES_HPP

#include <string>


// system parameters
struct SysPara {
    int N_CH;               // number of chains in system
    int N_AA;               // number of amino acids in chains
    std::string AA_seq;
    double L;               // side length of cubic simulation cell

    int NBin;               // number of energy bins for histogram
    double EMin;
    double EMax;
    double BinW;
    double EStart;          // energy below which the sim ends the preSAMC moves
    int tStart;             // minimum number of preSAMC moves

    int stepit;              //const int nstep = 8*N_AA*N_CH;
    unsigned long int T_0;
    unsigned long int T_MAX;
    unsigned long int T_WRITE;
    unsigned long int T_BC_RESET;

    double GAMMA_0;

    int NBOX;
    double LBOX;            // has to be larger than the biggest interaction radius (SW_HUGE = 7.4)

    int WT_WIGGLE;          // weight wiggle
    int WT_PHI;             // weight rotPhi
    int WT_PSI;             // weight rotPsi
    int WT_TRANS;           // weight translation
    int WT_ROSENC;
    int WT_ROSENN;
    // movement restraints
    double DISP_MAX;
    double DPHI_MAX;
    double DPSI_MAX;

    bool EBIN_TRUNC_UP;     // sorting of integer energy state into upper or lower bin (nessessary to differentiate unambiguously in order to reproduce lngE)
    bool FIX_lngE;
        bool HB_CONTMAT;
        bool WRITE_CONFIG;
            std::vector<double> CONFIG_E;
            double CONFIG_V;
};

struct Header {
    std::string confnm;     // input file name initial configuration
    std::string paranm;     // input file name system parameters
    std::string lngEnm;     // input file name ln g(E)
    std::string dbposi;     // output file name debug position
    std::string hbmatr;     // output file name HB matrices
    std::string lognm;      // output file system log

    std::ofstream os_log;   // ofstream of log file
};


#endif