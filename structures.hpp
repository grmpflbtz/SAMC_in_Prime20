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
        bool HB_ContMat;
        bool tGyr;
        bool wConfig;
            std::vector<double> ConfigE;
            double ConfigV;
};

struct Header {
    std::string confnm;     // input file name initial configuration
    std::string paranm;     // input file name system parameters
    std::string lngEnm;     // input file name ln g(E)
    std::string rrunnm;     // input file name rerun file input
    std::string dbposi;     // output file name debug position
    std::string hbmatr;     // output file name HB matrices
    std::string tGyrnm;     // output file name tensor of gyration
    std::string lognm;      // output file system log

    std::ofstream os_log;   // ofstream of log file
};

struct Output {
    long unsigned int *H;           // energy bin visitation histogram
    double *lngE;                   // lng(E)
    double *contHB;                 // hydrogen bond contact matrices
    double **rGyr;                  // sum radius of gyration for each chain and energy bin
    double *rGyrCur;                // current radius of gyration
    double ***tGyrEig;              // sum of principal moments of tensor of gyration for each chain and energy bin
    double **tGyrEigCur;            // current principal moments of tensor of gyration
    int *conf_n;                    // # of configs written
    int *conf_wt;                   // last time writing config for this energy
    long unsigned int nattempt[4];  // # of attempted moves
    long unsigned int naccept[4];   // # of accepted moves
};


#endif