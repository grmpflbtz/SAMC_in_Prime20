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
    int neighUpdate;        // update step for neighbor list
    int t_NLUpdate;         // steps passed since last neighbor list update
    int NeighListTest;      // tests integrity of LinkList and NeighHead every step [on/off]
    int BondLengthTest;     // check bond length in whole chain at T_WRITE [on/off]

    int WT_WIGGLE;          // weight wiggle
    int WT_PIVOT;           // weight pivot moves
    int WT_TRANS;           // weight translation
    int WT_ROT;             // weight rotation 
    int WT_ROSENC;
    int WT_ROSENN;
    // movement restraints
    double DISP_MAX;        // maximum local displacement
    double DPIV_MAX;        // maximum angle pivot rotation
    double DTRN_MAX;        // maximum chain translation
    double DROT_MAX;        // maximum angle chain rotation

    int cluster_opt;        // cluster optimization [on/off]

    bool FIX_lngE;          // fixed DOS run - no SAMC, only production of observables
    bool HB_ContMat;        // Hydrogen bond contact matrices
    bool Ree;               // end-to-end distance distribution
    bool tGyr;              // tensor of gyration
    bool intrainterMol;     // intra- vs inter-molecular energies
    bool Et;                // energy time development
    bool dihedral;          // dihedral angles
    bool wConfig;           // write configurations for energies specified by vector ConfigE
        double conf_EMax;
        double conf_EMin;
        double conf_dE;
        int conf_dt;
        int conf_Nmax;

    long int Seed; // seed for random number generator
    int add_Seed;           // integer added to Seed (manipulate seed with program arguments)

    time_t starttime;
};

struct Header {
    std::string confnm;         // input file name initial configuration
    std::string paranm;         // input file name system parameters
    std::string lngEnm;         // input file name ln g(E)
    std::string rrunnm;         // input file name rerun file input
    std::string dbposi;         // output file name debug position
    std::string iniconf;        // output file name config created by newChain()
    std::string hbmatr;         // output file name HB matrices
    std::string reenm;          // output file name end-to-end distance distribution
    std::string tGyrnm;         // output file name tensor of gyration
    std::string intrainterMol;  // output file name inter-vs-intra-molecular energies
    std::string grdcnm;         // output file name ground state configuration
    std::string enertm;         // output file name energy time development
    std::string snapshots;      // output file name snapshots
    std::string dihedPhinm;     // output file name dihedral angle Phi
    std::string dihedPsinm;     // output file name dihedral angle Psi
    std::string lognm;          // output file system log

    std::ofstream os_log;   // ofstream of log file
};

struct Output {
    long unsigned int *H;           // energy bin visitation histogram
    double *lngE;                   // lng(E)
    double *contHB;                 // hydrogen bond contact matrices
    double *Ree2;                   // squared end-to-end distance distribution
    double **rGyr;                  // sum radius of gyration for each chain and energy bin
    double *rGyrCur;                // current radius of gyration
    double ***tGyrEig;              // sum of principal moments of tensor of gyration for each chain and energy bin
    double **tGyrEigCur;            // current principal moments of tensor of gyration
        int tGyr_freq;              // frequency of gyration radius observable calculation
    double vdWener[2];              // intra- and inter-molecular van-der-Waals energy
    double HBener[2];               // intra- and inter-molecular HB energy
    double *intrainterE;            // inter- and intra-molecular energies
        int intrainterE_freq;       // frequency of interintraE calculations
    double *Et;                     // energy time development
    long unsigned int ***dihePhi;
    long unsigned int ***dihePsi;
    int conf_Ntot;                  // # of total configs recorded
    double conf_Eprev;              // previously recorded energy
    unsigned long conf_tprev;        // previously step recorded at
    long unsigned int nattempt[5];  // # of attempted moves
    long unsigned int naccept[5];   // # of accepted moves
};


#endif