/* FILE: prime20_samc_main.cpp
 * 
 *                          X         X     XXXX     XXXX  XX     XX
 *                          XX       XX    XX  XX     XX   XXX    XX
 *                          XXX     XXX   XX    XX    XX   XXXX   XX
 *                          XXXX   XXXX  XXX    XXX   XX   XX XX  XX
 *                          XX XX XX XX  XXXXXXXXXX   XX   XX  XX XX
 *                          XX  XXX  XX  XXX    XXX   XX   XX   XXXX
 *                          XX   X   XX  XX      XX   XX   XX    XXX
 *                          XX       XX  XX      XX  XXXX  XX     XX
 * 
 *              :'''''''''''''''''''''''''''''''''''''''''''''''''''''''''':
 *              : SAMC simulation for A-beta peptides in the Prime20 model :
 *              :,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,:
 */
/* changelog:
 * 2019-02-08   Start. initialize bead & Amino acid class.
 * 2019-02-22   newChain function set up. trying to construct first amino acid - method via angles or distances?
 * 2019-03-01   first amino acid constructed - geometric solution - positions correct
 * 2019-03-08   whole chain created (flat). angles and distances correct (test program "config_check.cpp")
 * 2019-03-22   EnOv_bead() finished and integrated into wiggle() - extensive tests missing
 * 2019-03-28   abandon hierarchical energy calculation -> change to neighbour lists, which should be more efficient for multiple chains (legacy version 20190328)
 * 2019-04-03   new neighbourlist concept -> linked lists (legacy version 20190403)
 * 2019-05-09   HBcheck created and implemented in wiggle. still needs checking in running simulation.
 * 2019-05-13   wiggle implemented in newChain. including HBList and HBDist updates.
 * 2019-05-16   periodic boundary condition and shortest distance calculation implemented in wiggle() and distVecBC()
 *              split parts in newChain(): only creation now. Movement to randomize starting configuration is now in main() (legacy version in 20190516)
 * 2019-05-29   EO_SegBead() is new function, for energy calculation and overlapp check. One can specify the segment the bead is checked against in the parameters (legacy version in 20190523)
 * 2019-06-06   start to implement tracker of applied periodic boundary conditions in order to restore real coordinates for rotation moves (legacy version in 20190606)
 * 2019-06-19   added rotPsi() -> all move functions are constructed 
 * 2019-06-26   tweaked EO_SegSeg() to count 2x energy, when moving along one segment to account for energy contribution in other segment
 * 2019-07-30   in wiggle(): start switch(j) cascade for re-considering brokenHB for new HB
 * 2019-08-07   noticed that in HB_check() the periodic boundary conditions were not considered → use distVecBC() instead of manual distance calculations (legacy version in 20190807)
 *              implementation of HB at chain ends → rotation of NCa in NCaC plane by angle 2.668
 * 2019-08-15   introduced "class Chain" to enable multiple chains in system (legacy version 20190815)
 * 2019-10-15   two chains working in system. All HB problems seem to be solved (further, more extensive test necessary).
 *              Write-function for backup files implemented.
 * 2019-10-25   read functions for coordinate input and lngE input implemented
 * 2019-11-11   added extra_lng() function: enables independent lngE[] preset from input file
 * 2019-12-17   backup in archive | change move selection to real_dist with trunc() | 
 * 2020-01-10   corrected EBin selection to not include energies < EMin
 * 2020-02-03   readCoord() function read BC coordinates → now: reads real coords (R script transforms backup_coord_output)
 *              tranlsation move implementation startet (legacy version 20200203) - done
 * 2020-02-05   transition to SWBB = 4.5 and 
 *              squeeze parameters. New chain creation function necessary (legacy version 20200205)
 * 2020-02-14   introduce production (approximate lng) and measure run (fixed lng - measure geometric observables). (legacy version 20200214)
 * 2020-08-04   changed nstep to number of beads in system
 * 
 * MOVED TO GITHUB DOCUMENTATION
 * 
 */
/* abbreviations:
 * SC   - side chain
 * BB   - backbone
 * SW   - square well
 * SQ   - squared
 * HC   - Hard core
 * HS   - hard sphere
 * SQZ  - squeeze factor
 */

#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <math.h>
#include <tuple>
#include <random>
#include <thread>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <iterator>

#include "structures.hpp"
#include "prime20_samc_chain.hpp"
#include "prime20_samc_para.hpp"
#include "timer.hpp"

using namespace std;

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXX  >>   SIMULATION PARAMETERS - GLOBAL VARIABLES   <<  XXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

int *neighHead;
int *neighList;

int **HBList;
int **HBLcpy;

double **NCDist;
double **NCDcpy;

mt19937 rng(time(NULL));                // constructor for random number generator
//mt19937 rng(44);                        // debug

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  >>   FUNCTION DECLARATIONS   <<  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double dotPro(double vecA[], double vecB[]);                                            // dot product
double absVec(double vec[]);                                                            // absolute of vector
tuple<double,double,double> crossPro(double vecA[], double vecB[]);                     // cross product
tuple<double,double,double> distVecBC(SysPara *sp, Bead vecA, Bead vecB);               // distance vector

int program_start_print(ostream &os);                                                   // prints simulation title to os
int command_print(Header *hd, ostream &os);                                             // prints input/output file names
int system_parameter_print(SysPara *sp, ostream &os);                                   // prints system parameters to os
int sim_parameter_print(SysPara *sp, ostream &os);                                      // prints simulation parameters to os
int cur_time_print(ostream &os);                                                        // prints current time

int CommandInitialize(int argc, char *argv[], Header *hd);                              // identify file names from command input

bool newChain(SysPara *sp, Chain Chn[], int chnNum);                                    // creates new chain
bool readParaInput(SysPara *sp, Header *hd);                                            // read system parameters from file
bool readCoord(SysPara *sp, Header *hd, Chain Chn[]);                                   // read chain config from file
bool readPrevRunInput(SysPara *sp, Header *hd, Output *ot, Chain Chn[], long unsigned int &tcont, double &gammasum);      // reads lngE, H, gammasum, and t from input file
bool read_lngE(SysPara *sp, Header *hd, Output *ot);                                    // reads lngE data from file

bool outputPositions(SysPara *sp, Header *hd, std::string fnm, Chain Chn[], int mode, double ener); // writes positions to file "fnm"
bool output_HBmat(SysPara *sp, Header *hd, Output *ot, int step);                       // write hydrogen bond matrix
bool output_Ree(SysPara *sp, Header *hd, Output *ot, int step);                         // write end-to-end distance
bool output_tGyr(SysPara *sp, Header *hd, Output *ot, int step);                        // write tensor of gyration
bool output_vdW(SysPara *sp, Header *hd, Output *ot, int step);                         // write van-der-Waals energy
bool output_Et(SysPara *sp, Header *hd, Output *ot, int step, int init);                // write energy time development
bool BackupSAMCrun(SysPara *sp, Header *hd, Output *ot, Chain Chn[], Timer &Timer, unsigned long int t, double gammasum, double gamma, double E);    // backup function in SAMC run
bool BackupProdRun(SysPara *sp, Header *hd, Output *ot, Timer &Timer, unsigned long int t);         // backup of observables for production run

bool HBcheck(SysPara *sp, Chain Chn[], int iN, int iC);                                 // check if HB exists and update HBList
double E_single(Chain Chn[], int h1, int i1, int h2, int i2, double d_sq);              // energy of single SC interaction
double EO_SegBead(SysPara *sypa, Chain Chn[], int h1, int i1, int j1, int sp, int ep, int EOswitch); // SC interaction energy of one SC Bead (j1 must be 3). Or overlapp check for any AmAc[i1].Bd[j1]. Versus chain segment [sp, ep).
double EO_SegSeg(SysPara *sp, Chain Chn[], int sp1, int ep1, int sp2, int ep2, int EOswitch);        // SC interaction energy of segment [sp1,ep1) versus segment [sp2,ep2). also overlapp check
double E_check(SysPara *sp, Chain Chn[]);                                               // recalculate energy from scratch
bool acceptance(double lngEold, double lngEnew);                                        // SAMC acceptance function

bool wiggle(SysPara *sp, Chain Chn[], int h, int i, int j, double &deltaE);             // small displacement of Chn[h].AmAc[i].Bd[j]
bool rotPhi(SysPara *sypa, Chain Chn[], int i1, int high, double &deltaE);              // rotation around N_AA-Ca axis of i1-th amino acid
bool rotPsi(SysPara *sypa, Chain Chn[], int i1, int high, double &deltaE);              // rotation around Ca-C axis of i1-th amino acid
bool translation(SysPara *sp, Chain Chn[], int iChn, double &deltaE);                   // translation move of the whole chain
bool rotation(SysPara *sp, Chain Chn[], int iChn, double &deltaE);                      // rotation move of the whole Chain

bool checkBndLngth(SysPara *sypa, Chain Chn[], int sp, int ep);                         // check all bond length from Chn[sp/N_AA].AmAc[sp%N_AA] to Chn[ep/N_AA].AmAc[ep%N_AA]
bool resetBCcouter(SysPara *sp, Chain Chn[]);                                           // resets the counter of boundary crossings so that the real coordinates move back to the simulation box

int calc_gyration_radius(SysPara *sp, Output *ot, Chain Chn[], int eBin);               // calculate radius of gyration
int calc_gyration_tensor(SysPara *sp, Output *ot, Chain Chn[], int eBin);               // calculate and sum up tensor of gyration
int calc_vanderWaals(SysPara *sp, Output *ot, Chain Chn[], int eBin, double &vdW_intra, double &vdW_inter, int recalc);                   // calculate van-der-Waals energy

int assignBox(SysPara *sp, Bead Bd);                                                    // assigns neighbour list box to bead
int LinkListInsert(SysPara *sp, Chain Chn[], int i1, int j1);                           // insert particle AmAc[i1].Bd[j1] into Linked List
int LinkListUpdate(SysPara *sp, Chain Chn[], int i1, int j1);                           // update Linked List position of AmAc[i1].Bd[j1]

int output_memory_deallocation(SysPara *sp, Output *ot);                               // deallocating momory of Output

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  >>   MAIN FUNCTION   <<  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

int main(int argc, char *argv[]) 
{
    SysPara *sp = new SysPara;                              // system parameters
    Header *hd = new Header;                                // file names for input/output
    Output *ot = new Output;                                // observables for output

    Timer Timer;                                            // timer for simulation
    Chain *Chn;                                             // chains in system
    Bead *BdCpy;                                            // copy of Beads
    ofstream backup;
    std::ostringstream oss;
    string filename;
    double Eold, Enew, Egrd, deltaE;                        // energy values
    double vdW_intra, vdW_inter;
    double dist[3];                                         // distance vector
    double distabs;                                         // absolute distance
    bool accept;                                            // acceptance value of move in SAMC
    int i_rand, ip, jp, moveselec, movetype;                // variables for performing moves
    int oldBox;                                             // variables for neighbour list calculations
    long unsigned int step, tcont;                          // time variable keeping track of # steps passed
    int it;
    int eBin_n, eBin_o;                                     // energy bin of new and old energy value
    double gamma, gammasum;                                 // gamma value and sum over gamma(t)
    double m_total;                                         // total mass of the system

    if(CommandInitialize(argc, argv, hd) == -1) {
        delete sp; delete hd; delete ot;
        return 0;
    }
    program_start_print(std::cout);
    command_print(hd, std::cout);
    program_start_print(hd->os_log);
    command_print(hd, hd->os_log);

    if( !readParaInput(sp, hd) ) {
        delete sp; delete hd; delete ot;
        return 0;
    }

    system_parameter_print(sp, std::cout);
    system_parameter_print(sp, hd->os_log);
    sim_parameter_print(sp, std::cout);
    sim_parameter_print(sp, hd->os_log);


    neighHead = new int[sp->NBOX*sp->NBOX*sp->NBOX];
    neighList = new int[4*sp->N_AA*sp->N_CH];
    HBList = new int*[sp->N_CH*sp->N_AA];
    HBLcpy = new int*[sp->N_CH*sp->N_AA];
    NCDist = new double*[sp->N_CH*sp->N_AA];
    NCDcpy = new double*[sp->N_CH*sp->N_AA];
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) { 
        HBList[i] = new int[2];
        HBLcpy[i] = new int[2];
        NCDist[i] = new double[sp->N_CH*sp->N_AA];
        NCDcpy[i] = new double[sp->N_CH*sp->N_AA];
    }

    Chn = new Chain[sp->N_CH];
    BdCpy = new Bead[4*sp->N_AA*sp->N_CH];

    ot->H = new long unsigned int[sp->NBin];
    ot->lngE = new double[sp->NBin];
    ot->contHB = new double[sp->NBin * sp->N_CH*sp->N_AA * sp->N_CH*sp->N_AA];
    ot->Ree2 = new double[sp->N_CH * sp->NBin];
    ot->rGyr = new double*[sp->N_CH];
    for( int i=0; i<sp->N_CH; i++ ) { ot->rGyr[i] = new double[sp->NBin]; }
    ot->rGyrCur = new double[sp->N_CH];
    ot->tGyrEig = new double**[sp->N_CH];
    for( int i=0; i<sp->N_CH; i++) {
        ot->tGyrEig[i] = new double*[sp->NBin]; 
        for( int j=0; j<sp->NBin; j++ ) {
            ot->tGyrEig[i][j] = new double[3]{};
        }
    }
    ot->tGyrEigCur = new double*[sp->N_CH];
    for( int i=0; i<sp->N_CH; i++ ) { ot->tGyrEigCur[i] = new double[3]; }
    ot->vdWener = new double[sp->NBin * 2];
    ot->Et = new double[sp->T_WRITE];
    ot->conf_n = new int[sp->ConfigE.size()];
    ot->conf_wt = new int[sp->ConfigE.size()];

    fill_n(neighHead, sp->NBOX*sp->NBOX*sp->NBOX, -1), fill_n(neighList, 4*sp->N_AA*sp->N_CH, -1);

    cur_time_print(std::cout);
    cur_time_print(hd->os_log);

    // initialization of some values
    for( int i=0; i<5; i++ ) {
        ot->nattempt[i] = 0;
        ot->naccept[i] = 0;
    }
    Egrd = 0.0;
    // lngE, H, gammasum, t: read from input file or start new
    if( readPrevRunInput(sp, hd, ot, Chn, tcont, gammasum) ) {
        gamma = sp->GAMMA_0*sp->T_0/tcont;
    }
    else {
        // initialize lngE and H and co.
        for( int i=0; i<sp->NBin; i++ ) {
            ot->lngE[i] = 0;
            ot->H[i] = 0;
        }
        gamma = sp->GAMMA_0;
        gammasum = 0.0;
        tcont = 0;
    }
    // check if there is an extra input file for lngE
    if( read_lngE(sp, hd, ot ) ) {
        tcont = 0;
        if(sp->FIX_lngE) {
            std::cout << "production run with fixed ln g(E)" << std::endl;
        }
    }
    if(sp->HB_ContMat) {
        for( int i=0; i<sp->NBin; i++ ) { 
            ot->H[i] = 0; 
            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    ot->contHB[i*sp->N_CH*sp->N_AA*sp->N_CH*sp->N_AA + j*sp->N_CH*sp->N_AA + k] = 0;
                }
            }
        }
    }
    if(sp->Ree) {
        for( int i=0; i<sp->N_CH*sp->NBin; i++ ) {
            ot->Ree2[i] = 0;
        }
    }
    if(sp->vdWener) {
        for( int i=0; i<sp->NBin * 2; i++ ) {
            ot->vdWener[i] = 0;
        }
    }
    if(sp->Et) {
        output_Et(sp, hd, ot, 0, 0);
    }
    if(sp->wConfig) {
        for( int i=0; i<sp->ConfigE.size(); i++ ) {
            ot->conf_n[i] = 0;
            ot->conf_wt[i] = 0;
        }
    }


    // geometry: read from input file or create new
    if( !readCoord(sp, hd, Chn) ) {
        std::cout << "building new chain(s) ... " << std::flush;
        for( int i=0; i<sp->N_CH; i++ ) {
            newChain(sp, Chn, i);
            for( int m=0; m<i+1; m++) {
                for( int j=0; j<sp->N_AA; j++ ) {
                    for( int k=0; k<4; k++ ) {
                        dist[2] = Chn[m].AmAc[j].Bd[k].getR(2) + (sp->L)/((double)sp->N_CH);
                        Chn[m].AmAc[j].Bd[k].setBC( 2, Chn[m].AmAc[j].Bd[k].getBC(2) + floor(dist[2]/sp->L) );
                        dist[2] = dist[2] - sp->L*floor(dist[2]/sp->L);
                        Chn[m].AmAc[j].Bd[k].setR( Chn[m].AmAc[j].Bd[k].getR(0), Chn[m].AmAc[j].Bd[k].getR(1), dist[2] );
                        LinkListUpdate(sp, Chn, (m*sp->N_AA)+j, k);
                    }
                }
            }
        }
        if(EO_SegSeg(sp, Chn, 0, sp->N_CH*sp->N_AA, 0, sp->N_CH*sp->N_AA, 0) == -1 ) { std::cout << "overlapp in chain !!!" << std::endl; }
        else { std::cout << "complete" << std::endl; }
    }


    /*std::cout << "chain masses: " << std::endl << "Chn[0].getM()=" << Chn[0].getM() << std::endl << "Chn[1].getM()=" << Chn[1].getM() << std::endl;*/


    // initialize HB list and N-C distance list
    for(int i = 0; i < sp->N_CH*sp->N_AA; i++) {
        HBList[i][0] = -1;  HBLcpy[i][0] = -1;
        HBList[i][1] = -1;  HBLcpy[i][1] = -1;
        for(int j = 0; j < sp->N_CH*sp->N_AA; j++) {
            if( (i/sp->N_AA != j/sp->N_AA) ||  ((i/sp->N_AA == j/sp->N_AA) && (abs(i-j) > 3)) ) {
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0], Chn[j/sp->N_AA].AmAc[j%sp->N_AA].Bd[2]);
                distabs = dotPro(dist, dist);
                if(distabs < SW2HUGE) {                // exclude unrelevant distances
                    NCDist[i][j] = distabs;
                    NCDcpy[i][j] = NCDist[i][j];
                    continue;
                }
            }
            NCDist[i][j] = -1;
            NCDcpy[i][j] = -1;
        }
    }
    // first energy calculation - no HB, yet
    Eold = 0.0;
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
            if( HBcheck(sp, Chn, i, j) ) {
                Eold -= 1.0;
            }
        }
        Eold += EO_SegBead(sp, Chn, i/sp->N_AA, i%sp->N_AA, 3, i, sp->N_CH*sp->N_AA, 1);
        // overlapp check
        for( int j=0; j<4; j++ ) {
            if( EO_SegBead(sp, Chn, i/sp->N_AA, i%sp->N_AA, j, i, sp->N_CH*sp->N_AA, 0) < 0.0 ) {
                std::cout << "--- ERROR ---\toverlapp in configuration detected" << std::endl;
            }
        }
    }

    // position check file output
    outputPositions(sp,hd, hd->dbposi, Chn, 0, Eold);

    // move new chains to randomize initial configuration. no energy-dependent acception criterion. all legal moves are accepted
    step = 0;
    while( true ) {
        // end pre-SAMC movement after tStart moves and within desired energy window
        if( (step >= sp->tStart) && (Eold >= sp->EMin) && (Eold < sp->EStart) ) {
            eBin_o = floor((Eold - sp->EMin)/sp->BinW);
            if( eBin_o == sp->NBin) { eBin_o = sp->NBin-1; }
            std::cout << "pre-SAMC move " << step << "\r" << std::flush;
            break;
        }

        if( step%1000 == 0 ) {
            std::cout << "pre-SAMC move " << step << "\r" << std::flush;
        }
        // HBList copy
        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
            HBLcpy[k][0] = HBList[k][0];
            HBLcpy[k][1] = HBList[k][1];
        }
        deltaE = 0.0;
        // select move type
        moveselec = trunc( ( (double)rng()/( (double)rng.max()+1 ) )*(sp->WT_WIGGLE + sp->WT_PHI + sp->WT_PSI + sp->WT_TRANS + sp->WT_ROT));
        if( moveselec < sp->WT_WIGGLE )                                                    { movetype = 0; }    // wiggle
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI )                                    { movetype = 1; }    // rotPhi
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI )                         { movetype = 2; }    // rotPsi
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI+sp->WT_TRANS )            { movetype = 3; }    // translation
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI+sp->WT_TRANS+sp->WT_ROT ) { movetype = 4; }    // rotation
        else { 
            hd->os_log<< "--- ERROR ---\tno movetype was selected" << endl; hd->os_log.close();
            std::cout << "--- ERROR ---\tno movetype was selected" << endl; return 0; 
        }
        switch( movetype ) {
            case 0:
                i_rand = trunc(((double)rng()/((double)rng.max()+1))*(sp->N_CH*sp->N_AA*4));
                ip = i_rand/4;                              // amino acid identifier
                jp = i_rand%4;                              // bead
                accept = wiggle(sp, Chn, ip/sp->N_AA, ip%sp->N_AA, jp, deltaE);
                break;
            case 1:
                ip = trunc(((double)rng()/((double)rng.max()+1))*(sp->N_CH*sp->N_AA));  // amino acid identifier of the rotation origin
                jp = trunc(((double)rng()/((double)rng.max()+1))*2);                                      // rotate lower or higher part
                accept = rotPhi(sp, Chn, ip, jp, deltaE);
                break;
            case 2:
                ip = trunc(((double)rng()/((double)rng.max()+1))*(sp->N_CH*sp->N_AA));  // amino acid identifier of the rotation origin
                jp = trunc(((double)rng()/((double)rng.max()+1))*2);                    // rotate lower or higher part
                accept = rotPsi(sp, Chn, ip, jp, deltaE);
                break;
            case 3:
                ip = trunc(((double)rng()/((double)rng.max()+1))*sp->N_CH);
                accept = translation(sp, Chn, ip, deltaE);
                break;
            case 4:
                ip = trunc(((double)rng()/((double)rng.max()+1))*sp->N_CH);
                accept = rotation(sp, Chn, ip, deltaE);
                break;
            default:
                std::cerr << "Unusable move type selected: movetype = " << movetype << endl;
        }

        if(accept) {        // legal move
            Eold += deltaE;
            // sync HBDist[] and HBDcpy[]
            switch( movetype ) {
                case 0:
                    if( jp == 0 ) {
                        for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                            NCDcpy[ip][m] = NCDist[ip][m];
                        }
                    }
                    if( jp == 2 ) {
                        for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                            NCDcpy[m][ip] = NCDist[m][ip];
                        }
                    }
                    break;
                default:
                    for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                        for( int n=0; n<sp->N_CH*sp->N_AA; n++ ) {
                            NCDcpy[m][n] = NCDist[m][n];
                            NCDcpy[n][m] = NCDist[n][m];
                        }
                    }
            }
            // neighbour list is updated in move function - only if it returns true            
        }
        else {              // illegal move
            // reset HBList and HBDist
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                HBList[k][0] = HBLcpy[k][0];
                HBList[k][1] = HBLcpy[k][1];
            }
            // all legal steps are accepted → why copy NCDist[][] ? 
            switch( movetype ) {                                // non-SAMC move: all legal moves are accepted -> why copy NCDist[][] ?
                case 0:
                    if( jp == 0 ) {
                        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) { NCDist[ip][k] = NCDcpy[ip][k]; }
                    }
                    if( jp == 2 ) {
                        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) { NCDist[k][ip] = NCDcpy[k][ip]; }
                    }
                    break;
                default:
                    for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                        for( int n=0; n<sp->N_CH*sp->N_AA; n++ ) {
                            NCDist[m][n] = NCDcpy[m][n];
                            NCDist[n][m] = NCDcpy[n][m];
                        }
                    }
            }
        }

        // check bond length after every move - if this failes: abort run
        if( !checkBndLngth(sp, Chn, 0, sp->N_CH*sp->N_AA) ) {
            hd->os_log<< std::endl << "bond length error: t=" << step << std::endl << "Energy = " << Eold << std::endl; hd->os_log.close();
            std::cerr << std::endl << "bond length error: t=" << step << std::endl << "Energy = " << Eold << std::endl;
            outputPositions(sp,hd, hd->dbposi, Chn, 1, Eold);
            this_thread::sleep_for(chrono::milliseconds(200));
            return 0;
        }
        if( abs(Eold-E_check(sp, Chn)) > 0.01 ) {
            hd->os_log<< std::endl << "--- ERROR ---\tenergies are not equal: Eold=" << std::fixed << std::setprecision(3) << Eold << " Ecur=" << E_check(sp, Chn) << std::endl; hd->os_log.close();
            std::cerr << std::endl << "--- ERROR ---\tenergies are not equal: Eold=" << std::fixed << std::setprecision(3) << Eold << " Ecur=" << E_check(sp, Chn) << std::endl;
            std::cerr.precision(3); std::cerr << std::fixed;
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                std::cerr << "k" << k << "\t";
                for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                    std::cerr << NCDist[k][m] << "\t";
                } std::cerr << std::endl;
            } std::cerr << std::endl;
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                std::cerr << k << "  " << HBList[k][0] << "\t" << HBList[k][1] << std::endl;
            }
            outputPositions(sp, hd, hd->dbposi, Chn, 1, Eold);
            this_thread::sleep_for(chrono::milliseconds(200));

            return 0;
        }
        step++;
    }
    std::cerr << std::endl << "starting with energy E=" << Eold << std::endl;
    // configuration when starting SAMC
    outputPositions(sp,hd, hd->dbposi, Chn, 1, Eold);

    /*              XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    XXXXXXXXXXXXXXXXXXXX    SAMC loop    XXXXXXXXXXXXXXXXXXXX
                    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */    

    for( step=tcont; step<sp->T_MAX; step++ ) {
        // lets go
        for( it=0; it<sp->stepit; it++ ) {

            if( step%((int)1e2) == 0 ) {
                Timer.PrintProgress(step-tcont, sp->T_MAX-tcont);
            }

            // select move type
            moveselec = trunc( ( (double)rng()/( (double)rng.max()+1 ) )*(sp->WT_WIGGLE + sp->WT_PHI + sp->WT_PSI + sp->WT_TRANS + sp->WT_ROT));
            //moveselec = trunc(realdist01(rng)*(sp->WT_WIGGLE + sp->WT_PHI + sp->WT_PSI + sp->WT_TRANS));
            if( moveselec < sp->WT_WIGGLE )                                                      { movetype = 0; }  // wiggle
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI) )                                    { movetype = 1; }  // rotPhi
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI) )                         { movetype = 2; }  // rotPsi
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI+sp->WT_TRANS) )            { movetype = 3; }  // translation
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI+sp->WT_TRANS+sp->WT_ROT) ) { movetype = 4; }  // rotation
            else { 
                hd->os_log<< std::endl << "--- ERROR ---\tno movetype was selected" << std::endl; hd->os_log.close();
                std::cerr << std::endl << "--- ERROR ---\tno movetype was selected" << std::endl; return 0; 
            }
            // chain copy setup
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                HBLcpy[k][0] = HBList[k][0];            // improve efficiency! → select copied part based on movetype
                HBLcpy[k][1] = HBList[k][1];
            }
            deltaE = 0.0;
            ot->nattempt[movetype]++;
            switch( movetype ) {
                case 0:
                    i_rand = trunc(((double)rng()/((double)rng.max()+1))*(sp->N_CH*sp->N_AA*4));
                    ip = i_rand/4;                              // amino acid identifier
                    jp = i_rand%4;                              // bead
                    BdCpy[ip*4+jp] = Chn[ip/sp->N_AA].AmAc[ip%sp->N_AA].Bd[jp];
                    accept = wiggle(sp, Chn, ip/sp->N_AA, ip%sp->N_AA, jp, deltaE);
                    break;
                case 1:
                    ip = trunc(((double)rng()/((double)rng.max()+1))*(sp->N_CH*sp->N_AA));  // amino acid identifier of the rotation origin
                    jp = trunc(((double)rng()/((double)rng.max()+1))*2);                    // rotate lower or higher part
                    switch( jp ) {
                        case 0:
                            for( int i=(ip/sp->N_AA)*sp->N_AA; i<ip+1; i++ ) {
                                for( int j=0; j<4; j++ ) {
                                    BdCpy[i*4+j] = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j];
                                }
                            }
                            break;
                        case 1:
                            for( int i=ip; i<((ip/sp->N_AA)+1)*sp->N_AA; i++ ) {
                                for( int j=0; j<4; j++ ) {
                                    BdCpy[i*4+j] = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j];
                                }
                            }
                    }
                    accept = rotPhi(sp, Chn, ip, jp, deltaE);
                    break;
                case 2:
                    ip = trunc(((double)rng()/((double)rng.max()+1))*(sp->N_CH*sp->N_AA));  // amino acid identifier of the rotation origin
                    jp = trunc(((double)rng()/((double)rng.max()+1))*2);                    // rotate lower or higher part
                    switch( jp ) {
                        case 0:
                            for( int i=(ip/sp->N_AA)*sp->N_AA; i<ip+1; i++ ) {
                                for( int j=0; j<4; j++ ) {
                                    BdCpy[i*4+j] = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j];
                                }
                            }
                            break;
                        case 1:
                            for( int i=ip; i<((ip/sp->N_AA)+1)*sp->N_AA; i++ ) {
                                for( int j=0; j<4; j++ ) {
                                    BdCpy[i*4+j] = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j];
                                }
                            }
                    }
                    accept = rotPsi(sp, Chn, ip, jp, deltaE);
                    break;
                case 3:
                    ip = trunc(((double)rng()/((double)rng.max()+1))*sp->N_CH);
                    for( int i=0; i<sp->N_AA; i++ ) {
                        for( int j=0; j<4; j++ ) {
                            BdCpy[ip*sp->N_AA*4+i*4+j] = Chn[ip].AmAc[i].Bd[j];
                        }
                    }
                    accept = translation(sp, Chn, ip, deltaE);
                    break;
                case 4:
                    ip = trunc(((double)rng()/((double)rng.max()+1))*sp->N_CH);
                    for( int i=0; i<sp->N_AA; i++ ) {
                        for( int j=0; j<4; j++ ) {
                            BdCpy[ip*sp->N_AA*4+i*4+j] = Chn[ip].AmAc[i].Bd[j];
                        }
                    }
                    accept = rotation(sp, Chn, ip, deltaE);
                    break;
                default:
                    std::cerr << "Unusable move type selected: movetype = " << movetype << endl;
            }

            if( accept ) {      // legal move → check energy window and acceptance criterion
                Enew = Eold + deltaE;
                // check if Enew is in energy window → if not: don't count attempt
                if( Enew<sp->EMin || Enew > sp->EMax+0.00001 ) {
                    switch( movetype ) {
                        case 0:     // wiggle
                            Chn[ip/sp->N_AA].AmAc[ip%sp->N_AA].Bd[jp] = BdCpy[4*ip+jp];                     // reset to old coordinates
                            if( jp == 0 ) {
                                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) { NCDist[ip][k] = NCDcpy[ip][k]; } // reset NCDist[][]
                            }
                            if( jp == 2 ) {
                                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) { NCDist[k][ip] = NCDcpy[k][ip]; } // reset NCDist[][]
                            }
                            break;
                        case 3:     // translation
                            for( int i=0; i<sp->N_AA; i++ ) {
                                for( int j=0; j<4; j++ ) {
                                    oldBox = Chn[ip].AmAc[i].Bd[j].getBox();
                                    Chn[ip].AmAc[i].Bd[j] = BdCpy[ip*4*sp->N_AA+4*i+j];
                                    Chn[ip].AmAc[i].Bd[j].setBox(oldBox);
                                    LinkListUpdate(sp, Chn, ip*sp->N_AA+i, j);
                                }
                                for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                                    if( j<ip*sp->N_AA || j>=(ip+1)*sp->N_AA ) {
                                        NCDist[i+ip*sp->N_AA][j] = NCDcpy[i+ip*sp->N_AA][j];
                                        NCDist[j][i+ip*sp->N_AA] = NCDcpy[j][i+ip*sp->N_AA];
                                    }
                                }
                            }
                            break;
                        case 4:     // chain rotation
                            for( int i=0; i<sp->N_AA; i++ ) {
                                for( int j=0; j<4; j++ ) {
                                    oldBox = Chn[ip].AmAc[i].Bd[j].getBox();
                                    Chn[ip].AmAc[i].Bd[j] = BdCpy[ip*4*sp->N_AA+4*i+j];
                                    Chn[ip].AmAc[i].Bd[j].setBox(oldBox);
                                    LinkListUpdate(sp, Chn, ip*sp->N_AA+i, j);
                                }
                                for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                                    if( j<ip*sp->N_AA || j>=(ip+1)*sp->N_AA ) {
                                        NCDist[i+ip*sp->N_AA][j] = NCDcpy[i+ip*sp->N_AA][j];
                                        NCDist[j][i+ip*sp->N_AA] = NCDcpy[j][i+ip*sp->N_AA];
                                    }
                                }
                            }
                            break;
                        default:    // pivot rotations
                            switch( jp ) {
                                case 0:
                                    for( int i=(ip/sp->N_AA)*sp->N_AA; i<ip+1; i++ ) {
                                        for( int j=0; j<4; j++ ) {
                                            oldBox = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBox();       // oldBox has to be transfered in order for LinListUpdate to update correctly
                                            Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j] = BdCpy[4*i+j];          // reset to old coordinates
                                            Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].setBox(oldBox);
                                            LinkListUpdate(sp, Chn, i, j);                              // reset LinkedList
                                        }
                                        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {                          // reset NCDist[][]
                                            if( j<(ip/sp->N_AA)*sp->N_AA || j>ip-1 ) {
                                                NCDist[i][j] = NCDcpy[i][j];
                                                NCDist[j][i] = NCDcpy[j][i];
                                            }
                                        }
                                    }
                                    break;
                                case 1:
                                    for( int i=ip; i<((ip/sp->N_AA)+1)*sp->N_AA; i++ ) {
                                        for( int j=0; j<4; j++ ) {
                                            oldBox = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBox();       // oldBox has to be transfered in order for LinListUpdate to update correctly
                                            Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j] = BdCpy[4*i+j];          // reset to old coordinates
                                            Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].setBox(oldBox);
                                            LinkListUpdate(sp, Chn, i, j);                              // reset LinkedList
                                        }
                                        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {                          // reset NCDist[][]
                                            if( j<ip+1 || j>=((ip/sp->N_AA)+1)*sp->N_AA ) {
                                                NCDist[i][j] = NCDcpy[i][j];
                                                NCDist[j][i] = NCDcpy[j][i];
                                            }
                                        }
                                    }
                            }
                    }

                    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
                        HBList[i][0] = HBLcpy[i][0];
                        HBList[i][1] = HBLcpy[i][1];
                    }
                    ot->nattempt[movetype]--;
                    /*Ecur = E_check(Chn);
                    if( abs(Ecur-Eold) > 0.01 ) {
                        outputPositions(Chn, "AnorLondo.xyz", 1);
                        std::cerr << "hahahahahahah versager t=" << t << endl;
                    }*/
                    continue;
                }
                // SAMC acceptance step
                if(sp->EBIN_TRUNC_UP) { eBin_n = floor(((Enew-sp->EMin)/sp->BinW)+0.00001); }
                else { eBin_n = floor(((Enew-sp->EMin)/sp->BinW)-0.00001); }
                if( eBin_n == sp->NBin) { eBin_n = sp->NBin-1; }            // E=0 lands in non-existing bin → belongs to highest bin
                if( eBin_n == -1 && Enew > sp->EMin-0.00001) { eBin_n = 0; }                  // E=EMin lands in non-existing bin if EBIN_TRUNK_UP = false → belongs to lowest bin

                accept = acceptance(ot->lngE[eBin_o], ot->lngE[eBin_n]);
            }
            if( accept ) {      // Move accepted → update lngE, H, ...
                ot->naccept[movetype]++;
                ot->H[eBin_n]++;
                if(!sp->FIX_lngE) ot->lngE[eBin_n] += gamma;
                Eold = Enew;
                eBin_o = eBin_n;
                // sync HBDist[] and HBDcpy[]
                switch( movetype ) {
                    case 0:
                        if( jp == 0 ) {
                            for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                                NCDcpy[ip][m] = NCDist[ip][m];
                            }
                        }
                        if( jp == 2 ) {
                            for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                                NCDcpy[m][ip] = NCDist[m][ip];
                            }
                        }
                        break;
                    default:
                        for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                            for( int n=0; n<sp->N_CH*sp->N_AA; n++ ) {
                                NCDcpy[m][n] = NCDist[m][n];
                                NCDcpy[n][m] = NCDist[n][m];
                            }
                        }
                }
                // neighbour list is updated in move functions - only if it returns true
            }
            else {
                ot->H[eBin_o]++;
                if(!sp->FIX_lngE) ot->lngE[eBin_o] += gamma;
                // Bead and neighbour list reset
                switch( movetype ) {
                    case 0:     // wiggle
                        Chn[ip/sp->N_AA].AmAc[ip%sp->N_AA].Bd[jp] = BdCpy[4*ip+jp];                         // reset old coordinates
                        if( jp == 0 ) {
                            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) { NCDist[ip][k] = NCDcpy[ip][k]; }          // reset NCDist[][]
                        }
                        if( jp == 2 ) {
                            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) { NCDist[k][ip] = NCDcpy[k][ip]; }          // reset NCDist[][]
                        }
                        break;
                    case 3:     // translation
                        for( int i=0; i<sp->N_AA; i++ ) {
                            for( int j=0; j<4; j++ ) {
                                oldBox = Chn[ip].AmAc[i].Bd[j].getBox();
                                Chn[ip].AmAc[i].Bd[j] = BdCpy[ip*4*sp->N_AA+4*i+j];
                                Chn[ip].AmAc[i].Bd[j].setBox(oldBox);
                                LinkListUpdate(sp, Chn, ip*sp->N_AA+i, j);
                            }
                            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                                if( j<ip*sp->N_AA || j>=(ip+1)*sp->N_AA ) {
                                    NCDist[i+ip*sp->N_AA][j] = NCDcpy[i+ip*sp->N_AA][j];
                                    NCDist[j][i+ip*sp->N_AA] = NCDcpy[j][i+ip*sp->N_AA];
                                }
                            }
                        }
                        break;
                    case 4:     // chain rotation
                        for( int i=0; i<sp->N_AA; i++ ) {
                            for( int j=0; j<4; j++ ) {
                                oldBox = Chn[ip].AmAc[i].Bd[j].getBox();
                                Chn[ip].AmAc[i].Bd[j] = BdCpy[ip*4*sp->N_AA+4*i+j];
                                Chn[ip].AmAc[i].Bd[j].setBox(oldBox);
                                LinkListUpdate(sp, Chn, ip*sp->N_AA+i, j);
                            }
                            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                                if( j<ip*sp->N_AA || j>=(ip+1)*sp->N_AA ) {
                                    NCDist[i+ip*sp->N_AA][j] = NCDcpy[i+ip*sp->N_AA][j];
                                    NCDist[j][i+ip*sp->N_AA] = NCDcpy[j][i+ip*sp->N_AA];
                                }
                            }
                        }
                        break;
                    default:    // pivot rotation
                        switch( jp ) {
                            case 0:
                                for( int i=(ip/sp->N_AA)*sp->N_AA; i<ip+1; i++ ) {
                                    for( int j=0; j<4; j++ ) {
                                        oldBox = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBox();           // oldBox has to be transfered in order for LinListUpdate to update correctly
                                        Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j] = BdCpy[4*i+j];              // reset old coordinates
                                        Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].setBox(oldBox);
                                        LinkListUpdate(sp, Chn, i, j);                                  // reset LinkedList
                                    }
                                    for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {                              // reset NCDist[][]
                                        if( j<(ip/sp->N_AA)*sp->N_AA || j>ip-1 ) {
                                            NCDist[i][j] = NCDcpy[i][j];
                                            NCDist[j][i] = NCDcpy[j][i];
                                        }
                                    }
                                }
                                break;
                            case 1:
                                for( int i=ip; i<((ip/sp->N_AA)+1)*sp->N_AA; i++ ) {
                                    for( int j=0; j<4; j++ ) {
                                        oldBox = Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBox();           // oldBox has to be transfered in order for LinListUpdate to update correctly
                                        Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j] = BdCpy[4*i+j];              // reset old coordinates
                                        Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].setBox(oldBox);
                                        LinkListUpdate(sp, Chn, i, j);                                  // reset LinkedList
                                    }
                                    for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {                              // reset NCDist[][]
                                        if( j<ip+1 || j>=((ip/sp->N_AA)+1)*sp->N_AA ) {
                                            NCDist[i][j] = NCDcpy[i][j];
                                            NCDist[j][i] = NCDcpy[j][i];
                                        }
                                    }
                                }
                        }
                }
                // reset HBList and HBDist
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    HBList[k][0] = HBLcpy[k][0];
                    HBList[k][1] = HBLcpy[k][1];
                }
            }
            
            if(sp->HB_ContMat) {
                for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
                    if( HBList[i][0] > -1 ) { ot->contHB[ eBin_o*sp->N_CH*sp->N_AA*sp->N_CH*sp->N_AA + i*sp->N_CH*sp->N_AA + HBList[i][0] ] += 1; }
                }
            }
            if(sp->Ree) {
                for( int i=0; i<sp->N_CH; i++ ) {
                    std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[i].AmAc[0].Bd[0], Chn[i].AmAc[sp->N_AA-1].Bd[3]);
                    ot->Ree2[sp->N_CH*eBin_o + i] += dotPro(dist, dist);
                }
            }
            if(sp->wConfig) {
                for( int i=0; i<sp->ConfigE.size(); i++ ) {
                    if( abs(Eold - sp->ConfigE.at(i)) <= sp->ConfigV && ot->conf_n[i]<10 ) {
                        if( ot->conf_wt[i] + 10000 < step ) {
                            oss.str("");
                            oss << "coordinates_E" << i << "_" << ot->conf_n[i] << ".xyz";
                            filename = oss.str();
                            outputPositions(sp,hd, filename, Chn, 1, Eold);
                            ot->conf_n[i]++;
                            ot->conf_wt[i]=step;
                        }
                    }
                }
            }
            if(sp->tGyr) {
                calc_gyration_tensor(sp, ot, Chn, eBin_o);
                calc_gyration_radius(sp, ot, Chn, eBin_o);
                
                for( int i=0; i<sp->N_CH; i++ ) {
                    if( ot->rGyrCur[i] - (ot->tGyrEigCur[i][0]+ot->tGyrEigCur[i][1]+ot->tGyrEigCur[i][2]) > 1e-10) {
                        std::cout << std::fixed << std::setprecision(5) << std::endl << "--- ERROR ---\tgyration radius does not match eigenvalues. Chn[" << i << "]" << std::endl << "             \trGyr=" << ot->rGyrCur[i] << "  tGyrX²+tGyrY²+tGyrZ²=" << (ot->tGyrEigCur[i][0]+ot->tGyrEigCur[i][1]+ot->tGyrEigCur[i][2]) << std::endl;
                        hd->os_log << std::endl << "--- ERROR ---\tgyration radius does not match eigenvalues. Chn[" << i << "]" << std::endl << "             \trGyr=" << ot->rGyrCur[i] << "  tGyrX²+tGyrY²+tGyrZ²=" << (ot->tGyrEigCur[i][0]+ot->tGyrEigCur[i][1]+ot->tGyrEigCur[i][2]) << std::endl;
                    }
                }
            }
            if(sp->vdWener){
                if( movetype==0 && jp!=3 && (step+it)!=0 ) { 
                    calc_vanderWaals(sp, ot, Chn, eBin_o, vdW_intra, vdW_inter, 0); }
                else {calc_vanderWaals(sp, ot, Chn, eBin_o, vdW_intra, vdW_inter, 1); }
            }
            if( Eold < Egrd ) {
                Egrd = Eold;
                outputPositions(sp, hd, hd->grdcnm, Chn, 0, Eold);
            }

            gammasum += gamma;

        }   // end of one MC step

        //gamma update
        if(!sp->FIX_lngE){
            gamma = sp->GAMMA_0*sp->T_0/((double)max(sp->T_0,step));       // from Liang et al.
        }

        //  energy time development observable
        if(sp->Et) {
            ot->Et[step%sp->T_WRITE] = Eold;
        }

        // write Backup-File
        if( (step+1)%sp->T_WRITE == 0 ) {
            if(!sp->FIX_lngE) {
                BackupSAMCrun(sp, hd, ot, Chn, Timer, step, gammasum, gamma, Eold);
            } else {
                BackupProdRun(sp, hd, ot, Timer, step);
            }
            if(sp->HB_ContMat) {
                output_HBmat(sp, hd, ot, step+1);
            }
            if(sp->Ree) {
                output_Ree(sp, hd, ot, step+1);
            }
            if(sp->tGyr) {
                output_tGyr(sp, hd, ot, step+1);
            }
            if(sp->vdWener) {
                output_vdW(sp, hd, ot, step+1);
            }
            if(sp->Et) {
                output_Et(sp, hd, ot, step+1, 1);
            }
        }

        // Energy check
        if( (step+1)%10000 == 0 ) {
            if( abs(Eold-E_check(sp, Chn)) > 0.01 ) {
                hd->os_log<< endl << "--- std:: ---\tenergies are not equal: Eold=" << std::fixed << std::setprecision(3) << Eold << " Ecur=" << E_check(sp, Chn) << " at step=" << step << endl;
                std::cerr << endl << "--- ERROR ---\tenergies are not equal: Eold=" << std::fixed << std::setprecision(3) << Eold << " Ecur=" << E_check(sp, Chn) << " at step=" << step << endl;
                std::cerr << endl << "HBList" << endl;
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    std::cerr << k << "  " << HBList[k][0] << "\t" << HBList[k][1] << std::endl << std::flush;
                }
                std::cerr << endl;

                outputPositions(sp,hd, hd->dbposi, Chn, 1, Eold);

                std::cerr << "eternal darkness awaits…" << endl << std::flush;
            }
        }

    }   // end of SAMC main loop

    // calculated positions after SAMC loop
    outputPositions(sp,hd, hd->dbposi, Chn, 1, Eold);

    Timer.PrintProgress(step-tcont, sp->T_MAX-tcont);
    this_thread::sleep_for(chrono::milliseconds(200));
    std::printf("\npraise the sun ☀\n");

    Timer.endProgram();


    hd->os_log.close();

    output_memory_deallocation(sp, ot);

    delete[] neighHead;
    delete[] neighList;
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) { 
        delete[] HBList[i];
        delete[] HBLcpy[i];
        delete[] NCDist[i];
        delete[] NCDcpy[i];
    }
    delete[] HBList;
    delete[] HBLcpy;
    delete[] NCDist;
    delete[] NCDcpy;

    delete[] Chn;
    delete[] BdCpy;

    delete sp;
    delete hd;
    delete ot;

    return 0;
}


//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  >>   FUNCTION DEFINTIONS   <<  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double dotPro(double vecA[], double vecB[])
{
    double result = vecA[0]*vecB[0] + vecA[1]*vecB[1] + vecA[2]*vecB[2];
    return result;
}

double absVec(double vec[])
{
    double result = sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    return result;
}

tuple<double,double,double> crossPro(double vecA[], double vecB[])
{
    double x, y, z;
    x = vecA[1]*vecB[2] - vecA[2]*vecB[1];
    y = vecA[2]*vecB[0] - vecA[0]*vecB[2];
    z = vecA[0]*vecB[1] - vecA[1]*vecB[0];
    return make_tuple(x,y,z);
}
// distance vector from vecA to vecB including periodic boundary conditions
tuple<double,double,double> distVecBC(SysPara *sp, Bead vecA, Bead vecB)
{
    // with periodic boundary conditions
    double x, y, z;
    x = vecB.getR(0)-vecA.getR(0) - sp->L*round( (vecB.getR(0)-vecA.getR(0)) / (double) sp->L );
    y = vecB.getR(1)-vecA.getR(1) - sp->L*round( (vecB.getR(1)-vecA.getR(1)) / (double) sp->L );
    z = vecB.getR(2)-vecA.getR(2) - sp->L*round( (vecB.getR(2)-vecA.getR(2)) / (double) sp->L );
    return make_tuple(x,y,z);
}

//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXXXX  CONSOLE HEAD PRINT  XXXXXXXXXXXX
int program_start_print(ostream &os) 
{
    os << "###########################################" << std::endl
       << "              SAMC in Prime20              " << std::endl
       << "###########################################" << std::endl;
    return 0;
}
int command_print(Header *hd, ostream &os)
{
    os << "Input Initial configuration: " << hd->confnm << std::endl
       << "Input System parameters:     " << hd->paranm << std::endl
       << "Input DOS:                   " << hd->lngEnm << std::endl
       << "Input rerun:                 " << hd->rrunnm << std::endl
       << "Output debug position:       " << hd->dbposi << std::endl
       << "Output HB matrix:            " << hd->hbmatr << std::endl
       << "Output inertia tensor:       " << hd->tGyrnm << std::endl
       << "Output van-der-Waals energy: " << hd->vdWnm  << std::endl
       << "Output ground state config.: " << hd->grdcnm << std::endl
       << "Output simulation log:       " << hd->lognm  << std::endl;
    return 0;
}
int system_parameter_print(SysPara *sp, ostream &os)
{
    os << ">> system information <<" << std::endl
       << "No. of chains:            " << sp->N_CH << std::endl
       << "Degree of polymerization: " << sp->N_AA << std::endl
       << "Amino acid sequence:      " << sp->AA_seq << std::endl
       << "Length of simulation box: " << sp->L << std::endl;

    return 0;
}
int sim_parameter_print(SysPara *sp, ostream &os)
{
    int obs = 0;
    os << ">> SAMC parameters <<" << std::endl
       << "Allowed energy range:     [" << sp->EMin << ";" << sp->EMax << "]" << std::endl
       << "# of energy bins:         " << sp ->NBin << std::endl
       << "# of SAMC steps:          " << sp->T_MAX << std::endl;
        if(sp->FIX_lngE == true ) {
            os << "Using fixed lng(E)" << std::endl; }
        else {
            os << "T_0:                      " << sp->T_0 << std::endl
               << "gamma_0:                  " << sp->GAMMA_0 << std::endl; }
    os << "frequency of DOS output   " << sp->T_WRITE << std::endl
       << ">> observables <<" << std::endl;
        if( sp->HB_ContMat == true ) { obs = 1;
            os << "- Hydrogen bond contact matrices" << std::endl; }
        if( sp->tGyr == true ) { obs = 1;
            os << "- Tensor of gyration" << std::endl; }
        if( sp->wConfig == true ) { obs = 1;
            os << "- Configurations at energies:";
            for( int i=0; i<sp->ConfigE.size(); i++ ) { os << " " << sp->ConfigE.at(i); }
            os << std::endl << "  → with deviation tolerance dE = " << sp->ConfigV << std::endl; }
        if( sp->Et == true ) { obs = 1;
            os << "- E(t) energy time development" << std::endl; }
        if( sp->vdWener == true ) { obs = 1;
            os << "- van-der-Waals energy (intra- & inter-chain)" << std::endl; }
        if(obs == 0) {
            os << "... no observables" << std::endl;
        }

    return 0;
}

//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXX  INPUT OUTPUT FUNCTIONS  XXXXXXXXXX

// identify file names from command input
int CommandInitialize(int argc, char *argv[], Header *hd)
{
    bool display_help = false;
    std::string opt1, opt2;
    //default values
    hd->confnm = "config_ini.xyz";
    hd->paranm = "Sys_param.dat";
    hd->lngEnm = "lngE_input.dat";
    hd->rrunnm = "rerun.dat";
    hd->dbposi = "AnorLondo.xyz";
    hd->hbmatr = "HBmat.dat";
    hd->reenm  = "ReeDAve.dat";
    hd->tGyrnm = "tGyr.dat";
    hd->vdWnm  = "vdWE.dat";
    hd->grdcnm = "grdConfig.xyz";
    hd->enertm = "Et.dat";
    hd->lognm  = "out.log";

    //reading arguments
    if((argc+1)%2 == 0) {
        for( int i=1; i<argc; i+=2 ) {
            opt1 = argv[i];
            opt2 = argv[i+1];
            if( opt1.compare("-c") == 0 ) {
                hd->confnm = opt2; }
            else if( opt1.compare("-p") == 0 ) {
                hd->paranm = opt2; }
            else if( opt1.compare("-l") == 0 ) {
                hd->lngEnm = opt2; }
            else if( opt1.compare("-r") == 0 ) {
                hd->rrunnm = opt2; }
            else if( opt1.compare("-d") == 0 ) {
                hd->dbposi = opt2; }
            else if( opt1.compare("-h") == 0 ) {
                hd->hbmatr = opt2; }
            else if( opt1.compare("-g") == 0 ) {
                hd->tGyrnm = opt2; }
            else {
                display_help = true;
            }
        }
    }
    else {display_help = true;}

    if(display_help == true) {
        std::cout << "Unable to process argument input." << std::endl
                  << "Usage: " << argv[0] << " [OPTION1] [FILE1] [OPTION2] [FILE2] ..." << std::endl
                  << "  -c, input initial configuration" << std::endl
                  << "  -p, input system parameters" << std::endl
                  << "  -l, input density of states (lng[E])" << std::endl
                  << "  -r, input rerun data" << std::endl
                  << "  -d, output positions for debugging" << std::endl
                  << "  -h, output hydrogen bond matrices" << std::endl
                  << "  -g, output tensor of gyration" << std::endl
                  << "If not specified default file names will be used" << std::endl;
        return -1;
    }

    hd->os_log.open(hd->lognm);

    return 0;
}
// prints simulation start time
int cur_time_print(ostream &os)
{
    time_t curtime;
    time(&curtime);
    os << std::endl << "Start time: " << ctime(&curtime) << endl;
    
    return 0;
}
// creates new chain Chn with N_AA amino acids specified in AA_seq. returns true if successful
bool newChain(SysPara *sp, Chain Chn[], int chnNum)
{
    AmiAc AAinsert;
    double rad, ZrotAng, XrotAng;                           // radius for y-z calc, rotation angle
    double variaMtrx[3];                                    // variable matrix for side chain
    double ZrotMtrx[3][3], XrotMtrx[3][3];                  // rotation matrix around Z and X axis
    double tranVec[3], newpos[3], preRot[3], postRotZ[3];   // translation vector, new position after translation, coord. before rotation, coord. after rotation

    Chn[chnNum].setChnNo(chnNum);
    Chn[chnNum].AmAc.clear();
    for (int i=0; i<sp->N_AA; i++) {
        AAinsert.Setup(sp->AA_seq, i);                                      // sets up bond lengths, angles, squeeze parameters
        Chn[chnNum].addM( MASS_N+MASS_C+MASS_O+MASS_R(sp->AA_seq.at(i)) );  // add mass of amino acid to chain mass
        Chn[chnNum].AmAc.push_back(AAinsert);                               // adds amino acid to chain
    }

        // rotation matrix around Z axis
    ZrotAng = ANGL_NCaC + ANGL_CNCa - ANGL_CaCN - M_PI;
    ZrotMtrx[0][0] = cos(ZrotAng);
    ZrotMtrx[0][1] = -sin(ZrotAng);
    ZrotMtrx[0][2] = 0.0;
    ZrotMtrx[1][0] = sin(ZrotAng);
    ZrotMtrx[1][1] = cos(ZrotAng);
    ZrotMtrx[1][2] = 0.0;
    ZrotMtrx[2][0] = 0.0;
    ZrotMtrx[2][1] = 0.0;
    ZrotMtrx[2][2] = 1.0;
        // rotation matrix around X axis
    XrotAng = -11.0*M_PI/10.0;
    XrotMtrx[0][0] = 1.0;
    XrotMtrx[0][1] = 0.0;
    XrotMtrx[0][2] = 0.0;
    XrotMtrx[1][0] = 0.0;
    XrotMtrx[1][1] = cos(XrotAng);
    XrotMtrx[1][2] = -sin(XrotAng);
    XrotMtrx[2][0] = 0.0;
    XrotMtrx[2][1] = sin(XrotAng);
    XrotMtrx[2][2] = cos(XrotAng);
        // translation vector
    tranVec[0] = - (BND_NCa-BND_CaC*cos(ANGL_NCaC) + BND_CN*cos(ANGL_CaCN-ANGL_NCaC));
    tranVec[1] = - (BND_CaC*sin(ANGL_NCaC) + BND_CN*sin(ANGL_CaCN-ANGL_NCaC));
    tranVec[2] = - 0.0;

    for(int i = 0; i < sp->N_AA; i++) {
        // relocation of existing chain
        if( i != 0 ) {
            // translation + double rotation
            for(int j = 0; j < i; j++) {            // amino acid
                for(int k = 0; k < 4; k++) {        // bead
                    for(int m = 0; m < 3; m++) {    // coordinate
                        preRot[m] = Chn[chnNum].AmAc[j].Bd[k].getR(m) + tranVec[m];            // translation
                    }
                    for(int n = 0; n < 3; n++) {    // coordinate
                        postRotZ[n] = preRot[0]*ZrotMtrx[n][0] + preRot[1]*ZrotMtrx[n][1] + preRot[2]*ZrotMtrx[n][2];       // rotation around Z axis
                    }
                    for(int n = 0; n < 3; n++) {
                        newpos[n] = postRotZ[0]*XrotMtrx[n][0] + postRotZ[1]*XrotMtrx[n][1] + postRotZ[2]*XrotMtrx[n][2];   // rotation around X axis
                    }
                    Chn[chnNum].AmAc[j].Bd[k].setR(newpos[0], newpos[1], newpos[2]);
                }
            }
        }
        // amino acid construction
            // backbone
        Chn[chnNum].AmAc[i].Bd[0].setR(0.0,0.0,0.0);
        Chn[chnNum].AmAc[i].Bd[1].setR(BND_NCa, 0.0, 0.0);
        Chn[chnNum].AmAc[i].Bd[2].setR(BND_NCa-BND_CaC*cos(ANGL_NCaC), BND_CaC*sin(ANGL_NCaC), 0.0);

            // side chain
        rad = Chn[chnNum].AmAc[i].getDisR(1)*sin(Chn[chnNum].AmAc[i].getAng(0));
        variaMtrx[0] = BND_NCa - Chn[chnNum].AmAc[i].getDisR(1)*cos(Chn[chnNum].AmAc[i].getAng(0));
        variaMtrx[1] = (    (variaMtrx[0]-Chn[chnNum].AmAc[i].Bd[2].getR(0))*(variaMtrx[0]-Chn[chnNum].AmAc[i].Bd[2].getR(0)) +         // (x - xC)^2
                            Chn[chnNum].AmAc[i].Bd[2].getR(1)*Chn[chnNum].AmAc[i].Bd[2].getR(1) + rad*rad -                             // + yC^2 + r^2
                            Chn[chnNum].AmAc[i].getDisR(2)*Chn[chnNum].AmAc[i].getDisR(2) )                                             // - CR^2
                            / ( 2*Chn[chnNum].AmAc[i].Bd[2].getR(1) );                                                          // / 2*yC
        variaMtrx[2] = -sqrt( rad*rad - variaMtrx[1]*variaMtrx[1] );

        Chn[chnNum].AmAc[i].Bd[3].setR(variaMtrx[0], variaMtrx[1], variaMtrx[2]);
    }
    // PBC, bead type and mass
    for(int i = 0; i < sp->N_AA; i++) {
        for(int j = 0; j < 4; j++) {
            for( int k=0; k<3; k++ ) {
                Chn[chnNum].AmAc[i].Bd[j].addBC(k, floor(Chn[chnNum].AmAc[i].Bd[j].getR(k)/sp->L ));
            }
            Chn[chnNum].AmAc[i].Bd[j].setR(Chn[chnNum].AmAc[i].Bd[j].getR(0) - sp->L*floor(Chn[chnNum].AmAc[i].Bd[j].getR(0)/sp->L), Chn[chnNum].AmAc[i].Bd[j].getR(1) - sp->L*floor(Chn[chnNum].AmAc[i].Bd[j].getR(1)/sp->L), Chn[chnNum].AmAc[i].Bd[j].getR(2) - sp->L*floor(Chn[chnNum].AmAc[i].Bd[j].getR(2)/sp->L));
            Chn[chnNum].AmAc[i].Bd[j].set_btype(j);
            switch(j) {
                case 0: Chn[chnNum].AmAc[i].Bd[j].setM(MASS_N); break;
                case 1: Chn[chnNum].AmAc[i].Bd[j].setM(MASS_C); break;
                case 2: Chn[chnNum].AmAc[i].Bd[j].setM(MASS_O); break;
                case 3: Chn[chnNum].AmAc[i].Bd[j].setM( MASS_R(Chn[chnNum].AmAc[i].get_AAalp()) ); break;
            }
                
        }
    }
    // create neighbour list
    for( int i=0; i<sp->N_AA; i++ ) {
        for ( int j=0; j<4; j++ ) {
            Chn[chnNum].AmAc[i].Bd[j].setBox(assignBox(sp, Chn[chnNum].AmAc[i].Bd[j]));
        }
    }
    for(int i = chnNum*sp->N_AA; i < (chnNum+1)*sp->N_AA; i++) {
        for(int j = 0; j < 4; j++) {
            LinkListInsert(sp, Chn, i, j);
        }
    }

    return true;
}
// read system arameters from file
bool readParaInput(SysPara *sp, Header *hd) 
{
    ifstream ifstr;
    stringstream ss_line, ss_ener;
    std::string s_line, option, value;
    double d;
    bool read_all = true;

    std::cout << "reading parameter input ... ";
    hd->os_log<< "reading parameter input ... ";

    int read_NCH = 0;
    int read_NAA = 0;
    int read_AAS = 0;
    int read_L   = 0;
    int read_NBin= 0;
    int read_EMin= 0;
    int read_EMax= 0;
    int read_ESt = 0;
    int read_tSt = 0;
    int read_T0  = 0;
    int read_TMax= 0;
    int read_TWrt= 0;
    int read_TBCR= 0;
    int read_Gam0= 0;
    int read_NBox= 0;
    int read_LBox= 0;
    int read_WTWi= 0;
    int read_WTPh= 0;
    int read_WTPs= 0;
    int read_WTTr= 0;
    int read_WTRo= 0;
    int read_Disp= 0;
    int read_DPhi= 0;
    int read_DPsi= 0;
    int read_DTrn= 0;
    int read_DRot= 0;
    int read_ETru= 0;
    int read_FixL= 0;
    int read_HBCM= 0;
    int read_Ree = 0;
    int read_tGyr= 0;
    int read_vdWe= 0;
    int read_Et  = 0;
    int read_WCon= 0;
    int read_ConE= 0;
    int read_ConV= 0;

    ifstr.open(hd->paranm);
    if( ifstr.is_open() ) {
        while( std::getline( ifstr, s_line ) ) {
            ss_line.clear();
            if( s_line.length() != 0 && s_line.at(0)!='#' && s_line.at(0)!='\n' && s_line.at(0)!='\0' ) {
                ss_line.str(s_line);
                std::getline(ss_line, option, '=');
                std::getline(ss_line, value, ';');
                if( option.compare("N_CH")==0 ) {
                    sp->N_CH = stoi(value, nullptr);     read_NCH = 1; }
                else if( option.compare("N_AA")==0 ) {
                    sp->N_AA = stoi(value, nullptr);     read_NAA = 1; }
                else if( option.compare("AA_seq")==0 ) {
                    sp->AA_seq = value;                  read_AAS = 1; }
                else if( option.compare("L")==0 ) {
                    sp->L = stod(value, nullptr);        read_L = 1; }
                else if( option.compare("NBin")==0 ) {
                    sp->NBin = stoi(value, nullptr);     read_NBin = 1; }
                else if( option.compare("EMin")==0 ) {
                    sp->EMin = stod(value, nullptr);     read_EMin = 1; }
                else if( option.compare("EMax")==0 ) {
                    sp->EMax = stod(value, nullptr);     read_EMax = 1; }
                else if( option.compare("EStart")==0 ) {
                    sp->EStart = stod(value, nullptr);   read_ESt = 1; }
                else if( option.compare("tStart")==0 ) {
                    sp->tStart = stoi(value, nullptr);   read_tSt = 1; }
                else if( option.compare("T_0")==0 ) {
                    istringstream constr(value); constr >> d;
                    sp->T_0 = round(d);                  read_T0 = 1; }
                else if( option.compare("T_MAX")==0 ) {
                    istringstream constr(value); constr >> d;
                    sp->T_MAX = round(d);                read_TMax = 1; }
                else if( option.compare("T_WRITE")==0 ) {
                    istringstream constr(value); constr >> d;
                    sp->T_WRITE = round(d);              read_TWrt = 1; }
                else if( option.compare("T_BC_RESET")==0 ) {
                    istringstream constr(value); constr >> d;
                    sp->T_BC_RESET = round(d);           read_TBCR = 1; }
                else if( option.compare("GAMMA_0")==0 ) {
                    sp->GAMMA_0 = stod(value, nullptr);  read_Gam0 = 1; }
                else if( option.compare("NBOX")==0 ) {
                    sp->NBOX = stoi(value, nullptr);     read_NBox = 1; }
                else if( option.compare("LBOX")==0 ) {
                    sp->LBOX = stod(value, nullptr);     read_LBox = 1; }
                else if( option.compare("WT_WIGGLE")==0 ) {
                    sp->WT_WIGGLE = stoi(value, nullptr); read_WTWi = 1; }
                else if( option.compare("WT_PHI")==0 ) {
                    sp->WT_PHI = stoi(value, nullptr);   read_WTPh = 1; }
                else if( option.compare("WT_PSI")==0 ) {
                    sp->WT_PSI = stoi(value, nullptr);   read_WTPs = 1; }
                else if( option.compare("WT_TRANS")==0 ) {
                    sp->WT_TRANS = stoi(value, nullptr); read_WTTr = 1; }
                else if( option.compare("WT_ROT")==0 ) {
                    sp->WT_ROT = stoi(value, nullptr); read_WTRo = 1; }
                else if( option.compare("DISP_MAX")==0 ) {
                    sp->DISP_MAX = stod(value, nullptr); read_Disp = 1; }
                else if( option.compare("DPHI_MAX")==0 ) {
                    sp->DPHI_MAX = stod(value, nullptr); read_DPhi = 1; }
                else if( option.compare("DPSI_MAX")==0 ) {
                    sp->DPSI_MAX = stod(value, nullptr); read_DPsi = 1; }
                else if( option.compare("DTRN_MAX")==0 ) {
                    sp->DTRN_MAX = stod(value, nullptr); read_DTrn = 1; }
                else if( option.compare("DROT_MAX")==0 ) {
                    sp->DROT_MAX = stod(value, nullptr); read_DRot = 1; }
                else if( option.compare("EBIN_TRUNC_UP")==0 ) {
                    if( value.compare("true")==0 ) { 
                        sp->EBIN_TRUNC_UP = true; read_ETru = 1; }
                    else if( value.compare("false")==0 ) { 
                        sp->EBIN_TRUNC_UP = false; read_ETru = 1; } }
                else if( option.compare("FIX_lngE")==0 ) {
                    if( value.compare("true")==0 ) { sp->FIX_lngE = true; read_FixL = 1; }
                    else if( value.compare("false")==0 ) { sp->FIX_lngE = false; read_FixL = 1; } }
                else if( option.compare("HB_ContMat")==0 ) {
                    if( value.compare("true")==0 ) { sp->HB_ContMat = true; read_HBCM = 1; }
                    else if( value.compare("false")==0 ) { sp->HB_ContMat = false; read_HBCM = 1; } }
                else if( option.compare("Ree")==0 ) {
                    if( value.compare("true")==0 ) { sp->Ree = true; read_Ree = 1; }
                    else if( value.compare("false")==0 ) { sp->Ree = false; read_Ree = 1; } }
                else if( option.compare("tGyr")==0 ) {
                    if( value.compare("true")==0 ) {sp->tGyr = true; read_tGyr = 1; }
                    else if( value.compare("false")==0 ) { sp->tGyr = false; read_tGyr = 1; } }
                else if( option.compare("vdWener")==0 ) {
                    if( value.compare("true")==0 ) {sp->vdWener = true; read_vdWe = 1; }
                    else if( value.compare("false")==0 ) { sp->vdWener = false; read_vdWe = 1; } }
                else if( option.compare("ener_t")==0 ) {
                    if( value.compare("true")==0 ) {sp->Et = true; read_Et = 1; }
                    else if( value.compare("false")==0 ) { sp->Et = false; read_Et = 1; } }
                else if( option.compare("wConfig")==0 ) {
                    if( value.compare("true")==0 ) { 
                        sp->wConfig = true;    read_WCon = 1;
                        for( int i=0; i<2; i++ ) {
                            std::getline(ifstr, s_line);
                            ss_line.clear();    ss_line.str(s_line);
                            std::getline(ss_line, option, '=');
                            std::getline(ss_line, value, ';');
                            if( option.compare("ConfigE")==0 ) {
                                ss_line.clear(); ss_line.str(value);
                                std::istream_iterator<std::string> begin(ss_line);
                                std::istream_iterator<std::string> end;
                                std::vector<std::string> v_ener(begin, end);
                                for( int j=0; j<v_ener.size(); j++ ) { sp->ConfigE.push_back(stod(v_ener.at(j))); }
                                read_ConE = 1;
                            }
                            else if( option.compare("ConfigV")==0 ) {
                                sp->ConfigV = stod(value, nullptr);
                                read_ConV = 1;
                            }
                        }
                    }
                    else if( value.compare("false")==0 ) { 
                        sp->wConfig = false; 
                        read_WCon = 1; read_ConE = 1; read_ConV = 1; } }
            }
        }
        sp->BinW = (sp->EMax - sp->EMin)/(double)sp->NBin;
        sp->stepit = 4*sp->N_AA*sp->N_CH;

        if( read_NCH  == 0 ){ read_all = false;
            std::cout  << "Warning! N_CH not read. " << std::endl; 
            hd->os_log << "Warning! N_CH not read. " << std::endl;}
        if( read_NAA  == 0 ){ read_all = false;
            std::cout  << "Warning! N_AA not read. " << std::endl; 
            hd->os_log << "Warning! N_AA not read. " << std::endl;}
        if( read_AAS  == 0 ){ read_all = false;
            std::cout  << "Warning! AAseq not read. " << std::endl; 
            hd->os_log << "Warning! AAseq not read. " << std::endl;}
        if( read_L    == 0 ){ read_all = false;
            std::cout  << "Warning! L not read. " << std::endl; 
            hd->os_log << "Warning! L not read. " << std::endl;}
        if( read_NBin == 0 ){ read_all = false;
            std::cout  << "Warning! NBin not read. " << std::endl; 
            hd->os_log << "Warning! NBin not read. " << std::endl;}
        if( read_EMin == 0 ){ read_all = false;
            std::cout  << "Warning! EMin not read. " << std::endl; 
            hd->os_log << "Warning! EMin not read. " << std::endl;}
        if( read_EMax == 0 ){ read_all = false;
            std::cout  << "Warning! EMax not read. " << std::endl; 
            hd->os_log << "Warning! EMax not read. " << std::endl;}
        if( read_ESt  == 0 ){ read_all = false;
            std::cout  << "Warning! EStart not read. " << std::endl; 
            hd->os_log << "Warning! EStart not read. " << std::endl;}
        if( read_tSt  == 0 ){ read_all = false;
            std::cout  << "Warning! tStart not read. " << std::endl; 
            hd->os_log << "Warning! tStart not read. " << std::endl;}
        if( read_T0   == 0 ){ read_all = false;
            std::cout  << "Warning! T_0 not read. " << std::endl; 
            hd->os_log << "Warning! T_0 not read. " << std::endl;}
        if( read_TMax == 0 ){ read_all = false;
            std::cout  << "Warning! T_Max not read. " << std::endl; 
            hd->os_log << "Warning! T_Max not read. " << std::endl;}
        if( read_TWrt == 0 ){ read_all = false;
            std::cout  << "Warning! T_Write not read. " << std::endl; 
            hd->os_log << "Warning! T_Write not read. " << std::endl;}
        if( read_TBCR == 0 ){ read_all = false;
            std::cout  << "Warning! T_BCReset not read. " << std::endl; 
            hd->os_log << "Warning! T_BCReset not read. " << std::endl;}
        if( read_Gam0 == 0 ){ read_all = false;
            std::cout  << "Warning! gamma_0 not read. " << std::endl; 
            hd->os_log << "Warning! gamma_0 not read. " << std::endl;}
        if( read_NBox == 0 ){ read_all = false;
            std::cout  << "Warning! NBox not read. " << std::endl; 
            hd->os_log << "Warning! NBox not read. " << std::endl;}
        if( read_LBox == 0 ){ read_all = false;
            std::cout  << "Warning! LBox not read. " << std::endl; 
            hd->os_log << "Warning! LBox not read. " << std::endl;}
        if( read_WTWi == 0 ){ read_all = false;
            std::cout  << "Warning! WT_wiggle not read. " << std::endl; 
            hd->os_log << "Warning! WT_wiggle not read. " << std::endl;}
        if( read_WTPh == 0 ){ read_all = false;
            std::cout  << "Warning! WT_phi not read. " << std::endl; 
            hd->os_log << "Warning! WT_phi not read. " << std::endl;}
        if( read_WTPs == 0 ){ read_all = false;
            std::cout  << "Warning! WT_psi not read. " << std::endl; 
            hd->os_log << "Warning! WT_psi not read. " << std::endl;}
        if( read_WTTr == 0 ){ read_all = false;
            std::cout  << "Warning! WT_trans not read. " << std::endl; 
            hd->os_log << "Warning! WT_trans not read. " << std::endl;}
        if( read_WTRo == 0 ){ read_all = false;
            std::cout  << "Warning! WT_rot not read. " << std::endl; 
            hd->os_log << "Warning! WT_rot not read. " << std::endl;}
        if( read_Disp == 0 ){ read_all = false;
            std::cout  << "Warning! DispMax not read. " << std::endl; 
            hd->os_log << "Warning! DispMax not read. " << std::endl;}
        if( read_DPhi == 0 ){ read_all = false;
            std::cout  << "Warning! dPhi_Max not read. " << std::endl; 
            hd->os_log << "Warning! dPhi_Max not read. " << std::endl;}
        if( read_DPsi == 0 ){ read_all = false;
            std::cout  << "Warning! dPsi_Max not read. " << std::endl; 
            hd->os_log << "Warning! dPsi_Max not read. " << std::endl;}
        if( read_DTrn == 0 ){ read_all = false;
            std::cout  << "Warning! dTrn_Max not read. " << std::endl; 
            hd->os_log << "Warning! dTrn_Max not read. " << std::endl;}
        if( read_DRot == 0 ){ read_all = false;
            std::cout  << "Warning! dRot_Max not read. " << std::endl; 
            hd->os_log << "Warning! dRot_Max not read. " << std::endl;}
        if( read_ETru == 0 ){ read_all = false;
            std::cout  << "Warning! EtruncUp not read. " << std::endl; 
            hd->os_log << "Warning! EtruncUp not read. " << std::endl;}
        if( read_FixL == 0 ){ read_all = false;
            std::cout  << "Warning! Fix_lngU not read. " << std::endl; 
            hd->os_log << "Warning! Fix_lngU not read. " << std::endl;}
        if( read_HBCM == 0 ){ read_all = false;
            std::cout  << "Warning! HB_ContMat not read. " << std::endl; 
            hd->os_log << "Warning! HB_ContMat not read. " << std::endl;}
        if( read_Ree  == 0 ){ read_all = false;
            std::cout  << "Warning! Ree not read. " << std::endl; 
            hd->os_log << "Warning! Ree not read. " << std::endl;}
        if( read_tGyr == 0 ){ read_all = false;
            std::cout  << "Warning! tGyr not read. " << std::endl; 
            hd->os_log << "Warning! tGyr not read. " << std::endl;}
        if( read_vdWe == 0 ){ read_all = false;
            std::cout  << "Warning! vdWener not read. " << std::endl; 
            hd->os_log << "Warning! vdWener not read. " << std::endl;}
        if( read_Et == 0 ){ read_all = false;
            std::cout  << "Warning! Et not read. " << std::endl; 
            hd->os_log << "Warning! Et not read. " << std::endl;}
        if( read_WCon == 0 ){ read_all = false;
            std::cout  << "Warning! wConfig not read. " << std::endl; 
            hd->os_log << "Warning! wConfig not read. " << std::endl;}
        if( read_ConE == 0 ){ read_all = false;
            std::cout  << "Warning! ConfigE not read. " << std::endl; 
            hd->os_log << "Warning! ConfigE not read. " << std::endl;}
        if( read_ConV == 0 ){ read_all = false;
            std::cout  << "Warning! ConfigV not read. " << std::endl; 
            hd->os_log << "Warning! ConfigV not read. " << std::endl;}

        ifstr.close();
        if(read_all) {
            std::cout << "complete" << std::endl;
            hd->os_log<< "complete" << std::endl;
            return true;  
        } 
        return false;
    }
    else {
        hd->os_log<< "failed" << std::endl << " -- ERROR --\t" << hd->paranm << " not found" << std::endl;
        std::cout << "failed" << std::endl << " -- ERROR --\t" << hd->paranm << " not found" << std::endl;
        return false;
    }
}
// read chain config from file
bool readCoord(SysPara *sp, Header *hd, Chain Chn[])
{
    AmiAc AAinsert;
    ifstream input;
    double distV[3], newposBC[3];
    double x, y, z;
    string ID;
    int i = 0;
    int NaaFile;

    input.open(hd->confnm);
    if( input.is_open() ) {
        std::cout << "reading coordiates from input file '" << hd->confnm << "'" << std::endl;
        for( int j=0; j<sp->N_CH; j++ ) {
            Chn[j].setChnNo(j);
            Chn[j].AmAc.clear();
            for( int k=0; k<sp->N_AA; k++ ) {
                AAinsert.Setup(sp->AA_seq, k);
                Chn[j].addM( MASS_N+MASS_C+MASS_O+MASS_R(sp->AA_seq.at(i)) );   // add mass of amino acid to chain mass
                Chn[j].AmAc.push_back(AAinsert);
            }
        }
        input >> NaaFile; input.ignore();
        while( input.good() && !input.eof() ) {
            input >> ID >> x >> y >> z;
            Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setR(x, y, z);
            Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].set_btype(i%4);

            switch(i%4) {
                case 0: Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setM(MASS_N); break;
                case 1: Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setM(MASS_C); break;
                case 2: Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setM(MASS_O); break;
                case 3: Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setM( MASS_R(Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].get_AAalp()) ); break;
            }
            input.ignore();
            i++;
            if(i == NaaFile) { break; }
        }
        if(i != sp->N_CH*sp->N_AA*4) {
            hd->os_log<< "--- ERROR ---\tincorrect number of beads: N_BB_real = " << i << "  N_BB_should = " << sp->N_AA*4 << std::endl;
            std::cout << "--- ERROR ---\tincorrect number of beads: N_BB_real = " << i << "  N_BB_should = " << sp->N_AA*4 << std::endl;
        }

        // PBC
        for( int i=0; i<sp->N_CH*sp->N_AA*4; i++ ) {
            for( int k=0; k<3; k++ ) {
                Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setBC(k, floor(Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].getR(k) / sp->L) );
                newposBC[k] = Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].getR(k) - sp->L*floor(Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].getR(k)/(double)sp->L);
            }
            Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setR(newposBC[0], newposBC[1], newposBC[2]);
        }

        // create neighbour list
        for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
            for( int j=0; j<4; j++ ) {
                Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].setBox( assignBox(sp, Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j]) );
            }
        }
        for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
            for( int j=0; j<4; j++ ) {
                LinkListInsert(sp, Chn, i, j);
            }
        }
        input.close();
        sp->tStart = 0;
        return true;
    }
    else {
        std::cerr << "unable to open coordinate input file '" << hd->confnm << "'" << endl;
        return false;
    }
}
// reads lngE, H, gammasum, and t from input file
bool readPrevRunInput(SysPara *sp, Header *hd, Output *ot, Chain Chn[], long unsigned int &tcont, double &gammasum)
{
    ifstream input;
    stringstream ss_line;
    std::string s_line, s_token;
    char delim = ' ';
    int num;
    double Ebin1, Ebin2;

    input.open(hd->rrunnm);
    if( input.is_open() ) {
        std::cout << "reading rerun data from input file '" << hd->rrunnm << "'" << std::endl;
        while( true ) {
            std::getline( input, s_line);
            if( s_line.length() == 0 ) {
                hd->os_log<< "--- ERROR ---\tencountered EOF inside input file head" << std::endl;
                std::cout << "--- ERROR ---\tencountered EOF inside input file head" << std::endl;
                input.close();
                return false;
            }
            if(s_line[0] == '#') {
                ss_line.str(s_line);
                std::getline(ss_line, s_token, delim);
                std::getline(ss_line, s_token, delim);
                if( s_token.length() == 0 ) continue;
                if( s_token.compare("number") == 0 ) {
                    std::getline(ss_line,s_token, delim);
                    std::getline(ss_line,s_token, delim);
                    std::getline(ss_line,s_token, delim);
                    std::getline(ss_line,s_token, delim);
                    tcont = stoi(s_token, nullptr);
                    std::cout << "  continue after step " << tcont << std::endl;
                }
                else if( s_token.compare("gammasum") == 0 ) {
                    std::getline(ss_line, s_token, delim);
                    std::getline(ss_line, s_token, delim);
                    gammasum = stod(s_token, nullptr);
                    std::cout << "  with gammasum = " << gammasum << std::endl;
                }
            }
            else break;
        }
        if( s_line[0] != 'b' ) {
            std::cout << "Aaaaaarg! What a terrible misfortune came upon us in these miserable times of pain and agony?" << std::endl << " → expected 'bin' at start of line. Instead '" << s_line << "' was read" << std::endl;
            input.close();
            return false;
        }
        for( int i=0; i<sp->NBin; i++ ) {
            if( input.good() ) {
                input >> num >> Ebin1 >> Ebin2 >> ot->lngE[i] >> ot->H[i];
                if( num != i ) {
                    hd->os_log<< "--- ERROR ---\tBad line number " << num << " != " << i << std::endl;
                    std::cout << "--- ERROR ---\tBad line number " << num << " != " << i << std::endl;
                    input.close();
                    return false;
                }
                //std::cout << "  " << i << " " << ot->lngE[i] << std::endl;
            }
            else {
                hd->os_log<< "--- ERROR ---\tencountered EOF inside lngE line << " << i << " expected " << sp->NBin-1 << std::endl;
                std::cout << "--- ERROR ---\tencountered EOF inside lngE line << " << i << " expected " << sp->NBin-1 << std::endl;
                input.close();
                return false;
            }
        }
        input.close();
        return true;
    }

    std::cout << "unable to open input file '" << hd->rrunnm << "'" << std::endl;
    return false;
}
// reads lngE data from file
bool read_lngE(SysPara *sp, Header *hd, Output *ot)
{
    ifstream input;
    std::string s_line;
    int num;
    long unsigned int H;
    double Ebin1, Ebin2;

    input.open(hd->lngEnm);
    if( input.is_open() ) {
        std::cout << "reading lngE from input file '" << hd->lngEnm << "'" << std::endl;
        std::getline(input, s_line);
        if(s_line[0] != 'b') {
            std::cout << "unexpected character at start of line in file '" << hd->lngEnm << "': expected 'b', instead '" << s_line[0] << std::endl;
            input.close();
            return false;
        }
        for( int i=0; i<sp->NBin; i++ ) {
            if( input.good() ) {
                input >> num >> Ebin1 >> Ebin2 >> ot->lngE[i] >> H;
                if( num != i ) {
                    hd->os_log<< "--- ERROR ---\tBad line number " << num << " != " << i << std::endl;
                    std::cout << "--- ERROR ---\tBad line number " << num << " != " << i << std::endl;
                    input.close();
                    return false;
                }
            }
            else {
                hd->os_log<< "--- ERROR ---\tencountered EOF inside lngE line << " << i << " expected " << sp->NBin-1 << std::endl;
                std::cout << "--- ERROR ---\tencountered EOF inside lngE line << " << i << " expected " << sp->NBin-1 << std::endl;
                input.close();
                return false;
            }
        }
        input.close();
        return true;
    }

    std::cout << "unable to open input file '" << hd->lngEnm << "'" << std::endl;
    return false;
}
// writes file called "AnorLondo.txt" with coordinates of beads (debug purpose)
bool outputPositions(SysPara *sp, Header *hd, std::string fnm, Chain Chn[], int mode, double ener)
{
    ofstream Checkpos;
    if(mode == 0) {
        Checkpos.open(fnm, ios::out);
    }
    else {
        Checkpos.open(fnm, ios::app);
    }
    if( Checkpos.is_open() ) {
        Checkpos << "# Config of " << sp->N_CH << " " << sp->N_AA << "-mer(s) with seq. " << sp->AA_seq << " at E=" << ener << std::endl;
        Checkpos << sp->N_CH*sp->N_AA*4 << std::endl;
        Checkpos << setprecision(3) << std::fixed;
        for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
            for( int j=0; j<4; j++ ) {
                Checkpos << Chn[i/sp->N_AA].getChnNo() << "-" << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].get_AAalp() << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].get_btype() << std::setw(3) << std::setfill('0') << i*4+j << " "
                         << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getR(0) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBC(0)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getR(1) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBC(1)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getR(2) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[j].getBC(2)*sp->L << endl;
            }
        }

        Checkpos << std::endl;
        Checkpos.close();
        return true;
    }
    else {
        hd->os_log<< "--- ERROR ---\tfailed to open " << fnm << endl;
        std::cerr << "--- ERROR ---\tfailed to open " << fnm << endl;
        return false;
    }
}
// write hydrogen bond matrix
bool output_HBmat(SysPara *sp, Header *hd, Output *ot, int step)
{
    ofstream ostr;
    ostr.open(hd->hbmatr, ios::out);
    if(ostr.is_open() ) {
        ostr << "# Hydrogen Bind contact matrices after " << step+1 << " steps" << std::endl;
        std::setprecision(3); std::fixed;
        for( int i=0; i<sp->NBin; i++ ) {
            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    ostr << ot->contHB[i*sp->N_CH*sp->N_AA*sp->N_CH*sp->N_AA + j*sp->N_CH*sp->N_AA + k]/ot->H[i] << " ";
                }
                ostr << std::endl;
            }
        }
        ostr.close();
        return true;
    }
    else {
        hd->os_log<< std::endl << "error opening " << hd->hbmatr << std::endl;
        std::cout << std::endl << "error opening " << hd->hbmatr << std::endl;
        return false;
    }
}
// write end-to-end-distance
bool output_Ree(SysPara *sp, Header *hd, Output *ot, int step)
{
    ofstream ostr;
    ostr.open(hd->reenm, ios::out);
    if( ostr.is_open() ) {
        ostr << "# Squared end-to-end distances after " << step << " steps" << std::endl;
        ostr << "bin H";
        for( int i=0; i<sp->N_CH; i++ ) { ostr << " chn" << i; } ostr << std::endl;
        for( int i=0; i<sp->NBin; i++ ) {
            ostr << i << " " << ot->H[i];
            for( int j=0; j<sp->N_CH; j++ ) { ostr << " " << std::setprecision(4) << std::fixed << ot->Ree2[2*i + j]/(double)ot->H[i]; }
            ostr << std::endl;
        }
        ostr.close();
        return true;
    }
    else {
        hd->os_log<< std::endl << "error opening " << hd->reenm << std::endl;
        std::cout << std::endl << "error opening " << hd->reenm << std::endl;
        return false;
    }
}
// write tensor of gyration
bool output_tGyr(SysPara *sp, Header *hd, Output *ot, int step)
{
    ofstream ostr;
    ostr.open(hd->tGyrnm, ios::out);
    if( ostr.is_open() ) {
        ostr << "# Eigenvalues of tensor of gyration after " << step << " steps" << std::endl
             << "bin H Gxx Gyy Gzz" << std::endl;
        for( int i=0; i<sp->N_CH; i++ ) {
            ostr << "chn" << i << std::endl;
            for( int j=0; j<sp->NBin; j++ ) {
                ostr << j << " " << ot->H[j];
                for( int k=0; k<3; k++ ) {
                    ostr << " " << std::setprecision(4) << std::fixed << (ot->tGyrEig[i][j][k])/(ot->H[j]);
                }   ostr << std::endl;
            }
        }
        ostr.close();
        return true;
    }
    else {
        hd->os_log<< std::endl << "--- ERROR ---\tcould not open file " << hd->tGyrnm << std::endl;
        std::cout << std::endl << "--- ERROR ---\tcould not open file " << hd->tGyrnm << std::endl;
        return false;
    }
}
// write van-der-Waals energy
bool output_vdW(SysPara *sp, Header *hd, Output *ot, int step)
{
    ofstream ostr;
    ostr.open(hd->vdWnm, ios::out);
    if( ostr.is_open() ) {
        ostr << "# van-der-Waals energy after " << step << " steps" << std::endl
             << "bin H intra inter" << std::endl;
        for( int i=0; i<sp->NBin; i++ ) {
            ostr << i << " " << ot->H[i];
            for( int j=0; j<2; j++ ) {
                ostr << " " << std::setprecision(4) << std::fixed << (ot->vdWener[i*2 + j])/(ot->H[i]);
            } ostr << std::endl;
        }
        ostr.close();
        return true;
    }
    else {
        hd->os_log<< std::endl << "--- ERROR ---\tcould not open file " << hd->vdWnm << std::endl;
        std::cout << std::endl << "--- ERROR ---\tcould not open file " << hd->vdWnm << std::endl;
        return false;
    }
}
// write energy time development
bool output_Et(SysPara *sp, Header *hd, Output *ot, int step, int init)
{
    ofstream ostr;
    double Eold, Enew;

    Eold = +2.0;

    if( init == 0 ) {
        ostr.open(hd->enertm, ios::out);
        if( ostr.is_open() ) {
            ostr << "# energy time development" << std::endl 
                 << "# step E" << std::endl;
            ostr.close();
            return true;
        }
        else {
            hd->os_log<< std::endl << "--- ERROR ---\tcould not open file " << hd->enertm << std::endl;
            std::cout << std::endl << "--- ERROR ---\tcould not open file " << hd->enertm << std::endl;
            return false;
        }
    }
    else if( init == 1 ) {
        ostr.open(hd->enertm, ios::app);
        if( ostr.is_open() ) {
            for( int i=0; i<sp->T_WRITE; i++ ) {
                Enew = ot->Et[i];
                if( abs(Enew-Eold) >= 1.0 ) {
                    ostr << ((step-1)/sp->T_WRITE)*sp->T_WRITE + i << " " << Enew << std::endl;
                    Eold = Enew;
                }
            }
            ostr.close();
            return true;
        }
        else {
            hd->os_log<< std::endl << "--- ERROR ---\tcould not open file " << hd->enertm << std::endl;
            std::cout << std::endl << "--- ERROR ---\tcould not open file " << hd->enertm << std::endl;
            return false;
        }
    }
    else {
        hd->os_log<< std::endl << "--- ERROR ---\tcould not write to file " << hd->enertm << std::endl;
        std::cout << std::endl << "--- ERROR ---\tcould not write to file " << hd->enertm << std::endl;
        return false;
    }

}
// writes backup file in SAMC run
bool BackupSAMCrun(SysPara *sp, Header *hd, Output *ot, Chain Chn[], Timer &Timer, unsigned long int t, double gammasum, double gamma, double E)
{
    std::ostringstream oss;
    oss << "SAMCbackup_" << t/sp->T_WRITE << ".dat";
    string name = oss.str();
    ofstream backup;
    backup.open(name);
    if( backup.is_open() ) {
        backup << "# SAMC simulation of " << sp->N_CH << " PRIME20 " << sp->N_AA << "-mer(s)" << std::endl;
        backup << "# length of simulation box L = " << sp->L << std::endl;
        backup << "# sequence " << sp->AA_seq << std::endl;
        backup << "# number of MC steps: " << t+1 << ", with " << sp->stepit << " moves per step" << std::endl;
        backup << "# gammasum = " << gammasum << " (current gamma = " << gamma << ")" << std::endl;
        backup << "# gamma_0 = " << sp->GAMMA_0 << ", T_0 = " << sp->T_0 << std::endl;
        backup << "# current runtime: " << Timer.curRunTime() << std::endl;
        backup << "# energy window: [" << sp->EMin << ";" << sp->EMax << "] in " << sp->NBin << " steps (bin width = " << sp->BinW << ")" << std::endl;
        backup << "# accepted " << ot->naccept[0] << " of " << ot->nattempt[0] << " (" << 100*(double)ot->naccept[0]/(double)ot->nattempt[0] << "%) of local moves" << std::endl;
        backup << "# accepted " << ot->naccept[1] << " of " << ot->nattempt[1] << " (" << 100*(double)ot->naccept[1]/(double)ot->nattempt[1] << "%) of pivot (Phi) moves" << std::endl;
        backup << "# accepted " << ot->naccept[2] << " of " << ot->nattempt[2] << " (" << 100*(double)ot->naccept[2]/(double)ot->nattempt[2] << "%) of pivot (Psi) moves" << std::endl;
        backup << "# accepted " << ot->naccept[3] << " of " << ot->nattempt[3] << " (" << 100*(double)ot->naccept[3]/(double)ot->nattempt[3] << "%) of translation moves" << std::endl;
        backup << "# accepted " << ot->naccept[4] << " of " << ot->nattempt[4] << " (" << 100*(double)ot->naccept[4]/(double)ot->nattempt[4] << "%) of chain rotation moves" << std::endl;
        backup << "bin  from  to  lng  H" << std::endl;
        for( int i=0; i<sp->NBin; i++ ) {
            backup << i << " " << sp->EMin+i*sp->BinW << " " << sp->EMin+(i+1)*sp->BinW << " " << ot->lngE[i] << " " << ot->H[i] << std::endl;
        }
        backup << "# current configuration ( E=" << E << ")" << std::endl;
        backup << "beadID  x  y  z  BCx  BCy  BCz" << std::endl;
        for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
            backup << std::setw(2) << std::setfill('0') << i << "N " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(2) << std::endl;
            backup << std::setw(2) << std::setfill('0') << i << "C " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(2) << std::endl;
            backup << std::setw(2) << std::setfill('0') << i << "O " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(2) << std::endl;
            backup << std::setw(2) << std::setfill('0') << i << "R " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(2) << std::endl;
        }
        backup.close();
        return true;
    }
    else {
        hd->os_log<< std::endl << "--- ERROR ---\tfailed to open " << name << std::endl;
        std::cout << std::endl << "--- ERROR ---\tfailed to open " << name << std::endl;
        return false;
    }
}
bool BackupProdRun(SysPara *sp, Header *hd, Output *ot, Timer &Timer, unsigned long int t)
{
    std::ostringstream oss;
    oss << "results_" << t/sp->T_WRITE << ".dat";
    string name = oss.str();
    ofstream results;
    results.open(name);
    if( results.is_open() ) {
        results << "# Observables for " << sp->N_CH << " PRIME20 " << sp->N_AA << "-mer(s)" << std::endl;
        results << "# length of simulation box L = " << sp->L << std::endl;
        results << "# sequence " << sp->AA_seq << std::endl;
        results << "# number of MC steps: " << t+1 << ", with " << sp->stepit << " moves per step" << std::endl;
        results << "# current runtime: " << Timer.curRunTime() << std::endl;
        results << "# energy window: [" << sp->EMin << ";" << sp->EMax << "] in " << sp->NBin << " steps (bin width = " << sp->BinW << ")" << std::endl;
        results << "# measured observables are:" << std::endl;
        results << "#\tvisit histogram H" << std::endl;
        if(sp->HB_ContMat)   { results << "#\tHB contact matrices (file HBmat.dat)" << std::endl; }
        if(sp->wConfig) { results << "#\tconfiguration snapshots (file .xyz)" << std::endl; }
        results << "Bin from to H" << std::endl;
        for( int i=0; i<sp->NBin; i++ ) {
            results << i << " " << sp->EMin+i*sp->BinW << " " <<sp->EMin+(i+1)*sp->BinW << " " << ot->H[i] << std::endl;
        }
        results.close();
        return true;
    }
    else {
        hd->os_log<< endl << "--- ERROR ---\tfailed to open " << name << std::endl;
        std::cout << endl << "--- ERROR ---\tfailed to open " << name << std::endl;
        return false;
    }
}
// returns true and updates HBList if H-Bond between N(iN) and C(iC) is formed. iN is index of h*N_C+i
bool HBcheck(SysPara *sp, Chain Chn[], int iN, int iC)
{
    double distV[3];                                        // distance vector
    double dist2;                                           // distance squared
    double nh[3], co[3], ho[3], h1[3], h2[3], rot_n[3];     // vector projections for angle calculation and helping vectors
    double nhl, col, hol, adl, nho, coh;                    // absolutes of vector projections
    double cos_a, sin_a, rot_nl;                            // angle/vector for chain end rotation matrix
    double rotMtrx[3][3];                                   // rotation matrix for chain end nh co calculation

    // condition 1:     N and C separated by at least 3 amino acids and within relevant distance
    if( NCDist[iN][iC] < 0 ) { return false; }
    if( (iN/sp->N_AA == iC/sp->N_AA) && (abs(iN-iC) < 4) ) { return false; }

    // condition 2:     distance < SW radius
    if( NCDist[iN][iC] < SW2_BB ) {
        // condition 3:     both partners are available
        if( HBList[iN][0] == -1 && HBList[iC][1] == -1 ) {
            // condition 4:     constraints of distances to neighbours of partner (HB stabilisation)
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0], Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1]);
            if( dotPro(distV, distV) < DIA2_HB_NCa ) {return false; }
            if( (iC+1)%sp->N_AA != 0 ) {
                std::tie(distV[0], distV[1], distV[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0], Chn[iC/sp->N_AA].AmAc[(iC+1)%sp->N_AA].Bd[0]);
                if( dotPro(distV, distV) < DIA2_HB_NN ) { return false; }
            }
            if( iN%sp->N_AA != 0 ) {
                std::tie(distV[0], distV[1], distV[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[(iN-1)%sp->N_AA].Bd[2], Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2]);
                if( dotPro(distV, distV) < DIA2_HB_CC ) { return false; }
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[1], Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2]);
            if( dotPro(distV, distV) < DIA2_HB_CCa ) { return false; }

            // condition 5:     calculate H and O positions and check if angle is within range (120°,180°)
            // H
            if( iN%sp->N_AA == 0 ) {     // NH vector is calculated by rotating NCa in NCaC plane by angle 2.668
                std::tie(h1[0], h1[1], h1[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0], Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[1]);
                std::tie(h2[0], h2[1], h2[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[1], Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[2]);
                cos_a = -0.456; sin_a = 0.89;
                std::tie(rot_n[0], rot_n[1], rot_n[2]) = crossPro( h2, h1);
                rot_nl = absVec(rot_n);
                rot_n[0] /= rot_nl; rot_n[1] /= rot_nl; rot_n[2] /= rot_nl;
                rotMtrx[0][0] = cos_a + rot_n[0]*rot_n[0]*( 1-cos_a );
                rotMtrx[0][1] = rot_n[0]*rot_n[1]*( 1-cos_a ) - rot_n[2]*sin_a;
                rotMtrx[0][2] = rot_n[0]*rot_n[2]*( 1-cos_a ) + rot_n[1]*sin_a;
                rotMtrx[1][0] = rot_n[1]*rot_n[0]*( 1-cos_a ) + rot_n[2]*sin_a;
                rotMtrx[1][1] = cos_a + rot_n[1]*rot_n[1]*( 1-cos_a );
                rotMtrx[1][2] = rot_n[1]*rot_n[2]*( 1-cos_a ) - rot_n[0]*sin_a;
                rotMtrx[2][0] = rot_n[2]*rot_n[0]*( 1-cos_a ) - rot_n[1]*sin_a;
                rotMtrx[2][1] = rot_n[2]*rot_n[1]*( 1-cos_a ) + rot_n[0]*sin_a;
                rotMtrx[2][2] = cos_a + rot_n[2]*rot_n[2]*( 1-cos_a );
                for( int i=0; i<3; i++ ) {
                    nh[i] = h1[0]*rotMtrx[0][i] + h1[1]*rotMtrx[1][i] + h1[2]*rotMtrx[2][i];
                }
                nhl = absVec(nh);
                nh[0] /= nhl; nh[1] /= nhl; nh[2] /= nhl;
            }
            else {
                std::tie(h1[0], h1[1], h1[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[(iN-1)%sp->N_AA].Bd[2], Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[1]);
                std::tie(h2[0], h2[1], h2[2]) = distVecBC(sp, Chn[iN/sp->N_AA].AmAc[(iN-1)%sp->N_AA].Bd[2], Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0]);
                adl = (dotPro(h1, h2)) / (dotPro(h1, h1));
                for( int i=0; i<3; i++ ) {
                    nh[i] = Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0].getR(i) - Chn[iN/sp->N_AA].AmAc[(iN-1)%sp->N_AA].Bd[2].getR(i)-h1[i]*adl - sp->L*round( (Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0].getR(i) - Chn[iN/sp->N_AA].AmAc[(iN-1)%sp->N_AA].Bd[2].getR(i)-h1[i]*adl) / (double)sp->L );
                }
                nhl = sqrt(dotPro(nh, nh));
                nh[0] /= nhl;   nh[1] /= nhl;   nh[2] /= nhl;
            }
            // O
            if( (iC+1)%sp->N_AA == 0 ) {       // CO vector is calculated by rotating CaC in NCaC plane by angle 2.622
                std::tie(h1[0], h1[1], h1[2]) = distVecBC(sp, Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1], Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2]);
                std::tie(h2[0], h2[1], h2[2]) = distVecBC(sp, Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[0], Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1]);
                cos_a = 0.497; sin_a = -0.868;
                std::tie(rot_n[0], rot_n[1], rot_n[2]) = crossPro( h2, h1);
                rot_nl = absVec(rot_n);
                rot_n[0] /= rot_nl; rot_n[1] /= rot_nl; rot_n[2] /= rot_nl;
                rotMtrx[0][0] = cos_a + rot_n[0]*rot_n[0]*( 1-cos_a );
                rotMtrx[0][1] = rot_n[0]*rot_n[1]*( 1-cos_a ) - rot_n[2]*sin_a;
                rotMtrx[0][2] = rot_n[0]*rot_n[2]*( 1-cos_a ) + rot_n[1]*sin_a;
                rotMtrx[1][0] = rot_n[1]*rot_n[0]*( 1-cos_a ) + rot_n[2]*sin_a;
                rotMtrx[1][1] = cos_a + rot_n[1]*rot_n[1]*( 1-cos_a );
                rotMtrx[1][2] = rot_n[1]*rot_n[2]*( 1-cos_a ) - rot_n[0]*sin_a;
                rotMtrx[2][0] = rot_n[2]*rot_n[0]*( 1-cos_a ) - rot_n[1]*sin_a;
                rotMtrx[2][1] = rot_n[2]*rot_n[1]*( 1-cos_a ) + rot_n[0]*sin_a;
                rotMtrx[2][2] = cos_a + rot_n[2]*rot_n[2]*( 1-cos_a );
                for( int i=0; i<3; i++ ) {
                    co[i] = h1[0]*rotMtrx[0][i] + h1[1]*rotMtrx[1][i] + h1[2]*rotMtrx[2][i];
                }
                col = absVec(co);
                co[0] *= ((1.2)/col); co[1] *= ((1.2)/col); co[2] *= ((1.2)/col);
            }
            else {
                std::tie(h1[0], h1[1], h1[2]) = distVecBC(sp, Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1], Chn[iC/sp->N_AA].AmAc[(iC+1)%sp->N_AA].Bd[0]);
                std::tie(h2[0], h2[1], h2[2]) = distVecBC(sp, Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1], Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2]);
                adl = (dotPro(h1, h2)) / (dotPro(h1, h1));
                for( int i=0; i<3; i++ ) {
                    co[i] = Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2].getR(i) - Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1].getR(i)-h1[i]*adl - sp->L*round( (Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2].getR(i) - Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[1].getR(i)-h1[i]*adl) / (double)sp->L );
                }
                col = sqrt(dotPro(co, co));
                co[0] *= (1.2)/col;   co[1] *= (1.2)/col;   co[2] *= (1.2)/col;
            }
            for( int i=0; i<3; i++ ) {
                ho[i] = Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2].getR(i)+co[i] - Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0].getR(i)-nh[i] - sp->L*round( (Chn[iC/sp->N_AA].AmAc[iC%sp->N_AA].Bd[2].getR(i)+co[i] - Chn[iN/sp->N_AA].AmAc[iN%sp->N_AA].Bd[0].getR(i)-nh[i]) / (double)sp->L );
            }
            hol = sqrt(dotPro(ho, ho));
            // cosinus calculation
            nho = -dotPro(nh, ho)/hol;              // NH distance is 1
            coh = dotPro(co, ho)/(1.2*hol);         // CO distance is 1.2

            if( nho <= -0.5 && coh <= -0.5 ) {      // both angles have to be between 120° and 180°
                HBList[iN][0] = iC;
                HBList[iC][1] = iN;
                return true;
            }
        }
    }

    return false;
}
//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXXXXX  ENERGY FUNCTIONS  XXXXXXXXXXXXX

// energy of a single SC interaction between AmAc[i1] and AmAc[i2]. distance d has to be squared
double E_single(Chain Chn[], int h1, int i1, int h2, int i2, double d_sq)
{
    if(h1==h2) {
        if( (h1 == h2) && (abs(i1-i2) < 2) ) return 0;
    }
    if(d_sq <= SWDia(Chn[h1].AmAc[i1], Chn[h2].AmAc[i2]) * SWDia(Chn[h1].AmAc[i1], Chn[h2].AmAc[i2])) {
        return SWDepth(Chn[h1].AmAc[i1], Chn[h2].AmAc[i2]);
    }
    return 0;
}
// SC interaction energy of one SC Bead (j1 must be 3). Or overlapp check for any AmAc[i1].Bd[j1]. Versus chain segment [sp, ep).
double EO_SegBead(SysPara *sypa, Chain Chn[], int h1, int i1, int j1, int sp, int ep, int EOswitch)
{
    int neighBox[3], centrBox[3];
    int searchBox, index;
    int Box, dbx, dby, dbz;
    double distV[3];
    double dist2;
    double energy = 0;

    if( EOswitch == 1 && j1 != 3 ) { std::cout << "--- ERROR ---\tEO_SegBead() energy calculation for j1=" << j1 << endl; return -1; }

    Box = assignBox(sypa, Chn[h1].AmAc[i1].Bd[j1]);
    centrBox[0] = Box % sypa->NBOX;
    centrBox[1] = (Box/sypa->NBOX) % sypa->NBOX;
    centrBox[2] = Box/(double)(sypa->NBOX*sypa->NBOX);
    for( int i=0; i<27; i++ ) {
        dbx = i%3; dby = (i/3)%3; dbz = i/9;
        neighBox[0] = (centrBox[0] + dbx-1 + sypa->NBOX) % sypa->NBOX;
        neighBox[1] = (centrBox[1] + dby-1 + sypa->NBOX) % sypa->NBOX;
        neighBox[2] = (centrBox[2] + dbz-1 + sypa->NBOX) % sypa->NBOX;
        searchBox = neighBox[0] + neighBox[1]*sypa->NBOX + neighBox[2]*sypa->NBOX*sypa->NBOX;
        index = neighHead[searchBox];
        while( index != -1 ) {
            if( EOswitch == 0 ) {       // overlapp check
                if( index/4 == h1*sypa->N_AA+i1 ) {
                    index = neighList[index];
                    continue;
                }
                if( index >= sp*4 && index < ep*4 ) {
                    std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[h1].AmAc[i1].Bd[j1], Chn[index/(4*sypa->N_AA)].AmAc[(index/4)%sypa->N_AA].Bd[index%4]);
                    dist2 = dotPro(distV, distV);
                    if( dist2 < DiaSQ(Chn, h1, i1, j1, index/(4*sypa->N_AA), (index/4)%sypa->N_AA, index%4) ) { 
                        //std::cerr << "overlapp: C" << h1 << "A" << i1 << "B" << j1 << " - C" << index/(4*sypa->N_AA) << "A" << (index/4)%sypa->N_AA << "B" << index%4 << std::endl;
                        return -1.0; 
                    }
                }
            }
            else {                      // energy calculation
                if( index%4 == 3 ) {
                    if( index >= sp*4 && index < ep*4 ) {
                        std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[h1].AmAc[i1].Bd[j1], Chn[index/(4*sypa->N_AA)].AmAc[(index/4)%sypa->N_AA].Bd[index%4]);
                        dist2 = dotPro(distV, distV);
                        energy += E_single(Chn, h1, i1, index/(4*sypa->N_AA), (index/4)%sypa->N_AA, dist2);

                        /*if( E_single(Chn, h1, i1, index/(4*N_AA), (index/4)%N_AA, dist2) != 0.0 ) {
                            std::cerr << endl << "energy contribution " << h1*N_AA+i1 << "  and " << index/4 << "   E = " << E_single(Chn, h1, i1, index/(4*N_AA), (index/4)%N_AA, dist2) << endl;
                        }*/

                    }
                }
            }
            index = neighList[index];
        }
    }
    return energy;
}
// SC interaction energy of segment [sp1,ep1) versus segment [sp2,ep2). also overlapp check
double EO_SegSeg(SysPara *sp, Chain Chn[], int sp1, int ep1, int sp2, int ep2, int EOswitch) 
{
    double dEnergy;
    double energy = 0;

    for( int i=sp1; i<ep1; i++ ) {
        for( int j=0; j<4; j++ ) {
            if( EOswitch == 1 && j != 3) continue;  // skip energy calculation for non-SC beads
            dEnergy = EO_SegBead(sp, Chn, i/sp->N_AA, i%sp->N_AA, j, sp2, ep2, EOswitch);
            if( dEnergy == -1 ) return -1;          // exit if overlapp check fails
            energy += dEnergy;
        }
    }

    return energy;
}
// calculates total energy from scratch
double E_check(SysPara *sp, Chain Chn[])
{
    int HBL_check[sp->N_CH*sp->N_AA][2], HBerrN[sp->N_CH*sp->N_AA][2], HBerrC[sp->N_CH*sp->N_AA][2];
    double NCD_check[sp->N_CH*sp->N_AA][sp->N_CH*sp->N_AA];
    double Energy, distabs;
    double distV[3];

    // create new HBList → HBL_check
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
        HBL_check[i][0] = HBList[i][0];
        HBL_check[i][1] = HBList[i][1];
        HBList[i][0] = -1;
        HBList[i][1] = -1;
        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
            if( (i/sp->N_AA == j/sp->N_AA) && (abs(i-j) < 4) ) {
                NCD_check[i][j] = -1.0;
                continue;
            }
            NCD_check[i][j] = NCDist[i][j];
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sp, Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0], Chn[j/sp->N_AA].AmAc[j%sp->N_AA].Bd[2]);
            distabs = dotPro(distV, distV);
            if( distabs < SW2HUGE ) {
                NCDist[i][j] = distabs;
            }
            else {
                NCDist[i][j] = -1.0;
            }
        }
    }
    // energy calculation
    Energy = 0.0;
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
            if( HBcheck(sp, Chn, i, j) ) { 
                Energy -= 1.0; 
            }
        }
        Energy += EO_SegBead(sp, Chn, i/sp->N_AA, i%sp->N_AA, 3, i, sp->N_CH*sp->N_AA, 1);
    }
    // disagreements between HBList and HBL_check
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
        if( HBL_check[i][0] != HBList[i][0] ) { 
            HBerrN[i][0] = HBList[i][0];
            HBerrN[i][1] = HBL_check[i][0];
        }
        else { HBerrN[i][0] = -1; HBerrN[i][1] = -1; }
        if( HBL_check[i][1] != HBList[i][1] ) {
            HBerrC[i][0] = HBList[i][1];
            HBerrC[i][1] = HBL_check[i][1];
        }
        else { HBerrC[i][0] = -1; HBerrC[i][1] = -1; }
    }
    // if multiple HB partners are possible, this does not affect Energy but is still listet as error
    // this should not lead to an error message, thus this is corrected below
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
        if( HBerrN[i][0] > -1 && HBerrN[i][1] > -1 ) {
            if( HBerrC[HBerrN[i][0]][1] == -1 && HBerrC[HBerrN[i][1]][0] == -1 ) {
                HBerrC[HBerrN[i][0]][0] = -1;
                HBerrC[HBerrN[i][1]][1] = -1;
                HBerrN[i][0] = -1;
                HBerrN[i][1] = -1;
            }
        }
        if( HBerrC[i][0] > -1 && HBerrC[i][1] > -1 ) {
            if( HBerrN[HBerrC[i][0]][1] == -1 && HBerrN[HBerrC[i][1]][0] == -1 ) {
                HBerrN[HBerrC[i][0]][0] = -1;
                HBerrN[HBerrC[i][1]][1] = -1;
                HBerrC[i][0] = -1;
                HBerrC[i][1] = -1;
            }
        }
    }
    // error messages and HBList reset
    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
        if( HBerrN[i][0] > -1 || HBerrN[i][1] > -1 ) {
            std::cout << "\n\tHB mismatch detected → N: check " << i << "-" << HBerrN[i][0] << " and running " << i << "-" << HBerrN[i][1] << endl;
        }
        if( HBerrC[i][0] > -1 || HBerrC[i][1] > -1 ) {
            std::cout << "\n\tHB mismatch detected → C: check " << i << "-" << HBerrC[i][0] << " and running " << i << "-" << HBerrC[i][1] << endl;
        }
        HBList[i][0] = HBL_check[i][0];
        HBList[i][1] = HBL_check[i][1];
        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
            if( abs(NCD_check[i][j] - NCDist[i][j]) > 0.01 ) {
                std::cout << endl << "\tN-" << i << " C-" << j << " dist² mismatch:  newly calculated d²=" << NCDist[i][j] << " running value d²=" << NCD_check[i][j] << endl;
            }
            NCDist[i][j] = NCD_check[i][j];
        }
    }

    return Energy;
}
// SAMC acceptance function
bool acceptance(double lngEold, double lngEnew)
{
    //uniform_real_distribution<double> distribution(0,1);
    //double MCrand = distribution(rng);
    double MCrand = (double)rng()/(double)rng.max();
    if( lngEnew <= lngEold ) return true;
    if( exp(lngEold-lngEnew) > MCrand ) return true;
    return false;
}
//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXXXXXX  MOVE FUNCTIONS  XXXXXXXXXXXXXX

// small displacement of Chn[h].AmAc[i].Bd[j]
bool wiggle(SysPara *sp, Chain Chn[], int h, int i, int j, double &deltaE)
{
    Bead cpy;                       // copy of bead to be wiggled, wiggled bead
    double disp[3], dist[3];        // random displacement vector, distance vector
    double newx, newy, newz;        // new coordinates for wiggled bead
    double Eold, distabs;
    int HBN_broken, HBC_broken, HBEnd_broken;

    // calculate old energy contribution of SC if j==3
    if( j == 3 ) {
        Eold = EO_SegBead(sp, Chn, h, i, j, 0, sp->N_AA*sp->N_CH, 1);
    }
    
    // calculate displacement vector (disp[]) and move bead
    cpy = Chn[h].AmAc[i].Bd[j];
    for(int k = 0; k < 3; k++) {
        disp[k] = ( ((double)rng()/(double)rng.max())*2 - 1. ) * sp->DISP_MAX;
    }
    newx = Chn[h].AmAc[i].Bd[j].getR(0) + disp[0];  Chn[h].AmAc[i].Bd[j].addBC(0, floor(newx/sp->L));   newx = newx - sp->L*floor(newx/sp->L);
    newy = Chn[h].AmAc[i].Bd[j].getR(1) + disp[1];  Chn[h].AmAc[i].Bd[j].addBC(1, floor(newy/sp->L));   newy = newy - sp->L*floor(newy/sp->L);
    newz = Chn[h].AmAc[i].Bd[j].getR(2) + disp[2];  Chn[h].AmAc[i].Bd[j].addBC(2, floor(newz/sp->L));   newz = newz - sp->L*floor(newz/sp->L);
    Chn[h].AmAc[i].Bd[j].setR( newx, newy, newz );
    // check all bond lengths: if outside allowed range => reject move.
    switch(j) {
        case 0:
            if (i>0) {
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i-1].Bd[1]);
                if( abs(absVec(dist)/PBND_CaN - 1) > BND_FLUCT )        { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i-1].Bd[2]);
                if( abs(absVec(dist)/BND_CN - 1) > BND_FLUCT )          { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[1]);
            if( abs(absVec(dist)/BND_NCa - 1) > BND_FLUCT )             { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[2]);
            if( abs(absVec(dist)/PBND_NC - 1) > BND_FLUCT )             { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[3]);
            if( abs(absVec(dist)/Chn[h].AmAc[i].getDisR(0) - 1) > BND_FLUCT )   { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            break;
        case 1:
            if (i>0) {
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i-1].Bd[1]);
                if( abs(absVec(dist)/PBND_CaCa - 1) > BND_FLUCT )       { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i-1].Bd[2]);
                if( abs(absVec(dist)/PBND_CCa - 1) > BND_FLUCT )        { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[0]);
            if( abs(absVec(dist)/BND_NCa - 1) > BND_FLUCT )             { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[2]);
            if( abs(absVec(dist)/BND_CaC - 1) > BND_FLUCT )             { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[3]);
            if( abs(absVec(dist)/Chn[h].AmAc[i].getDisR(1) - 1) > BND_FLUCT )   { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            if (i<sp->N_AA-1) {
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i+1].Bd[0]);
                if( abs(absVec(dist)/PBND_CaN - 1) > BND_FLUCT )        { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i+1].Bd[1]);
                if( abs(absVec(dist)/PBND_CaCa - 1) > BND_FLUCT )       { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            }
            break;
        case 2:
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[0]);
            if( abs(absVec(dist)/PBND_NC - 1) > BND_FLUCT )             { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[1]);
            if( abs(absVec(dist)/BND_CaC - 1) > BND_FLUCT )             { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[3]);
            if( abs(absVec(dist)/Chn[h].AmAc[i].getDisR(2) - 1) > BND_FLUCT )   { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            if (i<sp->N_AA-1) {
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i+1].Bd[0]);
                if( abs(absVec(dist)/BND_CN - 1) > BND_FLUCT )          { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
                std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i+1].Bd[1]);
                if( abs(absVec(dist)/PBND_CCa - 1) > BND_FLUCT )        { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            }
            break;
        case 3:
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[0]);
            if( abs(absVec(dist)/Chn[h].AmAc[i].getDisR(0) - 1) > BND_FLUCT )   { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[1]);
            if( abs(absVec(dist)/Chn[h].AmAc[i].getDisR(1) - 1) > BND_FLUCT )   { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[j], Chn[h].AmAc[i].Bd[2]);
            if( abs(absVec(dist)/Chn[h].AmAc[i].getDisR(2) - 1) > BND_FLUCT )   { Chn[h].AmAc[i].Bd[j] = cpy; return false; }
            break;
    }
    // perform overlap check
    if(EO_SegBead(sp, Chn, h, i, j, 0, sp->N_AA*sp->N_CH, 0) == -1) {
        Chn[h].AmAc[i].Bd[j] = cpy;
        return false;
    }
    // Break all H-Bonds affected by the move
    HBN_broken = -1; HBC_broken = -1; HBEnd_broken = -1;
    switch(j) {
        case 0:
            if( HBList[h*sp->N_AA+i][0] > -1 ) {
                HBN_broken = HBList[h*sp->N_AA+i][0];
                HBList[HBList[h*sp->N_AA+i][0]][1] = -1;
                HBList[h*sp->N_AA+i][0] = -1;
                deltaE += 1.0;
            }
            if( i%sp->N_AA != 0 ) {
                if( HBList[h*sp->N_AA+i-1][1] > -1 ) {
                    HBC_broken = HBList[h*sp->N_AA+i-1][1];
                    HBList[HBList[h*sp->N_AA+i-1][1]][0] = -1;
                    HBList[h*sp->N_AA+i-1][1] = -1;
                    deltaE += 1.0;
                }
            }
            if( (i+1)%sp->N_AA == 0 ) {
                if( HBList[h*sp->N_AA+i][1] > -1) {
                    HBEnd_broken = HBList[h*sp->N_AA+i][1];
                    HBList[HBList[h*sp->N_AA+i][1]][0] = -1;
                    HBList[h*sp->N_AA+i][1] = -1;
                    deltaE += 1.0;
                }
            }
            break;
        case 1:
            if( HBList[h*sp->N_AA+i][0] > -1 ) {
                HBN_broken = HBList[h*sp->N_AA+i][0];
                HBList[HBList[h*sp->N_AA+i][0]][1] = -1;
                HBList[h*sp->N_AA+i][0] = -1;
                deltaE += 1.0;
            }
            if( HBList[h*sp->N_AA+i][1] > -1 ) {
                HBC_broken = HBList[h*sp->N_AA+i][1];
                HBList[HBList[h*sp->N_AA+i][1]][0] = -1;
                HBList[h*sp->N_AA+i][1] = -1;
                deltaE += 1.0;
            }
            break;
        case 2:
            if( HBList[h*sp->N_AA+i][1] > -1 ) {
                HBC_broken = HBList[h*sp->N_AA+i][1];
                HBList[HBList[h*sp->N_AA+i][1]][0] = -1;
                HBList[h*sp->N_AA+i][1] = -1;
                deltaE += 1.0;
            }
            if( (i+1)%sp->N_AA != 0 ) {
                if( HBList[h*sp->N_AA+i+1][0] > -1 ) {
                    HBN_broken = HBList[h*sp->N_AA+i+1][0];
                    HBList[HBList[h*sp->N_AA+i+1][0]][1] = -1;
                    HBList[h*sp->N_AA+i+1][0] = -1;
                    deltaE += 1.0;
                }
            }
            if( i%sp->N_AA == 0 ) {
                if( HBList[h*sp->N_AA+i][0] > -1 ) {
                    HBEnd_broken = HBList[h*sp->N_AA+i][0];
                    HBList[HBList[h*sp->N_AA+i][0]][1] = -1;
                    HBList[h*sp->N_AA+i][0] = -1;
                    deltaE += 1.0;
                }
            }
    }
    // energy calculation
    switch(j) {
        case 0:
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                if( (h != k/sp->N_AA) || ((h == k/sp->N_AA) && (abs(h*sp->N_AA+i-k) > 3)) ) {
                    std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[0], Chn[k/sp->N_AA].AmAc[k%sp->N_AA].Bd[2]);
                    distabs = dotPro(dist, dist);
                    if( distabs < SW2HUGE )    { NCDist[h*sp->N_AA+i][k] = distabs; }
                    else                        { NCDist[h*sp->N_AA+i][k] = -1; }
                }
                if( HBcheck(sp, Chn, h*sp->N_AA+i, k) ) deltaE -= 1.0;
                if(i>0) { if( HBcheck(sp, Chn, k, h*sp->N_AA+i-1) ) deltaE -= 1.0; }
                if((i+1)%sp->N_AA == 0) { if(HBcheck(sp, Chn, k, h*sp->N_AA+i)) deltaE -= 1.0; }
            }
            break;
        case 1:
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                if( HBcheck(sp, Chn, h*sp->N_AA+i, k) ) deltaE -= 1.0;
                if( HBcheck(sp, Chn, k, h*sp->N_AA+i) ) deltaE -= 1.0;
            }
            break;
        case 2:
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                if( (h != k/sp->N_AA) || ((h == k/sp->N_AA) && (abs(h*sp->N_AA+i-k) > 3)) ) {
                    std::tie(dist[0], dist[1], dist[2]) = distVecBC(sp, Chn[h].AmAc[i].Bd[2], Chn[k/sp->N_AA].AmAc[k%sp->N_AA].Bd[0]);
                    distabs = dotPro(dist, dist);
                    if( distabs < SW2HUGE )    { NCDist[k][h*sp->N_AA+i] = distabs; }
                    else                        { NCDist[k][h*sp->N_AA+i] = -1; }
                }
                if(i==0) { if(HBcheck(sp, Chn, h*sp->N_AA+i, k)) deltaE -= 1.0; }
                if( HBcheck(sp, Chn, k, h*sp->N_AA+i) ) deltaE -= 1.0;
                if(i<sp->N_AA-1) { if( HBcheck(sp, Chn, h*sp->N_AA+i+1, k) ) deltaE -= 1.0; }
            }
            break;
        case 3:
            deltaE += EO_SegBead(sp, Chn, h, i, j, 0, sp->N_CH*sp->N_AA, 1) - Eold;
    }
    // previously broken HB can build new HB
    switch(j) {
        case 0:
            if( HBList[h*sp->N_AA+i][0] == HBN_broken ) HBN_broken = -1;
            if( i == 0 ) break;
            if( HBList[h*sp->N_AA+i-1][1] == HBC_broken ) HBC_broken = -1;
            if( i == sp->N_AA-1 ) { if( HBList[h*sp->N_AA+i][1] == HBEnd_broken ) HBEnd_broken = -1; }
            break;
        case 1:
            if( HBList[h*sp->N_AA+i][0] == HBN_broken ) HBN_broken = -1;
            if( HBList[h*sp->N_AA+i][1] == HBC_broken ) HBC_broken = -1;
            break;
        case 2:
            if( HBList[h*sp->N_AA+i][1] == HBC_broken ) HBC_broken = -1;
            if( i == sp->N_AA-1 ) break;
            if( HBList[h*sp->N_AA+i+1][0] == HBN_broken ) HBN_broken = -1;
            if( i == 0 ) { if( HBList[h*sp->N_AA+i][0] == HBEnd_broken ) HBEnd_broken = -1; }
    }
    if( j != 3 ) {
        if( HBN_broken > -1 ) {
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                if( HBcheck(sp, Chn, k, HBN_broken) ) deltaE -= 1.0;
            }
        }
        if( HBC_broken > -1 ) {
            for( int k=0; k<sp->N_AA; k++ ) {
                if( HBcheck(sp, Chn, HBC_broken, k) ) deltaE -= 1.0;
            }
        }
    }
    // HBEnd_broken
    if( i==0 && j==2 ) {
        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
            if( HBcheck(sp, Chn, 0, k) ) deltaE -= 1.0;
        }
        if( HBEnd_broken > -1 ) {
            if( HBEnd_broken != HBList[h*sp->N_CH][0] ) {
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    if( HBcheck(sp, Chn, k, HBEnd_broken) ) deltaE -= 1.0;
                }
            }
        }
    }
    if( i==sp->N_AA-1 && j==0 ) {
        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
            if( HBcheck(sp, Chn, k, sp->N_AA-1) ) deltaE -= 1.0;
        }
        if( HBEnd_broken > -1 ) {
            if( HBEnd_broken != HBList[h*sp->N_CH+sp->N_AA-1][1] ) {
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    if( HBcheck(sp, Chn, HBEnd_broken, k) ) deltaE -= 1.0;
                }
            }
        }
    }

    return true;
}
// rotation around N-Ca axis of i1-th amino acid. high=0 rotates lower part, high=1 rotates upper part
bool rotPhi(SysPara *sypa, Chain Chn[], int i1, int high, double &deltaE)
{
    Bead Bdcpy[4*sypa->N_AA*sypa->N_CH];                                // copy of rotated beads
    double angle, cos_a, sin_a, nsqrt;                      // variables for rotation matrix
    double Eold, dEhb, distabs;                             // assisting energy values for dE calculation
    double rotMtrx[3][3];                                   // rotation matrix
    double n[3], dVec[3], newpos[3], com[3], shift[3];
    int BrokenHB[sypa->N_AA*sypa->N_CH][2], nBHB;                       // list of broken HB
    int sp, ep, spHB, epHB;                                 // start and end of rotated chain segment

    // calculate Eold
    sp = (high == 0) ? ((i1/sypa->N_AA)*sypa->N_AA):i1;                 // identify start of rotated chain segment
    ep = (high == 0) ? (i1):(((i1/sypa->N_AA)+1)*sypa->N_AA);           // identify end of rotated chain segment
    Eold = EO_SegSeg(sypa, Chn, sp, ep, 0, sp, 1) + EO_SegSeg(sypa, Chn, sp, ep, ep, sypa->N_CH*sypa->N_AA, 1);
    // copy relevant part of the chain
    for( int i=sp; i<ep; i++ ) {
        for( int j=0; j<4; j++ ) {
            Bdcpy[i*4+j] = Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j];
        }
    }

    angle = sypa->DPHI_MAX*( ((double)rng()/(double)rng.max())*2 - 1. );
    cos_a = cos(angle); sin_a = sin(angle);
    std::tie(n[0], n[1], n[2]) = distVecBC(sypa, Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[0], Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1]);  nsqrt = absVec(n);
    n[0] /= nsqrt;  n[1] /= nsqrt;  n[2] /= nsqrt;
    if( high==0 ) { n[0] *= -1; n[1] *= -1; n[2] *= -1; }     // -n for low? is this nessessary to keep rotation angle distribution uniform?
    rotMtrx[0][0] = cos_a + n[0]*n[0]*(1-cos_a);
    rotMtrx[0][1] = n[0]*n[1]*(1-cos_a) - n[2]*sin_a;
    rotMtrx[0][2] = n[0]*n[2]*(1-cos_a) + n[1]*sin_a;
    rotMtrx[1][0] = n[1]*n[0]*(1-cos_a) + n[2]*sin_a;
    rotMtrx[1][1] = cos_a + n[1]*n[1]*(1-cos_a);
    rotMtrx[1][2] = n[1]*n[2]*(1-cos_a) - n[0]*sin_a;
    rotMtrx[2][0] = n[2]*n[0]*(1-cos_a) - n[1]*sin_a;
    rotMtrx[2][1] = n[2]*n[1]*(1-cos_a) + n[0]*sin_a;
    rotMtrx[2][2] = cos_a + n[2]*n[2]*(1-cos_a);

    if( high == 0 ) {
        for( int i=(i1/sypa->N_AA)*sypa->N_AA; i<i1; i++ ) {
            for( int j=0; j<4; j++ ) {
                for( int k=0; k<3; k++ ) {
                    dVec[k] = Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getR(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[0].getR(k) + sypa->L*(Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getBC(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[0].getBC(k) );
                }
                for( int k=0; k<3; k++ ) {
                    // periodic boundary conditions here and in high.
                    newpos[k] = dVec[0]*rotMtrx[0][k] + dVec[1]*rotMtrx[1][k] + dVec[2]*rotMtrx[2][k] + Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[0].getR(k);
                    Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setBC( k, Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[0].getBC(k) + floor(newpos[k]/sypa->L) );
                    newpos[k] = newpos[k] - sypa->L*floor(newpos[k]/sypa->L);
                }
                Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setR(newpos[0], newpos[1], newpos[2]);
                
                // overlapp check
                if( EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, 0, (i1/sypa->N_AA)*sypa->N_AA, 0) == -1 || EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, i1, sypa->N_CH*sypa->N_AA, 0) == -1 ) {
                    for( int k=(i/sypa->N_AA)*sypa->N_AA; k<i+1; k++ ) {
                        for( int m=0; m<4; m++ ) {
                            Chn[k/sypa->N_AA].AmAc[k%sypa->N_AA].Bd[m] = Bdcpy[k*4+m];        // reset chain & exit
                        }
                    }
                    return false;
                }
            }
        }
    }
    else if( high == 1 ) {
        for( int i=i1; i<((i1/sypa->N_AA)+1)*sypa->N_AA; i++ ) {
            for( int j=0; j<4; j++ ) {
                for( int k=0; k<3; k++ ) {
                    dVec[k] = Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getR(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getR(k) + sypa->L*( Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getBC(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getBC(k) );
                }
                for( int k=0; k<3; k++ ) {
                    newpos[k] = dVec[0]*rotMtrx[0][k] + dVec[1]*rotMtrx[1][k] + dVec[2]*rotMtrx[2][k] + Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getR(k);
                    Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setBC( k, Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getBC(k) + floor(newpos[k]/sypa->L) );
                    newpos[k] = newpos[k] - sypa->L*floor(newpos[k]/sypa->L);
                }
                Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setR(newpos[0], newpos[1], newpos[2]);
                
                // overlapp check
                if( EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, 0, i1, 0) == -1 || EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, ((i1/sypa->N_AA)+1)*sypa->N_AA, sypa->N_CH*sypa->N_AA, 0) == -1 ) {
                    for( int k=i1; k<i+1; k++ ) {
                        for( int m=0; m<4; m++ ) {
                            Chn[k/sypa->N_AA].AmAc[k%sypa->N_AA].Bd[m] = Bdcpy[k*4+m];        // reset chain & exit
                        }
                    }
                    return false;
                }
            }
        }
    }
    else {
        printf("--- ERROR ---\trotPhi(): int high has unacceptable value!\n");
    }
    // energy calculation
    // break old HB
    nBHB = -1;
    dEhb = 0.0;
    spHB = (high == 0) ? sp : (sp+1);
    epHB = (high == 0) ? (ep+1) : ep;
    for( int i=sp; i<epHB; i++ ) {
        if( (HBList[i][0] < spHB || HBList[i][0] >= ep) && HBList[i][0] > -1 ) {
            BrokenHB[++nBHB][0] = i;
            BrokenHB[nBHB][1] = HBList[i][0];
            HBList[HBList[i][0]][1] = -1;
            HBList[i][0] = -1;
            dEhb += 1.0;
        }
        if( (HBList[i][1] < spHB || HBList[i][1] >= ep) && HBList[i][1] > -1 ) {
            BrokenHB[++nBHB][1] = i;
            BrokenHB[nBHB][0] = HBList[i][1];
            HBList[HBList[i][1]][0] = -1;
            HBList[i][1] = -1;
            dEhb += 1.0;
        }
        // ATTENTION: maybe extra conditions have to be checked like HBList[i1][0] or HBList[i1-1][1] with all possible partners
    }
    // close new HB
    for( int m=sp; m<epHB; m++ ) {
        for( int n=0; n<sypa->N_CH*sypa->N_AA; n++ ) {
            if( n >= spHB && n < ep ) continue;
            if( (m/sypa->N_AA != n/sypa->N_AA) || (( m/sypa->N_AA == n/sypa->N_AA) && (abs(m-n) > 3)) ) {
                std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sypa, Chn[m/sypa->N_AA].AmAc[m%sypa->N_AA].Bd[0], Chn[n/sypa->N_AA].AmAc[n%sypa->N_AA].Bd[2]);
                distabs = dotPro(dVec, dVec);
                if( distabs < SW2HUGE )     { NCDist[m][n] = distabs; }
                else                        { NCDist[m][n] = -1; }
                std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sypa, Chn[m/sypa->N_AA].AmAc[m%sypa->N_AA].Bd[2], Chn[n/sypa->N_AA].AmAc[n%sypa->N_AA].Bd[0]);
                distabs = dotPro(dVec, dVec);
                if( distabs < SW2HUGE )     { NCDist[n][m] = distabs; }
                else                        { NCDist[n][m] = -1; }
            }
            else {
                NCDist[m][n] = -1;
                NCDist[n][m] = -1;
            }
            if(HBcheck(sypa, Chn, m, n)) dEhb -= 1.0;
            if(HBcheck(sypa, Chn, n, m)) dEhb -= 1.0;
        }
    }
    // ATTENTION: extra conditions and situations would have to be applied here as well
    if( nBHB >= 0 ) {       // previously broken HB can rebond
        for( int i=0; i<nBHB+1; i++ ) {
            for( int j=0; j<sypa->N_CH*sypa->N_AA; j++ ) {
                if(HBcheck(sypa, Chn, BrokenHB[i][0], j)) dEhb -= 1.0;
                if(HBcheck(sypa, Chn, j, BrokenHB[i][1])) dEhb -= 1.0;
            }
        }
    }
    for( int i=sp; i<ep; i++ ) {
        for( int j=0; j<4; j++ ) {
            LinkListUpdate(sypa, Chn, i, j);
        }
    }
    // SC interactions
    deltaE = EO_SegSeg(sypa, Chn, sp, ep, 0, sp, 1) + EO_SegSeg(sypa, Chn, sp, ep, ep, sypa->N_CH*sypa->N_AA, 1) + dEhb - Eold;

    return true;
}
// rotation around Ca-C axis of i1-th amino acid. high=0 rotates lower part, high=1 rotates upper part
bool rotPsi(SysPara *sypa, Chain Chn[], int i1, int high, double &deltaE)
{
    Bead Bdcpy[4*sypa->N_AA*sypa->N_CH];                                // copy of rotated beads
    double angle, cos_a, sin_a, nsqrt;                      // variables for rotation matrix
    double Eold, dEhb, distabs;                             // assisting energy values for dE calculation
    double rotMtrx[3][3];                                   // rotation matrix
    double n[3], dVec[3], newpos[3], com[3], shift[3];
    int BrokenHB[sypa->N_AA*sypa->N_CH][2], nBHB;                       // list of broken HB
    int sp, ep, spHB, epHB;                                 // start and end of rotated chain segment

    // calculate Eold
    sp = (high == 0) ? ((i1/sypa->N_AA)*sypa->N_AA):(i1+1);             // identify start of rotated chain segment
    ep = (high == 0) ? (i1+1):(((i1/sypa->N_AA)+1)*sypa->N_AA);         // identify end of rotated chain segment
    Eold = EO_SegSeg(sypa, Chn, sp, ep, 0, sp, 1) + EO_SegSeg(sypa, Chn, sp, ep, ep, sypa->N_CH*sypa->N_AA, 1);
    // copy relevant part of the chain
    for( int i=sp; i<ep; i++ ) {
        for( int j=0; j<4; j++ ) {
            Bdcpy[i*4+j] = Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j];
        }
    }

    angle = sypa->DPSI_MAX*( ((double)rng()/(double)rng.max())*2 - 1. );
    cos_a = cos(angle); sin_a = sin(angle);
    std::tie(n[0], n[1], n[2]) = distVecBC(sypa, Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1], Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[2]);  nsqrt = absVec(n);
    n[0] /= nsqrt;  n[1] /= nsqrt;  n[2] /= nsqrt;
    if( high==0 ) { n[0] = -n[0]; n[1] = -n[1]; n[2] = -n[2]; }     // -n for low? is this nessessary to keep rotation angle distribution uniform?
    rotMtrx[0][0] = cos_a + n[0]*n[0]*(1-cos_a);
    rotMtrx[0][1] = n[0]*n[1]*(1-cos_a) - n[2]*sin_a;
    rotMtrx[0][2] = n[0]*n[2]*(1-cos_a) + n[1]*sin_a;
    rotMtrx[1][0] = n[1]*n[0]*(1-cos_a) + n[2]*sin_a;
    rotMtrx[1][1] = cos_a + n[1]*n[1]*(1-cos_a);
    rotMtrx[1][2] = n[1]*n[2]*(1-cos_a) - n[0]*sin_a;
    rotMtrx[2][0] = n[2]*n[0]*(1-cos_a) - n[1]*sin_a;
    rotMtrx[2][1] = n[2]*n[1]*(1-cos_a) + n[0]*sin_a;
    rotMtrx[2][2] = cos_a + n[2]*n[2]*(1-cos_a);

    if( high == 0 ) {
        for( int i=(i1/sypa->N_AA)*sypa->N_AA; i<i1+1; i++ ) {
            for( int j=0; j<4; j++ ) {
                for( int k=0; k<3; k++ ) {
                    dVec[k] = Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getR(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getR(k) + sypa->L*(Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getBC(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getBC(k) );
                }
                for( int k=0; k<3; k++ ) {
                    // periodic boundary conditions here and in high.
                    newpos[k] = dVec[0]*rotMtrx[0][k] + dVec[1]*rotMtrx[1][k] + dVec[2]*rotMtrx[2][k] + Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getR(k);
                    Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setBC( k, Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[1].getBC(k) + floor(newpos[k]/sypa->L) );
                    newpos[k] = newpos[k] - sypa->L*floor(newpos[k]/sypa->L);
                }
                Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setR(newpos[0], newpos[1], newpos[2]);
                
                // overlapp check
                if( EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, 0, (i1/sypa->N_AA)*sypa->N_AA, 0) == -1 || EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, i1, sypa->N_CH*sypa->N_AA, 0) == -1 ) {
                    for( int k=(i/sypa->N_AA)*sypa->N_AA; k<i+1; k++ ) {
                        for( int m=0; m<4; m++ ) {
                            Chn[k/sypa->N_AA].AmAc[k%sypa->N_AA].Bd[m] = Bdcpy[k*4+m];        // reset chain & exit
                        }
                    }
                    return false;
                }
            }
        }
    }
    else if( high == 1 ) {
        for( int i=i1+1; i<((i1/sypa->N_AA)+1)*sypa->N_AA; i++ ) {
            for( int j=0; j<4; j++ ) {
                for( int k=0; k<3; k++ ) {
                    dVec[k] = Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getR(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[2].getR(k) + sypa->L*(Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].getBC(k) - Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[2].getBC(k) );
                }
                for( int k=0; k<3; k++ ) {
                    newpos[k] = dVec[0]*rotMtrx[0][k] + dVec[1]*rotMtrx[1][k] + dVec[2]*rotMtrx[2][k] + Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[2].getR(k);
                    Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setBC( k, Chn[i1/sypa->N_AA].AmAc[i1%sypa->N_AA].Bd[2].getBC(k) + floor(newpos[k]/sypa->L) );
                    newpos[k] = newpos[k] - sypa->L*floor(newpos[k]/sypa->L);
                }
                Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[j].setR(newpos[0], newpos[1], newpos[2]);
                
                // overlapp check
                if( EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, 0, i1+1, 0) == -1 || EO_SegBead(sypa, Chn, i/sypa->N_AA, i%sypa->N_AA, j, ((i/sypa->N_AA)+1)*sypa->N_AA, sypa->N_CH*sypa->N_AA, 0 ) ) {
                    for( int k=i1+1; k<i+1; k++ ) {
                        for( int m=0; m<4; m++ ) {
                            Chn[k/sypa->N_AA].AmAc[k%sypa->N_AA].Bd[m] = Bdcpy[k*4+m];        // reset chain & exit
                        }
                    }
                    return false;
                }
            }
        }
    }
    else {
        printf("--- ERROR ---\trotPsi(): int high has unacceptable value!\n");
    }
    // energy calculation
    // break old HB
    nBHB = -1;          // no. broken HB
    dEhb = 0.0;         // change in energy due to HB
    spHB = (high == 0) ? sp : (sp-1);
    epHB = (high == 0) ? (ep-1) : ep;
    for( int i=spHB; i<ep; i++ ) {
        if( (HBList[i][0] <= sp || HBList[i][0] >= epHB) && HBList[i][0] > -1 ) {
            BrokenHB[++nBHB][0] = i;
            BrokenHB[nBHB][1] = HBList[i][0];
            HBList[HBList[i][0]][1] = -1;
            HBList[i][0] = -1;
            dEhb += 1.0;
        }
        if( (HBList[i][1] <= sp || HBList[i][1] >= ep) && HBList[i][1] > -1 ) {
            BrokenHB[++nBHB][1] = i;
            BrokenHB[nBHB][0] = HBList[i][1];
            HBList[HBList[i][1]][0] = -1;
            HBList[i][1] = -1;
            dEhb += 1.0;
        }
        // ATTENTION: maybe extra conditions have to be checked like HBList[i1][0] or HBList[i1-1][1] with all possible partners
    }
    // close new HB
    for( int m=spHB; m<ep; m++ ) {
        for( int n=0; n<sypa->N_CH*sypa->N_AA; n++ ) {
            if( n >= sp && n < epHB) continue;
            if( (m/sypa->N_AA != n/sypa->N_AA) || (( m/sypa->N_AA == n/sypa->N_AA) && (abs(m-n) > 3)) ) {
                std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sypa, Chn[m/sypa->N_AA].AmAc[m%sypa->N_AA].Bd[0], Chn[n/sypa->N_AA].AmAc[n%sypa->N_AA].Bd[2]);
                distabs = dotPro(dVec, dVec);
                if( distabs < SW2HUGE )     { NCDist[m][n] = distabs; }
                else                        { NCDist[m][n] = -1; }
                std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sypa, Chn[m/sypa->N_AA].AmAc[m%sypa->N_AA].Bd[2], Chn[n/sypa->N_AA].AmAc[n%sypa->N_AA].Bd[0]);
                distabs = dotPro(dVec, dVec);
                if( distabs < SW2HUGE )     { NCDist[n][m] = distabs; }
                else                        { NCDist[n][m] = -1; }
            }
            else {
                NCDist[m][n] = -1;
                NCDist[n][m] = -1;
            }
            if(HBcheck(sypa, Chn, m, n)) dEhb -= 1.0;
            if(HBcheck(sypa, Chn, n, m)) dEhb -= 1.0;
        }
    }
    // ATTENTION: extra conditions and situations would have to be applied here as well
    if( nBHB >= 0 ) {       // previously broken HB can rebond
        for( int i=0; i<nBHB+1; i++ ) {
            for( int j=0; j<sypa->N_CH*sypa->N_AA; j++ ) {
                if(HBcheck(sypa, Chn, BrokenHB[i][0], j)) dEhb -= 1.0;
                if(HBcheck(sypa, Chn, j, BrokenHB[i][1])) dEhb -= 1.0;
            }
        }
    }
    for( int i=sp; i<ep; i++ ) {
        for( int j=0; j<4; j++ ) {
            LinkListUpdate(sypa, Chn, i, j);
        }
    }
    // SC interactions
    deltaE = EO_SegSeg(sypa, Chn, sp, ep, 0, sp, 1) + EO_SegSeg(sypa, Chn, sp, ep, ep, sypa->N_CH*sypa->N_AA, 1) + dEhb - Eold;

    return true;
}
// translation move of a whole chain
bool translation(SysPara *sp, Chain Chn[], int iChn, double &deltaE)
{
    Bead BdCpy[4*sp->N_AA];                 // copy of beads in moved chain
    double dVec[3], newpos[3];          // displacement vector, new position (before PBC)
    double Eold, dEhb, distabs;
    int BrokenHB[sp->N_AA][2], nBHB;

    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            BdCpy[i*4+j] = Chn[iChn].AmAc[i].Bd[j];
        }
    }
    Eold = EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, 0, iChn*sp->N_AA, 1) + EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 1);
    for( int i=0; i<3; i++ ) {
        dVec[i] = ( ((double)rng()/(double)rng.max())*2 - 1. )*sp->DTRN_MAX;
    }

    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            for( int k=0; k<3; k++ ) {
                newpos[k] = Chn[iChn].AmAc[i].Bd[j].getR(k) + dVec[k];
                Chn[iChn].AmAc[i].Bd[j].setBC(k, Chn[iChn].AmAc[i].Bd[j].getBC(k) + floor(newpos[k]/sp->L) );
                newpos[k] = newpos[k] - sp->L*floor(newpos[k]/sp->L);
            }
            Chn[iChn].AmAc[i].Bd[j].setR(newpos[0], newpos[1], newpos[2]);
            // overlapp check
            if( EO_SegBead(sp, Chn, iChn, i, j, 0, iChn*sp->N_AA, 0) == -1 || EO_SegBead(sp, Chn, iChn, i, j, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 0) == -1 ) {
                for( int m=0; m<i+1; m++ ) {
                    for( int n=0; n<4; n++ ) {
                        Chn[iChn].AmAc[m].Bd[n] = BdCpy[m*4+n];     // reset chain & exit
                    }
                }
                return false;
            }
        }
    }
    // energy calculation
    //break old HB
    nBHB = -1;
    dEhb = 0.0;
    for( int i=iChn*sp->N_AA; i<(iChn+1)*sp->N_AA; i++ ) {
        if( (HBList[i][0] < iChn*sp->N_AA || HBList[i][0] >= (iChn+1)*sp->N_AA ) && HBList[i][0] > -1 ) {
            BrokenHB[++nBHB][0] = i;
            BrokenHB[nBHB][1] = HBList[i][0];
            HBList[HBList[i][0]][1] = -1;
            HBList[i][0] = -1;
            dEhb += 1.0;
        }
        if( (HBList[i][1] < iChn*sp->N_AA || HBList[i][1] >= (iChn+1)*sp->N_AA ) && HBList[i][1] > -1 ) {
            BrokenHB[++nBHB][1] = i;
            BrokenHB[nBHB][0] = HBList[i][1];
            HBList[HBList[i][1]][0] = -1;
            HBList[i][1] = -1;
            dEhb += 1.0;
        }
    }
    // close new HB
    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
            if( j/sp->N_AA == iChn ) continue;
            std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sp, Chn[iChn].AmAc[i].Bd[0], Chn[j/sp->N_AA].AmAc[j%sp->N_AA].Bd[2]);
            distabs = dotPro(dVec, dVec);
            if( distabs < SW2HUGE ) { NCDist[iChn*sp->N_AA+i][j] = distabs; }
            else                    { NCDist[iChn*sp->N_AA+i][j] = -1; }
            std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sp, Chn[iChn].AmAc[i].Bd[2], Chn[j/sp->N_AA].AmAc[j%sp->N_AA].Bd[0]);
            distabs = dotPro(dVec, dVec);
            if( distabs < SW2HUGE ) { NCDist[j][iChn*sp->N_AA+i] = distabs; }
            else                    { NCDist[j][iChn*sp->N_AA+i] = -1; }
            if( HBcheck(sp, Chn, iChn*sp->N_AA+i, j) )  dEhb -= 1.0;
            if( HBcheck(sp, Chn, j, iChn*sp->N_AA+i) )  dEhb -= 1.0;
        }
    }
    // previously broken HB can rebond
    if( nBHB > -1 ) {
        for( int i=0; i<nBHB+1; i++) {
            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                if( HBcheck(sp, Chn, BrokenHB[i][0], j) )   dEhb -= 1.0;
                if( HBcheck(sp, Chn, j, BrokenHB[i][1]) )   dEhb -= 1.0;
            }
        }
    }
    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            LinkListUpdate(sp, Chn, iChn*sp->N_AA+i, j);
        }
    }

    deltaE = EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, 0, iChn*sp->N_AA, 1) + EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 1) + dEhb - Eold;

    return true;
}
// rotation move of a whole chain
bool rotation(SysPara *sp, Chain Chn[], int iChn, double &deltaE)
{
    Bead BdCpy[4*sp->N_AA];
    double Eold, dEhb, dist2;
    double angle, axisPolar, axisAzim, cos_a, sin_a;
    double rotAxis[3], com[3], dVec[3], newpos[3];
    double rotMtrx[3][3];
    int BrokenHB[sp->N_AA*sp->N_CH][2], nBHB;

    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            BdCpy[i*4+j] = Chn[iChn].AmAc[i].Bd[j];
        }
    }
    Eold = EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, 0, iChn*sp->N_AA, 1) + EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 1);

    angle = ((double)rng()/(double)rng.max()) * sp->DROT_MAX;       // angle of rotation
    cos_a = cos(angle); sin_a = sin(angle);
    axisPolar = ((double)rng()/(double)rng.max()) * M_PI;           // polar angle of rotation axis
    axisAzim = ((double)rng()/(double)rng.max()) * M_PI*2.0;        // azimuthal angle of rotation axis
    rotAxis[0] = cos(axisAzim)*sin(axisPolar);  rotAxis[1] = sin(axisAzim)*sin(axisPolar);  rotAxis[2] = cos(axisPolar);
    // rotation matrix
    rotMtrx[0][0] = cos_a + rotAxis[0]*rotAxis[0]*(1-cos_a);
    rotMtrx[0][1] = rotAxis[0]*rotAxis[1]*(1-cos_a) - rotAxis[2]*sin_a;
    rotMtrx[0][2] = rotAxis[0]*rotAxis[2]*(1-cos_a) + rotAxis[1]*sin_a;
    rotMtrx[1][0] = rotAxis[1]*rotAxis[0]*(1-cos_a) + rotAxis[2]*sin_a;
    rotMtrx[1][1] = cos_a + rotAxis[1]*rotAxis[1]*(1-cos_a);
    rotMtrx[1][2] = rotAxis[1]*rotAxis[2]*(1-cos_a) - rotAxis[0]*sin_a;
    rotMtrx[2][0] = rotAxis[2]*rotAxis[0]*(1-cos_a) - rotAxis[1]*sin_a;
    rotMtrx[2][1] = rotAxis[2]*rotAxis[1]*(1-cos_a) + rotAxis[0]*sin_a;
    rotMtrx[2][2] = cos_a + rotAxis[2]*rotAxis[2]*(1-cos_a);
    // center of mass
    com[0]=0.0; com[1]=0.0; com[2]=0.0;
    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            for( int k=0; k<3; k++ ) {
                com[k] += (Chn[iChn].AmAc[i].Bd[j].getR( k ) + Chn[iChn].AmAc[i].Bd[j].getBC( k )*sp->L) * Chn[iChn].AmAc[i].Bd[j].getM();
            }
        }
    }
    com[0] /= Chn[iChn].getM();    com[1] /= Chn[iChn].getM();    com[2] /= Chn[iChn].getM();
    // distance calc - rotated vector
    for(int i=0; i<sp->N_AA; i++){
        for(int j=0; j<4; j++) {
            for(int k=0; k<3; k++ ){
                dVec[k] = Chn[iChn].AmAc[i].Bd[j].getR(k) + Chn[iChn].AmAc[i].Bd[j].getBC(k)*sp->L - com[k];
            }
            for( int k=0; k<3; k++ ) {
                newpos[k] = dVec[0]*rotMtrx[0][k] + dVec[1]*rotMtrx[1][k] + dVec[2]*rotMtrx[2][k] + com[k];
                Chn[iChn].AmAc[i].Bd[j].setBC( k, floor(newpos[k]/sp->L) );
                newpos[k] = newpos[k] - sp->L*floor(newpos[k]/sp->L);
            }
            Chn[iChn].AmAc[i].Bd[j].setR(newpos[0], newpos[1], newpos[2]);
            // overlapp check
            if( EO_SegBead(sp, Chn, iChn, i, j, 0, iChn*sp->N_AA, 0) == -1 || EO_SegBead(sp, Chn, iChn, i, j, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 0 ) ) {
                for( int k=0; k<i+1; k++ ) {
                    for( int m=0; m<4; m++ ) {
                        Chn[iChn].AmAc[k].Bd[m] = BdCpy[k*4+m];        // reset chain & exit
                    }
                }
                return false;
            }
        }
    }
    // energy calculation
    // break old HB
    nBHB = -1;          // no. broken HB
    dEhb = 0.0;         // change in energy due to HB
    for( int i=iChn*sp->N_AA; i<(iChn+1)*sp->N_AA; i++ ) {
        if( (HBList[i][0] <= iChn*sp->N_AA || HBList[i][0] >= (iChn+1)*sp->N_AA) && HBList[i][0] > -1 ) {
            BrokenHB[++nBHB][0] = i;
            BrokenHB[nBHB][1] = HBList[i][0];
            HBList[HBList[i][0]][1] = -1;
            HBList[i][0] = -1;
            dEhb += 1.0;
        }
        if( (HBList[i][1] <= iChn*sp->N_AA || HBList[i][1] >= (iChn+1)*sp->N_AA) && HBList[i][1] > -1 ) {
            BrokenHB[++nBHB][1] = i;
            BrokenHB[nBHB][0] = HBList[i][1];
            HBList[HBList[i][1]][0] = -1;
            HBList[i][1] = -1;
            dEhb += 1.0;
        }
    }
    // close new HB
    for( int m=iChn*sp->N_AA; m<(iChn+1)*sp->N_AA; m++ ) {
        for( int n=0; n<sp->N_CH*sp->N_AA; n++ ) {
            if( n >= iChn*sp->N_AA && n < (iChn+1)*sp->N_AA) continue;
            if( (m/sp->N_AA != n/sp->N_AA) || (( m/sp->N_AA == n/sp->N_AA) && (abs(m-n) > 3)) ) {
                std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sp, Chn[m/sp->N_AA].AmAc[m%sp->N_AA].Bd[0], Chn[n/sp->N_AA].AmAc[n%sp->N_AA].Bd[2]);
                dist2 = dotPro(dVec, dVec);
                if( dist2 < SW2HUGE )   { NCDist[m][n] = dist2; }
                else                    { NCDist[m][n] = -1; }
                std::tie(dVec[0], dVec[1], dVec[2]) = distVecBC(sp, Chn[m/sp->N_AA].AmAc[m%sp->N_AA].Bd[2], Chn[n/sp->N_AA].AmAc[n%sp->N_AA].Bd[0]);
                dist2 = dotPro(dVec, dVec);
                if( dist2 < SW2HUGE )   { NCDist[n][m] = dist2; }
                else                    { NCDist[n][m] = -1; }
            }
            else {
                NCDist[m][n] = -1;
                NCDist[n][m] = -1;
            }
            if(HBcheck(sp, Chn, m, n)) dEhb -= 1.0;
            if(HBcheck(sp, Chn, n, m)) dEhb -= 1.0;
        }
    }
    if( nBHB >= 0 ) {       // previously broken HB can rebond
        for( int i=0; i<nBHB+1; i++ ) {
            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                if(HBcheck(sp, Chn, BrokenHB[i][0], j)) dEhb -= 1.0;
                if(HBcheck(sp, Chn, j, BrokenHB[i][1])) dEhb -= 1.0;
            }
        }
    }
    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            LinkListUpdate(sp, Chn, iChn*sp->N_AA+i, j);
        }
    }

    deltaE = EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, 0, iChn*sp->N_AA, 1) + EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 1) + dEhb - Eold;

    return true;
}
// check all bond lenght from Chn[sp/N_AA].AmAc[sp%N_AA] to Chn[ep/N_AA].AmAc[ep%N_AA]
bool checkBndLngth(SysPara *sypa, Chain Chn[], int sp, int ep)
{
    double distV[3];
    double dist2;
    bool res;

    res = true;
    for( int i=sp; i<ep; i++ ) {
        if( (i+1)%sypa->N_AA != 0 ) {
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[0], Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[1]);
            if( abs( absVec(distV) / BND_NCa - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[0]] and AmAc[%d].Bd[1]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i, BND_NCa, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[0], Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[2]);
            if( abs( absVec(distV) / PBND_NC - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[0]] and AmAc[%d].Bd[2]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i, PBND_NC, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[0], Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[3]);
            if( abs( absVec(distV) / Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].getDisR(0) - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[0]] and AmAc[%d].Bd[3]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].getDisR(0), BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[1], Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[2]);
            if( abs( absVec(distV) / BND_CaC - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[1]] and AmAc[%d].Bd[2]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i, BND_CaC, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[1], Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[3]);
            if( abs( absVec(distV) / Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].getDisR(1) - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[1]] and AmAc[%d].Bd[3]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].getDisR(1), BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[1], Chn[i/sypa->N_AA].AmAc[(i+1)%sypa->N_AA].Bd[0]);
            if( abs( absVec(distV) / PBND_CaN - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[1]] and AmAc[%d].Bd[0]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i+1, PBND_CaN, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[1], Chn[i/sypa->N_AA].AmAc[(i+1)%sypa->N_AA].Bd[1]);
            if( abs( absVec(distV) / PBND_CaCa - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[1]] and AmAc[%d].Bd[1]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i+1, PBND_CaCa, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[2], Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[3]);
            if( abs( absVec(distV) / Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].getDisR(2) - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[2]] and AmAc[%d].Bd[3]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].getDisR(2), BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[2], Chn[i/sypa->N_AA].AmAc[(i+1)%sypa->N_AA].Bd[0]);
            if( abs( absVec(distV) / BND_CN - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[2]] and AmAc[%d].Bd[0]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i+1, BND_CN, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[i%sypa->N_AA].Bd[2], Chn[i/sypa->N_AA].AmAc[(i+1)%sypa->N_AA].Bd[1]);
            if( abs( absVec(distV) / PBND_CCa - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[2]] and AmAc[%d].Bd[1]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        i, i+1, PBND_CCa, BND_FLUCT, absVec(distV) );
                        res = false;
            }
        }
        else {
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[0], Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[1]);
            if( abs( absVec(distV) / BND_NCa - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[0]] and AmAc[%d].Bd[1]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        sypa->N_AA-1, sypa->N_AA-1, BND_NCa, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[0], Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[2]);
            if( abs( absVec(distV) / PBND_NC - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[0]] and AmAc[%d].Bd[2]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        sypa->N_AA-1, sypa->N_AA-1, PBND_NC, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[0], Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[3]);
            if( abs( absVec(distV) / Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].getDisR(0) - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[0]] and AmAc[%d].Bd[3]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        sypa->N_AA-1, sypa->N_AA-1, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].getDisR(0), BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[1], Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[2]);
            if( abs( absVec(distV) / BND_CaC - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[1]] and AmAc[%d].Bd[2]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        sypa->N_AA-1, sypa->N_AA-1, BND_CaC, BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[1], Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[3]);
            if( abs( absVec(distV) / Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].getDisR(1) - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[1]] and AmAc[%d].Bd[3]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        sypa->N_AA-1, sypa->N_AA-1, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].getDisR(1), BND_FLUCT, absVec(distV) );
                        res = false;
            }
            std::tie(distV[0], distV[1], distV[2]) = distVecBC(sypa, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[2], Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].Bd[3]);
            if( abs( absVec(distV) / Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].getDisR(2) - 1.0 ) > BND_FLUCT+0.01 ) {
                printf("Invalid bond length between AmAc[%d].Bd[2]] and AmAc[%d].Bd[3]:\t should be %4.2f (+- %.3f)\tis %f\n", 
                        sypa->N_AA-1, sypa->N_AA-1, Chn[i/sypa->N_AA].AmAc[sypa->N_AA-1].getDisR(2), BND_FLUCT, absVec(distV) );
                        res = false;
            }
        }
    }

    return res;
}
//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXXX  MAINTENANCE FUNCTIONS  XXXXXXXXXX

// resets the counter of boundary crossings so that the real coordinates move back to the simulation box
bool resetBCcouter(SysPara *sp, Chain Chn[]) 
{
    double highestBC[3];

    for( int i=0; i<sp->N_CH; i++ ) {
        highestBC[0] = 0;  highestBC[1] = 0; highestBC[2] = 0;
        for( int j=0; j<sp->N_AA*4; j++ ) {
            for( int k=0; k<3; k++ ) {
                if( abs(Chn[i].AmAc[j/4].Bd[j%4].getBC(k)) > abs(highestBC[k]) ) {
                    highestBC[k] = Chn[i].AmAc[j/4].Bd[j%4].getBC(k);
                }
            }
        }
        for( int j=0; j<sp->N_AA*4; j++ ) {
            //Chn[i].AmAc[j/4].Bd[j%4].setR( Chn[i].AmAc[j/4].Bd[j%4].getR(0) - highestBC[0]*L, Chn[i].AmAc[j/4].Bd[j%4].getR(1) - highestBC[1]*L, Chn[i].AmAc[j/4].Bd[j%4].getR(2) - highestBC[2]*L );
            for( int k=0; k<3; k++ ) {
                Chn[i].AmAc[j/4].Bd[j%4].addBC(k, highestBC[k]);
            }
        }
    }

}

//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXXX  OBSERVABLE FUNCTIONS  XXXXXXXXXXX

// calculate radius of gyration
int calc_gyration_radius(SysPara *sp, Output *ot, Chain Chn[], int eBin)
{
    double com[3], dist[3];
    double rG2;

    for( int i=0; i<sp->N_CH; i++ ) { // chain
        // center of mass
        com[0]=0.0; com[1]=0.0; com[2]=0.0;
        for( int j=0; j<sp->N_AA; j++ ) {
            for( int k=0; k<4; k++ ) {
                com[ 0 ] += (Chn[i].AmAc[j].Bd[k].getR( 0 ) + Chn[i].AmAc[j].Bd[k].getBC( 0 )*sp->L) * Chn[i].AmAc[j].Bd[k].getM();
                com[ 1 ] += (Chn[i].AmAc[j].Bd[k].getR( 1 ) + Chn[i].AmAc[j].Bd[k].getBC( 1 )*sp->L) * Chn[i].AmAc[j].Bd[k].getM();
                com[ 2 ] += (Chn[i].AmAc[j].Bd[k].getR( 2 ) + Chn[i].AmAc[j].Bd[k].getBC( 2 )*sp->L) * Chn[i].AmAc[j].Bd[k].getM();
            }
        }
        com[0] /= Chn[i].getM();    com[1] /= Chn[i].getM();    com[2] /= Chn[i].getM();
        // radius of gyration squared
        rG2 = 0.0;
        for( int j=0; j<sp->N_AA; j++ ) { // amino acid
            for( int k=0; k<4; k++ ) { // bead
                for( int l=0; l<3; l++ ) { // coordinate
                    dist[l] = (Chn[i].AmAc[j].Bd[k].getR(l) + Chn[i].AmAc[j].Bd[k].getBC( l )*sp->L) - com[l];
                }
                rG2 += dotPro(dist, dist) * Chn[i].AmAc[j].Bd[k].getM();
            }
        }
        rG2 /= Chn[i].getM();
        ot->rGyrCur[i] = rG2;
        ot->rGyr[i][eBin] += rG2;
    }
    return 0;
}
// calculate tensor of gyration
int calc_gyration_tensor(SysPara *sp, Output *ot, Chain Chn[], int eBin)
{
    double com[3];                                      // center of mass 
    double rGxx, rGyy, rGzz, rGtemp;                    // principal moments of gyration tensor
    double Mat11, Mat12, Mat13, Mat22, Mat23, Mat33;    // gyration tensor elements
    double dx, dy, dz;                                  // distance to center of mass
    double A, B, C, D;                                  // coefficients of characteristic polynomial
    double p,q, discrim;                                // factors of reduced form and discriminante

    for( int i=0; i<sp->N_CH; i++ ) { // chain
        //center of mass
        com[0]=0.0; com[1]=0.0; com[2]=0.0;
        for( int j=0; j<sp->N_AA; j++ ) { // amino acid
            for( int k=0; k<4; k++ ) { // bead
                com[ 0 ] += (Chn[i].AmAc[j].Bd[k].getR( 0 ) + Chn[i].AmAc[j].Bd[k].getBC( 0 )*sp->L) * Chn[i].AmAc[j].Bd[k].getM();
                com[ 1 ] += (Chn[i].AmAc[j].Bd[k].getR( 1 ) + Chn[i].AmAc[j].Bd[k].getBC( 1 )*sp->L) * Chn[i].AmAc[j].Bd[k].getM();
                com[ 2 ] += (Chn[i].AmAc[j].Bd[k].getR( 2 ) + Chn[i].AmAc[j].Bd[k].getBC( 2 )*sp->L) * Chn[i].AmAc[j].Bd[k].getM();
            }
        }
        com[0] /= Chn[i].getM();    com[1] /= Chn[i].getM();    com[2] /= Chn[i].getM();


        // gyration tensor
        Mat11 = Mat12= Mat13 = Mat22 = Mat23 = Mat33 = 0.0;
        for( int j=0; j<sp->N_AA; j++ ) { // amino acid
            for( int k=0; k<4; k++ ) { // bead
                dx = (Chn[i].AmAc[j].Bd[k].getR(0) + Chn[i].AmAc[j].Bd[k].getBC( 0 )*sp->L) - com[0];
                dy = (Chn[i].AmAc[j].Bd[k].getR(1) + Chn[i].AmAc[j].Bd[k].getBC( 1 )*sp->L) - com[1];
                dz = (Chn[i].AmAc[j].Bd[k].getR(2) + Chn[i].AmAc[j].Bd[k].getBC( 2 )*sp->L) - com[2];
                Mat11 += dx*dx*Chn[i].AmAc[j].Bd[k].getM();
                Mat12 += dx*dy*Chn[i].AmAc[j].Bd[k].getM();
                Mat13 += dx*dz*Chn[i].AmAc[j].Bd[k].getM();
                Mat22 += dy*dy*Chn[i].AmAc[j].Bd[k].getM();
                Mat23 += dy*dz*Chn[i].AmAc[j].Bd[k].getM();
                Mat33 += dz*dz*Chn[i].AmAc[j].Bd[k].getM();
            }
        }
        Mat11/=Chn[i].getM(); Mat12/=Chn[i].getM(); Mat13/=Chn[i].getM();
        Mat22/=Chn[i].getM(); Mat23/=Chn[i].getM(); 
        Mat33/=Chn[i].getM();

        // characteristic polynomial
        A = 1.0;
        B = -Mat11-Mat22-Mat33;
        C = Mat11*(Mat22+Mat33) + Mat22*Mat33 - Mat12*Mat12 - Mat23*Mat23 - Mat13*Mat13;
        //D = Mat33*(Mat12*Mat12 - Mat11*Mat22) + Mat23*(2*Mat12*Mat13 + Mat11*Mat23) - Mat13*Mat13*Mat22;
        D = -Mat11*Mat22*Mat33 + Mat12*Mat12*Mat33 - 2*Mat12*Mat13*Mat23 + Mat13*Mat13*Mat22 + Mat11*Mat23*Mat23;

        p = (3.0*C - B*B) / (3.0);
        q = (2.0*B*B*B - 9.0*B*C + 27.0*D) / (27.0);
        discrim = (27.0*D*D + 4.0*B*B*B*D - 18.0*B*C*D + 4.0*C*C*C - B*B*C*C) / (108.0);

        // evaluate discriminate
        if( discrim > 1e-10 ) {
            std::cout << "\rinertia tensor: discriminante > 0  -->  complex solution in Chn[" << i << "]            " << std::endl << "--- Aborting calculation ---" << std::endl;
            return -1;
        }
        else if( discrim < 1e-10 && discrim > -1e-10 ) {
            if( (p != 0.0) && (q != 0.0) ) {
                rGxx = (B*B*B - 4.0*B*C + 9.0*D) / (3.0*C - B*B);
                rGyy = rGzz = (B*C - 9.0*D) / (6.0*C - 2.0*B*B);
            }
            else if( (p == 0.0) && (q == 0.0) ) {
                rGxx = rGyy = rGzz = -B / 3.0;
            }
            else {
                std::cout << "inertia tensor: discriminante = 0  -->  unknown case: p=" << p << " q=" << q << std::endl << "--- Aborting calculation ---" << std::endl;
                return -1;
            }
        }
        else if( discrim < 1e-10 ) {
            rGxx =  sqrt(-4.0*p/3.0) * cos( 1.0/3.0 * acos(-q/2.0 * sqrt(-27.0/(p*p*p)) )            ) - B/3.0;
            rGyy = -sqrt(-4.0*p/3.0) * cos( 1.0/3.0 * acos(-q/2.0 * sqrt(-27.0/(p*p*p)) ) + M_PI/3.0 ) - B/3.0;
            rGzz = -sqrt(-4.0*p/3.0) * cos( 1.0/3.0 * acos(-q/2.0 * sqrt(-27.0/(p*p*p)) ) - M_PI/3.0 ) - B/3.0;
        }
        // sorting the eigenvalues
        if (rGyy > rGxx) {
            rGtemp = rGxx;
            rGxx   = rGyy;
            rGyy   = rGtemp;
        }
        if (rGzz > rGxx) {
            rGtemp = rGxx;
            rGxx   = rGzz;
            rGzz   = rGtemp;
        }
        if (rGzz > rGyy) {
            rGtemp = rGyy;
            rGyy   = rGzz;
            rGzz   = rGtemp;
        }
        // current principal moments
        ot->tGyrEigCur[i][0] = rGxx;
        ot->tGyrEigCur[i][1] = rGyy;
        ot->tGyrEigCur[i][2] = rGzz;
        // summation of principal moments
        ot->tGyrEig[i][eBin][0] += rGxx;
        ot->tGyrEig[i][eBin][1] += rGyy;
        ot->tGyrEig[i][eBin][2] += rGzz;

    }

    return 0;
}
// calculate van-der-Waals energy
int calc_vanderWaals(SysPara *sp, Output *ot, Chain Chn[], int eBin, double &vdW_intra, double &vdW_inter, int recalc)
{
    int neighBox[3], centrBox[3];
    int searchBox, neighBead;
    int Box, dbx, dby, dbz;
    double distV[3];
    double dist2;

    if(recalc == 1) {
        vdW_intra = 0;
        vdW_inter = 0;
        for( int i=0; i< sp->N_CH; i++) {
            for( int j=0; j<sp->N_AA; j++ ) {
                Box = Chn[i].AmAc[j].Bd[3].getBox();
                centrBox[0] = Box % sp->NBOX;
                centrBox[1] = (Box/sp->NBOX) % sp->NBOX;
                centrBox[2] = Box/(double)(sp->NBOX*sp->NBOX);
                for( int m=0; m<27; m++ ) {
                    dbx = m%3; dby = (m/3)%3; dbz = m/9;
                    neighBox[0] = (centrBox[0] + dbx-1 + sp->NBOX) % sp->NBOX;
                    neighBox[1] = (centrBox[1] + dby-1 + sp->NBOX) % sp->NBOX;
                    neighBox[2] = (centrBox[2] + dbz-1 + sp->NBOX) % sp->NBOX;
                    searchBox = neighBox[0] + neighBox[1]*sp->NBOX + neighBox[2]*sp->NBOX*sp->NBOX;
                    neighBead = neighHead[searchBox];
                    while(neighBead!=-1) {
                        if( ( (i*sp->N_AA+j)*4 ) < neighBead ) {
                            if( neighBead%4 == 3 ) {
                                std::tie(distV[0], distV[1], distV[2]) = distVecBC(sp, Chn[i].AmAc[j].Bd[3], Chn[neighBead/(4*sp->N_AA)].AmAc[(neighBead/4)%sp->N_AA].Bd[neighBead%4]);
                                dist2 = dotPro(distV, distV);
                                // intra-chain SC contact
                                if( neighBead >= i*sp->N_AA*4 && neighBead < (i+1)*sp->N_AA*4 ) {
                                    vdW_intra += E_single(Chn, i, j, neighBead/(4*sp->N_AA), (neighBead/4)%sp->N_AA, dist2);
                                }
                                // inter-chain SC contact
                                else {
                                    vdW_inter += E_single(Chn, i, j, neighBead/(4*sp->N_AA), (neighBead/4)%sp->N_AA, dist2);
                                }
                            }
                        }
                        neighBead = neighList[neighBead];
                    }
                }
            }
        }
    }

    ot->vdWener[eBin*2 + 0] += vdW_intra;
    ot->vdWener[eBin*2 + 1] += vdW_inter;

    return 0;
}

//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXXXX  LINK LIST FUNCTIONS  XXXXXXXXXXX

// returns neighbour list box of bead
int assignBox(SysPara *sp, Bead Bd) 
{
    int xBox, yBox, zBox, Box;
    xBox = floor(Bd.getR(0) / sp->LBOX);
    yBox = floor(Bd.getR(1) / sp->LBOX);
    zBox = floor(Bd.getR(2) / sp->LBOX);
    Box = zBox*sp->NBOX*sp->NBOX + yBox*sp->NBOX + xBox;
    return Box;
}
// update linked list (neighbour list) -> insert AmAc[i1].Bd[j1]
int LinkListInsert(SysPara *sp, Chain Chn[], int i1, int j1)
{
    int start, i, prev, next;
    next = neighHead[Chn[i1/sp->N_AA].AmAc[i1%sp->N_AA].Bd[j1].getBox()];
    prev = -1;
    i = 0;
    while(i != -1){
        i = next;
        if( next < i1*4+j1) {
            if(prev == -1) {
                neighHead[Chn[i1/sp->N_AA].AmAc[i1%sp->N_AA].Bd[j1].getBox()] = i1*4+j1;
                neighList[i1*4+j1] = next;
                return 0;
            }
            neighList[prev] = i1*4+j1;
            neighList[i1*4+j1] = next;
            return 0;
        }
        prev = i;
        next = neighList[i];
    }
    std::cout << "LinkListUpdate() failed on AmAc[" << i1 << "].Bd[" << j1 << "]" << endl;
    return -1;
}
// update Linked List position of Chn[i1/N_AA].AmAc[i1%N_AA].Bd[j1]
int LinkListUpdate(SysPara *sp, Chain Chn[], int i1, int j1)
{
    int oldBox;
    int prevNei, nextNei;

    oldBox = Chn[i1/sp->N_AA].AmAc[i1%sp->N_AA].Bd[j1].getBox();
    Chn[i1/sp->N_AA].AmAc[i1%sp->N_AA].Bd[j1].setBox(assignBox( sp, Chn[i1/sp->N_AA].AmAc[i1%sp->N_AA].Bd[j1] ));
    if( oldBox != Chn[i1/sp->N_AA].AmAc[i1%sp->N_AA].Bd[j1].getBox() ) {
        prevNei = neighHead[ oldBox ];
        if( prevNei == i1*4+j1 ) {
            neighHead[ oldBox ] = neighList[ prevNei ];
            neighList[ prevNei ] = -1;
            prevNei = -1;
        }
        while(prevNei != -1) {
            nextNei = neighList[ prevNei ];
            if(nextNei == i1*4+j1) {
                neighList[ prevNei ] = neighList[ nextNei ];
                neighList[ nextNei ] = -1;
                prevNei = -1;
            }
            else {
                prevNei = nextNei;
            }
        }
        LinkListInsert(sp, Chn, i1, j1);
        //printf("\t\tneighlist was updated:\t\tchn[%d].Bd[%d]\tindex=%3d\toldBox=%3d\tnewBox=%3d\n", i1, j1, i1*4+j1, oldBox, AmAc[i1].Bd[j1].getBox());
    }
    else
        //printf("\t\tneighlist wasn't updated\tchn[%d].Bd[%d]\tindex=%3d\toldBox=%3d\n", i1, j1, i1*4+j1, oldBox);
    return 0;
}

// deallocate memory
int output_memory_deallocation(SysPara *sp, Output *ot)
{
    delete[] ot->H;
    delete[] ot->lngE;
    delete[] ot->contHB;
    delete[] ot->Ree2;
    for( int i=0; i<sp->N_CH; i++ ) { delete[] ot->rGyr[i]; }
    delete[] ot->rGyr;
    delete[] ot->rGyrCur;
    for( int i=0; i<sp->N_CH; i++ ) { 
        for( int j=0; j<sp->NBin; j++ ) { delete[] ot->tGyrEig[i][j]; }
        delete[] ot->tGyrEig[i]; }
    delete[] ot->tGyrEig;
    for( int i=0; i<sp->N_CH; i++ ) {
        delete[] ot->tGyrEigCur[i]; }
    delete[] ot->tGyrEigCur;
    delete[] ot->vdWener;
    delete[] ot->Et;
    delete[] ot->conf_n;
    delete[] ot->conf_wt;

    return 0;
}