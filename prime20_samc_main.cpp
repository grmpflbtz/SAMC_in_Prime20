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
//mt19937 rng(2);                        // debug

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  >>   FUNCTION DECLARATIONS   <<  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

double dotPro(double vecA[], double vecB[]);                                            // dot product
double absVec(double vec[]);                                                            // absolute of vector
tuple<double,double,double> crossPro(double vecA[], double vecB[]);                     // cross product
tuple<double,double,double> distVecBC(SysPara *sp, Bead vecA, Bead vecB);               // distance vector
double RND();                                                                           // RNG generating double (-1,1)
int CommandInitialize(int argc, char *argv[], Header *hd);                              // identify file names from command input
void ConsoleOutputHead(SysPara *sp);                                                    // writes simulation head to console
bool newChain(SysPara *sp, Chain Chn[], int chnNum);                                    // creates new chain
bool readParaInput(SysPara *sp, Header *hd);                                            // read system arameters from file
bool readCoord(SysPara *sp, Header *hd, Chain Chn[]);                                   // read chain config from file
bool readPrevRunInput(SysPara *sp, Chain Chn[], string inputFile, double lngE[], long unsigned int H[], long unsigned int &tcont, double &gammasum);      // reads lngE, H, gammasum, and t from input file
bool extra_lngE(SysPara *sp, Header *hd, double lngE[]);                          // reads lngE data from file
bool outputPositions(SysPara *sp, Chain Chn[], string name, int mode);                  // writes positions to file "name"
bool BackupSAMCrun(SysPara *sp, Chain Chn[], Timer &Timer, unsigned long int t, double gammasum, double gamma, unsigned long naccept[], unsigned long nattempt[], double lngE[], unsigned long H[], double E);    // backup function in SAMC run
bool BackupProdRun(SysPara *sp, Timer &Timer, unsigned long int t, unsigned long int H[]);           // backup of observables for production run
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
bool checkBndLngth(SysPara *sypa, Chain Chn[], int sp, int ep);                         // check all bond length from Chn[sp/N_AA].AmAc[sp%N_AA] to Chn[ep/N_AA].AmAc[ep%N_AA]
bool resetBCcouter(SysPara *sp, Chain Chn[]);                                           // resets the counter of boundary crossings so that the real coordinates move back to the simulation box
int assignBox(SysPara *sp, Bead Bd);                                                    // assigns neighbour list box to bead
int LinkListInsert(SysPara *sp, Chain Chn[], int i1, int j1);                           // insert particle AmAc[i1].Bd[j1] into Linked List
int LinkListUpdate(SysPara *sp, Chain Chn[], int i1, int j1);                           // update Linked List position of AmAc[i1].Bd[j1]

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX  >>   MAIN FUNCTION   <<  XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

int main(int argc, char *argv[]) 
{
    SysPara *sp = new SysPara;                              // system parameters
    Header *hd = new Header;                                // file names for input/output

    Timer Timer;                                            // timer for simulation
    Chain *Chn;                                             // chains in system
    Bead *BdCpy;                                            // copy of Beads
    uniform_real_distribution<double> realdist01(0.0,1.0);  // uniform distribution of real numbers
    ofstream backup;
    std::ostringstream oss;
    string filename;
    double Eold, Enew, Ecur, deltaE;                        // energy values
    double dist[3];                                         // distance vector
    double distabs;                                         // absolute distance
    bool accept;                                            // acceptance value of move in SAMC
    int i_rand, ip, jp, moveselec, movetype;                // variables for performing moves
    int oldBox, prevNei, nextNei;                           // variables for neighbour list calculations
    long unsigned int t, tcont, tTimer;                     // time variable keeping track of # steps passed
    int tstep, tE, twrite, tBCreset;
    long unsigned int *H;                                   // energy histogram
    int eBin_n, eBin_o;                                     // energy bin of new and old energy value
    long unsigned int nattempt[4], naccept[4];              // no. attempted moves and no. accepted moves
    double *lngE;                                           // ln g(E): logarithm of density of states g(E) → current approximation
    double *contHB;                                         // contact matrix of hydrogen bonds
    int *conf_n, *conf_write_t;                             // number of configurations written for the selcted energy and last time writing a config for this energy
    double gamma, gammasum, lngE_old, lngE_new;             // gamma value and sum over gamma(t), ln g(E_old) and ln g(E_new)

    CommandInitialize(argc, argv, hd);

    readParaInput(sp, hd);

    /*      sp->CONFIG_ENER = new double[sp->CONFIG_N]; sp->CONFIG_ENER[0] = -1.32; sp->CONFIG_ENER[1] = 0.4;
            sp->CONFIG_VAR = 0.02;*/

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
    H = new long unsigned int[sp->NBin];
    lngE = new double[sp->NBin];
    contHB = new double[sp->NBin * sp->N_CH*sp->N_AA * sp->N_CH*sp->N_AA];
    conf_n = new int[sp->CONFIG_N];
    conf_write_t = new int[sp->CONFIG_N];


    fill_n(neighHead, sp->NBOX*sp->NBOX*sp->NBOX, -1), fill_n(neighList, 4*sp->N_AA*sp->N_CH, -1);

    ConsoleOutputHead(sp);

    // initialization of some values
    for( int i=0; i<4; i++ ) {
        nattempt[i] = 0;
        naccept[i] = 0;
    }

    // lngE, H, gammasum, t: read from input file or start new
    if( readPrevRunInput(sp, Chn, "input.dat", lngE, H, tcont, gammasum) ) {
        gamma = sp->GAMMA_0*sp->T_0/tcont;
    }
    else {
        // initialize lngE and H and co.
        for( int i=0; i<sp->NBin; i++ ) {
            lngE[i] = 0;
            H[i] = 0;
        }
        gamma = sp->GAMMA_0;
        gammasum = 0.0;
        tcont = 0;
    }
    // check if there is an extra input file for lngE
    if( extra_lngE(sp, hd, lngE ) ) {
        tcont = 0;
        if(sp->MEASURE) {
            if(sp->HB_CONTMAT) {
                std::cout << "production run with fixed ln g(E)" << std::endl;
                for( int i=0; i<sp->NBin; i++ ) { 
                    H[i] = 0; 
                    for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                            contHB[i*sp->N_CH*sp->N_AA*sp->N_CH*sp->N_AA + j*sp->N_CH*sp->N_AA + k] = 0;
                        }
                    }
                }
            }
            if(sp->WRITE_CONFIG) {
                for( int i=0; i<sp->CONFIG_N; i++ ) {
                    conf_n[i] = 0;
                    conf_write_t[i] = 0;
                }
            }
        }
    }
    // geometry: read from input file or create new
    if( !readCoord(sp, hd, Chn) ) {
        std::cout << "creating new Chain" << std::endl;
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
    }

    // position check file output
    outputPositions(sp, Chn, hd->dbposi, 0);

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
    }
    // move new chains to randomize initial configuration. no energy-dependent acception criterion. all legal moves are accepted
    t = 0;
    while( true ) {
        // end pre-SAMC movement after tStart moves and within desired energy window
        if( (t >= sp->tStart) && (Eold >= sp->EMin) && (Eold < sp->EStart) ) {
            eBin_o = floor((Eold - sp->EMin)/sp->BinW);
            if( eBin_o == sp->NBin) { eBin_o = sp->NBin-1; }
            std::cout << "pre-SAMC move " << t << "\r" << std::flush;
            break;
        }

        if( t%1000 == 0 ) {
            std::cout << "pre-SAMC move " << t << "\r" << std::flush;
        }
        // HBList copy
        for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
            HBLcpy[k][0] = HBList[k][0];
            HBLcpy[k][1] = HBList[k][1];
        }
        deltaE = 0.0;
        // select move type
        moveselec = trunc(realdist01(rng)*(sp->WT_WIGGLE + sp->WT_PHI + sp->WT_PSI + sp->WT_TRANS));
        if( moveselec < sp->WT_WIGGLE )                                         { movetype = 0; }   // wiggle
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI )                         { movetype = 1; }   // rotPhi
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI )              { movetype = 2; }   // rotPsi
        else if( moveselec < sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI+sp->WT_TRANS ) { movetype = 3; }   // translation
        else { std::cout << "ERROR\tno movetype was selected" << endl; return 0; }

        switch( movetype ) {
            case 0:
                i_rand = trunc(realdist01(rng)*(sp->N_CH*sp->N_AA*4));
                ip = i_rand/4;                              // amino acid identifier
                jp = i_rand%4;                              // bead
                accept = wiggle(sp, Chn, ip/sp->N_AA, ip%sp->N_AA, jp, deltaE);
                break;
            case 1:
                ip = trunc(realdist01(rng)*(sp->N_CH*sp->N_AA));    // amino acid identifier of the rotation origin
                jp = trunc(realdist01(rng)*2);              // rotate lower or higher part
                accept = rotPhi(sp, Chn, ip, jp, deltaE);
                break;
            case 2:
                ip = trunc(realdist01(rng)*(sp->N_CH*sp->N_AA));    // amino acid identifier of the rotation origin
                jp = trunc(realdist01(rng)*2);              // rotate lower or higher part
                accept = rotPsi(sp, Chn, ip, jp, deltaE);
                break;
            case 3:
                ip = trunc(realdist01(rng)*sp->N_CH);
                accept = translation(sp, Chn, ip, deltaE);
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

        if( !checkBndLngth(sp, Chn, 0, sp->N_CH*sp->N_AA) ) {               // check bond length after every move - if this failes: abort run
            std::cerr << endl << "bond length error: t=" << t << std::endl;
            std::cerr << "Energy = " << Eold << std::endl;
            outputPositions(sp, Chn, hd->dbposi, 1);
            this_thread::sleep_for(chrono::milliseconds(200));
            return 0;
        }
        Ecur = E_check(sp, Chn);
        if( abs(Eold-Ecur) > 0.01 ) {
            std::cerr << endl << "ERROR → energies are not equal: Eold=" << std::fixed << std::setprecision(3) << Eold << " Ecur=" << Ecur << endl;
            std::cerr.precision(3); std::cerr << std::fixed;
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                std::cerr << "k" << k << "\t";
                for( int m=0; m<sp->N_CH*sp->N_AA; m++ ) {
                    std::cerr << NCDist[k][m] << "\t";
                } std::cerr << endl;
            } std::cerr << endl;
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                std::cerr << k << "  " << HBList[k][0] << "\t" << HBList[k][1] << endl;
            }
            outputPositions(sp, Chn, hd->dbposi, 1);
            this_thread::sleep_for(chrono::milliseconds(200));

            return 0;
        }
        t++;
    }
    std::cerr << std::endl << "starting with energy E=" << Eold << std::endl;
    // configuration when starting SAMC
    outputPositions(sp, Chn, hd->dbposi, 1);

    /*              XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    XXXXXXXXXXXXXXXXXXXX    SAMC loop    XXXXXXXXXXXXXXXXXXXX
                    XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX */    

    t = tcont, tE = 1, twrite = 1, tBCreset = 1, tTimer = 0;
    while( t<sp->T_MAX) {
        // lets go
        tstep = 0;          // nstep steps are merged into one SAMC loop
        while( tstep < sp->nstep ) {

            if( t%((int)1e3) == 0 ) {
                Timer.PrintProgress(tTimer, sp->T_MAX-tcont);
            }

            // select move type
            moveselec = trunc(realdist01(rng)*(sp->WT_WIGGLE + sp->WT_PHI + sp->WT_PSI + sp->WT_TRANS));
            if( moveselec < sp->WT_WIGGLE )                                             { movetype = 0; }   // wiggle
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI) )                           { movetype = 1; }   // rotPhi
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI) )                { movetype = 2; }   // rotPsi
            else if( moveselec < (sp->WT_WIGGLE+sp->WT_PHI+sp->WT_PSI+sp->WT_TRANS) )   { movetype = 3; }   // translation
            else { std::cerr << endl << "ERROR\tno movetype was selected" << endl; return 0; }
            // chain copy setup
            for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                HBLcpy[k][0] = HBList[k][0];                                                                // improve efficiency! → select copied part based on movetype
                HBLcpy[k][1] = HBList[k][1];
            }
            deltaE = 0.0;
            nattempt[movetype]++;
            switch( movetype ) {
                case 0:
                    i_rand = trunc(realdist01(rng)*(sp->N_CH*sp->N_AA*4));
                    ip = i_rand/4;                              // amino acid identifier
                    jp = i_rand%4;                              // bead
                    BdCpy[ip*4+jp] = Chn[ip/sp->N_AA].AmAc[ip%sp->N_AA].Bd[jp];
                    accept = wiggle(sp, Chn, ip/sp->N_AA, ip%sp->N_AA, jp, deltaE);
                    break;
                case 1:
                    ip = trunc(realdist01(rng)*(sp->N_CH*sp->N_AA));    // amino acid identifier of the rotation origin
                    jp = trunc(realdist01(rng)*2);              // rotate lower or higher part
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
                    ip = trunc(realdist01(rng)*(sp->N_CH*sp->N_AA));    // amino acid identifier of the rotation origin
                    jp = trunc(realdist01(rng)*2);              // rotate lower or higher part
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
                    ip = trunc(realdist01(rng)*sp->N_CH);
                    for( int i=0; i<sp->N_AA; i++ ) {
                        for( int j=0; j<4; j++ ) {
                            BdCpy[ip*sp->N_AA*4+i*4+j] = Chn[ip].AmAc[i].Bd[j];
                        }
                    }
                    accept = translation(sp, Chn, ip, deltaE);
                    break;
                default:
                    std::cerr << "Unusable move type selected: movetype = " << movetype << endl;
            }

            if( accept ) {      // legal move → check energy window and acceptance criterion
                Enew = Eold + deltaE;
                // check if Enew is in energy window → if not: don't count attempt
                if( Enew<sp->EMin || Enew > sp->EMax ) {
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
                        default:    // rotations
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
                    nattempt[movetype]--;
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

                accept = acceptance(lngE[eBin_o], lngE[eBin_n]);
            }
            if( accept ) {      // Move accepted → update lngE, H, ...
                naccept[movetype]++;
                H[eBin_n]++;
                if(!sp->MEASURE) lngE[eBin_n] += gamma;
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
                H[eBin_o]++;
                if(!sp->MEASURE) lngE[eBin_o] += gamma;
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
                    default:    // rotation
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
            if(sp->MEASURE) {
                if(sp->HB_CONTMAT) {
                    for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
                        if( HBList[i][0] > -1 ) { contHB[ eBin_o*sp->N_CH*sp->N_AA*sp->N_CH*sp->N_AA + i*sp->N_CH*sp->N_AA + HBList[i][0] ] += 1; }
                    }
                }
                if(sp->WRITE_CONFIG) {
                    for( int i=0; i<sp->CONFIG_N; i++ ) {
                        if( abs(Eold - sp->CONFIG_ENER[i]) <= sp->CONFIG_VAR && conf_n[i]<10 ) {
                            if( conf_write_t[i] + 10000 < t ) {
                                oss.str("");
                                oss << "coordinates_E" << i << "_" << conf_n[i] << ".xyz";
                                filename = oss.str();
                                backup.open(filename);
                                if(backup.is_open()) {
                                    backup << "# Coordinates of " << sp->N_CH << " " << sp->N_AA << "-mer(s) with sequence " << sp->AA_seq << " at E=" << Eold << std::endl;
                                    backup.close();
                                    outputPositions(sp, Chn, filename, 1);
                                    std::cout << std::endl << "wrote coordinate file " << filename << " at E=" << Eold << std::endl;
                                } else {
                                    std::cout << std::endl << "error opening " << filename << std::endl;
                                }
                                conf_n[i]++;
                                conf_write_t[i]=t;
                            }
                        }
                    }
                }
            }

            gammasum += gamma;
            tstep++;

            //check histogram anomalies
            for( int i=0; i<sp->NBin; i++ ) {
                if( H[i] > sp->nstep*sp->T_MAX) {
                    std::cout << std::endl << "What the actual fuck is happening in here?" << std::endl;
                    std::cout << "Eierkucheeeeeeeen" << std::endl;
                }
            }

        }   // end of one MC step
        //gamma update
        if(!sp->MEASURE){
            gamma = sp->GAMMA_0*sp->T_0/((double)max(sp->T_0,t));       // from Liang et al.
        }

        if( twrite == sp->T_WRITE ) {           // write Backup-File
            if(!sp->MEASURE) {
                BackupSAMCrun(sp, Chn, Timer, t, gammasum, gamma, naccept, nattempt, lngE, H, Eold);
            } else {
                BackupProdRun(sp, Timer, t, H);
                if(sp->HB_CONTMAT) {
                    filename = "HBmat.dat";
                    backup.open(filename, ios::out);
                    if(backup.is_open() ) {
                        backup << "# Hydrogen Bind contact matrices after " << t+1 << " steps" << std::endl;
                        std::setprecision(3); std::fixed;
                        for( int i=0; i<sp->NBin; i++ ) {
                            for( int j=0; j<sp->N_CH*sp->N_AA; j++ ) {
                                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                                    backup << contHB[i*sp->N_CH*sp->N_AA*sp->N_CH*sp->N_AA + j*sp->N_CH*sp->N_AA + k]/H[i] << " ";
                                }
                                backup << std::endl;
                            }
                        }
                        backup.close();
                    }
                    else {
                        std::cout << endl << "error opening " << filename << endl;
                    }
                }
            }
            twrite = 0;
        }

        // check bond length after every move - if this fails: abort run
        /*if( !checkBndLngth(Chn, 0, N_CH*N_AA) ) {
            std::cerr << endl << "bond length error" << endl;
            std::cerr << "Energy = " << Eold << endl;
            outputPositions(Chn, "AnorLondo.xyz", 1);
            this_thread::sleep_for(chrono::milliseconds(200));
            return 0;
        }*/

        // Energy check
        if( tE == 10000 ) {
            Ecur = E_check(sp, Chn);
            if( abs(Eold-Ecur) > 0.01 ) {
                std::cerr << endl << "ERROR → energies are not equal: Eold=" << std::fixed << std::setprecision(3) << Eold << " Ecur=" << Ecur << endl;
                std::cerr << endl << "HBList" << endl;
                for( int k=0; k<sp->N_CH*sp->N_AA; k++ ) {
                    std::cerr << k << "  " << HBList[k][0] << "\t" << HBList[k][1] << std::endl << std::flush;
                }
                std::cerr << endl;

                outputPositions(sp, Chn, hd->dbposi, 1);

                std::cerr << "eternal darkness awaits…" << endl << std::flush;
            }
            tE = 0;
        }

        t++; tE++; twrite++, tBCreset++, tTimer++;
    }   // end of SAMC main loop

    // calculated positions after SAMC loop
    outputPositions(sp, Chn, hd->dbposi, 1);

    Timer.PrintProgress(t, sp->T_MAX);
    this_thread::sleep_for(chrono::milliseconds(200));
    std::printf("\npraise the sun ☀\n");

    Timer.endProgram();

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
    delete[] contHB;
    delete[] conf_n;
    delete[] conf_write_t;

    delete sp;

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
// return random real number between -1 and 1
double RND()
{
    uniform_real_distribution<double> distribution(-1.0,1.0);
    return distribution(rng);
}
//          XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
//          XXXXXXXXXX  INPUT OUTPUT FUNCTIONS  XXXXXXXXXX

// identify file names from command input
int CommandInitialize(int argc, char *argv[], Header *hd)
{
    std::string opt1, opt2;
    //default values
    hd->confnm = "config_ini.xyz";
    hd->paranm = "Sys_param.dat";
    hd->dbposi = "AnorLondo.xyz";
    hd->lngEnm = "lngE_input.dat";

    //reading arguments
    for( int i=1; i<argc; i+=2 ) {
        opt1 = argv[i];
        opt2 = argv[i+1];
        if( opt1.compare("-c") == 0 ) {
            hd->confnm = opt2; }
        if( opt1.compare("-p") == 0 ) {
            hd->paranm = opt2; }
        if( opt1.compare("-d") == 0 ) {
            hd->dbposi = opt2; }
        if( opt1.compare("-l") == 0 ) {
            hd->lngEnm = opt2; }
    }

    return 0;
}
// Head of colsole output
void ConsoleOutputHead(SysPara *sp)
{
    time_t curtime;
    std::cout << "SAMC PRIME20 simulation" << endl;
    std::cout << "No. of chains " << sp->N_CH << ", chain length " << sp->N_AA << ", sequence " << sp->AA_seq << endl;
    std::cout << "duration T_MAX = " << sp->T_MAX << endl;
    if(!sp->MEASURE)
        std::cout << "SAMC run to apprximate ln[g(E)]." << endl;
    else
        std::cout << "production run to measure geometric observables." << endl;
    time(&curtime);
    std::cout << "start time: " << ctime(&curtime) << endl;
    
}
// creates new chain Chn with N_AA amino acids specified in AA_seq. returns true if successful
bool newChain(SysPara *sp, Chain Chn[], int chnNum)
{
    AmiAc AAinsert;
    double rad, ZrotAng, XrotAng;                           // radius for y-z calc, rotation angle
    double variaMtrx[3];                                    // variable matrix for side chain
    double ZrotMtrx[3][3], XrotMtrx[3][3];                  // rotation matrix around Z and X axis
    double tranVec[3], newpos[3], preRot[3], postRotZ[3];    // translation vector, new position after translation, coord. before rotation, coord. after rotation

    Chn[chnNum].setChnNo(chnNum);
    Chn[chnNum].AmAc.clear();
    for (int i=0; i<sp->N_AA; i++) {                            // converting char sequence to int sequence
        AAinsert.Setup(sp->AA_seq, i);
        Chn[chnNum].AmAc.push_back(AAinsert);
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
    XrotAng = 11.0*M_PI/10.0;
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
        variaMtrx[2] = sqrt( rad*rad - variaMtrx[1]*variaMtrx[1] );

        Chn[chnNum].AmAc[i].Bd[3].setR(variaMtrx[0], variaMtrx[1], variaMtrx[2]);
    }
    // move polymer to positive coordinates = apply periodic boundary conditions
    for(int i = 0; i < sp->N_AA; i++) {
        for(int j = 0; j < 4; j++) {
            for( int k=0; k<3; k++ ) {
                Chn[chnNum].AmAc[i].Bd[j].addBC(k, floor(Chn[chnNum].AmAc[i].Bd[j].getR(k)/sp->L ));
            }
            Chn[chnNum].AmAc[i].Bd[j].setR(Chn[chnNum].AmAc[i].Bd[j].getR(0) - sp->L*floor(Chn[chnNum].AmAc[i].Bd[j].getR(0)/sp->L), Chn[chnNum].AmAc[i].Bd[j].getR(1) - sp->L*floor(Chn[chnNum].AmAc[i].Bd[j].getR(1)/sp->L), Chn[chnNum].AmAc[i].Bd[j].getR(2) - sp->L*floor(Chn[chnNum].AmAc[i].Bd[j].getR(2)/sp->L));
                
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
    stringstream ss_line;
    std::string s_line, option, value;
    double d;

    std::cout << "reading parameter input ..." << std::endl;

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
    int read_Disp= 0;
    int read_DPhi= 0;
    int read_DPsi= 0;
    int read_ETru= 0;
    int read_Meas= 0;
    int read_HBCM= 0;
    int read_WCon= 0;

    ifstr.open(hd->paranm);
    if( ifstr.is_open() ) {
        while( getline( ifstr, s_line ) ) {
            ss_line.clear();
            if( s_line.length() != 0 && s_line.at(0)!='#' && s_line.at(0)!='\n' && s_line.at(0)!='\0' ) {
                ss_line.str(s_line);
                getline(ss_line, option, '=');
                getline(ss_line, value, ';');
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
                else if( option.compare("DISP_MAX")==0 ) {
                    sp->DISP_MAX = stod(value, nullptr); read_Disp = 1; }
                else if( option.compare("DPHI_MAX")==0 ) {
                    sp->DPHI_MAX = stod(value, nullptr); read_DPhi = 1; }
                else if( option.compare("DPSI_MAX")==0 ) {
                    sp->DPSI_MAX = stod(value, nullptr); read_DPsi = 1; }
                else if( option.compare("EBIN_TRUNC_UP")==0 ) {
                    if( value.compare("true")==0 ) { 
                        sp->EBIN_TRUNC_UP = true; read_ETru = 1; }
                    else if( value.compare("false")==0 ) { 
                        sp->EBIN_TRUNC_UP = false; read_ETru = 1; } }
                else if( option.compare("MEASURE")==0 ) {
                    if( value.compare("true")==0 ) { sp->MEASURE = true; read_Meas = 1; }
                    else if( value.compare("false")==0 ) { sp->MEASURE = false; read_Meas = 1; } }
                else if( option.compare("HB_CONTMAT")==0 ) {
                    if( value.compare("true")==0 ) { sp->HB_CONTMAT = true; read_HBCM = 1; }
                    else if( value.compare("false")==0 ) { sp->HB_CONTMAT = false; read_HBCM = 1; } }
                else if( option.compare("WRITE_CONFIG")==0 ) {
                    if( value.compare("true")==0 ) { sp->WRITE_CONFIG = true; read_WCon = 1; }
                    else if( value.compare("false")==0 ) { sp->WRITE_CONFIG = false; read_WCon = 1; } }
            }
        }
        sp->BinW = (sp->EMax - sp->EMin)/(double)sp->NBin;
        sp->nstep = 4*sp->N_AA*sp->N_CH;

        ifstr.close();
        return true;   
    }
    else {
        std::cout << " -- ERROR --    " << hd->paranm << " not found" << std::endl;
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
                Chn[j].AmAc.push_back(AAinsert);
            }
        }
        //getline( input, line );
        input >> NaaFile; input.ignore();
        while( input.good() && !input.eof() ) {
            input.ignore(3);
            input >> x >> y >> z;

            //std::cout << "Chn[" << (i/4)/N_AA << "].AmAc[" << (i/4)%N_AA << "].Bd[" << i%4 << "]\t" << x << "\t" << y << "\t" << z << endl;

            //if(i%N_AA == 0) 
            Chn[(i/4)/sp->N_AA].AmAc[(i/4)%sp->N_AA].Bd[i%4].setR(x, y, z);
            input.ignore();
            i++;
            if(i == NaaFile) { break; }
        }
        if(i != sp->N_CH*sp->N_AA*4) {
            std::cout << "ERROR\tincorrect number of beads: N_BB_real = " << i << "  N_BB_should = " << sp->N_AA*4 << std::endl;
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
bool readPrevRunInput(SysPara *sp, Chain Chn[], string inputFile, double lngE[], long unsigned int H[], long unsigned int &tcont, double &gammasum)
{
    ifstream input;
    stringstream ss_line;
    std::string s_line, s_token;
    char delim = ' ';
    int num;
    double Ebin1, Ebin2;

    input.open(inputFile);
    if( input.is_open() ) {
        std::cout << "reading lngE, H, tcont and gammasum from input file '" << inputFile << "'" << std::endl;
        while( true ) {
            getline( input, s_line);
            if( s_line.length() == 0 ) {
                std::cout << "ERROR: encountered EOF inside input file head" << std::endl;
                return false;
            }
            if(s_line[0] == '#') {
                ss_line.str(s_line);
                getline(ss_line, s_token, delim);
                getline(ss_line, s_token, delim);
                if( s_token.length() == 0 ) continue;
                if( s_token.compare("number") == 0 ) {
                    getline(ss_line,s_token, delim);
                    getline(ss_line,s_token, delim);
                    getline(ss_line,s_token, delim);
                    getline(ss_line,s_token, delim);
                    tcont = stoi(s_token, nullptr);
                    std::cout << "  continue after step " << tcont << std::endl;
                }
                else if( s_token.compare("gammasum") == 0 ) {
                    getline(ss_line, s_token, delim);
                    getline(ss_line, s_token, delim);
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
                input >> num >> Ebin1 >> Ebin2 >> lngE[i] >> H[i];
                if( num != i ) {
                    std::cout << "ERROR: Bad line number " << num << " != " << i << std::endl;
                    input.close();
                    return false;
                }
                std::cout << "  " << i << " " << lngE[i] << std::endl;
            }
            else {
                std::cout << "ERROR: encountered EOF inside lngE line << " << i << " expected " << sp->NBin-1 << std::endl;
                input.close();
                return false;
            }
        }
        input.close();
        return true;
    }

    std::cout << "unable to open input file '" << inputFile << "'" << std::endl;
    return false;
}
// reads lngE data from file
bool extra_lngE(SysPara *sp, Header *hd, double lngE[])
{
    ifstream input;
    std::string s_line;
    int num;
    long unsigned int H;
    double Ebin1, Ebin2;

    input.open(hd->lngEnm);
    if( input.is_open() ) {
        std::cout << "reading lngE from input file '" << hd->lngEnm << "'" << std::endl;
        getline(input, s_line);
        if(s_line[0] != 'b') {
            std::cout << "unexpected character at start of line in file '" << hd->lngEnm << "': expected 'b', instead '" << s_line[0] << std::endl;
            input.close();
            return false;
        }
        for( int i=0; i<sp->NBin; i++ ) {
            if( input.good() ) {
                input >> num >> Ebin1 >> Ebin2 >> lngE[i] >> H;
                if( num != i ) {
                    std::cout << "ERROR: Bad line number " << num << " != " << i << std::endl;
                    input.close();
                    return false;
                }
            }
            else {
                std::cout << "ERROR: encountered EOF inside lngE line << " << i << " expected " << sp->NBin-1 << std::endl;
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
bool outputPositions(SysPara *sp, Chain Chn[], string name, int mode)
{
    ofstream Checkpos;
    if(mode == 0) {
        Checkpos.open(name, ios::out);
    }
    else {
        Checkpos.open(name, ios::app);
    }
    if( Checkpos.is_open() ) {
        Checkpos << sp->N_CH*sp->N_AA*4 << std::endl;

        for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
            Checkpos << std::setw(2) << std::setfill('0') << i << "N " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(0) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(0)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(1) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(1)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(2) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(2)*sp->L << endl;
            Checkpos << std::setw(2) << std::setfill('0') << i << "C " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(0) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(0)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(1) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(1)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(2) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(2)*sp->L << endl;
            Checkpos << std::setw(2) << std::setfill('0') << i << "O " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(0) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(0)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(1) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(1)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(2) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(2)*sp->L << endl;
            Checkpos << std::setw(2) << std::setfill('0') << i << "R " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(0) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(0)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(1) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(1)*sp->L << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(2) + Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(2)*sp->L << endl;
        }

        Checkpos << std::endl;
        Checkpos.close();
        return true;
    }
    else {
        std::cerr << "error opening " << name << endl;
        return false;
    }
}
// writes backup file in SAMC run
bool BackupSAMCrun(SysPara *sp, Chain Chn[], Timer &Timer, unsigned long int t, double gammasum, double gamma, unsigned long naccept[], unsigned long nattempt[], double lngE[], unsigned long H[], double E)
{
    std::ostringstream oss;
    oss << "SAMCbackup_" << t/sp->T_WRITE << ".dat";
    string name = oss.str();
    ofstream backup;
    backup.open(name);
    if( backup.is_open() ) {
        backup << "# SAMC DoS of " << sp->N_CH << " PRIME20 " << sp->N_AA << "-mer(s)" << endl;
        backup << "# sequence " << sp->AA_seq << endl;
        backup << "# number of MC steps: " << t+1 << ", with " << sp->nstep << " moves per step" << endl;
        backup << "# gammasum = " << gammasum << " (current gamma = " << gamma << ")" << endl;
        backup << "# gamma_0 = " << sp->GAMMA_0 << ", constant until T_0 = " << sp->T_0 << endl;
        backup << "# current runtime: " << Timer.curRunTime() << endl;
        backup << "# energy window: [" << sp->EMin << ";" << sp->EMax << "] in " << sp->NBin << " steps (bin width = " << sp->BinW << ")" << endl;
        backup << "# accepted " << naccept[0] << " of " << nattempt[0] << " (" << 100*(double)naccept[0]/(double)nattempt[0] << "%) local moves" << endl;
        backup << "# accepted " << naccept[1] << " of " << nattempt[1] << " (" << 100*(double)naccept[1]/(double)nattempt[1] << "%) pivot (Phi) moves" << endl;
        backup << "# accepted " << naccept[2] << " of " << nattempt[2] << " (" << 100*(double)naccept[2]/(double)nattempt[2] << "%) pivot (Psi) moves" << endl;
        backup << "# accepted " << naccept[3] << " of " << nattempt[3] << " (" << 100*(double)naccept[3]/(double)nattempt[3] << "%) translation moves" << endl;
        backup << "bin  from  to  lng  H" << endl;
        for( int i=0; i<sp->NBin; i++ ) {
            backup << i << " " << sp->EMin+i*sp->BinW << " " << sp->EMin+(i+1)*sp->BinW << " " << lngE[i] << " " << H[i] << endl;
        }
        backup << "# current configuration ( E=" << E << ")" << endl;
        backup << "beadID  x  y  z  BCx  BCy  BCz" << std::endl;
        for( int i=0; i<sp->N_CH*sp->N_AA; i++ ) {
            backup << std::setw(2) << std::setfill('0') << i << "N " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[0].getBC(2) << endl;
            backup << std::setw(2) << std::setfill('0') << i << "C " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[1].getBC(2) << endl;
            backup << std::setw(2) << std::setfill('0') << i << "O " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[2].getBC(2) << endl;
            backup << std::setw(2) << std::setfill('0') << i << "R " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getR(2) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(0) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(1) << " " << Chn[i/sp->N_AA].AmAc[i%sp->N_AA].Bd[3].getBC(2) << endl;
        }
        backup.close();
        return true;
    }
    else {
        std::cout << endl << "error opening " << name << std::endl;
        return false;
    }
}
bool BackupProdRun(SysPara *sp, Timer &Timer, unsigned long int t, unsigned long int H[])
{
    std::ostringstream oss;
    oss << "results_" << t/sp->T_WRITE << ".dat";
    string name = oss.str();
    ofstream results;
    results.open(name);
    if( results.is_open() ) {
        results << "# Observables for " << sp->N_CH << " PRIME20 " << sp->N_AA << "-mer(s)" << std::endl;
        results << "# sequence " << sp->AA_seq << std::endl;
        results << "# number of MC steps: " << t+1 << ", with " << sp->nstep << " moves per step" << std::endl;
        results << "# current runtime: " << Timer.curRunTime() << std::endl;
        results << "# energy window: [" << sp->EMin << ";" << sp->EMax << "] in " << sp->NBin << " steps (bin width = " << sp->BinW << ")" << std::endl;
        results << "# measured observables are:" << std::endl;
        results << "\tvisit histogram H" << std::endl;
        if(sp->HB_CONTMAT)   { results << "#\tHB contact matrices (file HBmat.dat)" << std::endl; }
        if(sp->WRITE_CONFIG) { results << "#\tconfiguration snapshots (file .xyz)" << std::endl; }
        results << "Bin from to H" << std::endl;
        for( int i=0; i<sp->NBin; i++ ) {
            results << i << " " << sp->EMin+i*sp->BinW << " " <<sp->EMin+(i+1)*sp->BinW << " " << H[i] << std::endl;
        }
        results.close();
        return true;
    }
    else {
        std::cout << endl << "error opening " << name << std::endl;
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

    if( EOswitch == 1 && j1 != 3 ) { std::cout << "ERROR - EO_SegBead() energy calculation for j1=" << j1 << endl; return -1; }

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
                        //std::cerr << "overlapp: C" << h1 << "A" << i1 << "B" << j1 << " - C" << index/(4*N_AA) << "A" << (index/4)%N_AA << "B" << index%4 << std::endl;
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
    uniform_real_distribution<double> distribution(0,1);

    // biased moves? logEo -= gammasum*pio → ask Arne and Timur
    double MCrand = distribution(rng);
    //std::cout << endl << "lngEold=" << lngEold << " lngEnew=" << lngEnew << " exp(lngEnew-lngEold)=" << exp(lngEold-lngEnew) << " MCrand=" << MCrand;
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
        disp[k] = sp->DISP_MAX * RND();
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

    angle = sypa->DPHI_MAX*RND();
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
        printf("!ERROR! \trotPhi(): int high has unacceptable value!\n");
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

    angle = sypa->DPSI_MAX*RND();
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
        printf("!ERROR! \trotPsi(): int high has unacceptable value!\n");
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
    Bead Bdcpy[4*sp->N_AA];                 // copy of beads in moved chain
    double dVec[3], newpos[3];          // displacement vector, new position (before PBC)
    double Eold, dEhb, distabs;
    int BrokenHB[sp->N_AA][2], nBHB;

    for( int i=0; i<sp->N_AA; i++ ) {
        for( int j=0; j<4; j++ ) {
            Bdcpy[i*4+j] = Chn[iChn].AmAc[i].Bd[j];
        }
    }
    Eold = EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, 0, iChn*sp->N_AA, 1) + EO_SegSeg(sp, Chn, iChn*sp->N_AA, (iChn+1)*sp->N_AA, (iChn+1)*sp->N_AA, sp->N_CH*sp->N_AA, 1);
    for( int i=0; i<3; i++ ) {
        dVec[i] = RND()*(double)sp->L/2.0;
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
                        Chn[iChn].AmAc[m].Bd[n] = Bdcpy[m*4+n];     // reset chain & exit
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