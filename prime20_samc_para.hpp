/* FILE: prime20_samc_para.h
 * 
 *        :'''''''''''''''''''''''''''''''''''''''''''':
 *        :  Topological Parameters of Prime20 HEADER  :
 *        :............................................:
 * 
 */
#ifndef PRIME20_SAMC_PARA_HPP
#define PRIME20_SAMC_PARA_HPP

#include "prime20_samc_chain.hpp"

//-------------------------------------------------------------
//               >>  TOPOLOGICAL PARAMETERS  <<

//Backbone Bond lengths
const double BND_NCa = 1.460;
const double BND_CaC = 1.510;
const double BND_CN = 1.330;
const double BND_FLUCT = 0.02375;	// Allowed bond length fluctuation
//Backbone Pseudobonds
const double PBND_NC = 2.45;
const double PBND_CaN = 2.41;
const double PBND_CCa = 2.45;
const double PBND_CaCa = 3.80;
//Bond angles 
const double ANGL_NCaC = 1.937;     // 111.0 degrees
const double ANGL_CaCN = 2.025;	    // 116.0 deg
const double ANGL_CNCa = 2.129;	    // 122.0 deg
//Backbone bead diameters
const double DIA_N = 3.3;				
const double DIA_Ca = 3.7;
const double DIA_C = 4.0;
const double DIA_HUGE = 4.00;       // largest possible bead diameter
const double DIA2HUGE = 16.00;      // largest possible squared bead diameter
//Squared backbone bead diameters during HB stabilization (Nguyen2004, p.2921)
const double DIA2_HB_NCa = 25.00;
const double DIA2_HB_NN = 22.4676;
const double DIA2_HB_CCa = 23.6196;
const double DIA2_HB_CC = 23.3289;
//Backbone bead square well diameters
const double SW_BB = 4.5;			// Backbone SW diameter (N and C)
const double SW2_BB = 20.25;	    // Squared backbone SW diameter (N and C)
const double SW_HUGE = 7.4;         // Largest possible SW diameter
const double SW2HUGE = 54.76;		// Largest possible squared SW diameter
//Square-well parameter values (19 parameter set used. Some values are equal)
const double SWPAR2 = -0.200;
const double SWPAR3 = -0.203;
const double SWPAR4 = -0.210;
const double SWPAR5 = -0.205;
const double SWPAR6 = -0.201;
const double SWPAR7 = -0.084;
const double SWPAR8 = -0.585;
const double SWPAR9 = 0.253;
const double SWPAR10 = -0.136;
const double SWPAR11 = 0.073;
const double SWPAR12 = -0.148;
const double SWPAR13 = -0.139;
const double SWPAR14 = -0.116;	//Parameter 15/19, but 14/23
const double SWPAR15 = 0.015;		//Parameter 14/19, but 15/23
const double SWPAR16 = 0.015;		//Equals parameter 15/23
const double SWPAR17 = -0.116;	//Equals parameter 14/23
const double SWPAR18 = -0.086;	//Parameter 16/19
const double SWPAR19 = -0.086;	//Equals parameter 18/23
const double SWPAR20 = 0.074;		//Parameter 17/19
const double SWPAR21 = 0.074;		//Equals parameter 20/23
const double SWPAR22 = -0.080;	//Parameter 18/19
//Double-well SW parameters
//.....
//Squeeze parameters
const double SQZ1 = 4.40286;	//1.1436*(Ca+CO)/2
const double SQZ2 = 3.08;	    //0.88*(Ca+N)/2
const double SQZ3 = 3.2057585;	//0.87829*(CO+N)/2
const double SQZ4 = 2.64;	    //0.8*N
const double SQZ5 = 3.0852;	    //0.7713*CO

//Masses of beads. Not an original PRIME20 parameter!
const double MASS_N = 15.0;		//Backbone NH
const double MASS_C = 13.0;		//Backbone C(alpha)H
const double MASS_O = 28.0;		//Backbone CO
const double MASS_R(char AAnm);
//-------------------------------------------------------------
//            >>  FUNCTIONS PRODUCING PARAMETERS  <<

// square well diameter for chn[i] <-> chn[j]
const double SWDia(AmiAc AmAci, AmiAc AmAcj);

// square well depth for chn[i] <-> chn[j]
const double SWDepth(AmiAc AmAci, AmiAc AmAcj);

// side chain bead diameter for chn[i] <-> chn[j]
const double DiaSC(AmiAc AmAci, AmiAc AmAcj);

/* squared bead pair HC diameters, including reduced diameters for neighbouring residues
 * comparing chn[ha].AmiAc[ia].Bd[ja] with chn[hb].AmiAc[ib].Bd[jb]
 */
double DiaSQ(Chain chn[], int ha, int ia, int ja, int hb, int ib, int jb);

#endif