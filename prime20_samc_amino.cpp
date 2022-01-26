/* FILE: prime20_samc_amino.cpp
 *
 *          :'''''''''''''''''''''''':
 *          :  Amino acid class CPP  :
 *          :........................:
 * 
 */

#include "prime20_samc_amino.hpp"
#include "prime20_samc_para.hpp"

using namespace std;
//-------------------------------------------------------
// convert Amino Acid char to integer
const int AA2INT(string chain, int i)
{
	char c = chain.at(i);
	if(c=='G') return 0;
	if(c=='A') return 1;
	if(c=='V') return 2;
	if(c=='P') return 3;
	if(c=='T') return 4;
	if(c=='S') return 5;
	if(c=='N') return 6;
	if(c=='D') return 7;
	if(c=='R') return 8;
	if(c=='K') return 9;
	if(c=='E') return 10;
	if(c=='Q') return 11;
	if(c=='L') return 12;
	if(c=='I') return 13;
	if(c=='F') return 14;
	if(c=='Y') return 15;
	if(c=='W') return 16;
	if(c=='M') return 17;
	if(c=='C') return 18;
	if(c=='H') return 19;
	printf("invalid char in amino acid sequence: %c\n", c); return -1;
}
// sets up amino acid type & bead properties
void AmiAc::Setup(string chain, int i)
{
    // converting char sequence to int sequence
	AA_alp = chain.at(i);
	AA_num = AA2INT(chain, i);
	// side chain (pseudo-)bond lengths (CaR. NR, CR)
    switch(AA_num) {
		case 0: 
			BND_CaR = 0.0;
			PBND_NR = BND_NCa;
			PBND_CR = BND_CaC;
			ANGL_RCaN = 0.0;
			ANGL_RCaC = 0.0;
			SQZ6 = 2.600;
			SQZ7 = 2.610;
			SQZ8 = 4.013;
			SQZ9 = 4.064;
			SQZ10= 4.904;
			break;
		case 1: 
			BND_CaR = 1.600;
			PBND_NR = 2.500;
			PBND_CR = 2.560;
			ANGL_RCaN = 1.911;
			ANGL_RCaC = 1.934;
			SQZ6 = 3.312;
			SQZ7 = 3.000;
			SQZ8 = 4.353;
			SQZ9 = 4.598;
			SQZ10= 4.997;
			break;
		case 2: 
			BND_CaR = 2.002;
			PBND_NR = 2.775;
			PBND_CR = 2.959;
			ANGL_RCaN = 1.841;
			ANGL_RCaC = 1.991;
			SQZ6 = 3.570;
			SQZ7 = 3.378;
			SQZ8 = 4.635;
			SQZ9 = 4.754;
			SQZ10= 5.002;
			break;
		case 3: 
			BND_CaR = 1.926;
			PBND_NR = 1.851;
			PBND_CR = 2.995;
			ANGL_RCaN = 1.127;
			ANGL_RCaC = 2.108;
			SQZ6 = 3.133;
			SQZ7 = 3.298;
			SQZ8 = 4.665;
			SQZ9 = 3.884;
			SQZ10= 4.773;
			break;
		case 4: 
			BND_CaR = 1.981;
			PBND_NR = 2.650;
			PBND_CR = 2.900;
			ANGL_RCaN = 1.739;
			ANGL_RCaC = 1.948;
			SQZ6 = 3.447;
			SQZ7 = 3.290;
			SQZ8 = 4.573;
			SQZ9 = 4.617;
			SQZ10= 5.007;
			break;
		case 5: 
			BND_CaR = 1.967;
			PBND_NR = 2.650;
			PBND_CR = 2.800;
			ANGL_RCaN = 1.749;
			ANGL_RCaC = 1.859;
			SQZ6 = 3.331;
			SQZ7 = 3.128;
			SQZ8 = 4.380;
			SQZ9 = 4.507;
			SQZ10= 4.944;
			break;
		case 6: 
			BND_CaR = 2.510;
			PBND_NR = 3.050;
			PBND_CR = 3.350;
			ANGL_RCaN = 1.690;
			ANGL_RCaC = 1.927;
			SQZ6 = 3.607;
			SQZ7 = 3.565;
			SQZ8 = 4.680;
			SQZ9 = 4.633;
			SQZ10= 4.791;
			break;
		case 7: 
			BND_CaR = 2.500;
			PBND_NR = 3.100;
			PBND_CR = 3.250;
			ANGL_RCaN = 1.740;
			ANGL_RCaC = 1.843;
			SQZ6 = 3.751;
			SQZ7 = 3.435;
			SQZ8 = 4.558;
			SQZ9 = 4.785;
			SQZ10= 4.860;
			break;
		case 8: 
			BND_CaR = 4.200;
			PBND_NR = 4.500;
			PBND_CR = 4.800;
			ANGL_RCaN = 1.610;
			ANGL_RCaC = 1.819;
			SQZ6 = 4.827;
			SQZ7 = 4.651;
			SQZ8 = 5.535;
			SQZ9 = 5.703;
			SQZ10= 4.978;
			break;
		case 9: 
			BND_CaR = 3.550;
			PBND_NR = 4.050;
			PBND_CR = 4.250;
			ANGL_RCaN = 1.732;
			ANGL_RCaC = 1.873;
			SQZ6 = 4.384;
			SQZ7 = 4.191;
			SQZ8 = 5.163;
			SQZ9 = 5.323;
			SQZ10= 4.974;
			break;
		case 10: 
			BND_CaR = 3.180;
			PBND_NR = 3.780;
			PBND_CR = 3.930;
			ANGL_RCaN = 1.792;
			ANGL_RCaC = 1.894;
			SQZ6 = 4.175;
			SQZ7 = 3.997;
			SQZ8 = 5.074;
			SQZ9 = 5.162;
			SQZ10= 4.996;
			break;
		case 11: 
			BND_CaR = 3.300;
			PBND_NR = 3.750;
			PBND_CR = 4.000;
			ANGL_RCaN = 1.679;
			ANGL_RCaC = 1.859;
			SQZ6 = 4.139;
			SQZ7 = 3.996;
			SQZ8 = 5.062;
			SQZ9 = 5.134;
			SQZ10= 5.000;
			break;
		case 12: 
			BND_CaR = 2.625;
			PBND_NR = 3.290;
			PBND_CR = 3.500;
			ANGL_RCaN = 1.808;
			ANGL_RCaC = 1.970;
			SQZ6 = 3.918;
			SQZ7 = 3.724;
			SQZ8 = 4.863;
			SQZ9 = 4.936;
			SQZ10= 5.001;
			break;
		case 13: 
			BND_CaR = 2.400;
			PBND_NR = 3.050;
			PBND_CR = 3.300;
			ANGL_RCaN = 1.773;
			ANGL_RCaC = 1.976;
			SQZ6 = 3.740;
			SQZ7 = 3.626;
			SQZ8 = 4.867;
			SQZ9 = 4.867;
			SQZ10= 4.994;
			break;
		case 14: 
			BND_CaR = 3.425;
			PBND_NR = 3.650;
			PBND_CR = 4.050;
			ANGL_RCaN = 1.517;
			ANGL_RCaC = 1.805;
			SQZ6 = 3.991;
			SQZ7 = 3.973;
			SQZ8 = 4.827;
			SQZ9 = 4.780;
			SQZ10= 5.040;
			break;
		case 15: 
			BND_CaR = 3.843;
			PBND_NR = 4.050;
			PBND_CR = 4.300;
			ANGL_RCaN = 1.526;
			ANGL_RCaC = 1.695;
			SQZ6 = 4.208;
			SQZ7 = 4.246;
			SQZ8 = 4.978;
			SQZ9 = 4.898;
			SQZ10= 5.042;
			break;
		case 16: 
			BND_CaR = 3.881;
			PBND_NR = 4.100;
			PBND_CR = 4.350;
			ANGL_RCaN = 1.537;
			ANGL_RCaC = 1.706;
			SQZ6 = 4.460;
			SQZ7 = 4.187;
			SQZ8 = 4.963;
			SQZ9 = 5.180;
			SQZ10= 4.986;
			break;
		case 17: 
			BND_CaR = 3.400;
			PBND_NR = 3.800;
			PBND_CR = 4.050;
			ANGL_RCaN = 1.646;
			ANGL_RCaC = 1.824;
			SQZ6 = 4.205;
			SQZ7 = 4.032;
			SQZ8 = 5.067;
			SQZ9 = 5.206;
			SQZ10= 5.017;
			break;
		case 18: 
			BND_CaR = 2.350;
			PBND_NR = 2.800;
			PBND_CR = 3.100;
			ANGL_RCaN = 1.598;
			ANGL_RCaC = 1.829;
			SQZ6 = 3.516;
			SQZ7 = 3.350;
			SQZ8 = 4.501;
			SQZ9 = 4.560;
			SQZ10= 4.913;
			break;
		case 19: 
			BND_CaR = 3.160;
			PBND_NR = 3.450;
			PBND_CR = 3.830;
			ANGL_RCaN = 1.548;
			ANGL_RCaC = 1.826;
			SQZ6 = 3.886;
			SQZ7 = 3.838;
			SQZ8 = 4.790;
			SQZ9 = 4.766;
			SQZ10= 4.945;
			break;
	}
}
// access AA_num
int AmiAc::getAAnum()
{ 
	return AA_num; 
}
// access AA_alp
char AmiAc::get_AAalp()
{
	return AA_alp;
}
// access distances to side chain
double AmiAc::getDisR(int i)
{
	if (i == 0) {
		return PBND_NR;
	}
	if (i == 1) {
		return BND_CaR;
	}
	if (i == 2) {
		return PBND_CR;
	}
}
// access angles
double AmiAc::getAng(int i)
{
	if (i == 0) {
		return ANGL_RCaN;
	}
	if (i == 1) {
		return ANGL_RCaC;
	}
}
// Squeeze parameters
double AmiAc::getSQZ(int i)
{
	switch(i) {
		case 1:  return SQZ1;
		case 2:  return SQZ2;
		case 3:  return SQZ3;
		case 4:  return SQZ4;
		case 5:  return SQZ5;
		case 6:  return SQZ6;
		case 7:  return SQZ7;
		case 8:  return SQZ8;
		case 9:  return SQZ9;
		case 10: return SQZ10;
	}

}
// constructor
AmiAc::AmiAc(void)
{
}
// destructor
AmiAc::~AmiAc(void)
{
}