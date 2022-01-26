/* FILE: prime20_samc_para.cpp
 * 
 *        :'''''''''''''''''''''''''''':
 *        :  Parameters of Prime20 CPP :
 *        :............................:
 * 
 */
#include "prime20_samc_para.hpp"

using namespace std;
//-------------------------------------------------------------
// mass of side chain beads
const double MASS_R(char AAnm)
{
	if(AAnm=='G') return 1.0;		//G side chain: H
	if(AAnm=='A') return 15.0;		//A side chain: CH3
	if(AAnm=='V') return 43.0;		//V side chain: C3H7
	if(AAnm=='P') return 42.0;		//P side chain: C3H6
	if(AAnm=='T') return 45.0;		//T side chain: C2OH5
	if(AAnm=='S') return 31.0;		//S side chain: COH3
	if(AAnm=='N') return 58.0;		//N side chain: C2ONH4
	if(AAnm=='D') return 58.0;		//D side chain: C2O2H2
	if(AAnm=='R') return 101.0;		//R side chain: C4N3H11
	if(AAnm=='K') return 73.0;		//K side chain: C4NH11
	if(AAnm=='E')  return 72.0;		//E side chain: C3O2H4
	if(AAnm=='Q')  return 72.0;		//Q side chain: C3ONH6
	if(AAnm=='L')  return 57.0;		//L side chain: C4H9
	if(AAnm=='I')  return 57.0;		//I side chain: C4H9
	if(AAnm=='F')  return 91.0;		//F side chain: C7H7
	if(AAnm=='Y')  return 107.0;	//Y side chain: C7OH7
	if(AAnm=='W')  return 130.0;	//W side chain: C9NH8
	if(AAnm=='M')  return 75.0;		//M side chain: C3SH7
	if(AAnm=='C')  return 47.0;		//C side chain: CSH3
	if(AAnm=='H')  return 81.0;		//H side chain: C4N2H5
	std::cout << "invalid char in amino acid sequence: " << AAnm << std::endl; return -1;
}
// square well diameter for chn[i] <-> chn[j]
const double SWDia(AmiAc AmAci, AmiAc AmAcj)
{
    int a,b;
    a = AmAci.getAAnum();
    b = AmAcj.getAAnum();
    switch((int)min(a,b)) {
		case 0: return 0.0;
		case 1: switch((int)max(a,b)) {
				case 1: return 5.4;
				case 2: return 6.1;
				case 3: return 6.2;
				case 4: return 6.2;
				case 5: return 5.9;
				case 6: return 5.6;
				case 7: return 5.6;
				case 8: return 6.1;
				case 9: return 6.0;
				case 10: return 5.9;
				case 11: return 5.8;
				case 12: return 5.6;
				case 13: return 5.7;
				case 14: return 5.9;
				case 15: return 5.7;
				case 16: return 5.5;
				case 17: return 5.8;
				case 18: return 5.9;
				case 19: return 5.5;
			}
		case 2: switch((int)max(a,b)) {
				case 2: return 6.3;
				case 3: return 6.3;
				case 4: return 6.4;
				case 5: return 6.2;
				case 6: return 6.3;
				case 7: return 6.3;
				case 8: return 6.8;
				case 9: return 6.6;
				case 10: return 6.5;
				case 11: return 6.5;
				case 12: return 6.2;
				case 13: return 6.4;
				case 14: return 6.5;
				case 15: return 6.5;
				case 16: return 6.6;
				case 17: return 6.4;
				case 18: return 6.0;
				case 19: return 6.2;
			}
		case 3: switch((int)max(a,b)) {
				case 3: return 6.5;
				case 4: return 6.6;
				case 5: return 6.1;
				case 6: return 6.2;
				case 7: return 6.3;
				case 8: return 6.8;
				case 9: return 6.7;
				case 10: return 6.4;
				case 11: return 6.5;
				case 12: return 6.3;
				case 13: return 6.4;
				case 14: return 6.5;
				case 15: return 6.4;
				case 16: return 6.3;
				case 17: return 6.2;
				case 18: return 6.0;
				case 19: return 6.3;
			}
		case 4: switch((int)max(a,b)) {
				case 4: return 6.5;
				case 5: return 6.0;
				case 6: return 6.3;
				case 7: return 6.2;
				case 8: return 6.8;
				case 9: return 6.5;
				case 10: return 6.4;
				case 11: return 6.4;
				case 12: return 6.2;
				case 13: return 6.4;
				case 14: return 6.6;
				case 15: return 6.4;
				case 16: return 6.5;
				case 17: return 6.4;
				case 18: return 6.1;
				case 19: return 6.3;
			}
		case 5: switch((int)max(a,b)) {
				case 5: return 6.4;
				case 6: return 6.2;
				case 7: return 6.1;
				case 8: return 6.3;
				case 9: return 6.1;
				case 10: return 6.0;
				case 11: return 6.0;
				case 12: return 6.3;
				case 13: return 6.4;
				case 14: return 6.2;
				case 15: return 6.5;
				case 16: return 6.3;
				case 17: return 6.4;
				case 18: return 6.3;
				case 19: return 6.3;
			}
		case 6: switch((int)max(a,b)) {
				case 6: return 6.3;
				case 7: return 6.5;
				case 8: return 6.6;
				case 9: return 6.5;
				case 10: return 6.4;
				case 11: return 6.4;
				case 12: return 6.4;
				case 13: return 6.6;
				case 14: return 6.5;
				case 15: return 6.7;
				case 16: return 6.9;
				case 17: return 6.4;
				case 18: return 6.2;
				case 19: return 6.5;
			}
		case 7: switch((int)max(a,b)) {
				case 7: return 6.5;
				case 8: return 6.5;
				case 9: return 6.3;
				case 10: return 6.6;
				case 11: return 6.3;
				case 12: return 6.5;
				case 13: return 6.5;
				case 14: return 6.7;
				case 15: return 6.9;
				case 16: return 6.9;
				case 17: return 6.7;
				case 18: return 6.2;
				case 19: return 6.6;
			}
		case 8: switch((int)max(a,b)) {
				case 8: return 7.2;
				case 9: return 6.8;
				case 10: return 6.6;
				case 11: return 6.9;
				case 12: return 6.8;
				case 13: return 6.7;
				case 14: return 6.9;
				case 15: return 7.0;
				case 16: return 6.9;
				case 17: return 6.6;
				case 18: return 6.3;
				case 19: return 6.9;
			}
		case 9: switch((int)max(a,b)) {
				case 9: return 6.9;
				case 10: return 6.4;
				case 11: return 6.7;
				case 12: return 6.5;
				case 13: return 6.7;
				case 14: return 6.9;
				case 15: return 6.7;
				case 16: return 6.5;
				case 17: return 6.4;
				case 18: return 6.4;
				case 19: return 6.6;
			}
		case 10: switch((int)max(a,b)) {
				case 10: return 6.7;
				case 11: return 6.6;
				case 12: return 6.4;
				case 13: return 6.6;
				case 14: return 6.8;
				case 15: return 6.8;
				case 16: return 6.9;
				case 17: return 6.4;
				case 18: return 6.1;
				case 19: return 6.4;
			}
		case 11: switch((int)max(a,b)) {
				case 11: return 6.6;
				case 12: return 6.3;
				case 13: return 6.6;
				case 14: return 6.6;
				case 15: return 6.7;
				case 16: return 6.7;
				case 17: return 6.4;
				case 18: return 6.1;
				case 19: return 6.6;
			}
		case 12: switch((int)max(a,b)) {
				case 12: return 6.4;
				case 13: return 6.5;
				case 14: return 6.6;
				case 15: return 6.7;
				case 16: return 6.9;
				case 17: return 6.5;
				case 18: return 6.1;
				case 19: return 6.5;
			}
		case 13: switch((int)max(a,b)) {
				case 13: return 6.6;
				case 14: return 6.6;
				case 15: return 6.8;
				case 16: return 6.8;
				case 17: return 6.7;
				case 18: return 6.4;
				case 19: return 6.6;
			}
		case 14: switch((int)max(a,b)) {
				case 14: return 6.8;
				case 15: return 6.8;
				case 16: return 7.0;
				case 17: return 6.5;
				case 18: return 6.4;
				case 19: return 6.5;
			}
		case 15: switch((int)max(a,b)) {
				case 15: return 7.0;
				case 16: return 7.0;
				case 17: return 6.6;
				case 18: return 6.5;
				case 19: return 6.9;
			}
		case 16: switch((int)max(a,b)) {
				case 16: return 7.4;
				case 17: return 7.0;
				case 18: return 6.4;
				case 19: return 7.1;
			}
		case 17: switch((int)max(a,b)) {
				case 17: return 6.7;
				case 18: return 6.3;
				case 19: return 6.5;
			}
		case 18: switch((int)max(a,b)) {
				case 18: return 6.2;
				case 19: return 6.2;
			}
		case 19: return 6.7;
	}
}
// square well depth for chn[i] <-> chn[j]
const double SWDepth(AmiAc AmAci, AmiAc AmAcj)
{
	int a,b;
    a = AmAci.getAAnum();
    b = AmAcj.getAAnum();
    switch ((int)min(a,b)) {
		case 0: return 0.0;
		case 1: switch ((int)max(a,b)) {
			case 1: return SWPAR7;
			case 2: return SWPAR12;
			case 3: return SWPAR20;
			case 4: return SWPAR20;
			case 5: return SWPAR20;
			case 6: return SWPAR20;
			case 7: return SWPAR20;
			case 8: return SWPAR20;
			case 9: return SWPAR20;
			case 10: return SWPAR20;
			case 11: return SWPAR20;
			case 12: return SWPAR12;
			case 13: return SWPAR12;
			case 14: return SWPAR12;
			case 15: return SWPAR12;
			case 16: return SWPAR12;
			case 17: return SWPAR12;
			case 18: return SWPAR13;
			case 19: return SWPAR20;
		}
		case 2: switch ((int)max(a,b)) {
			case 2: return SWPAR2;
			case 3: return SWPAR15;
			case 4: return SWPAR15;
			case 5: return SWPAR15;
			case 6: return SWPAR16;
			case 7: return SWPAR15;
			case 8: return SWPAR16;
			case 9: return SWPAR16;
			case 10: return SWPAR15;
			case 11: return SWPAR16;
			case 12: return SWPAR2;
			case 13: return SWPAR2;
			case 14: return SWPAR3;
			case 15: return SWPAR3;
			case 16: return SWPAR3;
			case 17: return SWPAR2;
			case 18: return SWPAR13;
			case 19: return SWPAR16;
		}
		case 3: switch((int)max(a,b)) {
			case 3: return SWPAR21;
			case 4: return SWPAR21;
			case 5: return SWPAR21;
			case 6: return SWPAR21;
			case 7: return SWPAR21;
			case 8: return SWPAR21;
			case 9: return SWPAR21;
			case 10: return SWPAR21;
			case 11: return SWPAR21;
			case 12: return SWPAR15;
			case 13: return SWPAR15;
			case 14: return SWPAR15;
			case 15: return SWPAR15;
			case 16: return SWPAR15;
			case 17: return SWPAR15;
			case 18: return SWPAR15;
			case 19: return SWPAR21;
		}
		case 4: switch((int)max(a,b)) {
			case 4: return SWPAR18;
			case 5: return SWPAR18;
			case 6: return SWPAR18;
			case 7: return SWPAR18;
			case 8: return SWPAR18;
			case 9: return SWPAR18;
			case 10: return SWPAR18;
			case 11: return SWPAR18;
			case 12: return SWPAR15;
			case 13: return SWPAR15;
			case 14: return SWPAR15;
			case 15: return SWPAR18;
			case 16: return SWPAR18;
			case 17: return SWPAR17;
			case 18: return SWPAR17;
			case 19: return SWPAR18;
		}
		case 5: switch((int)max(a,b)) {
			case 5: return SWPAR18;
			case 6: return SWPAR18;
			case 7: return SWPAR18;
			case 8: return SWPAR18;
			case 9: return SWPAR18;
			case 10: return SWPAR18;
			case 11: return SWPAR18;
			case 12: return SWPAR15;
			case 13: return SWPAR15;
			case 14: return SWPAR15;
			case 15: return SWPAR18;
			case 16: return SWPAR18;
			case 17: return SWPAR17;
			case 18: return SWPAR17;
			case 19: return SWPAR18;
		}
		case 6: switch((int)max(a,b)) {
			case 6: return SWPAR22;
			case 7: return SWPAR18;
			case 8: return SWPAR18;
			case 9: return SWPAR18;
			case 10: return SWPAR18;
			case 11: return SWPAR22;
			case 12: return SWPAR16;
			case 13: return SWPAR16;
			case 14: return SWPAR16;
			case 15: return SWPAR19;
			case 16: return SWPAR19;
			case 17: return SWPAR17;
			case 18: return SWPAR17;
			case 19: return SWPAR22;
		}
		case 7: switch((int)max(a,b)) {
			case 7: return SWPAR9;
			case 8: return SWPAR10;
			case 9: return SWPAR10;
			case 10: return SWPAR9;
			case 11: return SWPAR18;
			case 12: return SWPAR15;
			case 13: return SWPAR15;
			case 14: return SWPAR15;
			case 15: return SWPAR18;
			case 16: return SWPAR18;
			case 17: return SWPAR15;
			case 18: return SWPAR17;
			case 19: return SWPAR18;
		}
		case 8: switch((int)max(a,b)) {
			case 8: return SWPAR11;
			case 9: return SWPAR11;
			case 10: return SWPAR10;
			case 11: return SWPAR18;
			case 12: return SWPAR16;
			case 13: return SWPAR16;
			case 14: return SWPAR16;
			case 15: return SWPAR19;
			case 16: return SWPAR16;
			case 17: return SWPAR17;
			case 18: return SWPAR17;
			case 19: return SWPAR18;
		}
		case 9: switch((int)max(a,b)) {
			case 9: return SWPAR11;
			case 10: return SWPAR10;
			case 11: return SWPAR18;
			case 12: return SWPAR16;
			case 13: return SWPAR16;
			case 14: return SWPAR16;
			case 15: return SWPAR19;
			case 16: return SWPAR16;
			case 17: return SWPAR17;
			case 18: return SWPAR17;
			case 19: return SWPAR18;
		}
		case 10: switch((int)max(a,b)) {
			case 10: return SWPAR9;
			case 11: return SWPAR18;
			case 12: return SWPAR15;
			case 13: return SWPAR15;
			case 14: return SWPAR15;
			case 15: return SWPAR18;
			case 16: return SWPAR18;
			case 17: return SWPAR15;
			case 18: return SWPAR17;
			case 19: return SWPAR18;
		}
		case 11: switch((int)max(a,b)) {
			case 11: return SWPAR22;
			case 12: return SWPAR16;
			case 13: return SWPAR16;
			case 14: return SWPAR16;
			case 15: return SWPAR19;
			case 16: return SWPAR19;
			case 17: return SWPAR17;
			case 18: return SWPAR17;
			case 19: return SWPAR22;
		}
		case 12: switch((int)max(a,b)) {
			case 12: return SWPAR2;
			case 13: return SWPAR2;
			case 14: return SWPAR3;
			case 15: return SWPAR3;
			case 16: return SWPAR3;
			case 17: return SWPAR2;
			case 18: return SWPAR13;
			case 19: return SWPAR16;
		}
		case 13: switch((int)max(a,b)) {
			case 13: return SWPAR2;
			case 14: return SWPAR3;
			case 15: return SWPAR3;
			case 16: return SWPAR3;
			case 17: return SWPAR2;
			case 18: return SWPAR13;
			case 19: return SWPAR16;
		}
		case 14: switch((int)max(a,b)) {
			case 14: return SWPAR5;
			case 15: return SWPAR5;
			case 16: return SWPAR5;
			case 17: return SWPAR3;
			case 18: return SWPAR13;
			case 19: return SWPAR16;
		}
		case 15: switch((int)max(a,b)) {
			case 15: return SWPAR6;
			case 16: return SWPAR6;
			case 17: return SWPAR4;
			case 18: return SWPAR14;
			case 19: return SWPAR19;
		}
		case 16: switch((int)max(a,b)) {
			case 16: return SWPAR5;
			case 17: return SWPAR4;
			case 18: return SWPAR14;
			case 19: return SWPAR19;
		}
		case 17: switch((int)max(a,b)) {
			case 17: return SWPAR2;
			case 18: return SWPAR13;
			case 19: return SWPAR17;
		}
		case 18: switch((int)max(a,b)) {
			case 18: return SWPAR8;
			case 19: return SWPAR17;
		}
		case 19: return SWPAR22;
	}
}
// side chain bead diameter for chn[i] <-> chn[j]
const double DiaSC(AmiAc AmAci, AmiAc AmAcj)
{
    int a,b;
    a = AmAci.getAAnum();
    b = AmAcj.getAAnum();
    switch ((int)min(a,b)) {
		case 0: 
            return 0.0;
		case 1: switch ((int)max(a,b)) {
			case 1: return 2.7;
			case 2: return 2.7;
			case 3: return 2.9;
			case 4: return 2.6;
			case 5: return 2.3;
			case 6: return 2.8;
			case 7: return 2.6;
			case 8: return 3.0;
			case 9: return 3.3;
			case 10: return 2.9;
			case 11: return 3.0;
			case 12: return 2.7;
			case 13: return 2.9;
			case 14: return 2.4;
			case 15: return 2.7;
			case 16: return 2.7;
			case 17: return 2.9;
			case 18: return 2.8;
			case 19: return 3.1;
		}
		case 2: switch ((int)max(a,b)) {
			case 2: return 3.3;
			case 3: return 3.3;
			case 4: return 2.8;
			case 5: return 2.8;
			case 6: return 3.1;
			case 7: return 3.0;
			case 8: return 3.1;
			case 9: return 3.1;
			case 10: return 3.1;
			case 11: return 3.3;
			case 12: return 3.0;
			case 13: return 3.3;
			case 14: return 3.2;
			case 15: return 3.0;
			case 16: return 2.9;
			case 17: return 3.0;
			case 18: return 2.9;
			case 19: return 3.1;
		}
		case 3: switch ((int)max(a,b)) {
			case 3: return 3.1;
			case 4: return 2.6;
			case 5: return 3.2;
			case 6: return 3.3;
			case 7: return 3.2;
			case 8: return 3.0;
			case 9: return 3.6;
			case 10: return 3.5;
			case 11: return 3.6;
			case 12: return 3.5;
			case 13: return 3.5;
			case 14: return 3.1;
			case 15: return 3.3;
			case 16: return 3.4;
			case 17: return 3.7;
			case 18: return 3.0;
			case 19: return 3.7;
		}
		case 4: switch ((int)max(a,b)) {
			case 4: return 2.9;
			case 5: return 2.9;
			case 6: return 3.1;
			case 7: return 3.1;
			case 8: return 3.2;
			case 9: return 3.1;
			case 10: return 3.1;
			case 11: return 3.3;
			case 12: return 3.2;
			case 13: return 3.0;
			case 14: return 2.8;
			case 15: return 3.2;
			case 16: return 3.3;
			case 17: return 3.6;
			case 18: return 2.7;
			case 19: return 2.9;
		}
		case 5: switch ((int)max(a,b)) {
			case 5: return 2.5;
			case 6: return 3.0;
			case 7: return 2.8;
			case 8: return 3.0;
			case 9: return 3.0;
			case 10: return 2.9;
			case 11: return 2.7;
			case 12: return 3.0;
			case 13: return 2.6;
			case 14: return 2.9;
			case 15: return 2.9;
			case 16: return 2.7;
			case 17: return 3.2;
			case 18: return 2.8;
			case 19: return 2.6;
		}
		case 6: switch ((int)max(a,b)) {
			case 6: return 3.3;
			case 7: return 3.2;
			case 8: return 2.9;
			case 9: return 3.2;
			case 10: return 3.1;
			case 11: return 3.5;
			case 12: return 3.4;
			case 13: return 2.8;
			case 14: return 2.7;
			case 15: return 3.3;
			case 16: return 2.8;
			case 17: return 3.5;
			case 18: return 3.1;
			case 19: return 3.4;
		}
		case 7: switch ((int)max(a,b)) {
			case 7: return 3.4;
			case 8: return 3.0;
			case 9: return 3.0;
			case 10: return 2.9;
			case 11: return 2.8;
			case 12: return 3.0;
			case 13: return 3.4;
			case 14: return 3.1;
			case 15: return 2.8;
			case 16: return 3.2;
			case 17: return 3.6;
			case 18: return 3.2;
			case 19: return 2.8;
		}
		case 8: switch ((int)max(a,b)) {
			case 8: return 3.2;
			case 9: return 3.9;
			case 10: return 3.1;
			case 11: return 3.6;
			case 12: return 3.4;
			case 13: return 3.6;
			case 14: return 3.3;
			case 15: return 3.1;
			case 16: return 3.0;
			case 17: return 3.7;
			case 18: return 3.3;
			case 19: return 3.5;
		}
		case 9: switch ((int)max(a,b)) {
			case 9: return 3.5;
			case 10: return 3.4;
			case 11: return 3.4;
			case 12: return 3.5;
			case 13: return 2.9;
			case 14: return 3.5;
			case 15: return 3.5;
			case 16: return 3.5;
			case 17: return 3.7;
			case 18: return 2.7;
			case 19: return 3.4;
		}
		case 10: switch ((int)max(a,b)) {
			case 10: return 3.2;
			case 11: return 2.9;
			case 12: return 3.3;
			case 13: return 3.2;
			case 14: return 3.3;
			case 15: return 3.3;
			case 16: return 3.5;
			case 17: return 3.3;
			case 18: return 2.7;
			case 19: return 3.3;
		}
		case 11: switch ((int)max(a,b)) {
			case 11: return 3.6;
			case 12: return 3.5;
			case 13: return 3.1;
			case 14: return 3.3;
			case 15: return 3.4;
			case 16: return 3.4;
			case 17: return 3.4;
			case 18: return 3.1;
			case 19: return 3.3;
		}
		case 12: switch ((int)max(a,b)) {
			case 12: return 3.4;
			case 13: return 3.4;
			case 14: return 3.4;
			case 15: return 3.2;
			case 16: return 3.4;
			case 17: return 3.6;
			case 18: return 3.4;
			case 19: return 3.2;
		}
		case 13: switch ((int)max(a,b)) {
			case 13: return 3.3;
			case 14: return 3.4;
			case 15: return 3.0;
			case 16: return 3.2;
			case 17: return 3.6;
			case 18: return 3.3;
			case 19: return 3.1;
		}
		case 14: switch ((int)max(a,b)) {
			case 14: return 3.3;
			case 15: return 3.2;
			case 16: return 3.4;
			case 17: return 3.2;
			case 18: return 3.2;
			case 19: return 2.9;
		}
		case 15: switch ((int)max(a,b)) {
			case 15: return 3.0;
			case 16: return 3.2;
			case 17: return 3.2;
			case 18: return 2.9;
			case 19: return 3.1;
		}
		case 16: switch ((int)max(a,b)) {
			case 16: return 3.7;
			case 17: return 3.2;
			case 18: return 3.3;
			case 19: return 3.2;
		}
		case 17: switch ((int)max(a,b)) {
			case 17: return 3.7;
			case 18: return 3.4;
			case 19: return 3.6;
		}
		case 18: switch ((int)max(a,b)) {
			case 18: return 2.1;
			case 19: return 2.8;
		}
		case 19: return 3.4;
	}
}
/* squared bead pair HC diameters, including reduced diameters for neighbouring residues
 * comparing chn[ha].AmiAc[ia].Bd[ja] with chn[hb].AmiAc[ib].Bd[jb]
 */
double DiaSQ(Chain chn[], int ha, int ia, int ja, int hb, int ib, int jb)
{
	int AA2_num;        // numerical value of other amino acid
	double d1, d2;		// HC diameters of AA1 and AA2
	int h1, h2, i1, i2, j1, j2;

    bool squeeze = true;

	if(chn[ha].AmAc[ia].getAAnum() == 0) return 0;
	if(chn[hb].AmAc[ib].getAAnum() == 0) return 0;

	// neighbouring residues => with SQZ, (pseudo-)bonded beads are given diameter 0.
	if(ha==hb) {
        if(ia == ib) return 0;
	    if( abs(ia-ib) == 1 ) {
            h1 = (ia<ib) ? ha:hb;   h2 = (ia<ib) ? hb:ha;
	    	i1 = (ia<ib) ? ia:ib;	i2 = (ia<ib) ? ib:ia;
	    	j1 = (ia<ib) ? ja:jb;   j2 = (ia<ib) ? jb:ja;
            if(squeeze) {
	    	    switch (j1) {
                    case 0: switch(j2) {	// N
                        case 0: return SQZ4 * SQZ4;
                        case 1: return SQZ2 * SQZ2;
                        case 2: d1 = DIA_N + DIA_C; return 0.25*d1*d1;
                        case 3: d1 = DIA_N + DiaSC(chn[h2].AmAc[i2], chn[h2].AmAc[i2]); return 0.25*d1*d1;
                    }
                    case 1: switch(j2) {	// Ca
                        case 0: return 0;
                        case 1: return 0;
                        case 2: return SQZ1 * SQZ1;
                        case 3: return chn[h2].AmAc[i2].getSQZ(9) * chn[h2].AmAc[i2].getSQZ(9);
                    }
                    case 2: switch(j2) {	// C
                        case 0: return 0;
                        case 1: return 0;
                        case 2: return SQZ5 * SQZ5;
                        case 3: return chn[h2].AmAc[i2].getSQZ(6) * chn[h2].AmAc[i2].getSQZ(6);
                    }
                    case 3: switch(j2) {	// R
                        case 0: return chn[h1].AmAc[i1].getSQZ(7) * chn[h1].AmAc[i1].getSQZ(7);
                        case 1: return chn[h1].AmAc[i1].getSQZ(8) * chn[h1].AmAc[i1].getSQZ(8);
                        case 2: d1 = DiaSC(chn[h1].AmAc[i1], chn[h1].AmAc[i1]) + DIA_C; return 0.25*d1*d1;
                        case 3: return DiaSC(chn[h1].AmAc[i1], chn[h2].AmAc[i2]) * DiaSC(chn[h1].AmAc[i1], chn[h2].AmAc[i2]);
                    }
                    default: return 0;
	    	    }
            }
            else {
                switch (j1) {
                    case 0: switch(j2) {
	    		    	case 0: return 0.5625*DIA_N*DIA_N;
	    		    	case 1: d1=DIA_N+DIA_Ca; return 0.25*d1*d1;
	    		    	case 2: d1=DIA_N+DIA_C; return 0.25*d1*d1;
	    		    	case 3: d1=DIA_N+DiaSC(chn[h2].AmAc[i2], chn[h2].AmAc[i1]); return 0.25*d1*d1;
	    		    }
	    		    case 1: switch(j2) {
	    		    	case 0: return 0;
	    		    	case 1: return 0;
	    		    	case 2: d1=DIA_Ca+DIA_C; return 0.25*d1*d1;
	    		    	case 3: d1=DIA_Ca+DiaSC(chn[h2].AmAc[i2], chn[h2].AmAc[i2]); return 0.25*d1*d1;
	    		    }
	    		    case 2: switch(j2) {
	    		    	case 0: return 0;
	    		    	case 1: return 0;
	    		    	case 2: return 0.5625*DIA_C*DIA_C;
	    		    	case 3: d1=DIA_C+DiaSC(chn[h2].AmAc[i2], chn[h2].AmAc[i2]); return 0.140625*d1*d1;     //0.140625 = 0.5625*0.25 = (3/4 * 1/2)^2
	    		    }
	    		    case 3: switch(j2) {
	    		    	case 0: d1=DIA_N+DiaSC(chn[h1].AmAc[i1], chn[h1].AmAc[i1]); return 0.140625*d1*d1;
	    		    	case 1: d1=DiaSC(chn[h1].AmAc[i1], chn[h1].AmAc[i1])+DIA_Ca; return 0.25*d1*d1;
	    		    	case 2: d1=DiaSC(chn[h1].AmAc[i1], chn[h1].AmAc[i1])+DIA_C; return 0.25*d1*d1;
	    		    	case 3: return DiaSC(chn[h1].AmAc[i1], chn[h2].AmAc[i2])*DiaSC(chn[h1].AmAc[i1], chn[h2].AmAc[i2]);
	    		    }
	    		    default: return 0;
                }
            }
	    }
    }
	if(squeeze) {
        if(ha==hb) {
            if(abs(ia-ib) == 2) {
                h1 = (ia<ib) ? ha:hb;   h2 = (ia<ib) ? hb:ha;
    	    	i1 = (ia<ib) ? ia:ib;	i2 = (ia<ib) ? ib:ia;
    	    	j1 = (ia<ib) ? ja:jb;   j2 = (ia<ib) ? jb:ja;
                if(j1==2 && j2==0) return SQZ3 * SQZ3;
                if(j1==2 && j2==3) return chn[h2].AmAc[i2].getSQZ(10) * chn[h2].AmAc[i2].getSQZ(10);
    	    }
        }
    }

	// non-neighbouring residues => normal case
	if (ja == jb) {
		switch(ja) {
			case 0: return DIA_N*DIA_N;
			case 1: return DIA_Ca*DIA_Ca;
			case 2: return DIA_C*DIA_C;
			case 3:  return DiaSC(chn[ha].AmAc[ia], chn[hb].AmAc[ib])*DiaSC(chn[ha].AmAc[ia], chn[hb].AmAc[ib]);
		}
	}
	switch (ja) {
		case 0: d1 = DIA_N; break;
		case 1: d1 = DIA_Ca; break;
		case 2: d1 = DIA_C; break;
		case 3: d1 = DiaSC(chn[ha].AmAc[ia], chn[ha].AmAc[ia]);
	}
	switch (jb) {
		case 0: d2 = DIA_N; break;
		case 1: d2 = DIA_Ca; break;
		case 2: d2 = DIA_C; break;
		case 3: d2 = DiaSC(chn[hb].AmAc[ib], chn[hb].AmAc[ib]);
	}
	return 0.25*(d1+d2)*(d1+d2);
}


std::vector<std::vector<double>> DiaSQValues;
// setup squared diameters values of matrix
int DiaSQValuesSetup(Chain chn[], int naa, int nch)
{
	double value;


	std::vector<double> dist;
	for( int i=0; i<4*naa*nch; i++ ) {
		for( int j=0; j<4*naa*nch; j++ ) {
			value = DiaSQ(chn, i/(4*naa), (i/4)%naa, i%4, j/(4*naa), (j/4)%naa, j%4);
			dist.push_back(value);
			//std::cout << i << "-" << j << ":" << value;
			//std::cout << "  " << std::flush;
		}
		//std::cout << std::endl;
		DiaSQValues.push_back(dist);
		dist.clear();
	}

	return 0;
}