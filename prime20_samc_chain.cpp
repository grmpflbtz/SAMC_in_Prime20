/* FILE: prime20_samc_chain.cpp
 *
 *          :''''''''''''''''''''''':
 *          :    Chain class CPP    :
 *          :.......................:
 * 
 */

#include "prime20_samc_chain.hpp"
//--------------------------------------------------
// set chain number
void Chain::setChnNo( int no )
{
    chnNo = no;
}
// get chain number
int Chain::getChnNo()
{
    return chnNo;
}
// add m to mass
void Chain::addM( double m )
{
    mass += m;
}
// get mass
double Chain::getM()
{
    return mass;
}
// constructor
Chain::Chain(void)
{
    mass = 0.0;
}
// destructor
Chain::~Chain(void)
{
}