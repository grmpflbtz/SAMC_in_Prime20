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
// constructor
Chain::Chain(void)
{
}
// destructor
Chain::~Chain(void)
{
}