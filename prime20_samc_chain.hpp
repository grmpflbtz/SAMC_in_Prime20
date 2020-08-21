/* FILE: prime20_samc_chain.h
 * 
 *        :'''''''''''''''''''''':
 *        :  Chain class HEADER  :
 *        :......................:
 * 
 */

#ifndef PRIME20_SAMC_CHAIN_HPP
#define PRIME20_SAMC_CHAIN_HPP

#include <vector>
#include "prime20_samc_amino.hpp"
//------------------------------------------------------
class Chain
{
    private:
        int chnNo;                              // chain number
        std::vector<char> sequence;             // Amino Acid sequence
        double mass;                            // total chain mass
    public:
        std::vector<AmiAc> AmAc;

        void setChnNo( int no );                // set chain number
        int getChnNo();                         // return chain number
        void addM( double m );                  // add m to mass of chain
        double getM();                          // get mass of chain

        Chain(void);
        ~Chain(void);
};

#endif