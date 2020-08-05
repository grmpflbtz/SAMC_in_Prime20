/* FILE: prime20_samc_amino.h
 *
 *       :''''''''''''''''''''''''''':
 *       :  Amino acid class HEADER  :
 *       :...........................:
 * 
 */
#ifndef PRIME20_SAMC_AMINO_HPP
#define PRIME20_SAMC_AMINO_HPP

#include <iostream>
#include <algorithm>
#include "prime20_samc_bead.hpp"

using namespace std;
//--------------------------------------------------
class AmiAc
{
    private:
        char AA_alp;        // type of amino acid (alphabetically)
        int AA_num;         // type of amino acid (numerically)
        int AA_pos;         // position in the polymer chain

        double BND_CaR;     // bond length between Ca and R
        double PBND_NR;     // pseudobond length between N and R
        double PBND_CR;     // pseudobond length between C and R
        double ANGL_RCaN;   // angle between R Ca N
        double ANGL_RCaC;   // angle between R Ca C

        double SQZ6;        // Squeeze parameter 6
        double SQZ7;        // Squeeze parameter 7
        double SQZ8;        // Squeeze parameter 8
        double SQZ9;        // Squeeze parameter 9
        double SQZ10;       // Squeeze parameter 10

    public:
        Bead Bd[4];       // Amino acid consists of 4 united atom beads

        void Setup(const char *chain, int i);                   // Sets up amino acid type & beads
        int getAAnum();                                         // access AA_num
        double getDisR(int i);                                  // access distance from side chain to bead i
        double getAng(int i);                                   // access angles
        double getSQZ(int i);                                   // access squeeze parameters
        double DiaSQ(const char *chain, int j, int ia, int ja); // squared bead pair HC diameters, including reduced diameters for neighbouring residues.

        AmiAc(void);
        ~AmiAc(void);
};

#endif