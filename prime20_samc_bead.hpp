/* FILE: prime20_samc_bead.h
 * 
 *        :''''''''''''''''''''':
 *        :  Bead class HEADER  :
 *        :.....................:
 * 
 */
#ifndef PRIME20_SAMC_BEAD_HPP
#define PRIME20_SAMC_BEAD_HPP

#include <iostream>

using namespace std;
//------------------------------------------------------
class Bead
{  
    friend class AmiAc;
    friend class Chain;

    private:
        double x;       // x coordinate
        double y;       // y coordinate
        double z;       // z coordinate
        double m;       // mass
        char btype;     // bead type
        int box;        // neighbour list box the particle is in
        int perBC[3];   // tracker of periodic boundary conditions applied to coordinates

    public:
        void setR(double rx, double ry, double rz);     // set position
        double getR(int i);                             // get position (0=x, 1=y, 2=z)
        void setM(double mass);                         // set mass
        double getM();                                  // get mass
        void set_btype(int i);                          // set bead type
        char get_btype();                               // get bead type
        void setBox(int i);
        int getBox();
        void setBC(int iCoord, int j);                  // set perBC value
        void addBC(int i, int j);                       // adds j to perBC[i]
        int getBC(int i);                               // get perBC value

        Bead(void);
        ~Bead(void);
};

#endif