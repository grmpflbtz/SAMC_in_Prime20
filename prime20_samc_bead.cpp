/* FILE: prime20_samc_bead.h
 * 
 *        :'''''''''''''''''':
 *        :  Bead class CPP  :
 *        :..................:
 * 
 */

#include "prime20_samc_bead.hpp"

using namespace std;
//--------------------------------------------------
// set position
void Bead::setR(double rx, double ry, double rz)
{
    x = rx;
    y = ry;
    z = rz;
}
// get position
double Bead::getR(int i)
{
    if (i==0) {
        return x;
    }
    if (i==1) {
        return y;
    }
    if (i==2) {
        return z;
    }
    else
    {
        cout << "error reading bead postiion" << endl;
        return 0;
    }
    
}
// set mass
void Bead::setM(double mass)
{
    m = mass;
}
// get mass
double Bead::getM()
{
    return m;
}
// set bead type
void Bead::set_btype(int i)
{
    switch(i) {
        case 0: btype = 'N'; break;
        case 1: btype = 'C'; break;
        case 2: btype = 'O'; break;
        case 3: btype = 'R'; break;
    }
}
// get bead type
char Bead::get_btype()
{
    return btype;
}
// set Box
void Bead::setBox(int i)
{
    box = i;
}
// get Box
int Bead::getBox()
{
    return box;
}
// set perBC
void Bead::setBC(int iCoord, int j)
{
    perBC[iCoord] = j;
}
// adds j to perBC[i]
void Bead::addBC(int i, int j)
{
    perBC[i] += j;
}
// get perBC
int Bead::getBC(int i)
{
    return perBC[i];
}
// constructor
Bead::Bead(void)
{
    perBC[0] = 0;
    perBC[1] = 0;
    perBC[2] = 0;
}
// destructor
Bead::~Bead(void)
{
}