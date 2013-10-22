#ifndef FM_H
#define FM_H

/////////////////////////
// OS-independent includes
/////////////////////////

#include <mex.h>
#include <engine.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>
#include <algorithm>
//#include "fheap/fib.h"
//#include "fheap/fibpriv.h"

/////////////////////////
// numerical macros
/////////////////////////
#undef		MIN
#undef		MAX
#define		MIN(a, b)       ((a) < (b) ? (a) : (b))			//!<	Returns the min value between a and b
#define		MINMIN(a,b,c)   ((a) < (b) ? MIN(a,c) : MIN(b,c))
#define		MAX(a, b)       ((a) > (b) ? (a) : (b))			//!<	Returns the max value between a and b
#define		MAXMAX(a,b,c)   ((a) > (b) ? MAX (a,c) : MAX (b,c))
#define         SCALE_01(x,rMin,rMax) ((x-rMin)/(rMax-rMin))
#define		ABS(a)       ((a) > 0 ? (a) : -(a))			//!<	Returns the absolute value a
#define		SIGN(a)       ((a) > 0 ? 1 : -1)				//!<	Returns the sign of a
#define		SQR(x)			((x)*(x))						//!<	Returns x square
#define		CUBE(x)			((x)*(x)*(x))					//!<	Returns x cube

#define         CLAMP_01(x)	if( (x)<0 ) x=0; if( (x)>1 ) x=1
#define         CLAMP(x, a,b)	if( (x)<a ) x=a; if( (x)>b ) x=b

#define         SWAP(x,y) x^=y; y^=x; x^=y
#define         ORDER(x,y) if(x>y){ GW_SWAP(x,y); }

/////////////////////////
// generic macros
/////////////////////////
/** a random number in [0-1] */
#define RAND ((double) (rand()%10000))/10000
/** a random number in [a,b] */
#define RAND_RANGE(a,b) (a)+((b)-(a))*((GW_Float) (rand()%10000))/10000
/** delete a single pointer */
#define DELETE(p) {if (p!=NULL) delete p; p=NULL;}
/** delete an array pointer */
#define DELETEARRAY(p) {if (p!=NULL) delete [] p; p=NULL;}

/////////////////////////
// some constants
/////////////////////////
#define TRUE  true
#define FALSE false
/** to make aproximate computations (derivation, float comparaisons ...) */
#define EPSILON 1e-6
/** very big number */
#define INFINITE 1e9

/////////////////////////
// numerical  constants
/////////////////////////

/** pi */
#define PI		3.1415926535897932384626433832795028841971693993751f
/** pi/2 */
#define HALFPI	1.57079632679489661923f
/** 2*pi */
#define TWOPI	6.28318530717958647692f
/** 1/pi */
#define GW_INVPI	0.31830988618379067154f
/** 180/pi */
#define RADTODEG(x)	(x)*57.2957795130823208768f
/** pi/180 */
#define DEGTORAD(x)	(x)*0.01745329251994329577f
/** e */
#define EXP		2.71828182845904523536f
/** 1/log10(2) */
#define ILOG2	3.32192809488736234787f
/** 1/3 */
#define INV3		0.33333333333333333333f
/** 1/6 */
#define INV6		0.16666666666666666666f
/** 1/9 */
#define INV7		0.14285714285714285714f
/** 1/9 */
#define INV9		0.11111111111111111111f
/** 1/255 */
#define INV255	0.00392156862745098039f
/** sqrt(2) */
#define SQRT2    1.41421356237f

//================================================================
typedef std::vector<int>        LISTofINT;
typedef std::vector<short>      LISTofSHORT;
typedef std::vector<float>      LISTofFLOAT;
typedef std::vector<bool>       LISTofBOOL;
typedef std::vector<LISTofINT>  LISTofLISTofINT;
typedef std::vector<LISTofBOOL> LISTofLISTofBOOL;

class MIN_PATH
{
  public:
  int*        V;
  float		  U;
  LISTofINT   points;
};
typedef std::vector<MIN_PATH> LISTofMIN_PATH;

class KEYPOINT
{
  public:
  int    point;
  int    V;
  float  U;
};
typedef std::vector<KEYPOINT> LISTofKEYPOINT;
//================================================================

inline int log2(int n)
{
	int x = 0;
	while (n > 1) 
	{
		x++;
		n /= 2;
	}
	return x;
}

inline int cmp(float a, float b)
{	
	if (a < b)
		return -1;
	if (a == b)
		return 0;
	return 1;
}

int sign(int v)
{
return v > 0 ? 1 : (v < 0 ? -1 : 0);
}



#endif // #ifndef FM_H
