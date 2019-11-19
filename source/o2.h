#ifndef O2_H
#define O2_H

#include <math.h>

#define C_TWO_PI (6.2831853071795864769252867665590058L)
#define C_PI (3.1415926535897932384626433832795029L)

typedef double _FPTYPE;

struct O2 {
  _FPTYPE x; /* X */
  _FPTYPE y; /* Y */
  _FPTYPE t; /* angle */
  _FPTYPE r; /* modulus */
};

#define O2sum(__sum,__s1,__s2) {			\
	(__sum).x=(__s1).x+(__s2).x;			\
	(__sum).y=(__s1).y+(__s2).y;			\
    }

#define O2comb(__sum,__c1,__s1,__c2,__s2) {	\
	(__sum).x=__c1*(__s1).x+__c2*(__s2).x;	\
	(__sum).y=__c1*(__s1).y+__c2*(__s2).y;	\
    }

#define O2prop(__sum,__c1,__s1) {		\
	(__sum).x=__c1*(__s1).x;		\
	(__sum).y=__c1*(__s1).y;		\
    }

#define O2accum(__acc,__s1) { \
	(__acc).x+=(__s1).x;  \
	(__acc).y+=(__s1).y;  \
    }

#define O2prod(__a,__b) ((__a).x*(__b).x+(__a).y*(__b).y)
#define O2vprod(__a,__b) ((__a).x*(__b).y-(__a).y*(__b).x)
#define O2norm2(__a) O2prod((__a),(__a))
#define O2norm(__a) (sqrt(O2norm2((__a))))

#define O2perp(__s, __sp) {  \
    (__sp).x= (__s).y;                        \
    (__sp).y= -(__s).x;                        \
}

#define O2phasediff(__phs,__a,__b) {				\
	(__phs).x=((__a).x*(__b).x+(__a).y*(__b).y);		\
	(__phs).y=((__a).y*(__b).x-(__a).x*(__b).y);		\
	(__phs).t=((__a).t-(__a).t);				\
    }

#define O2phasesum(__phs,__a,__b) {			\
	(__phs).x=((__a).x*(__b).x-(__a).y*(__b).y);	\
	(__phs).y=((__a).y*(__b).x+(__a).x*(__b).y);	\
	(__phs).t=((__a).t+(__a).t);			\
    }

#define polar_to_cartesian(__s) {			\
	(__s).x=((__s).r)*cos((__s).t);			\
	(__s).y=((__s).r)*sin((__s).t);			\
    }

#define cartesian_to_polar(__s) {			        \
    (__s).t=atan2((__s).y,(__s).x);                     \
    (__s).t=((__s).t>=0?(__s).t:(__s).t+C_TWO_PI);	\
    (__s).r=sqrt( ((__s).y*(__s).y) + ((__s).x*(__s).x));	\
}

#endif
