
//
// Created by ilaria on 2019-11-18.
//

#ifndef O3_H
#define O3_H
#include <math.h>
typedef double _FPTYPE;

struct O3 {
    _FPTYPE x; /* X */
    _FPTYPE y; /* Y */
    _FPTYPE z; /* Z */
    _FPTYPE p; /*angle x-y plane arctan(y/x)*/
    _FPTYPE t; /* angle z-xy arccos(z/r) */
    _FPTYPE r; /* modulus */
};

#define O3prod(__a,__b) ((__a).x*(__b).x+(__a).y*(__b).y +(__a).z*(__b).z)
#define O3norm2(__a) O3prod((__a),(__a))
#define O3norm(__a) (sqrt(O3norm2((__a))))

#define O3polar_to_cartesian(__s) {			\
	(__s).x=((__s).r)*sin((__s).t)*cos((__s).p);			\
	(__s).y=((__s).r)*sin((__s).t)*sin((__s).p);			\
	(__s).t=((__s).r)*sin((__s).t);			\
    }

#define O3cartesian_to_polar(__s) {			        \
    (__s).r=sqrt( ((__s).y*(__s).y) + ((__s).x*(__s).x) + ((__s).z*(__s).z));	\
    (__s).p=atan2((__s).y,(__s).x);                     \
    (__s).p=((__s).p>=0?(__s).p:(__s).p+C_TWO_PI);	\
    (__s).t=acos((__s).z,(__s).r);                     \
}

#endif //DEFAULT_NAME_O3_H
