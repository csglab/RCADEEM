// Copyright 2011 Hamed Shateri Najafabadi

/********************************************************************

This file is part of SEAMWU.

SEAMWU is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

SEAMWU is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with SEAMWU.  If not, see <http://www.gnu.org/licenses/>.

********************************************************************/


#include <math.h>


//******************************** These functions are provided by Hani Goodarzi (Princeton University) ******************

double gammln(double xx) {
 double x,y,tmp,ser;
 static double cof[6]={76.18009172947146,-86.50532032941677,
                       24.01409824083091,-1.231739572450155,
                       0.1208650973866179e-2,-0.5395239384953e-5};
 int j;
 y=x=xx;
 tmp=x+5.5;
 tmp -= (x+0.5)*log(tmp);
 ser=1.000000000190015;
 for (j=0;j<=5;j++) ser += cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}

double factln(int n) {
 static double a[101];
 if (n < 0) return 0.0; // this is in fact an error
 if (n <= 1) return 0.0;
 if (n <= 100) return a[n] ? a[n] : (a[n]=gammln(n+1.0));
 else return gammln(n+1.0);
}

double hypergeom(int i, int s1, int s2, int N) {
 return exp(factln(s1) + factln(N - s1) + factln(s2) + factln(N - s2) -
            factln(i) - factln(N) - factln(s1 - i) - factln(s2 - i)
            - factln(N - s1 - s2 + i));

}

double cumhyper(int i, int s1, int s2, int N) {

 int min = (s1<s2?s1:s2);
 double prod = 0.0;
 int a;
 double tmp = 0.0;

 for (a=i; a<=min; a++) {
   tmp = (double)hypergeom(a, s1, s2, N);
   prod += (double)hypergeom(a, s1, s2, N);
 }

 return prod;
}
