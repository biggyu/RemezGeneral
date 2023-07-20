#include <iostream>
#include <NTL/RR.h>
#include "RemezParam.h"
#include "RemezGeneralCos.h"
#include "Polynomial.h"
#include "func.h"
using namespace std;
using namespace NTL;
int main()
{
    long section_num = 6;
    long deg = 16;
    double sections[] = {-.75, -.625, -.5, -.25, -.125, 0, .125, .25, .5, .625, .75, .875};
    // long section_num = 2;
    // double sections[] = {-.75, -.625, -.5, -.25};
    RemezParam rmparm;
    Polynomial sin_cos_polynomial;
    RemezGeneralCos* poly_generator = new RemezGeneralCos(rmparm, section_num, sections, deg, (1 << 2));
    poly_generator->generate_optimal_poly(sin_cos_polynomial);
    // // sin_cos_polynomial.generate_poly_heap();
    sin_cos_polynomial.showcoeff();
    cout << "Hello World!" << endl;
    cout << endl;
}