#include <iostream>
#include <NTL/RR.h>
#include "RemezParam.h"
#include "RemezGeneralCos.h"
#include "RemezGeneralTanh.h"
#include "RemezGeneralGelu.h"
#include "Polynomial.h"
#include "func.h"
using namespace std;
using namespace NTL;
int main()
{
    long deg = 9;
    long section_num = 6;
    double sections[] = {-.75, -.6875, -.5, -.25, -.125, 0, .125, .25, .5, .625, .75, .875};
    // double sections[] = {-4.5, -4.4375, -2.5, -2.25, -.125, 0, .125, .25, 3.5, 3.625, 5.75, 5.875};
    // long section_num = 2;
    // double sections[] = {-.75, -.625, -.5, -.25};

    RemezParam rmparm;

    // Polynomial sin_cos_polynomial;
    // RemezGeneralCos* poly_generator = new RemezGeneralCos(rmparm, section_num, sections, deg, (1 << 2));
    // poly_generator->generate_optimal_poly(sin_cos_polynomial);
    // sin_cos_polynomial.showcoeff();

    // Polynomial tanh_polynomial;
    // RemezGeneralTanh* poly_generator = new RemezGeneralTanh(rmparm, section_num, sections, deg);
    // poly_generator->generate_optimal_poly(tanh_polynomial);
    // tanh_polynomial.showcoeff();
    
    Polynomial gelu_polynomial;
    RemezGeneralGelu* poly_generator = new RemezGeneralGelu(rmparm, section_num, sections, deg);
    poly_generator->generate_optimal_poly(gelu_polynomial);
    gelu_polynomial.showcoeff();

    // // sin_cos_polynomial.generate_poly_heap();
}