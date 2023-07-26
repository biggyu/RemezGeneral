#pragma once
#include"func.h"
#include"Point.h"
#include"Polynomial.h"
#include"RemezParam.h"
#include<string>
#include<vector>
#include<future>
#include<NTL/mat_RR.h>
#define _USE_MATH_DEFINES

using namespace std;
using namespace NTL;

class RemezGeneral {
private:
    RR* width;
    RR* sc;
    RR approx, max_err, min_err, current_err;
    Point *sample_point, *extreme_point;
    long extreme_count;
    RR* coeff;

public:
    RemezParam params;
//    long boundary_K;
    long section_num;
    double* sections;
    double chebeval_k;
    long deg;

    RemezGeneral() {}
    RemezGeneral(RemezParam _params, long _section_num, double* _sections, long _deg);

    virtual RR function_value(RR x) = 0;
    void better_initialize();
    void initialize();
    void getcoeffwitherr();
    void getextreme_local(Point* local_extreme_point, long &local_extreme_count, long section_ind);
    void getextreme();
    void choosemaxs();

    void generate_optimal_poly(Polynomial &poly);
    // void showcoeff();
};
