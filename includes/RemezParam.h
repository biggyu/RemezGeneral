#pragma once

#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

class RemezParam {
public:
	// double log_scan_step_diff = 9.5;
	double log_scan_step_diff = 11; //epsilon
	long binary_prec = 10; //극값 찾을 때 사용. sc/2^n까지 들어가서 확인->정확한 극값 찾을 수 있음.
	long RR_prec = 1500;
	long log_approx_degree = 120;
	long log_round_prec = 100;

	RemezParam() {}
	RemezParam(double scd, long bi, long rr, long app, long round)
		: log_scan_step_diff(scd), binary_prec(bi), RR_prec(rr), log_approx_degree(app), log_round_prec(round) {}

	RemezParam(const RemezParam& copy)
		: log_scan_step_diff(copy.log_scan_step_diff), binary_prec(copy.binary_prec), RR_prec(copy.RR_prec), log_approx_degree(copy.log_approx_degree), log_round_prec(copy.log_round_prec) {}
};
