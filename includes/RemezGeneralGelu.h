#pragma once

#include "RemezGeneral.h"

class RemezGeneralArcsin: public RemezGeneral {
public:
	RemezGeneralArcsin(RemezParam _param, long _section_num, double* _sections, long _deg)
		: RemezGeneral(_param, _section_num, _sections, _deg) {}

	RR function_value(RR x) {
		// return 1/(2*ComputePi_RR())*to_RR(asin(to_double(x)));
		// return 0.5 * x * (1 + tanh(sqrt(2 / ComputePi_RR()) * (x + 0.044715 * x * x * x)));
		// double tmp = to_double(x + 0.044715 * x * x * x);
		return 0.5 * x * (1 + tanh(sqrt(2 / M_PI) * to_double(x + 0.044715 * x * x * x)));
	}
};
