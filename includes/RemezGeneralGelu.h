#pragma once

#include "RemezGeneral.h"

class RemezGeneralGelu: public RemezGeneral {
public:
	RemezGeneralGelu(RemezParam _param, long _section_num, double* _sections, long _deg)
		: RemezGeneral(_param, _section_num, _sections, _deg) {}

	RR function_value(RR x) {
		RemezGeneralTanh tmp;
		// cout << tmp.function_value(x) << endl;
		// return tmp.function_value(x);
		return .5 * x * (1 + tmp.function_value(sqrt(2 / ComputePi_RR()) * (x + 0.044715 * x * x * x)));
	}
};
