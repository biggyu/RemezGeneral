#pragma once

#include "RemezGeneral.h"

class RemezGeneralArcsin: public RemezGeneral {
public:
	RemezGeneralArcsin(RemezParam _param, long _section_num, double* _sections, long _deg)
		: RemezGeneral(_param, 1, _sections, _deg) {}

	RR function_value(RR x) {
		// return 1/(2*ComputePi_RR())*to_RR(asin(to_double(x)));
		return 1/(2*ComputePi_RR())*arcsin(x);
	}
};
