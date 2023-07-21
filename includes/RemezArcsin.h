#pragma once

#include "Remez.h"

class RemezArcsin: public Remez {
public:
	RemezArcsin(RemezParam _param, double _log_width, long _deg)
		: Remez(_param, 1, _log_width, _deg) {}

	RR function_value(RR x) {
		// return 1/(2*ComputePi_RR())*to_RR(asin(to_double(x)));
		return 1/(2*ComputePi_RR())*arcsin(x);
	}
};
