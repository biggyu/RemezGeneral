#pragma once

#include "RemezGeneral.h"

class RemezGeneralCos : public RemezGeneral {
public:
	long scale_factor;

	RemezGeneralCos(RemezParam _param, long _section_num, double* _sections, long _deg, long _scale_factor)
		: RemezGeneral(_param, _section_num, _sections, _deg), scale_factor(_scale_factor) {}

	RR function_value(RR x) {
		if(scale_factor % 2 == 0) return cos(2 * ComputePi_RR() * (x - 0.25) / scale_factor);
		else return sin(2 * ComputePi_RR() * x / scale_factor);
	}
};
