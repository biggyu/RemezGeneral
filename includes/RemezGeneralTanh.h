#pragma once

#include "RemezGeneral.h"

class RemezGeneralTanh: public RemezGeneral {
public:
	RemezGeneralTanh(RemezParam _param, long _section_num, double* _sections, long _deg)
		: RemezGeneral(_param, _section_num, _sections, _deg) {}

	RR function_value(RR x) {
		return (exp(2 * x) - 1) / (exp(2 * x) + 1);
	}
};
