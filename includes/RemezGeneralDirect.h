#pragma once

#include "RemezGeneral.h"

class RemezGeneralDirect : public RemezGeneral {
public:
	RemezGeneralDirect(RemezParam _param, long _section_num, double* _sections, long _deg)
		: RemezGeneral(_param, _section_num, _sections, _deg) {}

	RR function_value(RR x) {
		return fracpart(x);
	}
};
