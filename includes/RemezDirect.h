#pragma once

#include "Remez.h"

class RemezDirect : public Remez {
public:
	RemezDirect(RemezParam _param, long _boundary_K, long _log_width, long _deg)
		: Remez(_param, _boundary_K, _log_width, _deg) {}

	RR function_value(RR x) {
		return fracpart(x);
	}
};
