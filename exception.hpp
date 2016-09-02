#pragma once
#include <stdexcept>

namespace frea {
	struct InvalidAxis : std::invalid_argument {
		using std::invalid_argument::invalid_argument;
	};
	struct NoValidAxis : std::runtime_error {
		using std::runtime_error::runtime_error;
	};
}
