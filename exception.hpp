#pragma once
#include <stdexcept>

namespace frea {
	struct InvalidAxis : std::invalid_argument {
		InvalidAxis(): std::invalid_argument("InvalidAxis") {}
	};
	struct NoValidAxis : std::runtime_error {
		NoValidAxis(): std::runtime_error("NoValidAxis") {}
	};
	struct InvalidFov : std::invalid_argument {
		InvalidFov(): std::invalid_argument("InvalidFov") {}
	};
	struct NoInverseMatrix : std::runtime_error {
		NoInverseMatrix(): std::runtime_error("NoInverseMatrix") {}
	};
}
