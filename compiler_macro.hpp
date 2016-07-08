#pragma once

//! 条件式の評価結果がfalseの場合、マクロ定義した箇所には到達しないことを保証
#define __assume(cond) do { if (!(cond)) __builtin_unreachable(); } while (0);
