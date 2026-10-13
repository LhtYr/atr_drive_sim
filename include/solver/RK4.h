#ifndef RK4_H
#define RK4_H

#include "Solver.h"

SystemState add_scaled(const SystemState state, const SystemDeriv deriv, double scale);

// RK4 合并公式：
// x_{k+1} = x_k + Ts/6 * (k1 + 2*k2 + 2*k3 + k4)
// 这里的 k1~k4 是各子步导数。
SystemState combine_rk4(const SystemState state,
                        const SystemDeriv k1,
                        const SystemDeriv k2,
                        const SystemDeriv k3,
                        const SystemDeriv k4,
                        double Ts);


#endif // RK4_H