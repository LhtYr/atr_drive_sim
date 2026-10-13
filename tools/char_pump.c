// 泵稳态特性扫描：转速 × 比流量 x = Q/n 满网格（覆盖试验包线）。
//
// 名义边界条件：
//   p_suction    = 1e5 Pa（吸入压力）
//   p_downstream = p_discharge（下游与出口同压 → Q_out = 0，纯泵本体特性）
//
// 泵模型为试验拟合特性（H-Q 二次式 + eta-Q/n 三次式，默认参数即标定值）。
// 每点由拟合扬程 H(Q,n) 定出 Δp = rho*g*H 再调 pump_evaluate 反解验证，
// x 取下降支顶点（≈8.07e-5）至效率拟合零点（≈3.05e-4）之间，保证满网格
// 全部落在包线内，供三维曲面绘图使用。
//
// 输出列：
//   omega_rpm   转速 rpm
//   dp_pump_Pa  泵压升 Pa
//   H_m         扬程 m
//   Q_pump      泵流量 m^3/s
//   T_pump      泵负载转矩 N·m
//   P_shaft_W   轴功率 = T_pump * omega
//   P_hyd_W     水力功率 = dp * Q_pump
//   efficiency  泵效率拟合值 eta(Q/n)

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Pump.h"
#include "char_sweep.h"

#define PUMP_P_SUCTION_PA  1.0e5	// 吸入压力 Pa
#define PUMP_RPM_MIN       1000.0	// 转速扫描下限 rpm
#define PUMP_RPM_MAX       20000.0	// 转速扫描上限 rpm
#define PUMP_X_MIN         8.1e-5	// Q/n 扫描下限 (L/s)/rpm（下降支顶点）
#define PUMP_X_MAX         3.0e-4	// Q/n 扫描上限 (L/s)/rpm（效率拟合零点内）
#define PUMP_G             9.80665	// 重力加速度 m/s^2
#define CHAR_PI            3.14159265358979323846

int char_sweep_pump(FILE* f, int n) {
	double* rpms = malloc((size_t)n * sizeof(double));
	double* xs = malloc((size_t)n * sizeof(double));
	if (rpms == NULL || xs == NULL) {
		free(rpms);
		free(xs);
		fprintf(stderr, "错误：内存分配失败\n");
		return -1;
	}
	char_linspace(rpms, PUMP_RPM_MIN, PUMP_RPM_MAX, n);
	char_linspace(xs, PUMP_X_MIN, PUMP_X_MAX, n);

	fprintf(f, "omega_rpm,dp_pump_Pa,H_m,Q_pump,T_pump,P_shaft_W,P_hyd_W,efficiency\n");

	const PumpParams* p = &pump_params;
	for (int i = 0; i < n; i++) {
		const double omega = rpms[i] * (2.0 * CHAR_PI / 60.0);
		for (int j = 0; j < n; j++) {
			// 由比流量定工况：Q_ls = x*n，Δp = rho*g*H(Q,n)
			const double q_ls = xs[j] * rpms[i];
			const double h_fit = p->h_qq * q_ls * q_ls
				+ p->h_qn * q_ls * rpms[i] + p->h_nn * rpms[i] * rpms[i];
			const double p_discharge = PUMP_P_SUCTION_PA
				+ p->rho * PUMP_G * h_fit;
			// p_downstream = p_discharge：出流为零，考察泵自身特性
			const PumpEval e = pump_evaluate(p,
				PUMP_P_SUCTION_PA, p_discharge, p_discharge, omega);

			// 满网格应全部落在包线内；若参数被改出包线则跳过，
			// 避免族曲线在 Q=0 处出现折回连线
			if (e.Q_pump <= 0.0) {
				continue;
			}

			const double p_shaft = e.T_pump * omega;
			const double p_hyd = e.dp_pump * e.Q_pump;

			fprintf(f, "%.1f,%.1f,%.6f,%.9g,%.6f,%.6f,%.6f,%.6f\n",
				rpms[i], e.dp_pump, e.H, e.Q_pump, e.T_pump,
				p_shaft, p_hyd, e.eta);
		}
	}

	free(rpms);
	free(xs);
	return 0;
}
