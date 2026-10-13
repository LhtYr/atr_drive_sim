// 涡轮稳态特性扫描：膨胀比 pi = p_plenum/p_downstream × 转速 网格。
//
// 名义边界条件（覆盖试验标定包线，见 参数.md：膨胀比 1.5–4.0、
// 入口温度 729–831 K、前腔压力 0.22–0.45 MPa）：
//   p_plenum = 4e5 Pa（前腔压力，取试验包线上沿）
//   p_in     = p_plenum（入口无节流，m_dot_in 不影响喷嘴流量特性）
//   T_in     = 800 K（试验温度区间中值），p_out = p_plenum / pi
//
// 输出列：
//   pi            膨胀比 p_plenum/p_downstream
//   omega_rpm     转速 rpm
//   m_dot_turb    通过涡轮质量流量 kg/s
//   mdot_reduced  换算流量 m_dot_turb*sqrt(T_in)/p_plenum (kg·K^0.5/(s·Pa))
//   T_turb        涡轮转矩 N·m
//   power_W       涡轮输出功率 W
//   efficiency    等熵效率 = power / (m_dot_turb*cp*T_in*(1-pi^(-(γ-1)/γ)))

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Turbine.h"
#include "char_sweep.h"

#define TURB_PLENUM_PA   4.0e5		// 前腔压力 Pa
#define TURB_T_IN_K      800.0		// 入口温度 K（试验区间 729–831 K 中值）
#define TURB_PI_MIN      1.0		// 膨胀比扫描下限（=1 无膨胀）
#define TURB_PI_MAX      4.0		// 膨胀比扫描上限（试验最大 4.0）
#define TURB_RPM_MIN     2000.0		// 转速扫描下限 rpm
#define TURB_RPM_MAX     60000.0	// 转速扫描上限 rpm
#define CHAR_PI          3.14159265358979323846	// 严格 -std=c11 下 glibc 不提供 M_PI

int char_sweep_turbine(FILE* f, int n) {
	double* pis = malloc((size_t)n * sizeof(double));
	double* rpms = malloc((size_t)n * sizeof(double));
	if (pis == NULL || rpms == NULL) {
		free(pis);
		free(rpms);
		fprintf(stderr, "错误：内存分配失败\n");
		return -1;
	}
	char_linspace(pis, TURB_PI_MIN, TURB_PI_MAX, n);
	char_linspace(rpms, TURB_RPM_MIN, TURB_RPM_MAX, n);

	fprintf(f, "pi,omega_rpm,m_dot_turb,mdot_reduced,T_turb,power_W,efficiency\n");

	const TurbineParams* p = &turbine_params;
	const double gamma_safe = fmax(p->gamma, 1.01);
	const double cp = gamma_safe * p->gas_R / (gamma_safe - 1.0);

	for (int i = 0; i < n; i++) {
		const double p_down = TURB_PLENUM_PA / pis[i];
		for (int j = 0; j < n; j++) {
			const double omega = rpms[j] * (2.0 * CHAR_PI / 60.0);
			const TurbineEval e = turbine_evaluate(p,
				TURB_PLENUM_PA, p_down, TURB_T_IN_K, 0.0,
				TURB_PLENUM_PA, p_down, omega);

			const double power = e.T_turb * omega;
			const double mdot_reduced = e.m_dot_turb * sqrt(TURB_T_IN_K) / TURB_PLENUM_PA;

			// 等熵可用焓降对应的理想功率；pi->1 时退化为 0，效率记 0
			const double dh_isen = cp * TURB_T_IN_K
				* (1.0 - pow(pis[i], -(gamma_safe - 1.0) / gamma_safe));
			const double power_ideal = e.m_dot_turb * dh_isen;
			const double efficiency = power_ideal > 1e-9 ? power / power_ideal : 0.0;

			fprintf(f, "%.6f,%.1f,%.9g,%.9g,%.9g,%.9g,%.6f\n",
				pis[i], rpms[j], e.m_dot_turb, mdot_reduced,
				e.T_turb, power, efficiency);
		}
	}

	free(pis);
	free(rpms);
	return 0;
}
