// 涡轮模型 7 点试验数据校验：用真实涡轮试验工况逐点比对模型流量，
// 并给出指定转速下的功率与 U/C 效率。流量误差超过 CHECK_TOL_PCT 判失败。
//
// 数据来源：涡轮试验报告（仓库根目录 参数.md），7 个工况点按列排列。
// 校验对象为全局 turbine_params（Turbine.c 默认值即试验标定值：
// gamma=1.18、gas_R=380 取固定圆整值，逐点物性差异计入误差；
// 若经 --config 覆盖，则校验覆盖后的参数）。
//
// 输出列：
//   point           工况点编号（1..7）
//   expansion_ratio 定压比（膨胀比 = p_in/p_out）
//   p_in_MPa        入口压力 MPa
//   T_in_K          入口温度 K
//   mdot_model      模型流量 kg/s
//   mdot_meas       实测流量 kg/s
//   err_pct         流量误差 %（(模型-实测)/实测）
//   u_over_c0       速比 u/c0（u = ω·r_mean，c0 = sqrt(2·cp·ΔT_ideal)）
//   eta_isen        U/C 多项式等熵效率
//   power_kW        指定转速下涡轮输出功率 kW（含机械效率）

#include <math.h>
#include <stdio.h>

#include "Turbine.h"
#include "char_sweep.h"

#define CHECK_RPM      20000.0	// 校验转速 rpm（与主仿真目标转速一致）
#define CHECK_TOL_PCT  1.0		// 流量误差判定阈值 %
#define CHECK_PI       3.14159265358979323846	// 严格 -std=c11 下 glibc 不提供 M_PI

// 试验工况点：温度 K / 定压比（膨胀比） / 入口压力 Pa / 实测流量 kg/s
typedef struct TurbCheckPoint {
	double T_in;
	double expansion_ratio;
	double p_in;
	double mdot_meas;
} TurbCheckPoint;

static const TurbCheckPoint CHECK_POINTS[] = {
	{ 729.0, 4.0, 0.424e6, 1.287 },
	{ 728.6, 2.0, 0.292e6, 0.878 },
	{ 749.7, 3.0, 0.402e6, 1.201 },
	{ 804.4, 2.3, 0.448e6, 1.291 },
	{ 808.7, 2.0, 0.374e6, 1.071 },
	{ 822.2, 2.0, 0.300e6, 0.848 },
	{ 830.9, 1.5, 0.221e6, 0.621 },
};
#define N_CHECK_POINTS ((int)(sizeof(CHECK_POINTS) / sizeof(CHECK_POINTS[0])))

int char_check_turbine(FILE* f, int n) {
	(void)n;	// 校验点固定为试验数据，与网格数无关

	const TurbineParams* p = &turbine_params;
	const double omega = CHECK_RPM * (2.0 * CHECK_PI / 60.0);
	const double cp = p->gamma * p->gas_R / (p->gamma - 1.0);

	fprintf(f, "point,expansion_ratio,p_in_MPa,T_in_K,mdot_model,mdot_meas,"
	           "err_pct,u_over_c0,eta_isen,power_kW\n");
	printf("涡轮 7 点试验数据校验（校验全局涡轮参数，默认即标定值，转速 %.0f rpm）\n",
	       CHECK_RPM);
	printf("%-4s %-6s %-8s %-7s %-10s %-10s %-8s %-7s %-8s %-9s\n",
	       "点", "膨胀比", "p_in/MPa", "T_in/K",
	       "模型流量", "实测流量", "误差%", "u/c0", "eta", "功率/kW");

	int failures = 0;
	for (int i = 0; i < N_CHECK_POINTS; i++) {
		const TurbCheckPoint* pt = &CHECK_POINTS[i];
		const double p_down = pt->p_in / pt->expansion_ratio;
		const TurbineEval e = turbine_evaluate(p,
			pt->p_in, p_down, pt->T_in, 0.0,
			pt->p_in, p_down, omega);

		const double err_pct = (e.m_dot_turb / pt->mdot_meas - 1.0) * 100.0;
		const double power = e.T_turb * omega;

		// 与 turbine_evaluate 内部一致的 U/C 速比与效率（用于报告）
		const double dT_ideal = pt->T_in
			* (1.0 - pow(p_down / pt->p_in, (p->gamma - 1.0) / p->gamma));
		const double c0 = sqrt(2.0 * cp * dT_ideal);
		const double uc = omega * p->mean_radius / c0;
		const double eta_isen = fmin(fmax(
			(p->eta_uc_a * uc + p->eta_uc_b) * uc + p->eta_uc_c, 0.0), 1.0);

		if (fabs(err_pct) > CHECK_TOL_PCT) {
			failures++;
		}
		fprintf(f, "%d,%.3f,%.3f,%.1f,%.6f,%.6f,%.4f,%.4f,%.4f,%.3f\n",
			i + 1, pt->expansion_ratio, pt->p_in / 1e6, pt->T_in,
			e.m_dot_turb, pt->mdot_meas, err_pct, uc, eta_isen, power / 1e3);
		printf("%-4d %-6.2f %-8.3f %-7.1f %-10.4f %-10.4f %+-8.3f %-7.3f %-8.4f %-9.1f%s\n",
			i + 1, pt->expansion_ratio, pt->p_in / 1e6, pt->T_in,
			e.m_dot_turb, pt->mdot_meas, err_pct, uc, eta_isen, power / 1e3,
			fabs(err_pct) > CHECK_TOL_PCT ? "  <-- FAIL" : "");
	}

	if (failures > 0) {
		printf("校验失败：%d/%d 个点流量误差超过 %.1f%%\n",
		       failures, N_CHECK_POINTS, CHECK_TOL_PCT);
		return -1;
	}
	printf("校验通过：%d 个点流量误差均 <= %.1f%%\n", N_CHECK_POINTS, CHECK_TOL_PCT);
	return 0;
}
