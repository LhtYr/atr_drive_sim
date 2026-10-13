// 涡轮模型动态稳态校验（附加诊断模式）：不绕开容腔动态，对每个试验工况点
// 给定引气边界（总压/总温/实测流量）积分两容腔动态方程至稳态，再检验
// 稳态前腔压力是否复现试验入口压力。
//
// 与代数校验（char_check_turbine）的关系：
//   - 代数校验固定腔压为试验值，比对喷嘴通流公式的流量误差；
//   - 本模式入口为流量源语义：引气流量 = 实测值直接注入前腔，稳态质量
//     平衡使通过流量恒等于实测值（构造上无信息量），全部模型残差转移到
//     前腔压力上——壅塞支下 m ∝ p_pl，故 p_pl/p_in 偏差 ≈ −代数流量误差。
//     两者数学等价、互为印证，且本模式额外检验了完整流路（出口节流会把
//     残差进一步放大到 p_pl 上）与时间积分链。
//
// 试验台接口语义：入口走流量源分支（A_duct <= 0 默认），注入实测流量；
// 出口为自由排气，全局默认 A_outlet = 0.01 即按试验台标定的通流面积
// （≈4× 喷嘴喉部），本模式直接使用全局参数、不再做试验台覆盖。
//
// 判定：稳态前腔压力相对试验入口压力偏差 > CHECK_TOL_PCT 或未收敛即失败。
// 数据来源：涡轮试验报告（仓库根目录 参数.md），7 个工况点按列排列。
// 校验对象为全局 turbine_params（默认值即试验标定值；--config 可覆盖）。
//
// 输出列：
//   point           工况点编号（1..7）
//   expansion_ratio 定压比（膨胀比 = p_in/p_out）
//   p_in_MPa        试验入口压力 MPa
//   T_in_K          引气总温 K
//   mdot_steady     稳态通过流量 kg/s（流量源语义下 = 实测值）
//   p_pl_steady     稳态前腔压力 MPa
//   p_err_pct       前腔压力偏差 %（(p_pl - p_in)/p_in，判定量）
//   u_over_c0       速比 u/c0（按稳态腔压计算）
//   eta_isen        U/C 多项式等熵效率
//   power_kW        指定转速下涡轮输出功率 kW（含机械效率）
//   t_conv_s        收敛耗时 s（流量漂移判据）
//   p_dn_p_back     稳态后腔/背压比（诊断出口节流，应 ≈1）

#include <math.h>
#include <stdio.h>

#include "Turbine.h"
#include "char_sweep.h"

#define CHECK_RPM      20000.0	// 校验转速 rpm（与主仿真目标转速一致）
#define CHECK_TOL_PCT  1.0		// 前腔压力偏差判定阈值 %
#define CHECK_PI       3.14159265358979323846	// 严格 -std=c11 下 glibc 不提供 M_PI

// 动态积分设置
#define DYN_DT         1e-4	// RK4 步长 s（容腔均压时间常数 ~9 ms，裕量充足）
#define DYN_T_MAX      20.0	// 积分上限 s
// 收敛判据：每 DYN_CHECK_STEPS 步比对三段流量，相对漂移 < DYN_TOL_DRIFT
// 即视为准稳态（流量源注入下稳态严格存在，漂移判据为一般性保护）
#define DYN_TOL_DRIFT  1e-6
#define DYN_CHECK_STEPS 2000	// 漂移检查窗口（0.2 s）

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

// 单工况点动态积分：两容腔压力经 RK4 演化至稳态，结果写入 dyn
typedef struct DynResult {
	double m_dot_turb;	// 稳态通过流量 kg/s
	double T_turb;		// 稳态涡轮转矩 N·m
	double t_conv;		// 收敛耗时 s（未收敛取 DYN_T_MAX）
	double p_plenum;	// 稳态前腔压力 Pa
	double p_downstream;	// 稳态后腔压力 Pa
	int converged;		// 是否在积分上限内收敛
} DynResult;

static void dyn_rhs(const TurbineParams* params, double p_in, double p_back,
                    double T_in, double m_dot_supply, double omega,
                    double p_pl, double p_dn, double* dpl, double* ddn) {
	const TurbineEval e = turbine_evaluate(params, p_in, p_back, T_in,
	                                       m_dot_supply, p_pl, p_dn, omega);
	*dpl = e.dp_plenum;
	*ddn = e.dp_downstream;
}

static DynResult run_dynamic(const TurbineParams* params,
                             const TurbCheckPoint* pt, double omega) {
	DynResult r = { 0 };
	// 引气边界：入口流量源注入实测流量，出口自由排气（全局 A_outlet 即试验台标定值）
	const double p_back = pt->p_in / pt->expansion_ratio;
	const double m_dot_supply = pt->mdot_meas;

	double p_pl = pt->p_in, p_dn = p_back;	// 初值：前腔=入口压力，后腔=背压
	double t = 0.0;
	int step = 0;
	double m_in_prev = -1.0, m_tb_prev = -1.0, m_out_prev = -1.0;
	while (t < DYN_T_MAX) {
		double k1p, k1d, k2p, k2d, k3p, k3d, k4p, k4d;
		dyn_rhs(params, pt->p_in, p_back, pt->T_in, m_dot_supply, omega,
		        p_pl, p_dn, &k1p, &k1d);
		dyn_rhs(params, pt->p_in, p_back, pt->T_in, m_dot_supply, omega,
		        p_pl + 0.5 * DYN_DT * k1p, p_dn + 0.5 * DYN_DT * k1d, &k2p, &k2d);
		dyn_rhs(params, pt->p_in, p_back, pt->T_in, m_dot_supply, omega,
		        p_pl + 0.5 * DYN_DT * k2p, p_dn + 0.5 * DYN_DT * k2d, &k3p, &k3d);
		dyn_rhs(params, pt->p_in, p_back, pt->T_in, m_dot_supply, omega,
		        p_pl + DYN_DT * k3p, p_dn + DYN_DT * k3d, &k4p, &k4d);
		p_pl += DYN_DT / 6.0 * (k1p + 2.0 * k2p + 2.0 * k3p + k4p);
		p_dn += DYN_DT / 6.0 * (k1d + 2.0 * k2d + 2.0 * k3d + k4d);
		t += DYN_DT;
		step++;

		// 漂移判据：窗口前后三段流量的相对变化均小于容差即视为准稳态
		if (step % DYN_CHECK_STEPS == 0) {
			const TurbineEval e = turbine_evaluate(params, pt->p_in, p_back, pt->T_in,
			                                       m_dot_supply, p_pl, p_dn, omega);
			if (m_in_prev >= 0.0 &&
			    fabs(e.m_dot_in - m_in_prev) <= DYN_TOL_DRIFT * e.m_dot_in &&
			    fabs(e.m_dot_turb - m_tb_prev) <= DYN_TOL_DRIFT * e.m_dot_turb &&
			    fabs(e.m_dot_out - m_out_prev) <= DYN_TOL_DRIFT * e.m_dot_out) {
				r.converged = 1;
				r.m_dot_turb = e.m_dot_turb;
				r.T_turb = e.T_turb;
				r.t_conv = t;
				r.p_plenum = p_pl;
				r.p_downstream = p_dn;
				return r;
			}
			m_in_prev = e.m_dot_in;
			m_tb_prev = e.m_dot_turb;
			m_out_prev = e.m_dot_out;
		}
	}
	// 未收敛：返回末态并标记
	const TurbineEval e = turbine_evaluate(params, pt->p_in, p_back, pt->T_in,
	                                       m_dot_supply, p_pl, p_dn, omega);
	r.m_dot_turb = e.m_dot_turb;
	r.T_turb = e.T_turb;
	r.t_conv = t;
	r.p_plenum = p_pl;
	r.p_downstream = p_dn;
	return r;
}

int char_check_turbine_dynamic(FILE* f, int n) {
	(void)n;	// 校验点固定为试验数据，与网格数无关

	const TurbineParams* p = &turbine_params;
	const double omega = CHECK_RPM * (2.0 * CHECK_PI / 60.0);
	const double cp = p->gamma * p->gas_R / (p->gamma - 1.0);

	fprintf(f, "point,expansion_ratio,p_in_MPa,T_in_K,mdot_steady,p_pl_steady_MPa,"
	           "p_err_pct,u_over_c0,eta_isen,power_kW,t_conv_s,p_dn_p_back\n");
	printf("涡轮 7 点动态稳态校验（引气流量源 + 容腔动态积分，转速 %.0f rpm）\n",
	       CHECK_RPM);
	printf("%-4s %-6s %-8s %-7s %-10s %-10s %-8s %-7s %-8s %-9s\n",
	       "点", "膨胀比", "p_in/MPa", "T_in/K",
	       "p_pl/MPa", "p偏差%", "u/c0", "eta", "功率/kW", "收敛/s");

	int failures = 0;
	for (int i = 0; i < N_CHECK_POINTS; i++) {
		const TurbCheckPoint* pt = &CHECK_POINTS[i];
		const DynResult dyn = run_dynamic(p, pt, omega);

		// 判定量：稳态前腔压力对试验入口压力的偏差（流量源语义下
		// 通过流量恒等于实测值，模型残差全部体现在前腔压力上）
		const double p_err_pct = (dyn.p_plenum / pt->p_in - 1.0) * 100.0;
		const double power = dyn.T_turb * omega;

		// U/C 速比与效率按稳态腔压计算（与 turbine_evaluate 内部公式一致）
		const double r_steady = dyn.p_downstream / dyn.p_plenum;
		const double dT_ideal = pt->T_in * (1.0 - pow(r_steady, (p->gamma - 1.0) / p->gamma));
		const double c0 = sqrt(2.0 * cp * fmax(dT_ideal, 0.0));
		const double uc = c0 > 1e-9 ? omega * p->mean_radius / c0 : 0.0;
		const double eta_isen = fmin(fmax(
			(p->eta_uc_a * uc + p->eta_uc_b) * uc + p->eta_uc_c, 0.0), 1.0);

		if (!dyn.converged || fabs(p_err_pct) > CHECK_TOL_PCT) {
			failures++;
		}
		fprintf(f, "%d,%.3f,%.3f,%.1f,%.6f,%.6f,%.4f,%.4f,%.4f,%.3f,%.4f,%.5f\n",
			i + 1, pt->expansion_ratio, pt->p_in / 1e6, pt->T_in,
			dyn.m_dot_turb, dyn.p_plenum / 1e6, p_err_pct, uc, eta_isen,
			power / 1e3, dyn.t_conv, dyn.p_downstream / (pt->p_in / pt->expansion_ratio));
		printf("%-4d %-6.2f %-8.3f %-7.1f %-10.4f %+-8.3f %-7.3f %-8.4f %-9.1f %-8.4f%s\n",
			i + 1, pt->expansion_ratio, pt->p_in / 1e6, pt->T_in,
			dyn.p_plenum / 1e6, p_err_pct, uc, eta_isen, power / 1e3, dyn.t_conv,
			!dyn.converged ? "  <-- 未收敛"
				: (fabs(p_err_pct) > CHECK_TOL_PCT ? "  <-- FAIL" : ""));
	}

	if (failures > 0) {
		printf("校验失败：%d/%d 个点未收敛或前腔压力偏差超过 %.1f%%\n",
		       failures, N_CHECK_POINTS, CHECK_TOL_PCT);
		return -1;
	}
	printf("校验通过：%d 个点均收敛且前腔压力偏差 <= %.1f%%\n",
	       N_CHECK_POINTS, CHECK_TOL_PCT);
	return 0;
}
