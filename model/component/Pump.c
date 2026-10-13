#include "Pump.h"

#include <math.h>
#include "utils/numeric_utils.h"

#define PUMP_PI 3.14159265358979323846	// 严格 -std=c11 下 glibc 不提供 M_PI
#define PUMP_G  9.80665					// 重力加速度 m/s^2

PumpParams pump_params = {
	.rho = 1000.0,			// 流体密度 kg/m^3
	.omega_floor = 1.0,		// 零速保护 rad/s

	.V_out = 1e-2,				// 出口等效容积 m^3
	.beta_eff = 1.0e9,			// 有效体积弹性模量 Pa
	.R_out = 1.0e9,			// 下游等效流阻 Pa/(m^3/s)
	.p_floor = 20000.0,		// 压力下限 Pa

	// 试验拟合特性（参数.md：H 单位 m，Q 单位 L/s，n 单位 r/min）
	.h_qq = -14.7994,		// 扬程 Q^2 项 m/(L/s)^2
	.h_qn = 2.3891e-3,		// 扬程 Q·n 项 m/((L/s)·rpm)
	.h_nn = 8.9369e-7,		// 扬程 n^2 项 m/rpm^2
	.e_x3 = -4.5796e9,		// 效率 (Q/n)^3 项
	.e_x2 = -2.3413e7,		// 效率 (Q/n)^2 项
	.e_x1 = 7580.5,			// 效率 (Q/n) 项
	.e_x0 = -0.0068608,		// 效率常数项
};

PumpView pump_view = {
	.p_discharge = 100000.0,	// 泵出口压力初值 Pa
};

PumpEval pump_evaluate(const PumpParams* params,
                       double p_suction_in,
                       double p_downstream_in,
                       double p_discharge_in,
                       double omega_shaft) {
	PumpEval out;
	const double p_suc = clamp_floor(p_suction_in, params->p_floor);
	const double p_dn = clamp_floor(p_downstream_in, params->p_floor);
	const double p_discharge_safe = clamp_floor(p_discharge_in, params->p_floor);
	const double omega_abs = fabs(omega_shaft);
	const double n_rpm = omega_abs * 30.0 / PUMP_PI;

	out.Q_out = fmax((p_discharge_safe - p_dn) / clamp_floor(params->R_out, 1.0), 0.0);

	// 试验拟合特性：H(Q,n) 为 Q 的下开口二次式。给定系统扬程 H_sys，
	// 解 a*Q^2 + b*Q + c = 0 取下降支根（泵稳定工作支）得运行流量；
	// 效率拟合 eta(Q/n) 须为正（试验包线内），否则视为无有效工况点。
	double Q_ls = 0.0;	// 泵流量 L/s
	double eta = 0.0;
	double H = params->h_nn * n_rpm * n_rpm;	// 关死扬程（无流量时的能力指示）
	if (n_rpm > 1e-6 && params->h_qq < -1e-12) {
		const double H_sys = fmax(p_discharge_safe - p_suc, 0.0)
			/ clamp_floor(params->rho * PUMP_G, 1.0);
		const double a = params->h_qq;
		const double b = params->h_qn * n_rpm;
		const double c = params->h_nn * n_rpm * n_rpm - H_sys;
		const double disc = b * b - 4.0 * a * c;
		if (disc >= 0.0) {
			const double q_desc = (-b - sqrt(disc)) / (2.0 * a);	// 下降支（较大根）
			if (q_desc > 0.0) {
				const double x = q_desc / n_rpm;
				const double e = ((params->e_x3 * x + params->e_x2) * x
					+ params->e_x1) * x + params->e_x0;
				if (e > 0.0) {
					Q_ls = q_desc;
					eta = fmin(e, 1.0);
					H = H_sys;
				}
			}
		}
	}

	out.Q_pump = Q_ls / 1000.0;	// L/s -> m^3/s
	out.eta = eta;
	out.H = H;
	out.dp_pump = params->rho * PUMP_G * H;

	const double V_safe = clamp_floor(params->V_out, 1e-6);
	out.dp_discharge = (params->beta_eff / V_safe) * (out.Q_pump - out.Q_out);

	// 轴功率 = 水力功率 / 效率；效率拟合无效时不产生转矩
	const double P_hyd = params->rho * PUMP_G * out.Q_pump * H;
	const double P_shaft = eta > 1e-9 ? P_hyd / eta : 0.0;
	const double omega_denom = clamp_floor(omega_abs, 1e-6);
	const double omega_knee = clamp_floor(params->omega_floor, 1.0);
	const double low_speed_scale = omega_abs / (omega_abs + omega_knee);
	out.T_pump = (P_shaft / omega_denom) * low_speed_scale;
	out.p_discharge = p_discharge_in;
	return out;
}

void pump_apply_state(PumpView* view,
                      const PumpParams* params,
                      double omega_mech_in,
                      double p_discharge_in,
                      const PumpEval* eval) {
	view->omega_mech = omega_mech_in;
	view->T_pump = eval->T_pump;
	view->Q_pump = eval->Q_pump;
	view->Q_out = eval->Q_out;
	view->dp_pump = eval->dp_pump;
	view->H = eval->H;
	view->p_discharge = clamp_floor(p_discharge_in, params->p_floor);
}
