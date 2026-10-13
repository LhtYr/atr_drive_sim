#include "Turbine.h"

#include <math.h>
#include "utils/numeric_utils.h"

// 默认参数即真实涡轮试验标定值（7 工况点，见 参数.md / turbine_params.txt）：
// 燃气物性 gamma=1.18、R=380 取 7 点均值圆整，喉部面积 0.00246 m^2 为试验给定，
// r_choke_nozzle=0.7 反映收扩喷嘴在高于经典临界压比时仍保持壅塞的特性。
TurbineParams turbine_params = {
	.gamma = 1.18,          // 比热比（燃气，标定值）
	.gas_R = 380.0,         // 比气体常数 J/(kg·K)（燃气，标定值）
	.eta_turb = 1.0,        // 涡轮效率（常数回退值，mean_radius <= 0 时生效）
	.omega_floor = 1.0,     // 角速度保护下限 rad/s

	.V_plenum = 1e-2,       // 前腔体积 m^3
	.V_downstream = 1e-2,   // 后腔体积 m^3
	.Cd_duct = 1.0,         // 引气管路流量系数
	.A_duct = 0.0,          // 引气管路等效导流面积 m^2（0 = 理想直连；按实际管路标定）
	.Cd_nozzle = 1.0,       // 前腔到后腔喷嘴流量系数（逐点物性反算 ≈1.001）
	.A_nozzle = 0.00246,    // 前腔到后腔喷嘴等效面积 m^2（试验喉部流通面积）
	.Cd_outlet = 1.0,       // 后腔出口流量系数
	.A_outlet = 0.01,       // 后腔出口等效面积 m^2（试验台自由排气标定值 ≈4× 喷嘴喉部）
	.p_floor = 20000.0,     // 压力下限 Pa

	.mean_radius = 0.055,   // 涡轮中径半径 m（中径 110 mm）
	.eta_mech = 0.98,       // 机械效率
	.eta_uc_a = -3.266,     // U/C 效率多项式二次项系数
	.eta_uc_b = 2.469,      // U/C 效率多项式一次项系数
	.eta_uc_c = 0.051,      // U/C 效率多项式常数项
	.r_choke_nozzle = 0.7,  // 喷嘴壅塞压比上限（收扩喷嘴标定值；<=0 取经典临界压比）
};

TurbineView turbine_view = {
	.m_dot_in = 0.05,       // 入口质量流量初值 kg/s
	.p_plenum = 400000.0,   // 前腔压力初值 Pa
	.p_downstream = 400000.0,   // 后腔压力初值 Pa
};

TurbineEval turbine_evaluate(const TurbineParams* params,
                             double p_in,
                             double p_out,
                             double T_in,
                             double m_dot,
                             double p_plenum_in,
                             double p_downstream_in,
                             double omega_shaft) {
	TurbineEval out;
	const double gamma_safe = clamp_floor(params->gamma, 1.01);
	const double T_gas = clamp_floor(T_in, 200.0);
	const double p_upstream = clamp_floor(p_in, params->p_floor);
	const double p_back = clamp_floor(p_out, params->p_floor);
	const double p_plenum_safe = clamp_floor(p_plenum_in, params->p_floor);
	const double p_downstream_safe = clamp_floor(p_downstream_in, params->p_floor);
	const double critical_ratio = pow(2.0 / (gamma_safe + 1.0), gamma_safe / (gamma_safe - 1.0));

	// 入口（流量源语义）：涡轮由引气驱动，引气过程给定总压/总温/流量上限
	// (p_in, T_in, m_dot)。A_duct <= 0 为理想直连，引气流量全部注入前腔，
	// 前腔压力由喷嘴通流质量平衡自定；A_duct > 0 时经孔板导流模型计入
	// 管路压降，并以引气流量上限封顶（引气能力不足时实际注入量回落）。
	if (params->A_duct <= 0.0) {
		out.m_dot_in = fmax(m_dot, 0.0);
	} else {
		const double duct_ratio = clamp_range(p_plenum_safe / p_upstream, 0.0, 1.0);
		double duct_coeff;
		if (duct_ratio <= critical_ratio) {
			duct_coeff = sqrt(gamma_safe / (params->gas_R * T_gas))
				* pow(2.0 / (gamma_safe + 1.0), (gamma_safe + 1.0) / (2.0 * (gamma_safe - 1.0)));
		} else {
			const double term = pow(duct_ratio, 2.0 / gamma_safe)
				- pow(duct_ratio, (gamma_safe + 1.0) / gamma_safe);
			duct_coeff = sqrt(fmax(2.0 * gamma_safe / (params->gas_R * T_gas * (gamma_safe - 1.0)) * term, 0.0));
		}
		const double m_dot_duct = params->Cd_duct * params->A_duct * p_upstream * duct_coeff;
		out.m_dot_in = fmin(fmax(m_dot, 0.0), m_dot_duct);
	}

	const double pressure_ratio = clamp_range(p_downstream_safe / p_plenum_safe, 0.0, 1.0);
	// 喷嘴壅塞压比上限：收扩喷嘴在高于经典临界压比时仍可保持喉部壅塞
	const double r_choke = params->r_choke_nozzle > 0.0 ? params->r_choke_nozzle : critical_ratio;

	if (pressure_ratio <= r_choke) {
		const double choked_coeff = sqrt(gamma_safe / (params->gas_R * T_gas))
			* pow(2.0 / (gamma_safe + 1.0), (gamma_safe + 1.0) / (2.0 * (gamma_safe - 1.0)));
		out.m_dot_turb = params->Cd_nozzle * params->A_nozzle * p_plenum_safe * choked_coeff;
	} else {
		const double term = pow(pressure_ratio, 2.0 / gamma_safe)
			- pow(pressure_ratio, (gamma_safe + 1.0) / gamma_safe);
		const double subsonic_coeff = sqrt(fmax(2.0 * gamma_safe / (params->gas_R * T_gas * (gamma_safe - 1.0)) * term, 0.0));
		out.m_dot_turb = params->Cd_nozzle * params->A_nozzle * p_plenum_safe * subsonic_coeff;
	}

	if (p_downstream_safe <= p_back) {
		out.m_dot_out = 0.0;
	} else {
		const double outlet_ratio = clamp_range(p_back / p_downstream_safe, 0.0, 1.0);
		if (outlet_ratio <= critical_ratio) {
			const double choked_coeff = sqrt(gamma_safe / (params->gas_R * T_gas))
				* pow(2.0 / (gamma_safe + 1.0), (gamma_safe + 1.0) / (2.0 * (gamma_safe - 1.0)));
			out.m_dot_out = params->Cd_outlet * params->A_outlet * p_downstream_safe * choked_coeff;
		} else {
			const double term = pow(outlet_ratio, 2.0 / gamma_safe)
				- pow(outlet_ratio, (gamma_safe + 1.0) / gamma_safe);
			const double subsonic_coeff = sqrt(fmax(2.0 * gamma_safe / (params->gas_R * T_gas * (gamma_safe - 1.0)) * term, 0.0));
			out.m_dot_out = params->Cd_outlet * params->A_outlet * p_downstream_safe * subsonic_coeff;
		}
	}

	const double V_plenum_safe = clamp_floor(params->V_plenum, 1e-6);
	const double V_downstream_safe = clamp_floor(params->V_downstream, 1e-6);
	out.dp_plenum = (params->gas_R * T_gas / V_plenum_safe) * (out.m_dot_in - out.m_dot_turb);
	out.dp_downstream = (params->gas_R * T_gas / V_downstream_safe) * (out.m_dot_turb - out.m_dot_out);

	out.p_plenum = p_plenum_in;
	out.p_downstream = p_downstream_in;

	if (p_plenum_safe <= p_downstream_safe) {
		out.T_turb = 0.0;
		return out;
	}

	const double cp = gamma_safe * params->gas_R / (gamma_safe - 1.0);
	const double T_out_ideal = T_gas * pow(p_downstream_safe / p_plenum_safe, (gamma_safe - 1.0) / gamma_safe);
	const double delta_T_ideal = fmax(T_gas - T_out_ideal, 0.0);

	// 等熵效率：mean_radius > 0 时用 U/C 多项式，否则回退常数 eta_turb
	double eta_isen = params->eta_turb;
	if (params->mean_radius > 0.0) {
		const double c0 = sqrt(2.0 * cp * delta_T_ideal);  // 等熵喷射速度 m/s
		if (c0 > 1e-9) {
			const double uc = fabs(omega_shaft) * params->mean_radius / c0;
			eta_isen = clamp_range(
				(params->eta_uc_a * uc + params->eta_uc_b) * uc + params->eta_uc_c, 0.0, 1.0);
		} else {
			eta_isen = 0.0;
		}
	}

	const double power_turbine = out.m_dot_turb * cp * delta_T_ideal * eta_isen;
	out.T_turb = params->eta_mech * power_turbine / magnitude_floor(omega_shaft, params->omega_floor);
	return out;
}

void turbine_apply_state(TurbineView* view,
                         const TurbineParams* params,
                         double p_plenum_in,
                         double p_downstream_in,
                         const TurbineEval* eval) {
	view->p_plenum = clamp_floor(p_plenum_in, params->p_floor);
	view->p_downstream = clamp_floor(p_downstream_in, params->p_floor);
	view->T_turb = eval->T_turb;
	view->m_dot_in = eval->m_dot_in;
	view->m_dot_turb = eval->m_dot_turb;
	view->m_dot_out = eval->m_dot_out;
}
