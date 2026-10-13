#include "SVPWM.h"

#include <math.h>
#include "utils/numeric_utils.h"
#include "utils/transforms.h"

SVPWM svpwm = {
	.linear_limit = 0.0,
};

SVPWMOutput svpwm_evaluate(const SVPWM* controller,
									  double v_dc_in,
									  double ud_ref,
									  double uq_ref,
									  double theta_e) {
	SVPWMOutput out;
	const double v_dc_safe = magnitude_floor(v_dc_in, 1.0);

	const double cos_theta = cos(theta_e);
	const double sin_theta = sin(theta_e);
	// 逆 Park + 逆 Clarke：dq 电压参考 -> 三相参考（共用实现）
	double u_alpha = 0.0;
	double u_beta = 0.0;
	double ua_ref = 0.0;
	double ub_ref = 0.0;
	double uc_ref = 0.0;
	dq_to_abc(ud_ref, uq_ref, cos_theta, sin_theta,
	          &u_alpha, &u_beta, &ua_ref, &ub_ref, &uc_ref);

	const double v_max_ref = fmax(fmax(ua_ref, ub_ref), uc_ref);
	const double v_min_ref = fmin(fmin(ua_ref, ub_ref), uc_ref);
	const double v_offset = -0.5 * (v_max_ref + v_min_ref);

	double ua_cmd = ua_ref + v_offset;
	double ub_cmd = ub_ref + v_offset;
	double uc_cmd = uc_ref + v_offset;

	const double v_phase_limit = controller->linear_limit * v_dc_safe;
	const double max_abs_cmd = clamp_floor(fmax(fmax(fabs(ua_cmd), fabs(ub_cmd)), fabs(uc_cmd)), 1e-9);
	out.modulation_scale = clamp_range(v_phase_limit / max_abs_cmd, 0.0, 1.0);

	ua_cmd *= out.modulation_scale;
	ub_cmd *= out.modulation_scale;
	uc_cmd *= out.modulation_scale;

	out.duty_a = clamp_range(0.5 + ua_cmd / v_dc_safe, 0.0, 1.0);
	out.duty_b = clamp_range(0.5 + ub_cmd / v_dc_safe, 0.0, 1.0);
	out.duty_c = clamp_range(0.5 + uc_cmd / v_dc_safe, 0.0, 1.0);

	out.gate_u_a = out.duty_a;
	out.gate_l_a = 1.0 - out.duty_a;
	out.gate_u_b = out.duty_b;
	out.gate_l_b = 1.0 - out.duty_b;
	out.gate_u_c = out.duty_c;
	out.gate_l_c = 1.0 - out.duty_c;
	return out;
}
