#include "Current.h"

#include <assert.h>
#include <math.h>
#include "utils/numeric_utils.h"
Current current = {
	.Ts = 5e-5,
	.kp_d = 0.692,
	.ki_d = 102.0,
	.kp_q = 0.692,
	.ki_q = 102.0,
	.uq_ref_max = 540.0,
	.uq_ref_min = -540.0,
	.ud_ref_max = 540.0,
	.ud_ref_min = -540.0,
	.anti_windup_gain = 0,
	.int_d = 0.0,	// 积分初始值
	.int_q = 0.0,	// 积分初始值
};

CurrentOutput current_evaluate(Current* controller,
								const MotorParams* motor_params,
								double id_ref,
								double iq_ref,
								double id_meas,
								double iq_meas,
								double omega_e,
								double v_dc,
								double linear_limit) {
	assert(controller);
	assert(motor_params);

	CurrentOutput out;
	out.ed = id_ref - id_meas;
	out.eq = iq_ref - iq_meas;

	const double vd_pi = controller->kp_d * out.ed + controller->int_d;
	const double vq_pi = controller->kp_q * out.eq + controller->int_q;

	const double Rs = motor_params->Rs;
	const double Ld = motor_params->Ld;
	const double Lq = motor_params->Lq;
	const double psi_f = motor_params->psi_f;

	const double vd_ff = Rs * id_meas - omega_e * Lq * iq_meas;
	const double vq_ff = Rs * iq_meas + omega_e * Ld * id_meas + omega_e * psi_f;

	const double ud_unsat = vd_pi + vd_ff;
	const double uq_unsat = vq_pi + vq_ff;

	const double v_phase_limit = magnitude_floor(v_dc, 1.0) * clamp_floor(linear_limit, 0.0);
	const double v_mag = sqrt(ud_unsat * ud_unsat + uq_unsat * uq_unsat);
	out.modulation_scale = (v_mag > 1e-9) ? clamp_range(v_phase_limit / v_mag, 0.0, 1.0) : 1.0;

	const double ud_sat = ud_unsat * out.modulation_scale;
	const double uq_sat = uq_unsat * out.modulation_scale;

	out.ud_ref = ud_sat;
	out.uq_ref = uq_sat;

	const double du_d = ud_sat - ud_unsat;
	const double du_q = uq_sat - uq_unsat;

	controller->int_d += (controller->ki_d * out.ed + controller->anti_windup_gain * du_d) * controller->Ts;
	controller->int_q += (controller->ki_q * out.eq + controller->anti_windup_gain * du_q) * controller->Ts;

	return out;
}
