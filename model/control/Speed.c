#include "Speed.h"

#include <math.h>
Speed speed = {
	.Ts = 5e-5,
	.kp = 0.15,
	.ki = 0.17,
	.iq_ref_min = -540.0,
	.iq_ref_max = 540.0,
	.anti_windup_gain = 0.0,
	.int_omega = 0.0,
};

SpeedOutput speed_evaluate(Speed* controller,
									  double omega_ref_mech,
									  double omega_mech) {
	SpeedOutput out;
	out.e_omega = omega_ref_mech - omega_mech;

	const double iq_unsat = controller->kp * out.e_omega + controller->int_omega;
	if (iq_unsat < controller->iq_ref_min) {
		out.iq_ref = controller->iq_ref_min;
	} else if (iq_unsat > controller->iq_ref_max) {
		out.iq_ref = controller->iq_ref_max;
	} else {
		out.iq_ref = iq_unsat;
	}

	const double diq_aw = out.iq_ref - iq_unsat;
	controller->int_omega += (controller->ki * out.e_omega + controller->anti_windup_gain * diq_aw) * controller->Ts;
	return out;
}
