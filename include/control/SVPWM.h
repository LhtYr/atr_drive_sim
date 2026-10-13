#ifndef SVPWM_H
#define SVPWM_H

typedef struct SVPWMOutput {
	double duty_a;
	double duty_b;
	double duty_c;
	double gate_u_a;
	double gate_l_a;
	double gate_u_b;
	double gate_l_b;
	double gate_u_c;
	double gate_l_c;
	double modulation_scale;
} SVPWMOutput;

typedef struct SVPWM {
	double linear_limit;
} SVPWM;
extern SVPWM svpwm;

SVPWMOutput svpwm_evaluate(const SVPWM* controller,
									  double v_dc_in,
									  double ud_ref,
									  double uq_ref,
									  double theta_e);

#endif // SVPWM_H
