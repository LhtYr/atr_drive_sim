#ifndef CURRENT_H
#define CURRENT_H

#include "component/Motor.h"

// 电流环控制器输出
typedef struct CurrentOutput {
	double ed;
	double eq;
	double ud_ref;
	double uq_ref;
	double modulation_scale;
} CurrentOutput;

// 电流环控制器：增益/限幅为参数，int_d/int_q 为运行时积分状态
typedef struct Current {
	double Ts;
	double kp_d;
	double ki_d;
	double kp_q;
	double ki_q;
	double uq_ref_max;
	double uq_ref_min;
	double ud_ref_max;
	double ud_ref_min;
	double anti_windup_gain;
	double int_d;
	double int_q;
} Current;
extern Current current;

CurrentOutput current_evaluate(Current* controller,
								const MotorParams* motor_params,
								double id_ref,
								double iq_ref,
								double id_meas,
								double iq_meas,
								double omega_e,
								double v_dc,
								double linear_limit);

#endif // CURRENT_H
