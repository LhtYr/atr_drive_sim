#ifndef SPEED_H
#define SPEED_H

typedef struct SpeedOutput {
	double e_omega;	// 速度误差
	double iq_ref;	// q轴电流参考值
} SpeedOutput;

typedef struct Speed {
	double Ts;
	double kp;
	double ki;
	double iq_ref_min;
	double iq_ref_max;
	double anti_windup_gain;
	double int_omega;
} Speed;
extern Speed speed;

SpeedOutput speed_evaluate(Speed* controller,
									  double omega_ref_mech,
									  double omega_mech);

#endif // SPEED_H
