#ifndef INVERTER_H
#define INVERTER_H

#include "SVPWM.h"

typedef struct InverterEval {
	double v_dc;
	double duty_a;
	double duty_b;
	double duty_c;
	double gate_u_a;
	double gate_l_a;
	double gate_u_b;
	double gate_l_b;
	double gate_u_c;
	double gate_l_c;
	double ua;
	double ub;
	double uc;
	double i_inv_dc;
	double p_inv;
} InverterEval;

// 逆变器参数（仿真期间保持不变）
typedef struct InverterParams {
	double f_sw;		// 开关频率 Hz（开关模型使用）
	int use_switching_model;	// 1 = 开关模型（载波比较），0 = 平均模型（默认，更快）
} InverterParams;

// 逆变器观测量（由求解器每步回写的运行时状态）
typedef struct InverterView {
	double v_dc;
	double ua;
	double ub;
	double uc;
	double duty_a;
	double duty_b;
	double duty_c;
	double gate_u_a;
	double gate_l_a;
	double gate_u_b;
	double gate_l_b;
	double gate_u_c;
	double gate_l_c;
	double i_inv_dc;
	double p_inv;
} InverterView;

extern InverterParams inverter_params;
extern InverterView inverter_view;

// 逆变器评估（纯计算）：v_dc 为当前直流母线电压，
// use_switching_model 选择开关模型/平均模型，sim_time_s 仅开关模型使用
InverterEval inverter_evaluate(double v_dc,
                               const SVPWMOutput* svpwm,
                               double ia,
                               double ib,
                               double ic,
                               double sim_time_s,
                               int use_switching_model,
                               double f_sw);

void inverter_apply_state(InverterView* view, const InverterEval* eval);

#endif // INVERTER_H
