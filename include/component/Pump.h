#ifndef PUMP_H
#define PUMP_H

typedef struct PumpEval {
	double dp_discharge; // 出口压力导数 Pa/s
	double T_pump;       // 泵负载转矩 N·m（符号约定：正值 = 泵吸收机械能）
	double Q_pump;       // 泵流量 m^3/s
	double Q_out;        // 出流量 m^3/s
	double dp_pump;      // 泵压升 Pa
	double H;            // 扬程 m（有流量时为系统扬程，否则为关死扬程）
	double eta;          // 效率拟合值（超出试验包线为 0）
	double p_discharge;  // 泵出口压力 Pa
} PumpEval;

// 泵参数（仿真期间保持不变）
typedef struct PumpParams {
	double rho;          // 流体密度 kg/m^3
	double omega_floor;  // 零速保护 rad/s

	double V_out;        // 出口等效容积 m^3
	double beta_eff;     // 有效体积弹性模量 Pa
	double R_out;        // 下游等效流阻 Pa/(m^3/s)
	double p_floor;      // 压力下限 Pa

	// 试验拟合特性（参数.md：H 单位 m，Q 单位 L/s，n 单位 r/min）
	// H = h_qq*Q^2 + h_qn*Q*n + h_nn*n^2
	double h_qq;
	double h_qn;
	double h_nn;
	// eta = e_x3*(Q/n)^3 + e_x2*(Q/n)^2 + e_x1*(Q/n) + e_x0
	double e_x3;
	double e_x2;
	double e_x1;
	double e_x0;
} PumpParams;

// 泵观测量（由求解器每步回写的运行时状态）
typedef struct PumpView {
	double omega_mech;   // 机械角速度 rad/s
	double T_pump;       // 泵负载转矩 N·m（正值 = 泵吸收机械能）
	double Q_pump;       // 泵流量 m^3/s
	double Q_out;        // 出流量 m^3/s
	double dp_pump;      // 泵压升 Pa
	double H;            // 扬程 m
	double p_discharge;  // 泵出口压力 Pa
} PumpView;

extern PumpParams pump_params;
extern PumpView pump_view;

// 泵水力评估（纯计算）：吸入/下游压力为外部边界输入，逐拍由调用方提供
PumpEval pump_evaluate(const PumpParams* params,
                       double p_suction_in,
                       double p_downstream_in,
                       double p_discharge_in,
                       double omega_shaft);

void pump_apply_state(PumpView* view,
                      const PumpParams* params,
                      double omega_mech_in,
                      double p_discharge_in,
                      const PumpEval* eval);

#endif // PUMP_H
