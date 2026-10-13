#ifndef TURBINE_H
#define TURBINE_H

typedef struct TurbineEval {
	double dp_plenum;       // 前腔压力导数 Pa/s
	double dp_downstream;   // 后腔压力导数 Pa/s
	double T_turb;         // 涡轮驱动转矩 N·m（符号约定：正值 = 涡轮输出机械能）
	double m_dot_in;       // 入口质量流量 kg/s
	double m_dot_turb;     // 通过涡轮质量流量 kg/s
	double m_dot_out;      // 出口质量流量 kg/s
	double p_plenum;       // 腔体压力 Pa
	double p_downstream;   // 下游压力 Pa
} TurbineEval;

// 涡轮参数（仿真期间保持不变）
typedef struct TurbineParams {
	double gamma;          // 比热比
	double gas_R;          // 比气体常数 J/(kg·K)
	double eta_turb;       // 涡轮效率（常数回退值，mean_radius <= 0 时生效）
	double omega_floor;    // 角速度保护下限 rad/s

	double V_plenum;       // 前腔体积 m^3
	double V_downstream;   // 后腔体积 m^3
	double Cd_duct;        // 引气管路流量系数
	double A_duct;         // 引气管路等效导流面积 m^2（<= 0 = 理想直连，无压降）
	double Cd_nozzle;      // 前腔到后腔喷嘴流量系数
	double A_nozzle;       // 前腔到后腔喷嘴等效面积 m^2
	double Cd_outlet;      // 后腔出口流量系数
	double A_outlet;       // 后腔出口等效面积 m^2
	double p_floor;        // 压力下限 Pa

	double mean_radius;    // 涡轮中径半径 m（<= 0 时回退常数效率 eta_turb）
	double eta_mech;       // 机械效率（乘在输出转矩上）
	double eta_uc_a;       // U/C 效率多项式二次项系数：eta = a·x² + b·x + c，x = u/c0
	double eta_uc_b;       // U/C 效率多项式一次项系数
	double eta_uc_c;       // U/C 效率多项式常数项
	double r_choke_nozzle; // 喷嘴维持壅塞的压比上限（<= 0 时取经典临界压比）
} TurbineParams;

// 涡轮观测量（由求解器每步回写的运行时状态）
typedef struct TurbineView {
	double T_turb;         // 涡轮驱动转矩 N·m（正值 = 涡轮输出机械能）
	double m_dot_in;       // 入口质量流量 kg/s
	double m_dot_turb;     // 通过涡轮质量流量 kg/s
	double m_dot_out;      // 出口质量流量 kg/s
	double p_plenum;       // 腔体压力 Pa
	double p_downstream;   // 下游压力 Pa
} TurbineView;

extern TurbineParams turbine_params;
extern TurbineView turbine_view;

// 涡轮气动评估（纯计算）：
// p_in/p_out/T_in/m_dot 为引气边界输入（引气总压/环境背压/引气总温/引气流量
// 上限），逐拍由调用方提供；入口为流量源语义，前腔压力由质量平衡自定
TurbineEval turbine_evaluate(const TurbineParams* params,
                             double p_in,
                             double p_out,
                             double T_in,
                             double m_dot,
                             double p_plenum_in,
                             double p_downstream_in,
                             double omega_shaft);

void turbine_apply_state(TurbineView* view,
                         const TurbineParams* params,
                         double p_plenum_in,
                         double p_downstream_in,
                         const TurbineEval* eval);

#endif // TURBINE_H
