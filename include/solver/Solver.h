#ifndef SOLVER_H
#define SOLVER_H

#include "component/Motor.h"
#include "component/Shaft.h"
#include "component/Pump.h"
#include "component/Turbine.h"
#include "component/Battery.h"

typedef struct SolverSwitches {
    int enable_inverter;
    int enable_motor;
    int enable_shaft;
    int enable_pump;
    int enable_turbine;
    int enable_battery;
    int enable_stage_iteration;
    double disabled_pump_torque;
    double disabled_turbine_torque;
    double disabled_battery_voltage;
} SolverSwitches;

// 系统“动态状态”向量：会被 RK4 积分推进的变量。
typedef struct SystemState {
    double id;  // dq 轴 d 分量电流 A
    double iq;  // dq 轴 q 分量电流 A
    double theta_mech;      // 机械角位置 rad
    double omega_mech;      // 机械角速度 rad/s
    double p_discharge;     // 泵出口压力 Pa
    double p_plenum;        // 涡轮前腔压力 Pa
    double p_downstream;    // 涡轮后腔压力 Pa
    double soc_batt;        // 电池 SOC
    double v_rc1_batt;      // 电池 RC 支路1极化电压 V
    double v_rc2_batt;      // 电池 RC 支路2极化电压 V
} SystemState;

// 包含系统的代数变量，用于 RK4 步进计算中间状态评估和模型状态同步
typedef struct SystemAlgebraic {
    ElectricalEval motor;
    MechanicalEval shaft;
    PumpEval pump;
    TurbineEval turbine;
    BatteryEval battery;
} SystemAlgebraic;

// 包含系统状态的时间导数，用于 RK4 步进计算
typedef struct SystemDeriv {
    double did;   // id 的时间导数
    double diq;   // iq 的时间导数
    double dtheta_mech;   // 机械角位置的时间导数（即机械角速度）
    double domega_mech;   // 机械角速度的时间导数（即机械角加速度）
    double dp_discharge;  // 泵出口压力时间导数
    double dp_plenum;     // 涡轮前腔压力时间导数
    double dp_downstream; // 涡轮后腔压力时间导数
    double dsoc_batt;     // 电池 SOC 时间导数
    double dv_rc1_batt;   // 电池 RC 支路1极化电压时间导数
    double dv_rc2_batt;   // 电池 RC 支路2极化电压时间导数
} SystemDeriv;

// RK4 步进器需要的系统评估结果：包含状态导数和代数变量
typedef struct EvalResult {
    SystemDeriv deriv;      // 状态导数
    SystemAlgebraic alg;    // 代数变量（组件评估结果）
 } EvalResult;

// 系统外部输入：控制参考 + 气路/液压边界条件
typedef struct SystemInput {
    double ud_ref;  // 电流环 d 轴电压参考 V
    double uq_ref;  // 电流环 q 轴电压参考 V

    // SVPWM 逆 Park 电压合成角通道（无感控制用）：
    // use_ctrl_theta=0（默认，有感）时 SVPWM 仍用真实电角 p*theta_mech，行为与改动前完全一致；
    // =1（无感）时用控制侧注入的冻结估计角 theta_ctrl_e（rad，数字控制一拍一更新语义）。
    int use_ctrl_theta;
    double theta_ctrl_e;

    double turb_p_in;   // 涡轮引气总压 Pa
    double turb_p_out;  // 涡轮出口背压 Pa
    double turb_T_in;   // 涡轮引气总温 K
    double turb_m_dot;  // 涡轮引气流量上限 kg/s（流量源语义，0 = 无供气）

    double pump_p_suction;      // 泵吸入口压力 Pa
    double pump_p_downstream;   // 泵下游压力 Pa

    double batt_p_load;         // 电池附加负载功率 W（使能电池时与逆变器功率叠加）
} SystemInput;

// 边界条件默认值（常量剖面的来源）：可经参数文件覆盖（键名边界 boundary.*，
// 见 src/param_file.c），须在构造 BoundaryProfile 之前完成覆盖
extern SystemInput boundary_defaults;

extern SolverSwitches solver_switches;

void solver_set_switches(SolverSwitches switches);

// 用各组件当前观测量初始化求解器内部状态向量（须在首次 step 之前调用）
void solver_init_state(void);

EvalResult evaluate(const SystemState state, const SystemInput input, double sim_time_s, double dt);
// 推进一个积分步：t_s 为步起点时间，dt 为积分步长（均由调用方显式提供）
void step(double t_s, double dt, const SystemInput input);

#endif // SOLVER_H