#ifndef SHAFT_H
#define SHAFT_H

typedef struct MechanicalEval{
    double domega_mech;     // 机械角速度的时间导数（即机械角加速度）
    double dtheta_mech;     // 机械角位置的时间导数（即机械角速度）
} MechanicalEval;

// 轴系参数（仿真期间保持不变）
typedef struct ShaftParams {
    double J_motor;     // 电机侧转动惯量 kg.m²
    double J_turb;      // 涡轮侧转动惯量 kg.m²
    double J_pump;      // 泵侧转动惯量 kg.m²

    double b_motor;     // 电机侧摩擦系数（N·m·s）
    double b_turb;      // 涡轮侧摩擦系数（N·m·s）
    double b_pump;      // 泵侧摩擦系数（N·m·s）
} ShaftParams;

// 轴系观测量（由求解器每步回写的运行时状态）
typedef struct ShaftView {
    double omega_mech;  // 机械角速度 rad/s
    double theta_mech;  // 机械角位置 rad
} ShaftView;

extern ShaftParams shaft_params;
extern ShaftView shaft_view;

// 轴系力学评估（纯计算）：转矩符号约定见各组件头文件
MechanicalEval shaft_evaluate_mechanics(const ShaftParams* params,
                                        double omega_mech_in,
                                        double T_em,
                                        double T_turb,
                                        double T_pump);

void shaft_apply_state(ShaftView* view, double theta_mech_in, double omega_mech_in);

#endif // SHAFT_H
