#ifndef MOTOR_H
#define MOTOR_H

// 电机电气侧评估输出（纯计算结果，不含任何状态）
typedef struct ElectricalEval{
    double did;
    double diq;
    double theta_e;
    double omega_e;
    double ud;
    double uq;
    double i_alpha;
    double i_beta;
    double ia;
    double ib;
    double ic;
    double T_em;
} ElectricalEval;

// 电机参数（仿真期间保持不变）
typedef struct MotorParams {
    double Rs;          // 定子电阻 Ohm
    double Ld;          // d 轴电感 H
    double Lq;          // q 轴电感 H
    double psi_f;       // 永磁体磁链 Wb
    int p;              // 极对数
} MotorParams;

// 电机观测量（由求解器每步回写的运行时状态，供控制与日志读取）
typedef struct MotorView {
    double ua, ub, uc;      // 三相电压 V
    double ia, ib, ic;      // 三相电流 A
    double id, iq;          // dq 轴电流 A
    double ud, uq;          // dq 轴电压 V
    double i_alpha, i_beta; // Clarke 变换后电流 A
    double u_alpha, u_beta; // Clarke 变换后电压 V
    double omega_e;         // 电角速度 rad/s
    double theta_e;         // 电角度 rad
    double T_em;            // 电磁转矩 N·m
} MotorView;

extern MotorParams motor_params;
extern MotorView motor_view;

// 电机电气侧评估（纯计算，不修改任何对象状态）
ElectricalEval motor_evaluate_electrical(const MotorParams* params,
                                         double id_in,
                                         double iq_in,
                                         double theta_mech,
                                         double omega_mech,
                                         double ua_in,
                                         double ub_in,
                                         double uc_in);

// 用积分后的系统状态回写电机观测量：
// ua_in/ub_in/uc_in 为逆变器输出电压；eval 提供 dq 电压与电流观测量。
void motor_apply_state(MotorView* view,
                       double ua_in,
                       double ub_in,
                       double uc_in,
                       double id_in,
                       double iq_in,
                       const ElectricalEval* eval);

#endif // MOTOR_H
