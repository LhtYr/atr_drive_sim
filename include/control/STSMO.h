#ifndef STSMO_H
#define STSMO_H

#include "component/Motor.h"

// 超螺旋滑模电流观测器（STSMO）输出：αβ 静止系。
// 观测电流 î_αβ 与实测电流收敛后，超螺旋切换项 z_αβ ≈ 反电动势 e_αβ，
// 经一阶低通后送 PLL 做位置估计。
typedef struct STSMOOutput {
    double emf_alpha;    // 低通滤波后的反电动势 α 分量 V
    double emf_beta;     // 低通滤波后的反电动势 β 分量 V
    double emf_mag;      // 反电动势幅值 V（归一化/诊断用）
    double i_hat_alpha;  // 观测电流 α A（诊断）
    double i_hat_beta;   // 观测电流 β A（诊断）
} STSMOOutput;

// 超螺旋滑模观测器：增益/低通为参数，其余为运行时状态（纯函数 + 全局单例）。
// 假设表贴式 PMSM（Ld≈Lq=L），αβ 系电感取 Ld。
typedef struct STSMO {
    double Ts;            // 离散步长 s（= control_dt）
    double k1;            // 超螺旋比例增益（V/A^0.5）
    double k2;            // 超螺旋积分增益（V/s）
    double emf_lpf_hz;    // 反电动势低通截止频率 Hz（抑制切换抖振）
    // 运行时状态
    double i_hat_alpha;   // 观测电流 α A
    double i_hat_beta;    // 观测电流 β A
    double v_alpha;       // 超螺旋积分项 α V（收敛后 ≈ 反电动势 α）
    double v_beta;        // 超螺旋积分项 β V
    double emf_alpha;     // 低通后反电动势 α V
    double emf_beta;      // 低通后反电动势 β V
} STSMO;

extern STSMO stsmo;

// 复位观测器状态（起动前调用）
void stsmo_reset(STSMO* obs);

// 单步更新：u_αβ 为上一控制拍施加的定子电压指令（合成自 ud/uq 与当时的控制角），
// i_αβ 为当前实测电流（Clarke 变换自三相传感器测量）。返回低通后的反电动势估计。
STSMOOutput stsmo_evaluate(STSMO* obs,
                           const MotorParams* motor_params,
                           double u_alpha,
                           double u_beta,
                           double i_alpha,
                           double i_beta);

#endif // STSMO_H
