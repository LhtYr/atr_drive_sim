#ifndef PLL_H
#define PLL_H

// 锁相环位置估计输出
typedef struct PLLOutput {
    double theta_e;   // 估计电角度 rad [0,2π)（已补偿反电动势低通相位滞后）
    double omega_e;   // 估计电角速度 rad/s（低通后，抑制切换纹波）
    double omega_m;   // 估计机械角速度 rad/s（= omega_e / p，供速度环反馈）
    double eps;       // 相位误差 ε = sin(θ_e − θ̂)（诊断，锁定后 →0）
} PLLOutput;

// 锁相环：增益/限幅为参数，theta_e/int_pll 为运行时状态（纯函数 + 全局单例）。
// 输入为 STSMO 低通后的反电动势 e_αβ，输出估计电角度/角速度。
typedef struct PLL {
    double Ts;            // 离散步长 s（= control_dt）
    double kp;            // 比例增益
    double ki;            // 积分增益
    double omega_e_max;   // 估计电角速度限幅 rad/s（未收敛时防发散）
    double emf_lpf_wc;    // STSMO 反电动势低通角频率 rad/s（相位滞后补偿，须与 stsmo.emf_lpf_hz 一致）
    double omega_lpf_hz;  // 转速估计输出低通截止 Hz（滤除超螺旋切换纹波，供速度环反馈/日志）
    int p;                // 极对数（换算机械角速度）
    // 运行时状态
    double theta_e;       // 内部锁定角 rad（锁定到低通后反电动势，未补偿）
    double int_pll;       // PI 积分项（≈ 估计电角速度）
    double omega_filt;    // 低通后的估计电角速度 rad/s（输出用，抑制纹波）
} PLL;

extern PLL pll;

// 复位 PLL（theta0 为初始角估计，通常取 0 或强制起始角）
void pll_reset(PLL* self, double theta0);

// 单步更新：输入 αβ 反电动势估计，输出估计角/角速度。
PLLOutput pll_evaluate(PLL* self, double emf_alpha, double emf_beta);

#endif // PLL_H
