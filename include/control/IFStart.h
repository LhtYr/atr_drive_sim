#ifndef IFSTART_H
#define IFSTART_H

// I/F 开环起动输出（供 SensorlessMgr 在阶段 0/1 使用）
typedef struct IFStartOutput {
    double theta_if_e;    // 强制电角度（含阻尼补偿）rad [0,2π)，注入 Park/SVPWM
    double omega_if_e;    // 强制电角速度（斜坡，未含补偿）rad/s，电流环前馈用
    double omega_if_m;    // 强制机械角速度 rad/s
    double id_ref;        // d 轴电流指令 A（恒 0）
    double iq_ref;        // q 轴电流指令 A（恒 i_start）
    int at_switch;        // 强制斜坡是否已到达 omega_switch（供切换状态机请求切换）
} IFStartOutput;

// I/F 开环起动器：参数 + 运行时状态（纯函数 + 全局单例，风格对齐现有控制器）。
// 生成随斜坡旋转的强制角 θ_if，并用观测器估计转速与强制转速的偏差注入阻尼（角度修正，
// 物理依据见 README 与 IFStart.c 顶部推导：角度修正才真正增大阻尼比，速率修正只增刚度）。
typedef struct IFStart {
    double Ts;              // 离散步长 s（= control_dt）
    int p;                  // 极对数
    double i_start;         // 强制 q 轴电流幅值 A（id=0, iq=i_start）
    double omega_switch_m;  // 切换机械角速度 rad/s（强制斜坡目标，达到后请求切换）
    double ramp_time;       // 强制转速斜坡时间 s（0 → omega_switch_m）
    double damp_k;          // 阻尼补偿时间常数 s（角度修正 θ_if -= damp_k·(ω̂_e−ω_if_e)）
    double damp_lpf_hz;     // 阻尼用 ω̂ 低通截止 Hz：低速反电动势小、PLL 估计含 ~266 Hz 纹波，
                            // 直接进 corr 会把纹波注入 θ_if（实测 δ/转速 ~330 Hz 抖动、δ 跌破 −90°）；
                            // 低通后只保留机械模态（ω_n≈26.5 Hz）附近的真实转速偏差供阻尼
    double theta_init_e;    // 强制角初值 rad（预定位假设：θ_true(0)≈0 时取 −π/2 使 δ(0)=−90° 平衡）
    // 运行时状态
    double omega_if_m;      // 当前强制机械角速度 rad/s
    double theta_ramp_e;    // 斜坡积分电角度（未补偿）rad
    double omega_est_filt_e;// 低通后的估计电角速度 rad/s（阻尼用运行时状态）
} IFStart;

extern IFStart ifstart;

// 复位 I/F 起动器（起动前调用）
void ifstart_reset(IFStart* self);

// 单步更新：omega_est_e 为 STSMO/PLL 给出的估计电角速度（阻尼用，非真值）。
// 返回强制角、强制转速与 I/F 阶段电流指令。
IFStartOutput ifstart_evaluate(IFStart* self, double omega_est_e);

#endif // IFSTART_H
