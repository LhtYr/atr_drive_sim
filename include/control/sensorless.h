#ifndef SENSORLESS_H
#define SENSORLESS_H

#include "component/Motor.h"

// 无感起动状态机阶段。
typedef enum SensorlessStage {
    SL_STAGE_IF = 0,      // I/F 开环起动：强制角 θ_if 斜坡 + 阻尼补偿，电流环在 θ_if 系闭环
    SL_STAGE_BLEND = 1,   // 切换过渡：θ_ctrl 由 θ_if 经最短角差插值到 θ̂，物理电流矢量保持不动
    SL_STAGE_CLOSED = 2,  // 闭环运行：速度环用 ω̂_m 反馈，电流环用 Park(θ̂) 估计 dq + ω̂_e 前馈
    SL_STAGE_GEN = 3,     // 发电稳压：速度环旁路，电压环为唯一 iq 源（iq<0 回馈母线），
                          // 转速解除钳制，由涡轮驱动转矩 = 发电制动转矩决定平衡点；
                          // STSMO/PLL/电流环/θ̂ 注入 SVPWM 原样续用，id_ref 保持 0（无弱磁）
} SensorlessStage;

// 单个控制拍的输出：交给 main.c 注入 SystemInput（ud/uq/theta_ctrl）与日志。
typedef struct SensorlessOutput {
    double ud_cmd;         // 电流环 d 轴电压指令 V（注入 SystemInput.ud_ref）
    double uq_cmd;         // 电流环 q 轴电压指令 V（注入 SystemInput.uq_ref）
    double theta_ctrl_e;   // 注入 SVPWM 的冻结控制电角度 rad（use_ctrl_theta=1）
    double id_ref_cmd;     // d 轴电流指令 A（日志）
    double iq_ref_cmd;     // q 轴电流指令 A（日志）
    double mod_scale;      // 电流环电压矢量缩放（1=未饱和，日志）
    double omega_est_m;    // 估计机械角速度 rad/s（速度环反馈 + 日志）
    double theta_est_e;    // PLL 估计电角度 θ̂ rad（诊断列 theta_err 的来源）
    int stage;             // 当前阶段 0/1/2/3
} SensorlessOutput;

// 起动状态机参数与运行时状态（纯函数 + 全局单例，风格对齐现有控制器）。
// 所有对外接口只允许使用：①实测三相电流 ia/ib/ic；②自身电压指令与估计角；③母线电压 v_dc。
// 严禁读取任何真值（theta_mech / omega_mech / motor_view.id/iq/theta_e/omega_e）。
typedef struct SensorlessMgr {
    double Ts;                 // 控制周期 s（= control_dt）
    // ---- 切换判据（阶段 0 -> 1）----
    double switch_omega_tol;   // |ω̂_e − ω_if_e| 容差 rad/s（电），确认估计转速已跟上强制转速
    double switch_angle_tol;   // PLL 锁定容差 rad：|ε| < sin(angle_tol) 视为已锁定（等价 |θ̂−θ_true| < angle_tol）
    int switch_lock_ticks;     // 上述判据需连续满足的控制拍数 N（滤除超螺旋纹波抖动）
    double lock_lpf_hz;        // 判据低通截止 Hz：对 |ε|、|ω̂_e−ω_if_e| 平滑后再比较，
                               // 使判据反映【平均锁定质量】而非被低速反电动势下的纹波瞬时值反复复位
    // ---- 切换过渡（阶段 1）----
    double blend_time;         // θ_if -> θ̂ 插值时长 s（建议 20~50 ms）
    // ---- 闭环运行（阶段 2）----
    double accel_time;         // 转速参考从切换转速斜坡升到目标的时间 s（限制 iq，避免阶跃冲击）
    // ---- 发电切换判据（阶段 2 -> 3，只用允许信号：参考斜坡到顶 + 低通后的速度环 iq_ref）----
    double gen_iq_th;          // 低通 iq_ref 切换阈值 A：≤ 此值判定电机扭矩需求消失（切入发电）。
                               // 默认 -1e9（不可能满足）= 关闭阶段 3，保证无 gen 键工况行为不变
    int gen_hold_ticks;        // 上述判据需连续满足的控制拍数（防供气爬坡期速度环输出纹波误触发）
    double gen_iq_lpf_hz;      // iq_ref 判据低通截止 Hz（平滑速度环输出纹波后再比较）
    // ---- 运行时状态 ----
    int stage;                 // 当前阶段
    int lock_count;            // 判据连续满足计数
    double blend_elapsed;      // 阶段 1 已历时 s
    double eps_abs_filt;       // |ε| 低通（判据用）
    double omega_err_filt;     // |ω̂_e − ω_if_e| 低通（判据用）
    double prev_u_alpha;       // 上一拍施加电压 α 分量 V（STSMO 的 ZOH 输入）
    double prev_u_beta;        // 上一拍施加电压 β 分量 V
    double ref_omega_m;        // 阶段 2 内部转速参考 rad/s（斜坡）
    double ref_ramp_rate;      // 阶段 2 参考斜坡速率 rad/s²（进入阶段 2 时锁定）
    double gen_iq_filt;        // 低通后的速度环 iq_ref A（发电切换判据用）
    int gen_count;             // 发电切换判据连续满足计数
} SensorlessMgr;

extern SensorlessMgr sensorless_mgr;

// 复位状态机（起动前调用）：清阶段/计数/blend 计时/ZOH 电压/参考斜坡。
void sensorless_reset(SensorlessMgr* mgr);

// 单拍推进。ia/ib/ic 为实测三相电流（视为电流传感器测量），v_dc 为母线电压测量，
// linear_limit 为 SVPWM 线性区上限，omega_ref_mech 为阶段 2 目标机械角速度 rad/s，
// v_ref 为阶段 3 电压环母线电压参考 V（与 sim.v_ref 键统一，仅阶段 3 生效）。
SensorlessOutput sensorless_evaluate(SensorlessMgr* mgr,
                                     const MotorParams* motor_params,
                                     double ia,
                                     double ib,
                                     double ic,
                                     double v_dc,
                                     double linear_limit,
                                     double omega_ref_mech,
                                     double v_ref);

#endif // SENSORLESS_H
