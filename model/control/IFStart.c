#include "IFStart.h"

#include <math.h>
#include "utils/numeric_utils.h"

// ============================ 默认参数与阻尼补偿物理依据 ============================
// 电机/轴系参数（同 Motor.c / Shaft.c）：Rs=0.0051Ω, Ld=Lq=34.6µH, ψ_f=0.0595Wb, p=2,
// J_total=7.74e-4 kg·m², b_total=4.8e-3 N·m·s。
//
// 【I/F 功角动态】id=0、iq=I_f 在【强制角 θ_if 坐标系】下给定，定子电流矢量落在强制 q 轴
// (θ_if+90°)。转子磁链 ψ_f 在 θ_true。电磁转矩 Tem=(3/2)p·ψ_f·iq_rotor=K·cosδ，
// 其中 K=(3/2)p·ψ_f·I_f，δ=θ_if−θ_true 为功角。机械方程 J·ω̇_m=K·cosδ−b·ω_m−T_load，
// 且 δ̇=p·(ω_if_m−ω_m)。空载稳态转子对齐定子磁场：λ=δ+90°≈0 → δ0≈−90°（λ 为转子滞后磁场角，
// 即经典"转子-定子磁场"功角；升速带载时 cosδ0=T_load/K>0 使 δ0 略向 0 侧移，实测约 −89°→−82°）。
// 小信号线性化：  Δδ̈ + (b/J)·Δδ̇ − (p·K·sinδ0/J)·Δδ = 0
// δ0≈−90° 时 sinδ0≈−1 → ω_n=√(pK/J)≈166 rad/s(26.5 Hz)，ζ=b/(2√(J·p·K))≈0.019（极欠阻尼，
// 升速时功角长时间振荡，是 I/F 失步主因）。
//
// 【为何用角度修正而非速率修正】
//  · 速率修正 ω_if_corr=ω_if−k(ω̂−ω_if)：δ̇=p(1+k)(ω_if−ω_m) → ω_n²→(1+k)ω_n²、阻尼项(b/J)不变，
//    ζ'=ζ/√(1+k) 反而【下降】（只增刚度，越补越振）。
//  · 角度修正 θ_if=θ_if_ramp−damp_k·(ω̂_e−ω_if_e)：ω_n 不变、阻尼项→(b+damp_k·p·K·|sinδ0|)，
//    ζ → ζ + (damp_k/2)·ω_n。取 damp_k=3.5e-3 s → ζ≈0.31（良阻尼，抑制功角振荡而不改自然频率）。
// 阻尼只用【可测量】的估计转速 ω̂（来自 STSMO/PLL 的反电动势估计），严禁使用真值或功角 δ。
//
// 【参数取值】
//  · i_start=60 A → K≈10.7 N·m，为空载升速需求（2000 rpm 时惯量+摩擦约 1.5 N·m）的 ~7×（带载校核留裕度），
//    且 < 2×稳态 iq(56.3 A)=112.6 A，满足切换过程电流冲击判据。I_f 越小则 cosδ0=T_load/K 越大、
//    δ0 越靠 0 侧（功角离 −90° 边界越远、字面 |δ|<90° 裕度越大），但拖载转矩 K 越小——二者权衡取 60 A。
//  · omega_switch_m=2000 rpm 机械角=209.44 rad/s → 反电动势 E=ω_e·ψ_f≈24.9 V，STSMO/PLL 信噪比可接受。
//  · ramp_time=0.3 s：0→2000 rpm 斜坡。三者均可经 sensorless.if.* 覆盖。
//
// 【强制角初值 theta_init_e —— 从平衡点无冲击起动】
//  若 θ_if(0)=θ_true(0)，则 δ(0)=0、Tem=K·cos0=K≈10.7 N·m（最大转矩），静止转子被阶跃转矩猛推，
//  激起欠阻尼功角模态大幅振荡甚至 δ 绕圈失步（叠加涡轮容腔吹除瞬态时更会飞车至数千 rpm；吹除瞬态
//  另由 config 将容腔初压设为背压消除，见 README）。稳定平衡功角为 δ0≈−90°（罗盘针对齐、Tem≈0）。
//  故假定起动前经【预定位】（真实驱动器用 DC 注入/高频注入建立已知初始位置）使 θ_true(0)≈0，
//  据此令 θ_if(0)=−π/2 → δ(0)=−90°、Tem(0)=0：定子电流矢量初始恰落在转子 d 轴（id_true≈I_f、
//  iq_true≈0），转子不受冲击；随后斜坡缓慢旋转 θ_if，δ 从 −90° 准静态爬到 −82°，转子平稳单调跟随。
//  【实测】I/F 全程 δ∈[−90°,−80°]（仅 t=0 初始平衡点取 −90°），经典功角 δ'=δ+90°∈[0°,9.6°]
//  远离失步边界 90°；转速单调跟随斜坡（最大回落 1.1 rpm≈0.05%）、无后仰。
//  注：θ_if(0) 若偏离 −π/2（如 −88°）会引入初始转矩阶跃，反而激起更大中期振荡（实测 δ 跌到 −90.3°），
//  故取精确平衡点 −π/2。theta_init_e 是【设计常量】而非运行时读取真值，符合无感约束。
IFStart ifstart = {
    .Ts = 5e-6,
    .p = 2,
    .i_start = 60.0,
    .omega_switch_m = 209.43951023931953,   // 2000 rpm 机械角
    .ramp_time = 0.3,
    .damp_k = 3.5e-3,
    .damp_lpf_hz = 70.0,                   // 阻尼用 ω̂ 低通：>ω_n(26.5Hz) 保留真实转速偏差，≪纹波(266Hz) 抑制注入
    .theta_init_e = -1.5707963267948966,    // −π/2：δ(0)=−90° 稳定平衡（预定位假设 θ_true(0)≈0）
    .omega_if_m = 0.0,
    .theta_ramp_e = -1.5707963267948966,
    .omega_est_filt_e = 0.0,
};

// 阻尼补偿内部约束（非对外参数）：
//  - 低速反电动势趋零、PLL 估计不可靠：零速时归一化把超螺旋抖噪放大为单位矢量，PLL 追逐噪声使
//    ω̂ 积分漂移（实测 t≈7 ms 时 ω̂_e≈716 rad/s，而转子实际≈−18 rpm）。该 garbage ω̂ 经 corr
//    反向拖拽 θ_if，造成起动初期转子"后仰"（ω_m 短暂负转、δ 跌破 −90°）。故设【淡入死区】：
//    ω_if 低于 IF_DAMP_GATE_ON_FRAC·ω_switch（反电动势尚不可测）时【不施加阻尼】，其后线性淡入，
//    达 IF_DAMP_GATE_FRAC·ω_switch 时全量——此时转子已实转、ω̂ 可信，阻尼才有效；
//  - 修正角限幅 IF_DAMP_CORR_MAX，防止未收敛时 ω̂ 残余误差污染强制角。
#define IF_DAMP_GATE_ON_FRAC 0.15          // 阻尼起始比例（300 rpm，反电动势 E≈3.7 V 起可测）
#define IF_DAMP_GATE_FRAC 0.5              // 阻尼全量比例（相对 omega_switch_m，1000 rpm 起全量）
#define IF_DAMP_CORR_MAX 0.349065850398864 // 修正角限幅 rad（≈20°）

void ifstart_reset(IFStart* self) {
    self->omega_if_m = 0.0;
    self->theta_ramp_e = wrap_2pi(self->theta_init_e);   // 从 δ(0)=−90° 平衡起步，消除飞车
    self->omega_est_filt_e = 0.0;
}

IFStartOutput ifstart_evaluate(IFStart* self, double omega_est_e) {
    // 强制机械角速度斜坡 0 → omega_switch_m（前向欧拉，达到后保持）
    const double ramp_rate = self->omega_switch_m / fmax(self->ramp_time, 1e-9);
    self->omega_if_m = clamp_range(self->omega_if_m + ramp_rate * self->Ts,
                                   0.0, self->omega_switch_m);
    const double omega_if_e = (double)self->p * self->omega_if_m;

    // 斜坡角积分（电角度），回绕 [0,2π)
    self->theta_ramp_e = wrap_2pi(self->theta_ramp_e + omega_if_e * self->Ts);

    // 阻尼补偿（角度修正）：θ_if = θ_ramp − damp_k·gate·(ω̂_e_filt − ω_if_e)
    // gate 带淡入死区：ω_if < ON_FRAC·ω_switch 时为 0（零速 ω̂ 不可信，避免后仰），其后线性升至全量。
    const double w_on = IF_DAMP_GATE_ON_FRAC * self->omega_switch_m;
    const double w_full = IF_DAMP_GATE_FRAC * self->omega_switch_m;
    const double gate = clamp_range((self->omega_if_m - w_on) / fmax(w_full - w_on, 1e-9),
                                    0.0, 1.0);
    // ω̂ 低通后再进 corr：低速反电动势小、PLL 估计含 ~266 Hz 纹波，若直接进 corr 会按 damp_k·gate
    // 把纹波注入 θ_if（实测 δ/转速出现 ~330 Hz 抖动、δ 跌破 −90°）。一阶低通截止 damp_lpf_hz
    // （>机械模态 ω_n≈34 Hz 保留真实转速偏差供阻尼，≪纹波 266 Hz 抑制其注入）。
    if (self->damp_lpf_hz > 0.0) {
        const double a_d = 1.0 - exp(-TWO_PI_RAD * self->damp_lpf_hz * self->Ts);
        self->omega_est_filt_e += a_d * (omega_est_e - self->omega_est_filt_e);
    } else {
        self->omega_est_filt_e = omega_est_e;
    }
    const double corr = clamp_range(self->damp_k * gate * (self->omega_est_filt_e - omega_if_e),
                                    -IF_DAMP_CORR_MAX, IF_DAMP_CORR_MAX);
    const double theta_if_e = wrap_2pi(self->theta_ramp_e - corr);

    IFStartOutput out;
    out.theta_if_e = theta_if_e;
    out.omega_if_e = omega_if_e;
    out.omega_if_m = self->omega_if_m;
    out.id_ref = 0.0;
    out.iq_ref = self->i_start;
    out.at_switch = (self->omega_if_m >= self->omega_switch_m - 1e-9);
    return out;
}
