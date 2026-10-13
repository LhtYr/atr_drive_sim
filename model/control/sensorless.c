#include "sensorless.h"

#include <math.h>

#include "utils/numeric_utils.h"
#include "utils/transforms.h"
#include "Current.h"
#include "Speed.h"
#include "Voltage.h"
#include "STSMO.h"
#include "PLL.h"
#include "IFStart.h"

// ================================ 切换过渡物理依据 ================================
// 【空载 I/F 稳态的功角】I/F 阶段电流矢量落在强制 q 轴（αβ 角 θ_if+90°），转子磁链在 θ_true。
// 电磁转矩 Tem=K·cosδ，K=(3/2)p·ψ_f·I_f，δ=θ_if−θ_true。永磁转子在定子磁场中趋于【对齐】
// （罗盘针效应）：空载稳态 θ_true → θ_if+90°，即 δ0≈−90°、转子【超前】强制角约 90°。
// 升速斜坡末端需少量转矩克服摩擦+惯量，δ0≈−85°（cosδ0≈0.09）。故切换时刻：
//   · 估计角 θ̂ ≈ θ_true ≈ θ_if+85°（|θ̂−θ_if|≈85°，并非小量！）；
//   · 物理电流矢量在 θ_if+90° ≈ θ̂+5°，即落在转子【d 轴附近】：id_rotor≈I_f·cos5°≈59.8 A，
//     iq_rotor≈I_f·sin5°≈5.2 A（转矩≈1.5 N·m，仅够克服摩擦）。
//
// 【为何不能照搬"iq_ref 恒为 I_f 交接给速度环"】若在过渡中令 θ_ctrl→θ̂ 而电流指令仍取
// id=0/iq=I_f（θ_ctrl 系），物理电流矢量将被拽到 θ̂+90°（转子 q 轴），δ→0、Tem→K≈10.7 N·m，
// 十几倍于摩擦转矩 → 转速冲击、失步。故本实现的过渡遵循两条不变量：
//   (1) 交换的是【换相角来源】θ_if→θ̂，而非把电流硬转到 q 轴；
//   (2) 过渡全程【物理电流矢量冻结】在 θ_if+90°（转矩≈摩擦，转速平稳）。
// 实现：θ_ctrl(s)=θ_if+s·Δ，Δ=wrap_pi(θ̂−θ_if)≈85°，s:0→1（T_blend 内）；同时在 θ_ctrl 系给
//   id_ref=I_f·sin(sΔ)、iq_ref=I_f·cos(sΔ)，
// 则物理电流矢量角 = θ_ctrl + atan2(iq_ref,id_ref) = θ_if+sΔ+(90°−sΔ) = θ_if+90°（恒定，与 s、Δ 符号无关）。
// s=1 时 θ_ctrl=θ̂，测量得 id≈59.8/iq≈5.2；随即进入阶段 2：id_ref→0（d 轴电流快速泄放，无转矩），
// 速度环积分器 int_omega 预置为该 iq（≈I_f·cosΔ），iq_ref 无跳变交接 → bumpless。
//
// 【阶段 2 转速参考斜坡】若沿用阶跃参考（0→目标），速度环比例项 kp·e_ω 会在交接瞬间把 iq_ref
// 顶到数百安（远超 2×稳态），故进入阶段 2 时以【当前 ω̂_m】为起点、accel_time 内斜坡升到目标，
// 使 e_ω 始终很小、iq 平稳爬升（默认 accel_time=1 s 时 iq≤~65 A < 2×稳态 56.3 A）。
//
// 【切换判据】空载下 |θ̂−θ_if|≈85° 恒不满足"小于容差"，故改用 PLL 锁定判据（等价于 |θ̂−θ_true| 小）：
//   ω_if 已到斜坡顶（at_switch）且 |ω̂_e−ω_if_e|<switch_omega_tol 且 |ε|<sin(switch_angle_tol)
//   连续 switch_lock_ticks 拍。ε=sin(θ_true−θ̂) 为 PLL 相位误差，锁定后 →0，直接对应 |theta_err| 小。
//
// ================================ 阶段 3：发电稳压切换 ================================
// 【任务剖面】阶段 2 稳速 20000 rpm 后，供气延迟斜坡启动、涡轮转矩渐增，速度环为守转速把 iq_ref
// 渐降（甚至转负=回馈制动，此乃正常物理）；当电机扭矩需求消失即切入阶段 3（发电）：速度环旁路、
// 电压环成为唯一 iq 源、转速解除钳制，由"涡轮驱动转矩 = 发电制动转矩"决定平衡转速。
//
// 【切换判据（只用允许信号）】① 参考斜坡已到顶 ref_omega_m == omega_ref_mech（已在稳速段）；
//   ② 低通后的速度环 iq_ref ≤ gen_iq_th（默认 +5 A，即扭矩需求消失）连续 gen_hold_ticks 拍。
//   用低通值（gen_iq_lpf_hz，τ≈3.2 ms）而非瞬时值：供气爬坡期速度环输出含纹波且可能已瞬时转负，
//   低通+连续计数防止纹波误触发。gen_iq_th 默认 -1e9 使判据【不可能满足】= 关闭阶段 3（回归不变）。
//
// 【电压环符号（继承，严禁改回）】发电（iq<0）把机械功率注入母线、抬升 v_dc，故稳定方向为
//   iq_ref ∝ (v_dc − v_ref)：v_dc 偏低→加大发电（iq 更负）抬升，v_dc 偏高→减小发电回落。见 Voltage.c。
//
// 【bumpless 交接】切换瞬间电压环积分器预置 voltage.int_v = iq_ref_now − voltage.kp·(v_dc − v_ref)。
//   电压环输出 iq = kp·(v_dc − v_ref) + int_v，代入预置值得 iq = iq_ref_now，即切入首拍 iq 指令
//   = 切换前速度环 iq_ref，连续无跳变（|Δiq| 仅来自相邻拍 v_dc 微变 × kp，远小于 10 A）。
//   本拍仍输出速度环 iq_ref、阶段 3 自下一拍接管，与阶段 1->2 速度环积分器预置交接完全同构。
//
// 【v_ref 可达性（已知坑）】发电只能把母线抬到电池 OCV(SOC) 之上；v_ref 须 ≥ 起点 OCV 才可达
//   （历史教训：v_ref=538 < OCV 不可达，积分器单向饱卷）。SOC≈0.8 起点 OCV≈540，v_ref=540 可达。
SensorlessMgr sensorless_mgr = {
    .Ts = 5e-6,
    .switch_omega_tol = 50.0,              // rad/s（电），约 ω_if_e(2000rpm)=419 的 12%
    .switch_angle_tol = 0.174532925199433, // 10° -> |ε|<sin(10°) 即判 PLL 锁定（|theta_err|<10°）
    .switch_lock_ticks = 200,              // 连续 1 ms
    .lock_lpf_hz = 100.0,                  // 判据低通：平滑 |ε|/|Δω̂|，抗低速反电动势下的纹波
    .blend_time = 0.03,                    // 30 ms 过渡
    .accel_time = 1.0,                     // 阶段 2 转速参考斜坡 2000->20000 rpm 历时 1 s
    // 发电切换判据（阶段 2 -> 3）：默认 gen_iq_th=-1e9 使判据【不可能满足】= 关闭阶段 3，
    // 保证 sensorless_start.txt 等无 gen 键工况逐拍不变；任务剖面经 sensorless.gen.* 键启用。
    .gen_iq_th = -1e9,                     // 低通 iq_ref ≤ 此值判电机扭矩需求消失（默认不可能满足）
    .gen_hold_ticks = 400,                 // 判据连续满足拍数（400×5µs=2 ms 去抖）
    .gen_iq_lpf_hz = 50.0,                 // iq_ref 判据低通截止 Hz（τ≈3.2 ms，平滑速度环输出纹波）
    .stage = SL_STAGE_IF,
    .lock_count = 0,
    .blend_elapsed = 0.0,
    .eps_abs_filt = 0.0,
    .omega_err_filt = 0.0,
    .prev_u_alpha = 0.0,
    .prev_u_beta = 0.0,
    .ref_omega_m = 0.0,
    .ref_ramp_rate = 0.0,
    .gen_iq_filt = 0.0,
    .gen_count = 0,
};

void sensorless_reset(SensorlessMgr* mgr) {
    mgr->stage = SL_STAGE_IF;
    mgr->lock_count = 0;
    mgr->blend_elapsed = 0.0;
    mgr->eps_abs_filt = 0.0;
    mgr->omega_err_filt = 0.0;
    mgr->prev_u_alpha = 0.0;
    mgr->prev_u_beta = 0.0;
    mgr->ref_omega_m = 0.0;
    mgr->ref_ramp_rate = 0.0;
    mgr->gen_iq_filt = 0.0;
    mgr->gen_count = 0;
}

SensorlessOutput sensorless_evaluate(SensorlessMgr* mgr,
                                     const MotorParams* motor_params,
                                     double ia,
                                     double ib,
                                     double ic,
                                     double v_dc,
                                     double linear_limit,
                                     double omega_ref_mech,
                                     double v_ref) {
    SensorlessOutput out;

    // ① 实测三相电流 -> αβ（Clarke，幅值不变）。这是控制侧唯一允许的电流信息。
    double i_alpha = 0.0;
    double i_beta = 0.0;
    clarke_abc_to_ab(ia, ib, ic, &i_alpha, &i_beta);

    // ② 观测：STSMO 用【上一拍施加电压】(ZOH) + 当前实测电流 -> 反电动势 -> PLL -> θ̂/ω̂。
    const STSMOOutput smo = stsmo_evaluate(&stsmo, motor_params,
                                           mgr->prev_u_alpha, mgr->prev_u_beta,
                                           i_alpha, i_beta);
    const PLLOutput pl = pll_evaluate(&pll, smo.emf_alpha, smo.emf_beta);
    const double theta_est = pl.theta_e;
    const double omega_est_e = pl.omega_e;
    const double omega_est_m = pl.omega_m;

    double theta_ctrl = 0.0;
    double id_ref = 0.0;
    double iq_ref = 0.0;
    double omega_ff_e = omega_est_e;   // 电流环前馈电角速度（各阶段按物理意义选取）

    if (mgr->stage == SL_STAGE_IF) {
        // 阶段 0：I/F 开环起动。强制角 θ_if（含阻尼补偿），电流环在 θ_if 系闭环 id=0/iq=I_f。
        const IFStartOutput ifo = ifstart_evaluate(&ifstart, omega_est_e);
        theta_ctrl = ifo.theta_if_e;
        id_ref = ifo.id_ref;      // 0
        iq_ref = ifo.iq_ref;      // I_f
        omega_ff_e = ifo.omega_if_e;

        // 切换判据：强制斜坡到顶 + PLL 锁定（转速跟上、相位误差小），连续 N 拍。
        // 低速（2000 rpm，E≈24.9 V）反电动势小、超螺旋纹波使 |ε|/|ω̂_e−ω_if_e| 瞬时抖动大，
        // 直接用瞬时值会反复复位计数；故先按 lock_lpf_hz 低通再比较，判据反映平均锁定质量。
        const double a_lock = 1.0 - exp(-TWO_PI_RAD * mgr->lock_lpf_hz * mgr->Ts);
        mgr->eps_abs_filt += a_lock * (fabs(pl.eps) - mgr->eps_abs_filt);
        mgr->omega_err_filt += a_lock * (fabs(omega_est_e - ifo.omega_if_e) - mgr->omega_err_filt);
        const int locked = ifo.at_switch &&
                           mgr->omega_err_filt < mgr->switch_omega_tol &&
                           mgr->eps_abs_filt < sin(mgr->switch_angle_tol);
        mgr->lock_count = locked ? (mgr->lock_count + 1) : 0;
        if (mgr->lock_count >= mgr->switch_lock_ticks) {
            mgr->stage = SL_STAGE_BLEND;
            mgr->blend_elapsed = 0.0;
            mgr->lock_count = 0;
        }
    } else if (mgr->stage == SL_STAGE_BLEND) {
        // 阶段 1：θ_if -> θ̂ 最短角差插值，物理电流矢量冻结在 θ_if+90°（转矩≈摩擦，无冲击）。
        const IFStartOutput ifo = ifstart_evaluate(&ifstart, omega_est_e);
        const double theta_if = ifo.theta_if_e;
        const double delta = wrap_pi(theta_est - theta_if);   // ≈ +85°（空载）
        mgr->blend_elapsed += mgr->Ts;
        double s = mgr->blend_elapsed / fmax(mgr->blend_time, 1e-9);
        if (s > 1.0) {
            s = 1.0;
        }
        theta_ctrl = wrap_2pi(theta_if + s * delta);
        id_ref = ifstart.i_start * sin(s * delta);
        iq_ref = ifstart.i_start * cos(s * delta);
        omega_ff_e = omega_est_e;

        if (s >= 1.0) {
            // 进入阶段 2：bumpless —— 速度环积分器预置为当前 iq_ref（≈I_f·cosΔ），
            // 转速参考以【当前 ω̂_m】为起点斜坡升到目标，避免阶跃导致 iq 冲击。
            mgr->stage = SL_STAGE_CLOSED;
            speed.int_omega = iq_ref;
            mgr->ref_omega_m = omega_est_m;
            mgr->ref_ramp_rate = (omega_ref_mech - omega_est_m) / fmax(mgr->accel_time, 1e-9);
        }
    } else if (mgr->stage == SL_STAGE_CLOSED) {
        // 阶段 2：闭环运行。转速参考斜坡 -> 速度环（ω̂_m 反馈）-> iq_ref；电流环在 θ̂ 系闭环。
        mgr->ref_omega_m += mgr->ref_ramp_rate * mgr->Ts;
        if (mgr->ref_ramp_rate >= 0.0) {
            mgr->ref_omega_m = fmin(mgr->ref_omega_m, omega_ref_mech);
        } else {
            mgr->ref_omega_m = fmax(mgr->ref_omega_m, omega_ref_mech);
        }
        const SpeedOutput spo = speed_evaluate(&speed, mgr->ref_omega_m, omega_est_m);
        theta_ctrl = theta_est;
        id_ref = 0.0;
        iq_ref = spo.iq_ref;
        omega_ff_e = omega_est_e;

        // 阶段 2 -> 3（发电）切换判据：只用允许信号——参考斜坡已到顶（ref_omega_m == omega_ref_mech）
        // 且【低通后】的速度环 iq_ref ≤ gen_iq_th 连续 gen_hold_ticks 拍。供气渐增使涡轮转矩接管，
        // 速度环为守转速把 iq_ref 拉低（甚至转负=回馈制动）；低通 iq_ref 降到阈值即判电机扭矩需求消失。
        // 用低通值（gen_iq_lpf_hz）而非瞬时值，防供气爬坡期速度环输出纹波误触发（正常物理：爬坡期
        // 速度环可能已输出负 iq 守转速）。gen_iq_th 默认 -1e9 使判据不可能满足（关闭阶段 3），
        // 既有 sensorless_start.txt 等无 gen 键工况逐拍不变。
        const double a_gen = 1.0 - exp(-TWO_PI_RAD * mgr->gen_iq_lpf_hz * mgr->Ts);
        mgr->gen_iq_filt += a_gen * (iq_ref - mgr->gen_iq_filt);
        const int ramp_done = (mgr->ref_omega_m == omega_ref_mech);
        const int gen_req = ramp_done && (mgr->gen_iq_filt <= mgr->gen_iq_th);
        mgr->gen_count = gen_req ? (mgr->gen_count + 1) : 0;
        if (mgr->gen_hold_ticks > 0 && mgr->gen_count >= mgr->gen_hold_ticks) {
            // bumpless 交接：电压环积分器预置 int_v = iq_ref_now − kp·(v_dc − v_ref)，使切入首拍
            // 电压环输出 iq = kp·(v_dc−v_ref) + int_v = iq_ref_now，iq 指令连续无跳变。本拍仍输出
            // 速度环 iq_ref，阶段 3 自下一拍起接管（与阶段 1->2 速度环积分器预置交接同构）。
            mgr->stage = SL_STAGE_GEN;
            voltage.int_v = iq_ref - voltage.kp * (v_dc - v_ref);
            mgr->gen_count = 0;
        }
    } else {
        // 阶段 3：发电稳压。速度环旁路，电压环为唯一 iq 源：iq_ref = voltage_evaluate(v_ref, v_dc).iq_ref
        // （稳定方向 iq_ref ∝ v_dc − v_ref，稳态为负=发电回馈母线抬升 v_dc）。转速解除钳制，由
        // "涡轮驱动转矩 = 发电制动转矩"决定平衡点。STSMO/PLL/电流环/θ̂ 注入 SVPWM 原样续用，
        // id_ref 保持 0（无弱磁）。ref_omega_m 冻结不再更新（速度环已旁路）。
        const VoltageOutput vo = voltage_evaluate(&voltage, v_ref, v_dc);
        theta_ctrl = theta_est;
        id_ref = 0.0;
        iq_ref = vo.iq_ref;
        omega_ff_e = omega_est_e;
    }

    // ③ 电流环：在 θ_ctrl 系闭合。测量电流经 Park(θ_ctrl)（估计/强制系，非真值系）。
    const double cos_t = cos(theta_ctrl);
    const double sin_t = sin(theta_ctrl);
    double id_meas = 0.0;
    double iq_meas = 0.0;
    park_ab_to_dq(i_alpha, i_beta, cos_t, sin_t, &id_meas, &iq_meas);
    const CurrentOutput co = current_evaluate(&current, motor_params,
                                              id_ref, iq_ref,
                                              id_meas, iq_meas,
                                              omega_ff_e, v_dc, linear_limit);
    out.ud_cmd = co.ud_ref;
    out.uq_cmd = co.uq_ref;
    out.mod_scale = co.modulation_scale;

    // ④ 由本拍电压指令 + θ_ctrl 合成 u_αβ，作为下一拍 STSMO 的 ZOH 输入（数字控制一拍一更新语义）。
    double u_alpha = 0.0;
    double u_beta = 0.0;
    double ua_unused = 0.0;
    double ub_unused = 0.0;
    double uc_unused = 0.0;
    dq_to_abc(out.ud_cmd, out.uq_cmd, cos_t, sin_t,
              &u_alpha, &u_beta, &ua_unused, &ub_unused, &uc_unused);
    mgr->prev_u_alpha = u_alpha;
    mgr->prev_u_beta = u_beta;

    out.theta_ctrl_e = theta_ctrl;
    out.id_ref_cmd = id_ref;
    out.iq_ref_cmd = iq_ref;
    out.omega_est_m = omega_est_m;
    out.theta_est_e = theta_est;
    out.stage = mgr->stage;
    return out;
}
