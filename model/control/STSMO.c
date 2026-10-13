#include "STSMO.h"

#include <math.h>
#include "utils/numeric_utils.h"

// 默认整定依据（Ts = control_dt = 5e-6 s，Ld=Lq=34.6 µH，Rs=0.0051 Ω，ψ_f=0.0595 Wb，p=2）：
//  - k2 须大于反电动势最大变化率 |de/dt|_max = E_max·ω_e_max；20000 rpm 时
//    E=ω_e·ψ_f≈249 V、ω_e≈4189 rad/s → 约 1.0e6 V/s，取 k2=3e6 留裕度（超螺旋积分项
//    v 才能跟上旋转反电动势的斜率）；
//  - k1 决定离散滑模到达速度与电流误差抖振带：稳态抖振带约 (k1·Ts/L)²，
//    k1=5 时约 0.5 A（远小于工作电流 ~100 A），既能快速到达又不致过冲发散；
//  - emf_lpf_hz 取 2 kHz：仍高于 20000 rpm 反电动势基波 666 Hz 约 3×（相位滞后由 PLL 补偿），
//    将 ~15 kHz 量级的超螺旋切换纹波衰减约 18 dB；截止越低，低速段（2000 rpm、E≈25 V）
//    反电动势幅值纹波相对越小、PLL 信噪比越好，但须留足 PLL 相位裕度（环路带宽 ~356 Hz）。
// 三者均可经参数文件 sensorless.stsmo.* 覆盖后重调。
STSMO stsmo = {
    .Ts = 5e-6,
    .k1 = 5.0,
    .k2 = 3.0e6,
    .emf_lpf_hz = 2000.0,
    .i_hat_alpha = 0.0,
    .i_hat_beta = 0.0,
    .v_alpha = 0.0,
    .v_beta = 0.0,
    .emf_alpha = 0.0,
    .emf_beta = 0.0,
};

// 防御性限幅：正常工况（电流 <540 A、反电动势 <249 V）远不触及，仅兜底参数失配导致的发散，
// 保证全程无 nan/inf。
#define STSMO_IHAT_CLAMP 1.0e4   // 观测电流限幅 A
#define STSMO_V_CLAMP 5.0e3      // 超螺旋积分项限幅 V

void stsmo_reset(STSMO* obs) {
    obs->i_hat_alpha = 0.0;
    obs->i_hat_beta = 0.0;
    obs->v_alpha = 0.0;
    obs->v_beta = 0.0;
    obs->emf_alpha = 0.0;
    obs->emf_beta = 0.0;
}

STSMOOutput stsmo_evaluate(STSMO* obs,
                           const MotorParams* motor_params,
                           double u_alpha,
                           double u_beta,
                           double i_alpha,
                           double i_beta) {
    const double Rs = motor_params->Rs;
    const double L = fmax(motor_params->Ld, 1e-12);   // 表贴式假设 Ld≈Lq，αβ 系电感取 Ld

    // 电流估计误差 e = î − i
    const double e_alpha = obs->i_hat_alpha - i_alpha;
    const double e_beta = obs->i_hat_beta - i_beta;

    // 超螺旋切换项 z = k1·|e|^0.5·sign(e) + v（收敛后 z ≈ 反电动势）
    const double z_alpha = obs->k1 * sqrt(fabs(e_alpha)) * (e_alpha > 0.0 ? 1.0 : (e_alpha < 0.0 ? -1.0 : 0.0)) + obs->v_alpha;
    const double z_beta = obs->k1 * sqrt(fabs(e_beta)) * (e_beta > 0.0 ? 1.0 : (e_beta < 0.0 ? -1.0 : 0.0)) + obs->v_beta;

    // 电流观测（前向欧拉）：dî/dt = (u − Rs·î − z)/L
    obs->i_hat_alpha += obs->Ts * (u_alpha - Rs * obs->i_hat_alpha - z_alpha) / L;
    obs->i_hat_beta += obs->Ts * (u_beta - Rs * obs->i_hat_beta - z_beta) / L;
    obs->i_hat_alpha = clamp_range(obs->i_hat_alpha, -STSMO_IHAT_CLAMP, STSMO_IHAT_CLAMP);
    obs->i_hat_beta = clamp_range(obs->i_hat_beta, -STSMO_IHAT_CLAMP, STSMO_IHAT_CLAMP);

    // 超螺旋积分项：dv/dt = k2·sign(e)
    obs->v_alpha += obs->Ts * obs->k2 * (e_alpha > 0.0 ? 1.0 : (e_alpha < 0.0 ? -1.0 : 0.0));
    obs->v_beta += obs->Ts * obs->k2 * (e_beta > 0.0 ? 1.0 : (e_beta < 0.0 ? -1.0 : 0.0));
    obs->v_alpha = clamp_range(obs->v_alpha, -STSMO_V_CLAMP, STSMO_V_CLAMP);
    obs->v_beta = clamp_range(obs->v_beta, -STSMO_V_CLAMP, STSMO_V_CLAMP);

    // 反电动势低通（一阶，精确 ZOH 离散）：抑制 z 的切换抖振，相位滞后由 PLL 补偿
    const double wc = TWO_PI_RAD * fmax(obs->emf_lpf_hz, 1e-6);
    const double a_lpf = 1.0 - exp(-wc * obs->Ts);
    obs->emf_alpha += a_lpf * (z_alpha - obs->emf_alpha);
    obs->emf_beta += a_lpf * (z_beta - obs->emf_beta);

    STSMOOutput out;
    out.emf_alpha = obs->emf_alpha;
    out.emf_beta = obs->emf_beta;
    out.emf_mag = sqrt(obs->emf_alpha * obs->emf_alpha + obs->emf_beta * obs->emf_beta);
    out.i_hat_alpha = obs->i_hat_alpha;
    out.i_hat_beta = obs->i_hat_beta;
    return out;
}
