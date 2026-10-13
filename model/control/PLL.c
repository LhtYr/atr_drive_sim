#include "PLL.h"

#include <math.h>
#include "utils/numeric_utils.h"

// 默认整定依据（Ts = control_dt = 5e-6 s）：
// 反电动势经幅值归一化后 ε=sin(θ_e−θ̂)∈[−1,1] 与转速无关，PLL 为二阶环：
// 自然频率 ω_n=√ki、阻尼比 ζ=kp/(2√ki)。取 ki=5e6 → ω_n≈2236 rad/s（带宽 ~356 Hz，
// 远高于起动升速速率），kp=3000 → ζ≈0.67（接近临界阻尼 0.707，快而不过冲）。
// 离散稳定性 Ts·ω_n≈0.011 ≪ 1。omega_e_max=8000 rad/s 覆盖 20000 rpm 电角速度
// ω_e=p·ω_m≈4189 rad/s 并留裕度，未收敛时防止积分发散。
PLL pll = {
    .Ts = 5e-6,
    .kp = 3000.0,
    .ki = 5.0e6,
    .omega_e_max = 8000.0,
    .emf_lpf_wc = 12566.370614359172,   // 2π·2000，与 stsmo.emf_lpf_hz 默认一致（manager 初始化时同步）
    .omega_lpf_hz = 120.0,
    .p = 2,
    .theta_e = 0.0,
    .int_pll = 0.0,
    .omega_filt = 0.0,
};

void pll_reset(PLL* self, double theta0) {
    self->theta_e = wrap_2pi(theta0);
    self->int_pll = 0.0;
    self->omega_filt = 0.0;
}

PLLOutput pll_evaluate(PLL* self, double emf_alpha, double emf_beta) {
    // 反电动势矢量指向 q 轴（θ_e+90°）：(e_α,e_β)=E(−sinθ_e, cosθ_e)。
    // 先旋转 −90° 对齐到转子 d 轴基准 r=(e_β,−e_α)=E(cosθ_e, sinθ_e)，再幅值归一化，
    // 使 ε=sin(θ_e−θ̂) 与 E=ω_e·ψ_f 无关（消除 PLL 增益随转速变化）。
    const double r_alpha = emf_beta;
    const double r_beta = -emf_alpha;
    const double mag = sqrt(r_alpha * r_alpha + r_beta * r_beta);
    const double inv = 1.0 / magnitude_floor(mag, 1e-6);
    const double e_a = r_alpha * inv;   // ≈ cos θ_e
    const double e_b = r_beta * inv;    // ≈ sin θ_e

    const double cos_t = cos(self->theta_e);
    const double sin_t = sin(self->theta_e);
    // ε = ê_β·cosθ̂ − ê_α·sinθ̂ = sin(θ_e − θ̂)：θ̂ 落后时 ε>0 → ω̂ 上升（负反馈锁定）
    const double eps = e_b * cos_t - e_a * sin_t;

    // PI → ω̂_e（积分项与总输出均限幅，防止未收敛时发散）
    self->int_pll += self->ki * eps * self->Ts;
    self->int_pll = clamp_range(self->int_pll, -self->omega_e_max, self->omega_e_max);
    const double omega_e_raw = clamp_range(self->kp * eps + self->int_pll,
                                           -self->omega_e_max, self->omega_e_max);

    // 内部锁定角积分用【原始】ω̂_e：保证快速跟踪（相位滞后由下方补偿）
    self->theta_e = wrap_2pi(self->theta_e + omega_e_raw * self->Ts);

    // 输出转速低通：超螺旋切换纹波 + 低速反电动势小（2000 rpm、E≈24.9 V）使 ε 抖动，
    // 经 kp·ε 放大后 ω̂_e 纹波可达十几个百分点（虽 θ̂ 跟踪真值仅数度，但该角度纹波按 2πf 微分
    // 进 ω̂ 后被显著放大）。对【输出/反馈用】转速做一阶低通（内部角跟踪用原始 ω̂_e，不受影响）：
    // 截止 120 Hz 远低于纹波主频（~266 Hz 及更高）、又远高于速度环带宽（数 Hz），
    // 纹波衰减 >2× 而对闭环滞后可忽略；随转速升高反电动势增大、SNR 改善，纹波进一步下降。
    if (self->omega_lpf_hz > 0.0) {
        const double a_om = 1.0 - exp(-TWO_PI_RAD * self->omega_lpf_hz * self->Ts);
        self->omega_filt += a_om * (omega_e_raw - self->omega_filt);
    } else {
        self->omega_filt = omega_e_raw;
    }

    // 反电动势低通相位滞后补偿：一阶 LPF 在频率 ω̂_e 处滞后 atan(ω̂_e/ω_c)，
    // 内部角锁定到滞后的反电动势（θ_true−φ），输出补上 φ 使 θ̂ 对齐真实转子角，
    // 否则稳态 theta_err 随转速升高而增大。补偿用低通后的 ω̂ 以免补偿角抖动。
    double theta_out = self->theta_e;
    if (self->emf_lpf_wc > 0.0) {
        theta_out = wrap_2pi(self->theta_e + atan2(self->omega_filt, self->emf_lpf_wc));
    }

    PLLOutput out;
    out.theta_e = theta_out;
    out.omega_e = self->omega_filt;
    out.omega_m = self->omega_filt / (double)self->p;
    out.eps = eps;
    return out;
}
