#include "BoundaryProfile.h"

#include "utils/numeric_utils.h"

SystemInput boundary_profile_constant(double t_s, void* ctx) {
    (void)t_s;
    return *(const SystemInput*)ctx;
}

BoundaryProfile boundary_profile_make_constant(const SystemInput* defaults) {
    BoundaryProfile profile = {
        .evaluate = boundary_profile_constant,
        .ctx = (void*)defaults,
    };
    return profile;
}

SystemInput boundary_profile_m_dot_ramp(double t_s, void* ctx) {
    const SystemInput* base = (const SystemInput*)ctx;
    SystemInput in = *base;
    // 延迟起坡：以 t_rel = t_s - ramp_t_start 为斜坡自变量。
    //   · t_rel < 0（t_s < ramp_t_start）：完全无供气，turb_m_dot = 0（任务剖面起动段涡轮零转矩）；
    //   · t_rel >= 0：自 ramp_t_start 起在 ramp_time 内由 m_dot_start 线性升至 base->turb_m_dot 后保持。
    // 默认 ramp_t_start = 0 时 t_rel = t_s，且 t_s >= 0 恒成立，与原斜坡逐拍一致（回归不变）。
    const double t_rel = t_s - boundary_ramp_config.ramp_t_start;
    if (t_rel < 0.0) {
        in.turb_m_dot = 0.0;
        return in;
    }
    const double frac = (boundary_ramp_config.ramp_time > 0.0)
        ? clamp_range(t_rel / boundary_ramp_config.ramp_time, 0.0, 1.0)
        : 1.0;
    in.turb_m_dot = boundary_ramp_config.m_dot_start +
                    (base->turb_m_dot - boundary_ramp_config.m_dot_start) * frac;
    return in;
}

BoundaryProfile boundary_profile_make_m_dot_ramp(const SystemInput* base) {
    BoundaryProfile profile = {
        .evaluate = boundary_profile_m_dot_ramp,
        .ctx = (void*)base,
    };
    return profile;
}
