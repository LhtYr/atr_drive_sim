#ifndef BOUNDARY_PROFILE_H
#define BOUNDARY_PROFILE_H

#include "Solver.h"

// 边界条件剖面：按仿真时间逐拍生成系统外部输入（控制参考除外）。
// 工况仿真（阶跃/斜坡/查表等）只需实现自定义 evaluate 并装配 BoundaryProfile。
typedef struct BoundaryProfile {
    SystemInput (*evaluate)(double t_s, void* ctx);   // 时间 -> 外部输入
    void* ctx;                                        // evaluate 的用户上下文
} BoundaryProfile;

// 常量剖面：始终返回 ctx 指向的默认输入（ctx = const SystemInput*）
SystemInput boundary_profile_constant(double t_s, void* ctx);

// 构造常量剖面（默认输入须在整个仿真期间存活）
BoundaryProfile boundary_profile_make_constant(const SystemInput* defaults);

// 供气斜坡剖面参数（经参数文件 boundary.m_dot_start / boundary.ramp_time /
// boundary.ramp_t_start 覆盖，默认值定义于 BoundaryDefaults.c）
typedef struct BoundaryRampConfig {
    double m_dot_start;   // 斜坡起始引气流量 kg/s
    double ramp_time;     // 斜坡时长 s（到达后保持 base->turb_m_dot 上限）
    double ramp_t_start;  // 延迟起坡时刻 s：t < ramp_t_start 时 turb_m_dot = 0（完全无供气），
                          // 自 ramp_t_start 起在 ramp_time 内由 m_dot_start 线性升至 base->turb_m_dot。
                          // 默认 0 时行为与原斜坡逐拍一致（t_rel = t_s - 0 = t_s）
} BoundaryRampConfig;

extern BoundaryRampConfig boundary_ramp_config;

// 构造引气流量（延迟）斜坡剖面：t < ramp_t_start 时 turb_m_dot = 0（无供气）；
// 自 ramp_t_start 起在 ramp_time 内从 m_dot_start 线性升至 base->turb_m_dot 后保持；
// 其余边界项沿用 base（base 须在仿真期间存活）。ramp_t_start = 0 时退化为原斜坡。
BoundaryProfile boundary_profile_make_m_dot_ramp(const SystemInput* base);

#endif // BOUNDARY_PROFILE_H
