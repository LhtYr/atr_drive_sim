#include "solver/Solver.h"
#include "solver/BoundaryProfile.h"

// 边界条件默认值（常量剖面的来源）：默认工况无供气，涡轮仅容腔吹除瞬态；
// 引气工况经参数文件覆盖（见 src/param_file.c 的 boundary.* 键）。
// 独立成文件是因为特性工具（build/char）不链接求解器，但经
// src/param_file.c 同样需要该符号。
SystemInput boundary_defaults = {
    .ud_ref = 0.0,           // 占位：实际控制参考在每拍由控制链回填
    .uq_ref = 0.0,

    .use_ctrl_theta = 0,     // 默认有感：SVPWM 用真实电角（行为与改动前一致）
    .theta_ctrl_e = 0.0,     // 占位：无感模式由控制拍注入冻结估计角

    .turb_p_in = 100000.0,   // 引气总压 Pa
    .turb_p_out = 100000.0,  // 出口背压 Pa
    .turb_T_in = 1000.0,     // 引气总温 K
    .turb_m_dot = 0.0,       // 引气流量上限 kg/s（0 = 无供气）

    .pump_p_suction = 100000.0,
    .pump_p_downstream = 100000.0,

    .batt_p_load = 0.0,
};

// 供气斜坡剖面默认参数（发电稳压能力边界测试用，
// 经参数文件 boundary.m_dot_start / boundary.ramp_time / boundary.ramp_t_start 覆盖）
BoundaryRampConfig boundary_ramp_config = {
    .m_dot_start = 0.3,   // kg/s
    .ramp_time = 10.0,    // s
    .ramp_t_start = 0.0,  // s（默认 0：t < 0 不存在，与原斜坡逐拍一致；
                          //      任务剖面经 boundary.ramp_t_start 延后起坡，起动段无供气）
};
