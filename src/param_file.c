#include "utils/param_file.h"

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "Motor.h"
#include "Shaft.h"
#include "Pump.h"
#include "Turbine.h"
#include "Battery.h"
#include "Inverter.h"
#include "Voltage.h"
#include "IFStart.h"
#include "STSMO.h"
#include "PLL.h"
#include "sensorless.h"
#include "solver/Solver.h"
#include "solver/BoundaryProfile.h"

// 仿真级参数 sim.* 键的绑定目标（声明见 utils/param_file.h）
double sim_v_ref = 540.0;

// 无感起动切换转速 shim：sensorless.if.omega_switch_rpm 以 rpm 给定，须换算为 rad/s 机械角
// 才能写入 ifstart.omega_switch_m，故先绑定到此中间变量（<0 = 参数文件未设置，沿用内置/CLI 值），
// 由 main.c 在无感初始化时换算（声明见 utils/param_file.h）。
double sensorless_if_omega_switch_rpm = -1.0;

// 参数绑定表：键名 -> 目标字段。is_int 标记整型字段（数值取整后写入）。
typedef struct ParamBinding {
    const char* name;
    void* ptr;
    int is_int;
} ParamBinding;

static const ParamBinding bindings[] = {
    {"motor.Rs", &motor_params.Rs, 0},
    {"motor.Ld", &motor_params.Ld, 0},
    {"motor.Lq", &motor_params.Lq, 0},
    {"motor.psi_f", &motor_params.psi_f, 0},
    {"motor.p", &motor_params.p, 1},

    {"shaft.J_motor", &shaft_params.J_motor, 0},
    {"shaft.J_turb", &shaft_params.J_turb, 0},
    {"shaft.J_pump", &shaft_params.J_pump, 0},
    {"shaft.b_motor", &shaft_params.b_motor, 0},
    {"shaft.b_turb", &shaft_params.b_turb, 0},
    {"shaft.b_pump", &shaft_params.b_pump, 0},

    {"battery.capacity", &battery_params.capacity, 0},
    {"battery.v_oc", &battery_params.v_oc, 0},
    {"battery.r_int", &battery_params.r_int, 0},
    {"battery.r1", &battery_params.r1, 0},
    {"battery.c1", &battery_params.c1, 0},
    {"battery.r2", &battery_params.r2, 0},
    {"battery.c2", &battery_params.c2, 0},
    {"battery.eta_coulomb", &battery_params.eta_coulomb, 0},
    {"battery.soc_init", &battery_view.soc, 0},

    {"pump.rho", &pump_params.rho, 0},
    {"pump.omega_floor", &pump_params.omega_floor, 0},
    {"pump.V_out", &pump_params.V_out, 0},
    {"pump.beta_eff", &pump_params.beta_eff, 0},
    {"pump.R_out", &pump_params.R_out, 0},
    {"pump.p_floor", &pump_params.p_floor, 0},
    {"pump.h_qq", &pump_params.h_qq, 0},
    {"pump.h_qn", &pump_params.h_qn, 0},
    {"pump.h_nn", &pump_params.h_nn, 0},
    {"pump.e_x3", &pump_params.e_x3, 0},
    {"pump.e_x2", &pump_params.e_x2, 0},
    {"pump.e_x1", &pump_params.e_x1, 0},
    {"pump.e_x0", &pump_params.e_x0, 0},
    {"pump.p_dis_init", &pump_view.p_discharge, 0},

    {"turbine.gamma", &turbine_params.gamma, 0},
    {"turbine.gas_R", &turbine_params.gas_R, 0},
    {"turbine.eta_turb", &turbine_params.eta_turb, 0},
    {"turbine.omega_floor", &turbine_params.omega_floor, 0},
    {"turbine.V_plenum", &turbine_params.V_plenum, 0},
    {"turbine.V_downstream", &turbine_params.V_downstream, 0},
    {"turbine.Cd_duct", &turbine_params.Cd_duct, 0},
    {"turbine.A_duct", &turbine_params.A_duct, 0},
    {"turbine.Cd_nozzle", &turbine_params.Cd_nozzle, 0},
    {"turbine.A_nozzle", &turbine_params.A_nozzle, 0},
    {"turbine.Cd_outlet", &turbine_params.Cd_outlet, 0},
    {"turbine.A_outlet", &turbine_params.A_outlet, 0},
    {"turbine.p_floor", &turbine_params.p_floor, 0},
    {"turbine.mean_radius", &turbine_params.mean_radius, 0},
    {"turbine.eta_mech", &turbine_params.eta_mech, 0},
    {"turbine.eta_uc_a", &turbine_params.eta_uc_a, 0},
    {"turbine.eta_uc_b", &turbine_params.eta_uc_b, 0},
    {"turbine.eta_uc_c", &turbine_params.eta_uc_c, 0},
    {"turbine.r_choke_nozzle", &turbine_params.r_choke_nozzle, 0},

    // 涡轮容腔初始压力（solver_init_state 复制入状态向量；设为背压可消除起动容腔吹除瞬态，
    // 建立"无供气空载"起动工况。默认 4 bar 预充会向转子倾倒 ~73 N·m 吹除转矩，见 README）
    {"turbine.p_plenum_init", &turbine_view.p_plenum, 0},
    {"turbine.p_downstream_init", &turbine_view.p_downstream, 0},

    // 气路边界条件默认值（常量剖面来源，仿真期间保持恒定）
    {"boundary.turb_p_in", &boundary_defaults.turb_p_in, 0},
    {"boundary.turb_p_out", &boundary_defaults.turb_p_out, 0},
    {"boundary.turb_T_in", &boundary_defaults.turb_T_in, 0},
    {"boundary.turb_m_dot", &boundary_defaults.turb_m_dot, 0},

    // 供气（延迟）斜坡剖面（电压模式转速扫掠 / 任务剖面渐增供气，见 BoundaryProfile.c）
    {"boundary.m_dot_start", &boundary_ramp_config.m_dot_start, 0},
    {"boundary.ramp_time", &boundary_ramp_config.ramp_time, 0},
    {"boundary.ramp_t_start", &boundary_ramp_config.ramp_t_start, 0},

    // 电压环控制器（发电稳压工况，见 model/control/Voltage.c）
    {"voltage.kp", &voltage.kp, 0},
    {"voltage.ki", &voltage.ki, 0},
    {"voltage.iq_ref_min", &voltage.iq_ref_min, 0},
    {"voltage.iq_ref_max", &voltage.iq_ref_max, 0},

    // 无感起动（I/F + STSMO + PLL + 切换状态机；见 model/control/{IFStart,STSMO,PLL,sensorless}.c）。
    // 这些键覆盖编译期默认值，且优先于命令行 --if-current/--switch-rpm（param_file_apply 在 parse_args 之后）。
    {"sensorless.if.i_start", &ifstart.i_start, 0},
    {"sensorless.if.omega_switch_rpm", &sensorless_if_omega_switch_rpm, 0},  // rpm，main.c 换算为 rad/s 机械角
    {"sensorless.if.ramp_time", &ifstart.ramp_time, 0},
    {"sensorless.if.damp_k", &ifstart.damp_k, 0},
    {"sensorless.if.damp_lpf_hz", &ifstart.damp_lpf_hz, 0},
    {"sensorless.stsmo.k1", &stsmo.k1, 0},
    {"sensorless.stsmo.k2", &stsmo.k2, 0},
    {"sensorless.stsmo.emf_lpf_hz", &stsmo.emf_lpf_hz, 0},  // main.c 无感初始化时同步 pll.emf_lpf_wc
    {"sensorless.pll.kp", &pll.kp, 0},
    {"sensorless.pll.ki", &pll.ki, 0},
    {"sensorless.switch.angle_tol_rad", &sensorless_mgr.switch_angle_tol, 0},
    {"sensorless.switch.omega_tol", &sensorless_mgr.switch_omega_tol, 0},
    {"sensorless.switch.blend_time", &sensorless_mgr.blend_time, 0},
    // 发电切换判据（阶段 2 -> 3，任务剖面用；默认 gen_iq_th=-1e9 关闭，无这些键则不进阶段 3）
    {"sensorless.gen.iq_th", &sensorless_mgr.gen_iq_th, 0},
    {"sensorless.gen.hold_ticks", &sensorless_mgr.gen_hold_ticks, 1},
    {"sensorless.gen.iq_lpf_hz", &sensorless_mgr.gen_iq_lpf_hz, 0},

    // 仿真级参数（非组件参数，由 src/main.c 读取）
    {"sim.v_ref", &sim_v_ref, 0},

    {"inverter.f_sw", &inverter_params.f_sw, 0},
    {"inverter.use_switching_model", &inverter_params.use_switching_model, 1},
};

// 去除字符串首尾空白（原地修改，返回首字符位置）
static char* trim(char* s) {
    while (isspace((unsigned char)*s)) {
        ++s;
    }
    char* end = s + strlen(s);
    while (end > s && isspace((unsigned char)end[-1])) {
        --end;
    }
    *end = '\0';
    return s;
}

int param_file_apply(const char* path) {
    FILE* file = fopen(path, "r");
    if (file == NULL) {
        fprintf(stderr, "无法打开参数文件：%s\n", path);
        return -1;
    }

    char line[256];
    int line_no = 0;
    int status = 0;
    while (fgets(line, sizeof(line), file) != NULL) {
        ++line_no;
        char* comment = strchr(line, '#');
        if (comment != NULL) {
            *comment = '\0';
        }
        char* entry = trim(line);
        if (*entry == '\0') {
            continue;
        }

        char* eq = strchr(entry, '=');
        if (eq == NULL) {
            fprintf(stderr, "%s:%d：缺少 '='：%s\n", path, line_no, entry);
            status = -1;
            break;
        }
        *eq = '\0';
        const char* key = trim(entry);
        char* value_text = trim(eq + 1);

        char* end = NULL;
        const double value = strtod(value_text, &end);
        if (end == value_text || *end != '\0') {
            fprintf(stderr, "%s:%d：无效数值：%s\n", path, line_no, value_text);
            status = -1;
            break;
        }

        const int count = (int)(sizeof(bindings) / sizeof(bindings[0]));
        int found = 0;
        for (int i = 0; i < count; ++i) {
            if (strcmp(key, bindings[i].name) == 0) {
                if (bindings[i].is_int) {
                    *(int*)bindings[i].ptr = (int)value;
                } else {
                    *(double*)bindings[i].ptr = value;
                }
                found = 1;
                break;
            }
        }
        if (!found) {
            fprintf(stderr, "%s:%d：未知参数键：%s\n", path, line_no, key);
            status = -1;
            break;
        }
    }

    fclose(file);
    return status;
}
