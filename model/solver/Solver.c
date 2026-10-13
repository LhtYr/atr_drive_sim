#include "solver/Solver.h"
#include "solver/RK4.h"
#include "SVPWM.h"
#include "Inverter.h"
#include "Motor.h"
#include "Shaft.h"
#include "Pump.h"
#include "Turbine.h"
#include "Battery.h"
#include "utils/numeric_utils.h"
#include "utils/transforms.h"

#include <math.h>

SystemState state;
SolverSwitches solver_switches = {
    .enable_inverter = 1,
    .enable_motor = 1,
    .enable_shaft = 1,
    .enable_pump = 1,
    .enable_turbine = 1,
    .enable_battery = 1,
    .enable_stage_iteration = 0,        // 是否启用子步迭代（启用后会比较慢），默认关闭
    .disabled_pump_torque = -0.0,       // 如果 pump 模块被禁用，使用该值作为泵的输出转矩（负值表示泵吸收机械能）
    .disabled_turbine_torque = 35,      // 如果 turbine 模块被禁用，使用该值作为涡轮的输出转矩（正值表示涡轮输出机械能）
    .disabled_battery_voltage = 540,    // 如果 battery 模块被禁用，使用该值作为电池的输出电压（V）
};

#define STAGE_ITER_MAX 80           // 每个仿真步内的最大子步迭代次数（如果启用子步迭代）
#define STAGE_ITER_RELAX 0.6        // 迭代变量更新的松弛因子，范围 (0, 1]，值越小越稳定但收敛更慢
#define STAGE_ITER_EPS_OMEGA 1e-4   // 机械角速度的迭代收敛阈值：如果启用子步迭代，要求每次迭代中机械角速度的变化小于该值才认为收敛
#define STAGE_ITER_EPS_PRESS 1e-2   // 压力的迭代收敛阈值：如果启用子步迭代，要求每次迭代中压力的变化小于该值才认为收敛

typedef struct CoupledIterVars {
    double v_dc;
    double omega_mech;
    double p_discharge;
    double p_plenum;
    double p_downstream;
} CoupledIterVars;

typedef struct StageEval {
    SystemAlgebraic alg;
    BatteryDeriv batt_deriv;
    InverterEval inverter;
    double battery_p_total;
} StageEval;

// 评估函数：评估系统的代数变量和电池状态导数
static StageEval evaluate_stage_algebraic_once(const SystemState* state,
                                               const SystemInput* input,
                                               double sim_time_s,
                                               const CoupledIterVars* vars) {
    StageEval out = {0};
    const double theta_e = motor_params.p * state->theta_mech;
    // SVPWM 逆 Park 电压合成角：无感模式（use_ctrl_theta=1）用控制侧注入的估计角，
    // 否则用真实电角（有感默认路径，与改动前完全一致）。
    // 注意：下方物理侧电流重构（dq_to_abc）始终用真实角 theta_e——state->id/iq 定义在真实转子系。
    const double theta_svpwm = input->use_ctrl_theta ? input->theta_ctrl_e : theta_e;
    double i_alpha_unused = 0.0;
    double i_beta_unused = 0.0;
    double ia_state = 0.0;
    double ib_state = 0.0;
    double ic_state = 0.0;
    dq_to_abc(state->id, state->iq, cos(theta_e), sin(theta_e),
              &i_alpha_unused, &i_beta_unused,
              &ia_state, &ib_state, &ic_state);

    const SVPWMOutput svpwm_out = solver_switches.enable_inverter
        ? svpwm_evaluate(&svpwm, vars->v_dc, input->ud_ref, input->uq_ref, theta_svpwm)
        : (SVPWMOutput){0};
    out.inverter = solver_switches.enable_inverter
        ? inverter_evaluate(vars->v_dc, &svpwm_out, ia_state, ib_state, ic_state, sim_time_s,
                            inverter_params.use_switching_model, inverter_params.f_sw)
        : (InverterEval){0};

    const BatteryState batt_state = {
        .soc = state->soc_batt,
        .v_rc1 = state->v_rc1_batt,
        .v_rc2 = state->v_rc2_batt,
    };
    if (solver_switches.enable_battery) {
        // 功率符号约定：正 = 放电，负 = 充电（制动回馈），允许双向流动
        out.battery_p_total = input->batt_p_load + out.inverter.p_inv;
        out.alg.battery = battery_evaluate_from_state(&battery_params, out.battery_p_total, &batt_state);
        out.batt_deriv = battery_evaluate_deriv_from_state(&battery_params, out.alg.battery.i_batt, &batt_state);
    } else {
        out.battery_p_total = input->batt_p_load;
        out.alg.battery.i_batt = 0.0;
        out.alg.battery.v_batt = solver_switches.disabled_battery_voltage;
        out.alg.battery.soc_next = state->soc_batt;
        out.alg.battery.v_rc1_next = state->v_rc1_batt;
        out.alg.battery.v_rc2_next = state->v_rc2_batt;
        out.batt_deriv = (BatteryDeriv){0};
    }

    out.alg.motor = solver_switches.enable_motor
        ? motor_evaluate_electrical(&motor_params,
                                    state->id,
                                    state->iq,
                                    state->theta_mech,
                                    vars->omega_mech,
                                    out.inverter.ua,
                                    out.inverter.ub,
                                    out.inverter.uc)
        : (ElectricalEval){0};

    out.alg.turbine = solver_switches.enable_turbine
        ? turbine_evaluate(&turbine_params,
                           input->turb_p_in,
                           input->turb_p_out,
                           input->turb_T_in,
                           input->turb_m_dot,
                           vars->p_plenum,
                           vars->p_downstream,
                           vars->omega_mech)
        : (TurbineEval){
            .dp_plenum = 0.0,
            .dp_downstream = 0.0,
            .T_turb = solver_switches.disabled_turbine_torque,
            .m_dot_in = 0.0,
            .m_dot_turb = 0.0,
            .m_dot_out = 0.0,
            .p_plenum = vars->p_plenum,
            .p_downstream = vars->p_downstream,
        };

    out.alg.pump = solver_switches.enable_pump
        ? pump_evaluate(&pump_params,
                        input->pump_p_suction,
                        input->pump_p_downstream,
                        vars->p_discharge,
                        vars->omega_mech)
        : (PumpEval){
            .dp_discharge = 0.0,
            .T_pump = solver_switches.disabled_pump_torque,
            .Q_pump = 0.0,
            .Q_out = 0.0,
            .dp_pump = 0.0,
            .H = 0.0,
            .p_discharge = vars->p_discharge,
        };

    out.alg.shaft = solver_switches.enable_shaft
        ? shaft_evaluate_mechanics(&shaft_params,
                                   vars->omega_mech,
                                   out.alg.motor.T_em,
                                   out.alg.turbine.T_turb,
                                   out.alg.pump.T_pump)
        : (MechanicalEval){0};

    return out;
}

void solver_set_switches(SolverSwitches switches) {
    solver_switches = switches;
}

// 用各组件当前观测量初始化求解器内部状态向量（须在首次 step 之前调用）
void solver_init_state(void) {
    state.id = motor_view.id;
    state.iq = motor_view.iq;
    state.theta_mech = shaft_view.theta_mech;
    state.omega_mech = shaft_view.omega_mech;
    state.p_discharge = pump_view.p_discharge;
    state.p_plenum = turbine_view.p_plenum;
    state.p_downstream = turbine_view.p_downstream;
    state.soc_batt = battery_view.soc;
    state.v_rc1_batt = battery_view.v_rc1;
    state.v_rc2_batt = battery_view.v_rc2;
}

// 系统评估函数：根据当前系统状态和输入参考，评估系统的状态导数和代数变量
EvalResult evaluate(const SystemState state, const SystemInput input, double sim_time_s, double dt){
    EvalResult result = {0};

    // 用于子步迭代的耦合变量初始值，直接使用当前系统状态和输入计算得到的相关变量作为初始值
    CoupledIterVars iter_vars = {   
        .v_dc = solver_switches.enable_battery ? battery_view.v_batt : solver_switches.disabled_battery_voltage,
        .omega_mech = state.omega_mech,
        .p_discharge = state.p_discharge,
        .p_plenum = state.p_plenum,
        .p_downstream = state.p_downstream,
    };

    StageEval stage_eval = {0};

    if (solver_switches.enable_stage_iteration) {
        for (int it = 0; it < STAGE_ITER_MAX; ++it) {
            stage_eval = evaluate_stage_algebraic_once(&state, &input, sim_time_s, &iter_vars);

            const double vdc_target = solver_switches.enable_battery
                ? stage_eval.alg.battery.v_batt
                : solver_switches.disabled_battery_voltage;
            const double omega_target = state.omega_mech + dt * stage_eval.alg.shaft.domega_mech;
            const double p_discharge_target = state.p_discharge + dt * stage_eval.alg.pump.dp_discharge;
            const double p_plenum_target = state.p_plenum + dt * stage_eval.alg.turbine.dp_plenum;
            const double p_downstream_target = state.p_downstream + dt * stage_eval.alg.turbine.dp_downstream;

            const double vdc_new = iter_vars.v_dc + STAGE_ITER_RELAX * (vdc_target - iter_vars.v_dc);
            const double omega_new = iter_vars.omega_mech + STAGE_ITER_RELAX * (omega_target - iter_vars.omega_mech);
            const double p_discharge_new = iter_vars.p_discharge + STAGE_ITER_RELAX * (p_discharge_target - iter_vars.p_discharge);
            const double p_plenum_new = iter_vars.p_plenum + STAGE_ITER_RELAX * (p_plenum_target - iter_vars.p_plenum);
            const double p_downstream_new = iter_vars.p_downstream + STAGE_ITER_RELAX * (p_downstream_target - iter_vars.p_downstream);

            const double d_vdc = fabs(vdc_new - iter_vars.v_dc);
            const double d_omega = fabs(omega_new - iter_vars.omega_mech);
            const double d_p_discharge = fabs(p_discharge_new - iter_vars.p_discharge);
            const double d_p_plenum = fabs(p_plenum_new - iter_vars.p_plenum);
            const double d_p_downstream = fabs(p_downstream_new - iter_vars.p_downstream);
            const double d_p_max = fmax(d_p_discharge, fmax(d_p_plenum, d_p_downstream));

            iter_vars.v_dc = vdc_new;
            iter_vars.omega_mech = omega_new;
            iter_vars.p_discharge = p_discharge_new;
            iter_vars.p_plenum = p_plenum_new;
            iter_vars.p_downstream = p_downstream_new;

            if (d_vdc < STAGE_ITER_EPS_PRESS && d_omega < STAGE_ITER_EPS_OMEGA && d_p_max < STAGE_ITER_EPS_PRESS) {
                break;
            }
        }
    }

    stage_eval = evaluate_stage_algebraic_once(&state, &input, sim_time_s, &iter_vars);
    result.alg = stage_eval.alg;

    // 从组件评估结果中提取状态导数
    result.deriv.did = result.alg.motor.did;
    result.deriv.diq = result.alg.motor.diq;
    result.deriv.dtheta_mech = result.alg.shaft.dtheta_mech;
    result.deriv.domega_mech = result.alg.shaft.domega_mech;
    result.deriv.dp_discharge = result.alg.pump.dp_discharge;
    result.deriv.dp_plenum = result.alg.turbine.dp_plenum;
    result.deriv.dp_downstream = result.alg.turbine.dp_downstream;
    result.deriv.dsoc_batt = stage_eval.batt_deriv.dsoc;
    result.deriv.dv_rc1_batt = stage_eval.batt_deriv.dv_rc1;
    result.deriv.dv_rc2_batt = stage_eval.batt_deriv.dv_rc2;
    return result;
}

void step(double t_s, double dt, const SystemInput input){
    EvalResult e1 = evaluate(state, input, t_s, dt);
    SystemState s2 = add_scaled(state, e1.deriv, 0.5 * dt);
    EvalResult e2 = evaluate(s2, input, t_s + 0.5 * dt, dt);
    SystemState s3 = add_scaled(state, e2.deriv, 0.5 * dt);
    EvalResult e3 = evaluate(s3, input, t_s + 0.5 * dt, dt);
    SystemState s4 = add_scaled(state, e3.deriv, dt);
    EvalResult e4 = evaluate(s4, input, t_s + dt, dt);

    state = combine_rk4(state, e1.deriv, e2.deriv, e3.deriv, e4.deriv, dt);
    state.soc_batt = clamp_range(state.soc_batt, 0.0, 1.0); // 确保 SOC 始终在合理范围内

    // 用新状态重新评估一次，保证观测量与积分后的状态一致
    EvalResult efinal = evaluate(state, input, t_s + dt, dt);

    const double vdc_final = solver_switches.enable_battery
        ? efinal.alg.battery.v_batt
        : solver_switches.disabled_battery_voltage;

    const SVPWMOutput svpwm_out = solver_switches.enable_inverter
        ? svpwm_evaluate(&svpwm,
                         vdc_final,
                         input.ud_ref,
                         input.uq_ref,
                         input.use_ctrl_theta ? input.theta_ctrl_e : motor_params.p * state.theta_mech)
        : (SVPWMOutput){0};
    const InverterEval inv_eval = solver_switches.enable_inverter
        ? inverter_evaluate(vdc_final,
                            &svpwm_out,
                            efinal.alg.motor.ia,
                            efinal.alg.motor.ib,
                            efinal.alg.motor.ic,
                            t_s + dt,
                            inverter_params.use_switching_model,
                            inverter_params.f_sw)
        : (InverterEval){.v_dc = vdc_final};
    inverter_apply_state(&inverter_view, &inv_eval);

    motor_apply_state(&motor_view, inv_eval.ua, inv_eval.ub, inv_eval.uc, state.id, state.iq, &efinal.alg.motor);
    shaft_apply_state(&shaft_view, state.theta_mech, state.omega_mech);
    turbine_apply_state(&turbine_view, &turbine_params, state.p_plenum, state.p_downstream, &efinal.alg.turbine);
    pump_apply_state(&pump_view, &pump_params, state.omega_mech, state.p_discharge, &efinal.alg.pump);
    BatteryEval batt_final = efinal.alg.battery;
    if (solver_switches.enable_battery) {
        const double p_total = clamp_floor(input.batt_p_load + inv_eval.p_inv, 0.0);
        batt_final.soc_next = state.soc_batt;
        batt_final.v_rc1_next = state.v_rc1_batt;
        batt_final.v_rc2_next = state.v_rc2_batt;
        battery_apply_state(&battery_view, p_total, &batt_final);
    } else {
        batt_final.i_batt = 0.0;
        batt_final.v_batt = solver_switches.disabled_battery_voltage;
        batt_final.soc_next = state.soc_batt;
        batt_final.v_rc1_next = state.v_rc1_batt;
        batt_final.v_rc2_next = state.v_rc2_batt;
        battery_apply_state(&battery_view, input.batt_p_load, &batt_final);
    }
}


