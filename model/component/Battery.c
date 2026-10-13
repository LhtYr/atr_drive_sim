#include "Battery.h"
#include "utils/numeric_utils.h"

#include <math.h>

BatteryParams battery_params = {
    .capacity = 40.0,       // 电池容量 Ah
    .v_oc = 540.0,          // 开路电压 V
    .r_int = 0.01,          // 欧姆内阻 R0 Ohm
    .r1 = 0.006,            // 极化支路1电阻 Ohm
    .c1 = 1200.0,           // 极化支路1电容 F
    .r2 = 0.004,            // 极化支路2电阻 Ohm
    .c2 = 8000.0,           // 极化支路2电容 F
    .eta_coulomb = 1.0,     // 库伦效率（放电/充电近似统一）
};

BatteryView battery_view = {
    .soc = 0.8,             // 初始 SOC
    .v_rc1 = 0.0,           // 初始极化电压 V
    .v_rc2 = 0.0,           // 初始极化电压 V
    .p_load = 0.0,          // 负载功率 W
    .i_batt = 0.0,          // 电池电流 A
    .v_batt = 540.0,        // 电池端电压 V（首步 v_dc 初值依赖此值）
};

BatteryDeriv battery_evaluate_deriv_from_state(const BatteryParams* params,
                                               double i_batt,
                                               const BatteryState* st) {
    BatteryDeriv deriv;
    const double r1_safe = clamp_floor(params->r1, 1e-9);
    const double c1_safe = clamp_floor(params->c1, 1e-9);
    const double r2_safe = clamp_floor(params->r2, 1e-9);
    const double c2_safe = clamp_floor(params->c2, 1e-9);
    const double q_coulomb = clamp_floor(params->capacity, 1e-9) * 3600.0;
    const double eta = clamp_floor(params->eta_coulomb, 0.0);

    deriv.dv_rc1 = -st->v_rc1 / (r1_safe * c1_safe) + i_batt / c1_safe;
    deriv.dv_rc2 = -st->v_rc2 / (r2_safe * c2_safe) + i_batt / c2_safe;
    deriv.dsoc = -eta * i_batt / q_coulomb;
    return deriv;
}

BatteryEval battery_evaluate_from_state(const BatteryParams* params,
                                        double p_load_in,
                                        const BatteryState* st) {
    BatteryEval out;
    const double v_oc_safe = clamp_floor(params->v_oc, 1.0);
    const double r0_safe = clamp_floor(params->r_int, 1e-9); 

    // 二阶 RC 等效端口：Vt = Voc - Vrc1 - Vrc2 - I*R0
    const double v_th = v_oc_safe - st->v_rc1 - st->v_rc2;  // 电压源部分的 Thevenin 电压 Vt

    // 由端口功率约束 P = (Vt - I*R0)*I 求 I，代入得：R0*I^2 - Vt*I + P = 0。
    // P > 0（放电）取小电流根；P < 0（充电）判别式恒大于 Vt^2，根自动为负。
    double i_sol = 0.0;
    const double p_max = (v_th * v_th) / (4.0 * r0_safe);
    const double p_req = fmin(p_load_in, p_max);
    const double disc = fmax(v_th * v_th - 4.0 * r0_safe * p_req, 0.0);
    i_sol = (v_th - sqrt(disc)) / (2.0 * r0_safe);
    out.i_batt = i_sol;
    out.v_batt = fmax(v_th - out.i_batt * r0_safe, 0.0);

    // 默认先回填当前状态，调用方可按导数更新 next。
    out.v_rc1_next = st->v_rc1;
    out.v_rc2_next = st->v_rc2;
    out.soc_next = clamp_range(st->soc, 0.0, 1.0);
    return out;
}

void battery_apply_state(BatteryView* view, double p_load_in, const BatteryEval* eval) {
    view->i_batt = eval->i_batt;
    view->v_batt = eval->v_batt;
    view->v_rc1 = eval->v_rc1_next;
    view->v_rc2 = eval->v_rc2_next;
    view->soc = eval->soc_next;
    view->p_load = p_load_in;
}
