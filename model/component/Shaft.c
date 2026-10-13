#include "Shaft.h"
#include <math.h>

ShaftParams shaft_params = {
    .J_motor = 0.000772,    // 电机侧转动惯量 kg.m²
    .J_turb  = 0.000001,    // 涡轮侧转动惯量 kg.m²
    .J_pump  = 0.000001,    // 泵侧转动惯量 kg.m²

    .b_motor = 0.0001,      // 电机侧摩擦系数（N·m·s）
    .b_turb  = 0.0047,         // 涡轮侧摩擦系数（N·m·s）
    .b_pump = 0.0,          // 泵侧摩擦系数（N·m·s）
};

ShaftView shaft_view = {0};

MechanicalEval shaft_evaluate_mechanics(const ShaftParams* params,
                                        double omega_mech_in,
                                        double T_em,
                                        double T_turb,
                                        double T_pump) {
    MechanicalEval out;
    double J_total = params->J_motor + params->J_turb + params->J_pump;
    double b_total = params->b_motor + params->b_turb + params->b_pump;

    out.domega_mech = (T_em + T_turb - T_pump - b_total * omega_mech_in) / J_total;
    out.dtheta_mech = omega_mech_in;
    return out;
}

void shaft_apply_state(ShaftView* view, double theta_mech_in, double omega_mech_in) {
    view->theta_mech = theta_mech_in;
    view->omega_mech = omega_mech_in;
}
