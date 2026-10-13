#include "RK4.h"

SystemState add_scaled(const SystemState state, const SystemDeriv deriv, double scale){
    SystemState out;
    out.id = state.id + scale * deriv.did;
    out.iq = state.iq + scale * deriv.diq;
    out.theta_mech = state.theta_mech + scale * deriv.dtheta_mech;
    out.omega_mech = state.omega_mech + scale * deriv.domega_mech;
    out.p_discharge = state.p_discharge + scale * deriv.dp_discharge;
    out.p_plenum = state.p_plenum + scale * deriv.dp_plenum;
    out.p_downstream = state.p_downstream + scale * deriv.dp_downstream;
    out.soc_batt = state.soc_batt + scale * deriv.dsoc_batt;
    out.v_rc1_batt = state.v_rc1_batt + scale * deriv.dv_rc1_batt;
    out.v_rc2_batt = state.v_rc2_batt + scale * deriv.dv_rc2_batt;
    return out;
}

SystemState combine_rk4(const SystemState state,
                        const SystemDeriv k1,
                        const SystemDeriv k2,
                        const SystemDeriv k3,
                        const SystemDeriv k4,
                        double Ts){
    SystemState out;
    out.id = state.id + (Ts / 6.0) * (k1.did + 2.0 * k2.did + 2.0 * k3.did + k4.did);
    out.iq = state.iq + (Ts / 6.0) * (k1.diq + 2.0 * k2.diq + 2.0 * k3.diq + k4.diq);
    out.theta_mech = state.theta_mech + (Ts / 6.0) * (k1.dtheta_mech + 2.0 * k2.dtheta_mech + 2.0 * k3.dtheta_mech + k4.dtheta_mech);
    out.omega_mech = state.omega_mech + (Ts / 6.0) * (k1.domega_mech + 2.0 * k2.domega_mech + 2.0 * k3.domega_mech + k4.domega_mech);
    out.p_discharge = state.p_discharge + (Ts / 6.0) * (k1.dp_discharge + 2.0 * k2.dp_discharge + 2.0 * k3.dp_discharge + k4.dp_discharge);
    out.p_plenum = state.p_plenum + (Ts / 6.0) * (k1.dp_plenum + 2.0 * k2.dp_plenum + 2.0 * k3.dp_plenum + k4.dp_plenum);
    out.p_downstream = state.p_downstream + (Ts / 6.0) * (k1.dp_downstream + 2.0 * k2.dp_downstream + 2.0 * k3.dp_downstream + k4.dp_downstream);
    out.soc_batt = state.soc_batt + (Ts / 6.0) * (k1.dsoc_batt + 2.0 * k2.dsoc_batt + 2.0 * k3.dsoc_batt + k4.dsoc_batt);
    out.v_rc1_batt = state.v_rc1_batt + (Ts / 6.0) * (k1.dv_rc1_batt + 2.0 * k2.dv_rc1_batt + 2.0 * k3.dv_rc1_batt + k4.dv_rc1_batt);
    out.v_rc2_batt = state.v_rc2_batt + (Ts / 6.0) * (k1.dv_rc2_batt + 2.0 * k2.dv_rc2_batt + 2.0 * k3.dv_rc2_batt + k4.dv_rc2_batt);
    return out;
}