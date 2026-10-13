#include "Motor.h"
#include "utils/transforms.h"
#include <math.h>

MotorParams motor_params = {
    .Rs = 0.0051,           // 定子电阻 Ohm
    .Ld = 0.0000346,        // d 轴电感 H
    .Lq = 0.0000346,        // q 轴电感 H
    // 磁链按弱磁设计点标定：id=0 时反电动势 ω_e·ψ_f 在约 25000 rpm 处
    // 达到 SVPWM 线性区上限 v_dc/sqrt(3) ≈ 311.8 V，高于该转速需弱磁控制
    .psi_f = 0.0595,        // 永磁体磁链 Wb
    .p = 2,                 // 极对数
};

MotorView motor_view = {0};

// 电机电气侧评估（纯计算，不修改对象状态）：
ElectricalEval motor_evaluate_electrical(const MotorParams* params,
                                         double id_in,
                                         double iq_in,
                                         double theta_mech,
                                         double omega_mech,
                                         double ua_in,
                                         double ub_in,
                                         double uc_in) {
    ElectricalEval out;
    // 机械到电气映射（极对数 p）
    out.omega_e = params->p * omega_mech;
    out.theta_e = params->p * theta_mech;

    // Clarke 变换：abc -> alpha-beta（共用实现，见 utils/transforms.h）
    double u_alpha = 0.0;
    double u_beta = 0.0;
    clarke_abc_to_ab(ua_in, ub_in, uc_in, &u_alpha, &u_beta);

    // Park 变换：alpha-beta -> dq
    const double cos_theta = cos(out.theta_e);
    const double sin_theta = sin(out.theta_e);
    out.ud = u_alpha * cos_theta + u_beta * sin_theta;
    out.uq = -u_alpha * sin_theta + u_beta * cos_theta;

    // dq 电流微分方程
    const double Ld_safe = fmax(params->Ld, 1e-9);
    const double Lq_safe = fmax(params->Lq, 1e-9);
    out.did = (out.ud - params->Rs * id_in + out.omega_e * params->Lq * iq_in) / Ld_safe;
    out.diq = (out.uq - params->Rs * iq_in - out.omega_e * params->Ld * id_in - out.omega_e * params->psi_f) / Lq_safe;

    // 逆 Park + 逆 Clarke：dq -> alpha-beta -> abc（输出观测，共用实现）
    dq_to_abc(id_in, iq_in, cos_theta, sin_theta,
              &out.i_alpha, &out.i_beta, &out.ia, &out.ib, &out.ic);

    // PMSM 转矩方程（幅值不变 Park 变换：T = 3/2 · p · (...)，p 为极对数）
    out.T_em = THREE_2 * params->p * (params->psi_f * iq_in + (params->Ld - params->Lq) * id_in * iq_in);
    return out;
}

void motor_apply_state(MotorView* view,
                       double ua_in,
                       double ub_in,
                       double uc_in,
                       double id_in,
                       double iq_in,
                       const ElectricalEval* eval) {
    view->ua = ua_in;
    view->ub = ub_in;
    view->uc = uc_in;
    view->id = id_in;
    view->iq = iq_in;
    view->ud = eval->ud;
    view->uq = eval->uq;
    view->i_alpha = eval->i_alpha;
    view->i_beta = eval->i_beta;
    view->ia = eval->ia;
    view->ib = eval->ib;
    view->ic = eval->ic;
    view->omega_e = eval->omega_e;
    view->theta_e = eval->theta_e;
    view->T_em = eval->T_em;
}
