#include "Voltage.h"
#include <assert.h>
#include <math.h>
#include "utils/numeric_utils.h"

Voltage voltage = {
    .Ts = 5e-5,
    .kp = 0.5,
    .ki = 50.0,
    .iq_ref_max = 540,
    .iq_ref_min = -540,
    .anti_windup_gain = 0,  // 反电动势前馈项的增益，0 表示不使用反电动势前馈
    .int_v = 0.0,
};

VoltageOutput voltage_evaluate(Voltage* controller,
                                double v_ref,
                                double v_meas)
{
    VoltageOutput out;
    // 发电（iq<0）把机械功率注入母线、抬升 v_dc，故稳定方向为
    // iq_ref ∝ (v_meas - v_ref)：v_dc 偏低 -> 加大发电抬升，v_dc 偏高 -> 减小发电回落。
    // （若取反号 v_ref - v_meas 则成正反馈：实测 v_ref=538 时 v_dc 反向漂至 543+）
    out.ev = v_meas - v_ref;
    const double iq_unsat = controller->kp * out.ev + controller->int_v;
    if (iq_unsat < controller->iq_ref_min) {
		out.iq_ref = controller->iq_ref_min;
	} else if (iq_unsat > controller->iq_ref_max) {
		out.iq_ref = controller->iq_ref_max;
	} else {
		out.iq_ref = iq_unsat;
	}
    const double diq_aw = out.iq_ref - iq_unsat;
    controller->int_v += (controller->ki * out.ev + controller->anti_windup_gain * diq_aw) * controller->Ts;
    return out;
}
