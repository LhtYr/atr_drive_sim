#include "Inverter.h"	// 逆变器模型：将 DC 母线电压转换为三相交流电压，供电机使用。
// 平均模型（默认，快）与开关模型（载波比较，慢但含开关纹波）按运行期开关切换。

#include <math.h>
#include "utils/numeric_utils.h"

InverterParams inverter_params = {
	.f_sw = 50000.0,			// 开关频率 Hz（开关模型使用）
	.use_switching_model = 0,	// 默认平均模型（结果与历史基线一致）
};

InverterView inverter_view = {
	.v_dc = 540.0,
};

// 三角载波：范围 [0, 1]，频率由 f_sw 指定
static double triangle_carrier(double sim_time_s, double f_sw) {
	const double f_sw_safe = clamp_floor(f_sw, 1.0);
	const double period = 1.0 / f_sw_safe;
	double phase = fmod(sim_time_s, period) / period;
	if (phase < 0.0) {
		phase += 1.0;
	}
	return (phase < 0.5) ? (2.0 * phase) : (2.0 * (1.0 - phase));
}

InverterEval inverter_evaluate(double v_dc,
                               const SVPWMOutput* svpwm,
                               double ia,
                               double ib,
                               double ic,
                               double sim_time_s,
                               int use_switching_model,
                               double f_sw) {
	InverterEval out;
	const double v_dc_safe = magnitude_floor(v_dc, 1.0);
	out.v_dc = v_dc_safe;
	
	out.duty_a = clamp_range(svpwm->gate_u_a, 0.0, 1.0);
	out.duty_b = clamp_range(svpwm->gate_u_b, 0.0, 1.0);
	out.duty_c = clamp_range(svpwm->gate_u_c, 0.0, 1.0);

	if (use_switching_model) {
		// 开关模型：占空比与三角载波比较得到门级信号，相电压取开关电平（相对母线中点）
		const double carrier = triangle_carrier(sim_time_s, f_sw);
		out.gate_u_a = (out.duty_a >= carrier) ? 1.0 : 0.0;
		out.gate_l_a = 1.0 - out.gate_u_a;
		out.gate_u_b = (out.duty_b >= carrier) ? 1.0 : 0.0;
		out.gate_l_b = 1.0 - out.gate_u_b;
		out.gate_u_c = (out.duty_c >= carrier) ? 1.0 : 0.0;
		out.gate_l_c = 1.0 - out.gate_u_c;

		out.ua = (out.gate_u_a - 0.5) * v_dc_safe;
		out.ub = (out.gate_u_b - 0.5) * v_dc_safe;
		out.uc = (out.gate_u_c - 0.5) * v_dc_safe;
	} else {
		// 平均模型：占空比直接映射为相电压（相对母线中点），门级信号直接透传
		out.gate_u_a = svpwm->gate_u_a;
		out.gate_l_a = svpwm->gate_l_a;
		out.gate_u_b = svpwm->gate_u_b;
		out.gate_l_b = svpwm->gate_l_b;
		out.gate_u_c = svpwm->gate_u_c;
		out.gate_l_c = svpwm->gate_l_c;

		out.ua = (out.duty_a - 0.5) * v_dc_safe;
		out.ub = (out.duty_b - 0.5) * v_dc_safe;
		out.uc = (out.duty_c - 0.5) * v_dc_safe;
	}
	
	const double p_ac = out.ua * ia + out.ub * ib + out.uc * ic;	// 交流侧功率
	out.i_inv_dc = p_ac / v_dc_safe;	// 直流侧电流（假设逆变器无损耗）
	out.p_inv = v_dc * out.i_inv_dc;	// 逆变器输入功率（假设逆变器无损耗，用原始 v_dc）
	return out;
}

void inverter_apply_state(InverterView* view, const InverterEval* eval) {
	view->v_dc = eval->v_dc;
	view->ua = eval->ua;
	view->ub = eval->ub;
	view->uc = eval->uc;
	view->duty_a = eval->duty_a;
	view->duty_b = eval->duty_b;
	view->duty_c = eval->duty_c;
	view->gate_u_a = eval->gate_u_a;
	view->gate_l_a = eval->gate_l_a;
	view->gate_u_b = eval->gate_u_b;
	view->gate_l_b = eval->gate_l_b;
	view->gate_u_c = eval->gate_u_c;
	view->gate_l_c = eval->gate_l_c;
	view->i_inv_dc = eval->i_inv_dc;
	view->p_inv = eval->p_inv;
}
