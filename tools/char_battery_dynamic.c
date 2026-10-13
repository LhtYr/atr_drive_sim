// 电池动态响应特性扫描：固定 SOC 下的阶跃负载响应。
//
// 与 char_battery.c 的准静态扫描对应（同锚点 SOC=0.8），但释放极化状态：
// t=0 施加阶跃负载功率，对两 RC 极化支路 RK4 积分，
// 考察“R0 瞬时跳变 + 极化指数漂移”的端口动态响应。
// SOC 在测试中保持恒定：库仑积分时间常数远大于极化时间常数
//（tau1 = r1*c1 = 7.2 s，tau2 = r2*c2 = 32 s），属更慢的时间尺度。
//
// 功率档固定 8 档（±0.5/±1/±1.5/±2 MW，整值便于族线图例），
// 时间格固定 1 s 间隔至 120 s（功率×时间满网格，三维曲面无空洞）。
//
// 输出列：
//   t_s        时间 s
//   p_load_W   阶跃负载功率 W（正 = 放电，负 = 充电）
//   i_batt_A   电池电流 A
//   v_batt_V   电池端电压 V
//   v_rc1_V    极化支路 1 电压 V
//   v_rc2_V    极化支路 2 电压 V

#include <stdio.h>

#include "Battery.h"
#include "char_sweep.h"

#define BATT_SOC       0.8		// 固定 SOC（与准静态扫描同）
#define BATT_T_END_S   120.0	// 阶跃响应时长 s（约 4 倍 tau2）
#define BATT_SUBSTEPS  20		// 每个输出间隔内的 RK4 积分步数
#define BATT_N_T       121		// 每档功率输出时间样点数（1 s 间隔）

static const double STEP_POWERS_W[] = {
	-2.0e6, -1.5e6, -1.0e6, -0.5e6, 0.5e6, 1.0e6, 1.5e6, 2.0e6,
};
#define BATT_N_STEPS ((int)(sizeof(STEP_POWERS_W) / sizeof(STEP_POWERS_W[0])))

// 极化状态导数：电流按当前状态由端口约束回解
static void polar_deriv(const BatteryParams* params, double p_load,
                        const BatteryState* st, double* dv_rc1, double* dv_rc2) {
	const BatteryEval e = battery_evaluate_from_state(params, p_load, st);
	const BatteryDeriv d = battery_evaluate_deriv_from_state(params, e.i_batt, st);
	*dv_rc1 = d.dv_rc1;
	*dv_rc2 = d.dv_rc2;
}

// 极化状态单步 RK4（不积分 SOC）
static void polar_rk4(const BatteryParams* params, double p_load, BatteryState* st, double dt) {
	double k1a, k1b, k2a, k2b, k3a, k3b, k4a, k4b;
	BatteryState s;

	polar_deriv(params, p_load, st, &k1a, &k1b);
	s = *st;
	s.v_rc1 += 0.5 * dt * k1a;
	s.v_rc2 += 0.5 * dt * k1b;
	polar_deriv(params, p_load, &s, &k2a, &k2b);
	s = *st;
	s.v_rc1 += 0.5 * dt * k2a;
	s.v_rc2 += 0.5 * dt * k2b;
	polar_deriv(params, p_load, &s, &k3a, &k3b);
	s = *st;
	s.v_rc1 += dt * k3a;
	s.v_rc2 += dt * k3b;
	polar_deriv(params, p_load, &s, &k4a, &k4b);

	st->v_rc1 += dt / 6.0 * (k1a + 2.0 * k2a + 2.0 * k3a + k4a);
	st->v_rc2 += dt / 6.0 * (k1b + 2.0 * k2b + 2.0 * k3b + k4b);
}

int char_sweep_battery_dynamic(FILE* f, int n) {
	(void)n;	// 时间格与功率档固定，与网格数无关
	fprintf(f, "t_s,p_load_W,i_batt_A,v_batt_V,v_rc1_V,v_rc2_V\n");

	const double dt_out = BATT_T_END_S / (double)(BATT_N_T - 1);
	const double dt = dt_out / (double)BATT_SUBSTEPS;
	for (int s = 0; s < BATT_N_STEPS; s++) {
		BatteryState st = { .soc = BATT_SOC, .v_rc1 = 0.0, .v_rc2 = 0.0 };
		for (int i = 0; i < BATT_N_T; i++) {
			const BatteryEval e = battery_evaluate_from_state(&battery_params,
				STEP_POWERS_W[s], &st);
			fprintf(f, "%.3f,%.1f,%.6f,%.6f,%.6f,%.6f\n",
				(double)i * dt_out, STEP_POWERS_W[s],
				e.i_batt, e.v_batt, st.v_rc1, st.v_rc2);
			for (int k = 0; k < BATT_SUBSTEPS; k++) {
				polar_rk4(&battery_params, STEP_POWERS_W[s], &st, dt);
			}
		}
	}
	return 0;
}
