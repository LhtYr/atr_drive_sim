// 电池准静态特性扫描：负载功率 p_load（正放电、负充电）一维网格。
//
// 状态取准静态工况点：SOC = 0.8、RC 极化电压均为 0
//（考察欧姆内阻主导的端口特性；极化动态特性需经主仿真观察）。
//
// 输出列：
//   p_load_W    负载功率 W（正 = 放电，负 = 充电）
//   i_batt_A    电池电流 A（正 = 放电，负 = 充电）
//   v_batt_V    电池端电压 V
//   p_term_W    端口电功率 = v_batt * i_batt
//   efficiency  电压口径效率（端口功率/内电源功率）：
//               放电 = v_batt/V_th，充电 = V_th/v_batt，
//               V_th = V_oc - v_rc1 - v_rc2（准静态下即 V_oc）

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Battery.h"
#include "char_sweep.h"

#define BATT_SOC       0.8		// 准静态 SOC
#define BATT_P_MAX_W   2.0e6	// 扫描功率幅值上限 W（±）

int char_sweep_battery(FILE* f, int n) {
	double* powers = malloc((size_t)n * sizeof(double));
	if (powers == NULL) {
		fprintf(stderr, "错误：内存分配失败\n");
		return -1;
	}
	char_linspace(powers, -BATT_P_MAX_W, BATT_P_MAX_W, n);

	fprintf(f, "p_load_W,i_batt_A,v_batt_V,p_term_W,efficiency\n");

	const BatteryState st = { .soc = BATT_SOC, .v_rc1 = 0.0, .v_rc2 = 0.0 };
	for (int i = 0; i < n; i++) {
		const BatteryEval e = battery_evaluate_from_state(&battery_params,
			powers[i], &st);

		const double p_term = e.v_batt * e.i_batt;
		// 端口功率由约束恒等于请求功率，往返口径效率恒为 1、无信息量；
		// 改电压口径：端口功率与内电源功率 (V_th*I) 之比，体现欧姆内阻损耗
		const double v_th = fmax(battery_params.v_oc - st.v_rc1 - st.v_rc2, 1.0);
		double efficiency = 1.0;
		if (e.i_batt > 1e-9) {					// 放电：V/V_th
			efficiency = e.v_batt / v_th;
		} else if (e.i_batt < -1e-9 && e.v_batt > 1e-9) {	// 充电：V_th/V
			efficiency = v_th / e.v_batt;
		}

		fprintf(f, "%.1f,%.6f,%.6f,%.1f,%.6f\n",
			powers[i], e.i_batt, e.v_batt, p_term, efficiency);
	}

	free(powers);
	return 0;
}
