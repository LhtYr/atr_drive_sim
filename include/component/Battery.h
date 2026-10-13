#ifndef BATTERY_H
#define BATTERY_H

// 二阶 RC 等效电路电池模型。
// 符号约定：
//   p_load > 0 放电，p_load < 0 充电（制动回馈）；
//   i_batt > 0 放电、< 0 充电；v_batt = V_th - i*R0；
//   dsoc = -eta*i_batt/Q（放电降、充电升）。

typedef struct BatteryEval {
    double i_batt;
    double v_batt;
    double v_rc1_next;
    double v_rc2_next;
    double soc_next;
} BatteryEval;

typedef struct BatteryState {
    double soc;
    double v_rc1;
    double v_rc2;
} BatteryState;

typedef struct BatteryDeriv {
    double dsoc;
    double dv_rc1;
    double dv_rc2;
} BatteryDeriv;

// 电池参数（仿真期间保持不变）
typedef struct BatteryParams {
    double capacity;	// 电池容量 Ah
    double v_oc;		// 开路电压 V
    double r_int;		// 欧姆内阻 R0 Ohm

    // 二阶 RC 极化网络参数
    double r1;		// 极化支路1电阻 Ohm
    double c1;		// 极化支路1电容 F
    double r2;		// 极化支路2电阻 Ohm
    double c2;		// 极化支路2电容 F

    double eta_coulomb;	// 库伦效率（放电/充电近似统一）
} BatteryParams;

// 电池观测量（由求解器每步回写的运行时状态）
typedef struct BatteryView {
    double soc;		// 电池荷电状态 0~1
    double v_rc1;	// RC 支路1极化电压 V
    double v_rc2;	// RC 支路2极化电压 V
    double p_load;	// 负载功率 W（正 = 放电，负 = 充电）
    double i_batt;	// 电池电流 A（正 = 放电，负 = 充电）
    double v_batt;	// 电池端电压 V
} BatteryView;

extern BatteryParams battery_params;
extern BatteryView battery_view;

void battery_apply_state(BatteryView* view, double p_load_in, const BatteryEval* eval);

BatteryEval battery_evaluate_from_state(const BatteryParams* params,
                                        double p_load_in,
                                        const BatteryState* st);
BatteryDeriv battery_evaluate_deriv_from_state(const BatteryParams* params,
                                               double i_batt,
                                               const BatteryState* st);

#endif // BATTERY_H