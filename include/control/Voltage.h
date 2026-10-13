#ifndef VOLTAGE_H
#define VOLTAGE_H

typedef struct VoltageOutput {
    double ev;       // 电压误差 v_meas - v_ref（稳定方向：发电抬升母线）
    double iq_ref;   // q 轴电流参考（<0 发电，>0 电动）
} VoltageOutput;

typedef struct Voltage {
    double Ts;
    double kp;
    double ki;
    double iq_ref_min;
    double iq_ref_max;
    double anti_windup_gain;
    double int_v;
} Voltage;

extern Voltage voltage;

VoltageOutput voltage_evaluate(Voltage* controller,
                                double v_ref,
                                double v_meas);

#endif // VOLTAGE_H
