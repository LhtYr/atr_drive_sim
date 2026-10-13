#ifndef LOGGER_H
#define LOGGER_H

#include <stdbool.h>
#include <stdio.h>

typedef struct MotorView MotorView;
typedef struct ShaftView ShaftView;
typedef struct PumpView PumpView;
typedef struct TurbineView TurbineView;
typedef struct InverterView InverterView;
typedef struct BatteryView BatteryView;

// 控制链观测量（不属于任何组件 View，由 main.c 每控制拍回填后供日志读取）
typedef struct ControlLogView {
    double mod_scale;    // 电流环电压矢量缩放（1 = 未饱和，<1 = 触及 SVPWM 线性区上限）
    double iq_ref_cmd;   // q 轴电流参考（外环输出）
    double ud_ref;       // 电流环 d 轴电压参考（饱和后）
    double uq_ref;       // 电流环 q 轴电压参考（饱和后）
    // ---- 无感起动诊断列（仅 --sensorless 选用 sensorless_columns 时输出；
    //      默认 25 列不含这三列，保证有感回归逐字节一致。读取真值仅供诊断，不回馈控制）----
    double theta_err_rad; // wrap(θ̂ − p·theta_mech) rad（位置估计误差，诊断列）
    double omega_est;     // 估计机械角速度 rad/s（= ω̂_e / p）
    double startup_stage; // 起动状态机阶段：0=I/F，1=切换过渡，2=闭环
} ControlLogView;

// 一次快照中所有可记录观测量的集合（列提取函数从中取值）
typedef struct LogViews {
    const MotorView* motor;
    const ShaftView* shaft;
    const PumpView* pump;
    const TurbineView* turbine;
    const InverterView* inverter;
    const BatteryView* battery;
    const ControlLogView* control;
} LogViews;

// 单列描述：表头名 + 取值函数。增删日志列只需改动列定义表。
typedef struct LoggerColumn {
    const char* name;
    double (*extract)(double t, const LogViews* views);
} LoggerColumn;

typedef struct SimulationLogger {
    FILE* stream;
    bool first_column;
    const LoggerColumn* columns;
    int column_count;
} SimulationLogger;

bool SimulationLogger_init(SimulationLogger* logger, const char* path);
// 覆盖默认列定义表（须在 write_header 之前调用）
void SimulationLogger_set_columns(SimulationLogger* logger,
                                  const LoggerColumn* columns,
                                  int column_count);
// 选用无感起动扩展列表（默认 25 列 + theta_err_rad/omega_est/startup_stage = 28 列）。
// 仅 --sensorless 时调用；默认路径仍用 25 列 default_columns，保证有感回归逐字节一致。
// 须在 SimulationLogger_write_header 之前调用。
void SimulationLogger_use_sensorless_columns(SimulationLogger* logger);
void SimulationLogger_write_header(SimulationLogger* logger);
void SimulationLogger_log_snapshot(SimulationLogger* logger,
                                   double t,
                                   const LogViews* views);
void SimulationLogger_close(SimulationLogger* logger);

#endif // LOGGER_H
