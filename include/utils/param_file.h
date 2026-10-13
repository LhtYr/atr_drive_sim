#ifndef PARAM_FILE_H
#define PARAM_FILE_H

// 组件参数文件加载：覆盖编译期内置的组件参数。
//
// 文件格式（每行一条）：
//   component.field = value
// 支持 '#' 行注释与空行；未知键或非法数值将报错并返回 -1。
//
// 示例：
//   # 电机磁链按弱磁设计点标定
//   motor.psi_f = 0.0595
//   battery.capacity = 40
//   battery.soc_init = 0.8        # 电池初始 SOC（观测量初值）
//   inverter.use_switching_model = 0

// 加载参数文件并覆盖组件参数，成功返回 0，失败返回 -1（错误信息已打印到 stderr）
int param_file_apply(const char* path);

// 仿真级参数（sim.* 键绑定，定义于 src/param_file.c）：
//   sim.v_ref -> sim_v_ref，电压环母线电压参考 V（命令行 --vref 优先）
extern double sim_v_ref;

// 无感起动切换转速 shim（sensorless.if.omega_switch_rpm 绑定，定义于 src/param_file.c）：
//   以 rpm 给定 I/F->闭环切换转速，<0 表示参数文件未设置（沿用内置/CLI 值）；
//   main.c 无感初始化时换算为 ifstart.omega_switch_m（rad/s 机械角）。
extern double sensorless_if_omega_switch_rpm;

#endif // PARAM_FILE_H
