#ifndef CHAR_SWEEP_H
#define CHAR_SWEEP_H

#include <stdio.h>

// 部件稳态特性扫描：对边界条件网格逐点调用组件的纯函数评估接口，
// 输出特性数据 CSV。不依赖求解器与控制链。
//
// 参数：
//   f     已打开的输出文件（表头由各扫描函数自行写入）
//   n     每个扫描维度的网格点数
// 返回值：成功 0，失败 -1
typedef int (*CharSweepFn)(FILE* f, int n);

// 生成 [x0, x1] 上 n 个等距点（n >= 2），写入 xs
void char_linspace(double* xs, double x0, double x1, int n);

// 各部件特性扫描
int char_sweep_turbine(FILE* f, int n);
int char_sweep_pump(FILE* f, int n);
int char_sweep_battery(FILE* f, int n);
// 电池动态响应特性（阶跃负载×时间，极化动态；与 char_sweep_battery 准静态对应）
int char_sweep_battery_dynamic(FILE* f, int n);

// 涡轮 7 点试验数据校验（n 不使用；流量误差超阈值返回 -1）
int char_check_turbine(FILE* f, int n);
// 涡轮 7 点动态稳态校验（不绕开容腔动态，积分至稳态后比流量）
int char_check_turbine_dynamic(FILE* f, int n);

#endif // CHAR_SWEEP_H
