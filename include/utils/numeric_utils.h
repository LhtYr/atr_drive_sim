#ifndef NUMERIC_UTILS_H
#define NUMERIC_UTILS_H

#include <math.h>

// 将一个值限制在一个下限以上
static inline double clamp_floor(double x, double floor_v) { 
    return x < floor_v ? floor_v : x;
}

// 将一个值限制在指定范围内
static inline double clamp_range(double x, double low, double high) {
    if (x < low) {
        return low;
    }
    if (x > high) {
        return high;
    }
    return x;
}

// 返回绝对值并保证不低于下限。
// 注意：结果恒为非负，符号会丢失（仅适用于已知非负或只关心幅值的量，
// 如母线电压幅值）；需要保符号时请勿使用。
static inline double magnitude_floor(double value, double floor_v) {
    const double abs_value = fabs(value);
    return abs_value < floor_v ? floor_v : abs_value;
}

// 2π 与 π 的字面量（工程统一避免依赖 M_PI，-std=c11 下 math.h 不保证暴露）
#define TWO_PI_RAD 6.283185307179586
#define PI_RAD 3.141592653589793

// 角度回绕到 [0, 2π)
static inline double wrap_2pi(double angle) {
    double r = fmod(angle, TWO_PI_RAD);
    return r < 0.0 ? r + TWO_PI_RAD : r;
}

// 角差回绕到 (-π, π]（最短角差路径，用于角度误差与插值）
static inline double wrap_pi(double angle) {
    double r = fmod(angle + PI_RAD, TWO_PI_RAD);
    if (r < 0.0) {
        r += TWO_PI_RAD;
    }
    return r - PI_RAD;
}

#endif // NUMERIC_UTILS_H