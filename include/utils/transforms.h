#ifndef TRANSFORMS_H
#define TRANSFORMS_H

// 三相/dq 坐标系公共常量与变换（幅值不变约定，Clarke 系数 2/3）。
// 全工程的坐标变换只允许通过本文件提供的函数实现，禁止在组件内手写。

// sqrt(3)/2，Clarke/Park 变换系数
#define SQRT3_2 0.8660254037844386
// 1/sqrt(3)，SVPWM 线性区最大电压矢量幅值（相对直流母线电压）
#define INV_SQRT3 0.5773502691896258
// 3/2，幅值不变约定下的功率/转矩系数
#define THREE_2 1.5

// Clarke 变换：abc -> alpha-beta（幅值不变，2/3 系数）
static inline void clarke_abc_to_ab(double a,
                                    double b,
                                    double c,
                                    double* alpha,
                                    double* beta) {
    *alpha = (2.0 / 3.0) * (a - 0.5 * b - 0.5 * c);
    *beta = (2.0 / 3.0) * (SQRT3_2 * b - SQRT3_2 * c);
}

// 正向 Park：alpha-beta -> dq（幅值不变，dq_to_abc 的正变换逆过程）。
// cos_theta/sin_theta 由调用方预先算好，避免同一角度重复求三角函数。
// theta 为 dq 系相对 alpha 轴的电角度（无感控制中传入估计角 θ̂）。
static inline void park_ab_to_dq(double alpha,
                                 double beta,
                                 double cos_theta,
                                 double sin_theta,
                                 double* d,
                                 double* q) {
    *d = alpha * cos_theta + beta * sin_theta;
    *q = -alpha * sin_theta + beta * cos_theta;
}

// 逆 Park + 逆 Clarke：dq -> abc。
// cos_theta/sin_theta 由调用方预先算好，避免同一角度重复求三角函数。
// alpha/beta 中间量也一并输出（调用方可作为观测量记录）。
static inline void dq_to_abc(double d,
                             double q,
                             double cos_theta,
                             double sin_theta,
                             double* alpha,
                             double* beta,
                             double* a,
                             double* b,
                             double* c) {
    *alpha = d * cos_theta - q * sin_theta;
    *beta = d * sin_theta + q * cos_theta;
    *a = *alpha;
    *b = -0.5 * (*alpha) + SQRT3_2 * (*beta);
    *c = -0.5 * (*alpha) - SQRT3_2 * (*beta);
}

#endif // TRANSFORMS_H
