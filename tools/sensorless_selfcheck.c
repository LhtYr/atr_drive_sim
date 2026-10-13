// 无感观测器离线自检：在合成 PMSM αβ 数据上验证 STSMO + PLL 收敛性。
//
// 用途：不依赖求解器/全系统仿真，毫秒级验证超螺旋观测器提取反电动势、PLL 锁定
// 转子角的正确性与增益整定，覆盖 I/F 起动段（0→2000 rpm）与闭环高速段（→20000 rpm）。
//
// 原理：给定已知转子角 θ_true(t)（斜坡升速）与真值系电流（id=0, iq=I_f），
// 由 PMSM αβ 电压方程反解出施加电压 u_αβ 与实测电流 i_αβ；把 (u_αβ, i_αβ) 喂给
// STSMO→PLL（观测器不知 θ_true），比对估计角/角速度与真值。
//
// 编译运行（工程根目录，一条命令）：
//   gcc -std=c11 -Wall -Wextra -O2 -Iinclude -Iinclude/component -Iinclude/control
//       -Iinclude/solver tools/sensorless_selfcheck.c model/control/STSMO.c
//       model/control/PLL.c -lm -o build/sl_selfcheck && ./build/sl_selfcheck
#include <math.h>
#include <stdio.h>

#include "STSMO.h"
#include "PLL.h"
#include "component/Motor.h"
#include "utils/numeric_utils.h"
#include "utils/transforms.h"

int main(int argc, char** argv) {
    // 电机参数（同 Motor.c 默认）
    const MotorParams mp = {.Rs = 0.0051, .Ld = 34.6e-6, .Lq = 34.6e-6, .psi_f = 0.0595, .p = 2};
    const double L = mp.Ld;
    const double Ts = 5e-6;

    // 可选：argv[1] 给出路径时导出逐步轨迹 CSV（诊断纹波/收敛用）
    FILE* trace = (argc > 1) ? fopen(argv[1], "w") : NULL;
    if (trace) {
        fprintf(trace, "t,theta_true,theta_est,theta_err_deg,wm_true,wm_est,wm_err_pct,emf_mag,eps\n");
    }

    stsmo.Ts = Ts;
    pll.Ts = Ts;
    pll.p = mp.p;
    pll.emf_lpf_wc = TWO_PI_RAD * stsmo.emf_lpf_hz;   // 相位滞后补偿与 STSMO 低通同步

    // 转速剖面：0→2000 rpm（0.3 s，I/F 段），保持 0.25 s，再 2000→20000 rpm（0.5 s，闭环段），保持 0.25 s
    const double rpm2rad = 2.0 * PI_RAD / 60.0;
    const double t1 = 0.3, t2 = 0.55, t3 = 1.05;
    const double w_lo = 2000.0 * rpm2rad, w_hi = 20000.0 * rpm2rad;
    const double I_f = 100.0;

    stsmo_reset(&stsmo);
    pll_reset(&pll, 0.0);

    double theta_true_e = 0.0;   // 真值电角度
    double max_err_low = 0.0;    // I/F 段 2000 rpm 稳态保持窗 [0.45,0.55] 最大角误差 rad
    double max_err_hi = 0.0;     // 高速段 20000 rpm 稳态保持窗 [1.15,1.30] 最大角误差 rad
    double w_err_hi = 0.0;       // 高速段稳态机械角速度估计误差 rad/s
    double peak_err_all = 0.0;   // 全程峰值角误差（含升速瞬态，参考）
    int nan_seen = 0;

    const int n_steps = (int)(1.3 / Ts);
    for (int k = 0; k < n_steps; ++k) {
        const double t = k * Ts;
        // 机械角速度剖面
        double w_mech;
        if (t < t1) {
            w_mech = w_lo * (t / t1);
        } else if (t < t2) {
            w_mech = w_lo;
        } else if (t < t3) {
            w_mech = w_lo + (w_hi - w_lo) * ((t - t2) / (t3 - t2));
        } else {
            w_mech = w_hi;
        }
        const double w_e = mp.p * w_mech;          // 电角速度
        theta_true_e = wrap_2pi(theta_true_e + w_e * Ts);

        const double ct = cos(theta_true_e), st = sin(theta_true_e);
        // 真值系电流 id=0, iq=I_f → αβ（逆 Park）：i_α=−I_f·sinθ, i_β=I_f·cosθ
        const double i_alpha = -I_f * st;
        const double i_beta = I_f * ct;
        // di_αβ/dt（解析）
        const double di_alpha = -I_f * ct * w_e;
        const double di_beta = -I_f * st * w_e;
        // 反电动势 e_α=−ω_e·ψ_f·sinθ, e_β=ω_e·ψ_f·cosθ
        const double emf_a_true = -w_e * mp.psi_f * st;
        const double emf_b_true = w_e * mp.psi_f * ct;
        // αβ 电压方程反解施加电压
        const double u_alpha = mp.Rs * i_alpha + L * di_alpha + emf_a_true;
        const double u_beta = mp.Rs * i_beta + L * di_beta + emf_b_true;

        // 观测器 + PLL（只吃 u_αβ 与 i_αβ，不知 θ_true）
        const STSMOOutput so = stsmo_evaluate(&stsmo, &mp, u_alpha, u_beta, i_alpha, i_beta);
        const PLLOutput po = pll_evaluate(&pll, so.emf_alpha, so.emf_beta);

        if (!isfinite(po.theta_e) || !isfinite(po.omega_e) || !isfinite(so.emf_alpha)) {
            nan_seen = 1;
        }
        const double err = fabs(wrap_pi(po.theta_e - theta_true_e));
        if (err > peak_err_all) {
            peak_err_all = err;
        }
        // I/F 段稳态保持窗 [0.45,0.55]（2000 rpm 升速+再整定已完成，PLL 应已锁定）
        if (t >= 0.45 && t <= 0.55 && err > max_err_low) {
            max_err_low = err;
        }
        // 高速段稳态保持窗 [1.15,1.30]（20000 rpm 升速已完成）
        if (t >= 1.15) {
            if (err > max_err_hi) {
                max_err_hi = err;
            }
            const double we = fabs(po.omega_m - w_mech);
            if (we > w_err_hi) {
                w_err_hi = we;
            }
        }
        if (trace) {
            fprintf(trace, "%.7f,%.6f,%.6f,%.5f,%.4f,%.4f,%.5f,%.5f,%.6f\n",
                    t, theta_true_e, po.theta_e, err * 180.0 / PI_RAD,
                    w_mech, po.omega_m, 100.0 * fabs(po.omega_m - w_mech) / w_hi,
                    so.emf_mag, po.eps);
        }
    }
    if (trace) {
        fclose(trace);
    }

    const double rad2deg = 180.0 / PI_RAD;
    printf("=== 无感观测器离线自检（STSMO + PLL） ===\n");
    printf("增益：stsmo.k1=%.3g k2=%.3g lpf=%.0fHz | pll.kp=%.3g ki=%.3g\n",
           stsmo.k1, stsmo.k2, stsmo.emf_lpf_hz, pll.kp, pll.ki);
    printf("I/F 段(2000 rpm) 最大角误差   : %.3f deg\n", max_err_low * rad2deg);
    printf("高速段(20000 rpm) 最大角误差  : %.3f deg\n", max_err_hi * rad2deg);
    printf("全程峰值角误差(含升速瞬态,参考): %.3f deg\n", peak_err_all * rad2deg);
    printf("高速段 最大机械角速度估计误差 : %.3f rad/s (%.3f%%)\n",
           w_err_hi, 100.0 * w_err_hi / w_hi);
    printf("nan/inf: %s\n", nan_seen ? "YES(FAIL)" : "no");

    // 判据：I/F 段 <10°，高速段 <5°，速度误差 <1%，无 nan
    const int ok = !nan_seen && max_err_low < 10.0 / rad2deg &&
                   max_err_hi < 5.0 / rad2deg && w_err_hi < 0.01 * w_hi;
    printf("结果：%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
