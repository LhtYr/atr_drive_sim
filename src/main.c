#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#ifdef _WIN32
#include <direct.h>
#define ATR_MKDIR(path, mode) _mkdir(path)
#else
#define ATR_MKDIR(path, mode) mkdir(path, mode)
#endif

#include "Logger.h"
#include "Current.h"
#include "Battery.h"
#include "Inverter.h"
#include "Solver.h"
#include "BoundaryProfile.h"
#include "Motor.h"
#include "Shaft.h"
#include "Pump.h"
#include "SVPWM.h"
#include "Speed.h"
#include "Voltage.h"
#include "Turbine.h"
#include "STSMO.h"
#include "PLL.h"
#include "IFStart.h"
#include "sensorless.h"
#include "utils/transforms.h"
#include "utils/numeric_utils.h"
#include "utils/param_file.h"

// 控制模式：speed = 速度环给出 iq 指令（默认）；voltage = 电压环为唯一 iq 指令源（发电稳压）
typedef enum ControlMode {
	CONTROL_MODE_SPEED = 0,
	CONTROL_MODE_VOLTAGE = 1,
} ControlMode;

// 仿真配置（命令行可覆盖，默认值与原硬编码一致）
typedef struct SimConfig {
	double duration_s;			// 仿真时长 s
	double model_dt;			// 模型积分步长 s
	double control_dt;			// 控制器更新周期 s
	double omega_ref_mech;		// 目标机械角速度 rad/s
	double v_ref;				// 电压环母线电压参考 V（--control-mode voltage 时生效）
	ControlMode control_mode;	// 外环模式（默认速度环）
	int sensorless;				// 无感起动开关（默认 0=有感，全路径与改动前一致；1=I/F+STSMO+PLL）
	int log_stride;				// 日志记录步长（每隔多少步记录一次）
	const char* output_path;	// CSV 输出路径
	const char* config_path;	// 组件参数文件路径（NULL = 仅用内置参数）
} SimConfig;

// --vref 是否由命令行显式指定（是则优先于参数文件 sim.v_ref）
static int g_vref_from_cli = 0;

static void print_usage(const char* prog) {
	printf("Usage: %s [options]\n"
	       "  --duration <s>       仿真时长（默认 20）\n"
	       "  --model-dt <s>       模型积分步长（默认 5e-7）\n"
	       "  --control-dt <s>     控制器更新周期（默认 5e-6）\n"
	       "  --omega-ref-rpm <n>  目标转速 rpm（默认 20000）\n"
	       "  --control-mode <m>   外环模式：speed（默认）或 voltage（电压环，速度环旁路）\n"
	       "  --vref <V>           电压环母线电压参考 V（默认 540，仅 voltage 模式生效）\n"
	       "  --sensorless         启用无感起动（I/F+STSMO+PLL 切换到速度闭环，默认关闭）\n"
	       "  --if-current <A>     I/F 起动电流幅值 A（默认 60，仅 --sensorless 生效）\n"
	       "  --switch-rpm <n>     I/F->闭环切换转速 rpm 机械角（默认 2000，仅 --sensorless 生效）\n"
	       "  --log-stride <n>     日志记录步长（默认 1000）\n"
	       "  --output <path>      CSV 输出路径（默认 results/simulation.csv）\n"
	       "  --config <path>      组件参数文件（component.field = value，见 utils/param_file.h）\n"
	       "  -h, --help           显示帮助\n",
	       prog);
}

// 解析命令行参数，成功返回 0，失败返回 -1（错误信息已打印）
static int parse_args(int argc, char** argv, SimConfig* cfg) {
	for (int i = 1; i < argc; ++i) {
		const char* opt = argv[i];
		if (strcmp(opt, "-h") == 0 || strcmp(opt, "--help") == 0) {
			print_usage(argv[0]);
			exit(0);
		}
		if (strcmp(opt, "--sensorless") == 0) {	// 布尔开关，无参数值
			cfg->sensorless = 1;
			continue;
		}
		const int known = strcmp(opt, "--duration") == 0 || strcmp(opt, "--model-dt") == 0 ||
		                  strcmp(opt, "--control-dt") == 0 || strcmp(opt, "--omega-ref-rpm") == 0 ||
		                  strcmp(opt, "--control-mode") == 0 || strcmp(opt, "--vref") == 0 ||
		                  strcmp(opt, "--if-current") == 0 || strcmp(opt, "--switch-rpm") == 0 ||
		                  strcmp(opt, "--log-stride") == 0 || strcmp(opt, "--output") == 0 ||
		                  strcmp(opt, "--config") == 0;
		if (!known) {
			fprintf(stderr, "未知选项：%s\n", opt);
			return -1;
		}
		if (i + 1 >= argc) {
			fprintf(stderr, "选项 %s 缺少参数值\n", opt);
			return -1;
		}
		const char* value = argv[++i];
		if (strcmp(opt, "--output") == 0) {
			cfg->output_path = value;
			continue;
		}
		if (strcmp(opt, "--config") == 0) {
			cfg->config_path = value;
			continue;
		}
		if (strcmp(opt, "--control-mode") == 0) {
			if (strcmp(value, "speed") == 0) {
				cfg->control_mode = CONTROL_MODE_SPEED;
			} else if (strcmp(value, "voltage") == 0) {
				cfg->control_mode = CONTROL_MODE_VOLTAGE;
			} else {
				fprintf(stderr, "--control-mode 仅支持 speed 或 voltage：%s\n", value);
				return -1;
			}
			continue;
		}
		char* end = NULL;
		const double num = strtod(value, &end);
		if (end == value || *end != '\0') {
			fprintf(stderr, "选项 %s 的参数无效：%s\n", opt, value);
			return -1;
		}
		if (strcmp(opt, "--duration") == 0) {
			cfg->duration_s = num;
		} else if (strcmp(opt, "--model-dt") == 0) {
			cfg->model_dt = num;
		} else if (strcmp(opt, "--control-dt") == 0) {
			cfg->control_dt = num;
		} else if (strcmp(opt, "--omega-ref-rpm") == 0) {
			cfg->omega_ref_mech = num * 0.10471975511965977;	// rpm -> rad/s（π/30）
		} else if (strcmp(opt, "--vref") == 0) {
			cfg->v_ref = num;
			g_vref_from_cli = 1;
		} else if (strcmp(opt, "--if-current") == 0) {
			ifstart.i_start = num;			// I/F 起动电流（参数文件 sensorless.if.i_start 可覆盖）
		} else if (strcmp(opt, "--switch-rpm") == 0) {
			ifstart.omega_switch_m = num * 0.10471975511965977;	// rpm -> rad/s 机械角
		} else if (strcmp(opt, "--log-stride") == 0) {
			cfg->log_stride = (int)num;
		} else {
			fprintf(stderr, "内部错误：未处理的选项 %s\n", opt);
			return -1;
		}
	}
	if (cfg->duration_s <= 0.0 || cfg->model_dt <= 0.0 || cfg->control_dt <= 0.0 ||
	    cfg->v_ref <= 0.0 || cfg->log_stride < 1) {
		fprintf(stderr, "duration/model-dt/control-dt/vref 必须为正，log-stride 必须 >= 1\n");
		return -1;
	}
	return 0;
}

int main(int argc, char** argv) {
	// 简单闭环示例：速度环 -> 电流环 -> step（SVPWM/逆变器由求解器内部处理）
	SimConfig cfg = {
		.duration_s = 20.0,
		.model_dt = 5e-7,				// 模型求解时间步长 s
		.control_dt = 5e-6,			// 控制器更新时间步长 s
		.omega_ref_mech = 2094.39510239,	// 20000 rpm
		.v_ref = 540.0,				// 电压环母线电压参考 V（电池箝位附近）
		.control_mode = CONTROL_MODE_SPEED,	// 默认速度环
		.sensorless = 0,			// 默认有感（关闭时全路径与改动前完全一致）
		.log_stride = 1000,
		.output_path = "results/simulation.csv",
		.config_path = NULL,
	};
	if (parse_args(argc, argv, &cfg) != 0) {
		print_usage(argv[0]);
		return 1;
	}
	if (cfg.config_path != NULL && param_file_apply(cfg.config_path) != 0) {
		return 1;
	}
	if (!g_vref_from_cli) {
		cfg.v_ref = sim_v_ref;	// 参数文件 sim.v_ref 覆盖默认值（命令行 --vref 优先）
	}

	// 确保输出文件所在目录存在（仅创建最后一级目录，忽略已存在的情况）
	{
		char dir_buf[512];
		const char* slash = strrchr(cfg.output_path, '/');
		if (slash != NULL && slash != cfg.output_path && (size_t)(slash - cfg.output_path) < sizeof(dir_buf)) {
			memcpy(dir_buf, cfg.output_path, (size_t)(slash - cfg.output_path));
			dir_buf[slash - cfg.output_path] = '\0';
			ATR_MKDIR(dir_buf, 0755);
		}
	}

	// 先写临时文件，仿真完成后原子改名，避免读者（如 notebook）在仿真进行中读到半截 CSV
	char tmp_path[576];
	const int tmp_len = snprintf(tmp_path, sizeof(tmp_path), "%s.tmp", cfg.output_path);
	if (tmp_len < 0 || (size_t)tmp_len >= sizeof(tmp_path)) {
		fprintf(stderr, "output path too long: %s\n", cfg.output_path);
		return 1;
	}

	SimulationLogger logger;
	if (!SimulationLogger_init(&logger, tmp_path)) {
		fprintf(stderr, "failed to open log file: %s\n", tmp_path);
		return 1;
	}
	// 无感模式扩展日志列（28 列：+theta_err_rad/omega_est/startup_stage）；
	// 默认有感仍为 25 列，保证回归逐字节一致。须在 write_header 之前选用。
	if (cfg.sensorless) {
		SimulationLogger_use_sensorless_columns(&logger);
	}
	SimulationLogger_write_header(&logger);

	// 仿真配置（来自命令行，见上方默认值）
	const double model_dt = cfg.model_dt;
	const double control_dt = cfg.control_dt;
	const int total_steps = (int)(cfg.duration_s / model_dt);	// 总仿真步数
	const int log_stride = cfg.log_stride;

	// 进度条配置
	const int progress_width = 30;		// 进度条宽度（字符数）
	const int progress_stride = (total_steps >= 100) ? (total_steps / 100) : 1;
	const char spinner[] = {'|', '/', '-', '\\'};

	// 控制周期与模型周期必须为整数倍关系（用于整数节拍同步）
	const int ctrl_ratio = (int)llround(control_dt / model_dt);
	if (fabs(ctrl_ratio * model_dt - control_dt) > 1e-15) {
		fprintf(stderr, "control_dt must be an integer multiple of model_dt\n");
		return 1;
	}

	// 求解器初始化：积分步长逐拍显式传入，这里只需用组件初值填充状态向量
	solver_init_state();
	speed.Ts = control_dt;
	current.Ts = control_dt;
	voltage.Ts = control_dt;
	svpwm.linear_limit = INV_SQRT3;	// 1/sqrt(3)，SVPWM 线性区最大电压矢量幅值（相对于 DC 母线电压的比值）

	// 无感起动初始化：复位状态机/观测器，绑定控制周期与极对数，
	// 并同步 PLL 的反电动势低通角频率（相位滞后补偿须与 stsmo.emf_lpf_hz 一致）
	if (cfg.sensorless) {
		// 参数文件 sensorless.if.omega_switch_rpm（rpm）覆盖切换转速：换算为 rad/s 机械角
		// （命令行 --switch-rpm 已在 parse_args 写入 ifstart.omega_switch_m，配置优先于 CLI）
		if (sensorless_if_omega_switch_rpm > 0.0) {
			ifstart.omega_switch_m = sensorless_if_omega_switch_rpm * 0.10471975511965977;
		}
		sensorless_reset(&sensorless_mgr);
		sensorless_mgr.Ts = control_dt;
		ifstart.Ts = control_dt;
		ifstart.p = motor_params.p;
		stsmo.Ts = control_dt;
		pll.Ts = control_dt;
		pll.p = motor_params.p;
		pll.emf_lpf_wc = TWO_PI_RAD * stsmo.emf_lpf_hz;
		stsmo_reset(&stsmo);
		pll_reset(&pll, 0.0);
		ifstart_reset(&ifstart);
	}

	// 控制目标
	const double omega_ref_mech = cfg.omega_ref_mech;
	// （id=0 时反电动势约 25000 rpm 到达母线电压极限，默认目标转速在线性调制区内）
	double id_ref_cmd = 0.0;
	double iq_ref_cmd = 50.0;

	// 运行时状态
	double ud_cmd = 0.0;			// 电流环输出：d轴电压参考
	double uq_cmd = 0.0;			// 电压环控制输出：q轴电压参考
	int last_percent = -1;			// 进度显示上次更新的百分比
	ControlLogView control_log = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};	// 控制链观测量（每控制拍回填，供日志读取；后 3 项为无感诊断列）

	// SVPWM 逆 Park 电压合成角通道：有感默认 use_ctrl_theta=0（用真实电角，行为与改动前一致）；
	// 无感时每控制拍置 1 并注入冻结的估计角 θ_ctrl（数字控制一拍一更新，整个控制周期内保持不变）
	int use_ctrl_theta = 0;
	double theta_ctrl_e = 0.0;

	// 无感起动诊断跟踪（仅日志/验收用，读取真值 theta_mech/omega_mech/iq 不参与控制）
	double sl_switch_t = -1.0;			// 阶段 0->1（开始过渡）的时刻 s
	double sl_switch_theta_err_deg = 0.0;		// 切换时刻 |theta_err| 度
	double sl_peak_err_blend_deg = 0.0;		// 过渡阶段 max|theta_err| 度
	double sl_peak_iq = 0.0;			// 起动全程（IF+过渡）max|iq| A
	double sl_closed_t = -1.0;			// 进入阶段 2 的时刻 s
	double sl_steady_err_deg = 0.0;			// 稳态段（末 1 s）max|theta_err| 度
	double sl_steady_werr_pct = 0.0;		// 稳态段 max|ω̂−ω_true|/ω_true %
	int sl_prev_stage = 0;
	int sl_bad = 0;					// nan/inf 标志

	// 任务剖面（阶段 3 发电稳压）诊断跟踪：切换时刻/转速、bumpless 交接跳变、稳压精度、
	// 观测器锁定、电压边界裕度。读取真值仅供验收/日志，绝不回馈控制。
	double sl_switch_omega = 0.0;			// 阶段 0->1 时机械角速度 rad/s
	double sl_closed_omega = 0.0;			// 阶段 1->2 时机械角速度 rad/s
	double iq_ref_prev = 0.0;			// 上一控制拍 iq_ref 指令 A（测交接跳变用）
	double sl_gen_t = -1.0;				// 阶段 2->3（切入发电）的时刻 s
	double sl_gen_omega = 0.0;			// 切入发电时机械角速度 rad/s
	double sl_gen_soc = 0.0;			// 切入发电时 SOC（对比终态验证发电段充电上升）
	double sl_gen_iq_jump = 0.0;			// 交接瞬间 max|Δiq_ref| A（bumpless 连续性）
	double sl_gen_vdc_transient = 0.0;		// 发电段瞬态 max|v_dc−v_ref| V（切入后 0.5 s 内）
	double sl_gen_vdc_steady = 0.0;			// 发电段稳态 max|v_dc−v_ref| V（末 1 s）
	double sl_gen_theta_err_max = 0.0;		// 发电段全程 max|theta_err| deg
	double sl_gen_theta_err_steady = 0.0;		// 发电段稳态 max|theta_err| deg（末 1 s）
	double sl_gen_mod_min = 1.0;			// 发电段 min mod_scale（=1 全程未触电压边界）
	const double soc_init_snapshot = battery_view.soc;	// 起点 SOC（param_file 已覆盖，摘要对比用）

	// 外部边界条件剖面选择：
	//   · 电压模式（gen_voltage_reg 能力边界测试）沿用供气斜坡（ramp_t_start 默认 0，t=0 起坡）；
	//   · 任务剖面（--sensorless + 延迟供气斜坡）经 boundary.ramp_t_start > 0 显式启用延迟斜坡，
	//     使起动段（t < ramp_t_start）完全无供气，涡轮供气渐增在阶段 2 之后才介入；
	//   · 其余（默认速度环、sensorless_start 无供气空载）用常量剖面，行为与改动前逐拍一致。
	// 默认值在 Solver.c 的 boundary_defaults，可经 --config 的 boundary.* 键覆盖；
	// 斜坡参数见 boundary_ramp_config（boundary.m_dot_start / ramp_time / ramp_t_start）。
	const int use_ramp_profile = (cfg.control_mode == CONTROL_MODE_VOLTAGE) ||
	                             (boundary_ramp_config.ramp_t_start > 0.0);
	const BoundaryProfile boundary_profile = use_ramp_profile
		? boundary_profile_make_m_dot_ramp(&boundary_defaults)
		: boundary_profile_make_constant(&boundary_defaults);

	// 能力边界实测跟踪（仅电压模式有效）
	double sat_t = -1.0;				// a) mod_scale 首次饱和（<1）的时刻 s
	double sat_omega = 0.0;				//    对应机械角速度 rad/s
	double max_iq_dev = 0.0;			// b) 饱和后 |iq - iq_ref| 最大偏差 A
	double vdev_t = -1.0;				// c) |v_dc - v_ref| 首次超 2 V 的时刻 s
	double vdev_omega = 0.0;			//    对应机械角速度 rad/s
	double vdev_v = 0.0;				//    对应 v_dc 实测值 V

	printf("Simulation beginning...\n");

	for (int k = 0; k < total_steps; ++k) {
		const double t = k * model_dt;
		if ((k % ctrl_ratio) == 0) {
			if (cfg.sensorless) {
				// 无感起动链路：仅用实测 ia/ib/ic + v_dc + 自身电压指令/估计角（禁止读取任何真值）。
				// I/F + STSMO + PLL + 切换状态机细节见 model/control/sensorless.c。
				const SensorlessOutput sl = sensorless_evaluate(&sensorless_mgr,
						&motor_params,
						motor_view.ia, motor_view.ib, motor_view.ic,
						inverter_view.v_dc, svpwm.linear_limit,
						omega_ref_mech, cfg.v_ref);
				ud_cmd = sl.ud_cmd;
				uq_cmd = sl.uq_cmd;
				id_ref_cmd = sl.id_ref_cmd;
				iq_ref_cmd = sl.iq_ref_cmd;
				use_ctrl_theta = 1;			// SVPWM 改用控制侧注入的冻结估计角
				theta_ctrl_e = sl.theta_ctrl_e;
				control_log.mod_scale = sl.mod_scale;
				control_log.iq_ref_cmd = iq_ref_cmd;
				control_log.ud_ref = ud_cmd;
				control_log.uq_ref = uq_cmd;

				// 诊断（读取真值仅供验收/日志，绝不回馈控制）：theta_err = wrap(θ̂ − p·theta_mech)
				const double theta_err_rad = wrap_pi(sl.theta_est_e -
						motor_params.p * shaft_view.theta_mech);
				const double theta_err_deg = fabs(theta_err_rad) * 57.29577951308232;
				// 回填无感诊断日志列（仅 sensorless_columns 选用时输出，默认 25 列不含）
				control_log.theta_err_rad = theta_err_rad;
				control_log.omega_est = sl.omega_est_m;
				control_log.startup_stage = (double)sl.stage;
				if (sl.stage != sl_prev_stage) {
					if (sl.stage == SL_STAGE_BLEND) {		// 0->1：PLL 锁定、开始过渡
						sl_switch_t = t;
						sl_switch_theta_err_deg = theta_err_deg;
						sl_switch_omega = shaft_view.omega_mech;
					} else if (sl.stage == SL_STAGE_CLOSED) {	// 1->2：过渡完成、进入闭环
						sl_closed_t = t;
						sl_closed_omega = shaft_view.omega_mech;
					} else if (sl.stage == SL_STAGE_GEN) {	// 2->3：扭矩需求消失、切入发电
						sl_gen_t = t;
						sl_gen_omega = shaft_view.omega_mech;
						sl_gen_soc = battery_view.soc;
					}
					sl_prev_stage = sl.stage;
				}
				if (sl.stage == SL_STAGE_BLEND && theta_err_deg > sl_peak_err_blend_deg) {
					sl_peak_err_blend_deg = theta_err_deg;
				}
				// 切换过程 iq 冲击：覆盖 IF + 过渡 + 进入闭环后 50 ms
				if (sl.stage != SL_STAGE_CLOSED || (sl_closed_t >= 0.0 && t < sl_closed_t + 0.05)) {
					const double iq_abs = fabs(motor_view.iq);
					if (iq_abs > sl_peak_iq) {
						sl_peak_iq = iq_abs;
					}
				}
				if (t >= cfg.duration_s - 1.0) {		// 稳态段（末 1 s）
					if (theta_err_deg > sl_steady_err_deg) {
						sl_steady_err_deg = theta_err_deg;
					}
					if (fabs(shaft_view.omega_mech) > 1.0) {
						const double werr_pct = fabs(sl.omega_est_m - shaft_view.omega_mech) /
								fabs(shaft_view.omega_mech) * 100.0;
						if (werr_pct > sl_steady_werr_pct) {
							sl_steady_werr_pct = werr_pct;
						}
					}
				}
				// 任务剖面阶段 3（发电）诊断：bumpless 交接跳变、v_dc 稳压精度、theta_err 锁定、
				// 电压边界裕度。仅在切入发电后累计，sensorless_start（不进阶段 3）不受影响。
				if (sl.stage == SL_STAGE_GEN) {
					if (t <= sl_gen_t + 0.002) {		// 交接窗（切入后 2 ms）相邻控制拍 iq_ref 跳变
						const double diq = fabs(iq_ref_cmd - iq_ref_prev);
						if (diq > sl_gen_iq_jump) {
							sl_gen_iq_jump = diq;
						}
					}
					const double vdev = fabs(inverter_view.v_dc - cfg.v_ref);
					if (t < sl_gen_t + 0.5 && vdev > sl_gen_vdc_transient) {
						sl_gen_vdc_transient = vdev;	// 切入瞬态窗（0.5 s）
					}
					if (theta_err_deg > sl_gen_theta_err_max) {
						sl_gen_theta_err_max = theta_err_deg;
					}
					if (sl.mod_scale < sl_gen_mod_min) {
						sl_gen_mod_min = sl.mod_scale;	// =1 全程未触电压边界
					}
					if (t >= cfg.duration_s - 1.0) {	// 稳态窗（末 1 s）
						if (vdev > sl_gen_vdc_steady) {
							sl_gen_vdc_steady = vdev;
						}
						if (theta_err_deg > sl_gen_theta_err_steady) {
							sl_gen_theta_err_steady = theta_err_deg;
						}
					}
				}
				iq_ref_prev = iq_ref_cmd;			// 每控制拍更新（交接跳变基准）
				if (!isfinite(ud_cmd) || !isfinite(uq_cmd) || !isfinite(theta_ctrl_e) ||
				    !isfinite(sl.omega_est_m) || !isfinite(sl.theta_est_e) ||
				    !isfinite(iq_ref_cmd)) {
					sl_bad = 1;
				}
			} else {
				if (cfg.control_mode == CONTROL_MODE_VOLTAGE) {
					// 电压模式：电压环为唯一 iq 指令源（发电稳压），速度环旁路；
					// 发电（iq<0）把机械功率注入母线抬升 v_dc，稳定方向为 iq_ref ∝ (v_meas - v_ref)
					const VoltageOutput voltage_out = voltage_evaluate(&voltage,
							cfg.v_ref,
							inverter_view.v_dc);
					iq_ref_cmd = voltage_out.iq_ref;
				} else {
					const SpeedOutput speed_out = speed_evaluate(&speed,
							omega_ref_mech,
							shaft_view.omega_mech);
					iq_ref_cmd = speed_out.iq_ref;	// 速度环输出作为 q 轴电流参考
				}
				id_ref_cmd = 0.0;				// d 轴电流参考保持 0（无弱磁）

				const CurrentOutput current_out = current_evaluate(&current,
						&motor_params,
						id_ref_cmd,
						iq_ref_cmd,
						motor_view.id,
						motor_view.iq,
						motor_view.omega_e,
						inverter_view.v_dc,
						svpwm.linear_limit);
				ud_cmd = current_out.ud_ref;
				uq_cmd = current_out.uq_ref;
				control_log.mod_scale = current_out.modulation_scale;
				control_log.iq_ref_cmd = iq_ref_cmd;
				control_log.ud_ref = current_out.ud_ref;
				control_log.uq_ref = current_out.uq_ref;

				if (cfg.control_mode == CONTROL_MODE_VOLTAGE) {
					// mod_scale = 1 为未饱和；首次 <1 即电压矢量触及 v_dc/√3 上限（能力边界）
					if (sat_t < 0.0 && current_out.modulation_scale < 1.0) {
						sat_t = t;
						sat_omega = shaft_view.omega_mech;
					}
					if (sat_t >= 0.0) {
						const double dev = fabs(motor_view.iq - iq_ref_cmd);
						if (dev > max_iq_dev) {
							max_iq_dev = dev;
						}
					}
					// 从边界点后起判，避免起机电气瞬态被误记为能力边界事件
					if (sat_t >= 0.0 && vdev_t < 0.0 && fabs(inverter_view.v_dc - cfg.v_ref) > 2.0) {
						vdev_t = t;
						vdev_omega = shaft_view.omega_mech;
						vdev_v = inverter_view.v_dc;
					}
				}
			}
		}

		SystemInput input = boundary_profile.evaluate(t, boundary_profile.ctx);
		input.ud_ref = ud_cmd;
		input.uq_ref = uq_cmd;
		input.use_ctrl_theta = use_ctrl_theta;	// 有感=0（真实电角），无感=1（下方冻结估计角）
		input.theta_ctrl_e = theta_ctrl_e;		// 无感：整个控制周期内保持不变（ZOH）

		step(t, model_dt, input);
		if ((k % log_stride) == 0) {
			const LogViews views = {
				.motor = &motor_view,
				.shaft = &shaft_view,
				.pump = &pump_view,
				.turbine = &turbine_view,
				.inverter = &inverter_view,
				.battery = &battery_view,
				.control = &control_log,
			};
			SimulationLogger_log_snapshot(&logger, t, &views);
		}

		// 进度显示，要极致性能可注释掉下面这段代码
		if ((k % progress_stride) == 0 || k == total_steps - 1) {
			const int percent = (int)(100.0 * (k + 1) / total_steps);
			if (percent != last_percent) {
				const int filled = (percent * progress_width) / 100;
				char bar[31];
				for (int i = 0; i < progress_width; ++i) {
					bar[i] = '.';
				}
				for (int i = 0; i < filled && i < progress_width; ++i) {
					bar[i] = '=';
				}
				if (percent < 100 && filled < progress_width) {
					bar[filled] = '>';
				}
				bar[progress_width] = '\0';
				printf("\r%c [%s] %3d%%  %5.2fs/%5.2fs",
				       spinner[(percent / 2) % 4],
				       bar,
				       percent,
				       t,
				       cfg.duration_s);
				fflush(stdout);	// 刷新输出缓冲区，确保进度信息及时显示
				last_percent = percent;
			}
		}// 进度显示结束

	}
	printf("\n");

	SimulationLogger_close(&logger);

	// 写入完毕，原子替换目标文件（POSIX 下 rename 为原子操作）
	if (rename(tmp_path, cfg.output_path) != 0) {
		fprintf(stderr, "failed to rename %s -> %s\n", tmp_path, cfg.output_path);
		return 1;
	}

	printf("Simulation completed\n");
	printf("Results logged to %s\n", cfg.output_path);

	if (cfg.control_mode == CONTROL_MODE_VOLTAGE) {
		const double rpm_per_rad_s = 9.549296585513721;	// rad/s -> rpm（30/π）
		printf("\n=== 发电稳压能力边界实测（电压环、无弱磁） ===\n");
		if (sat_t >= 0.0) {
			printf("a) mod_scale 首次饱和（电压矢量触及 v_dc/√3）：t = %.4f s，转速 = %.1f rpm\n",
			       sat_t, sat_omega * rpm_per_rad_s);
			printf("b) 该时刻之后 |iq - iq_ref| 最大偏差：%.2f A\n", max_iq_dev);
		} else {
			printf("a) mod_scale 全程未饱和（未达到能力边界）\n");
		}
		if (vdev_t >= 0.0) {
			printf("c) |v_dc - v_ref| 首次超过 2 V：t = %.4f s，转速 = %.1f rpm，v_dc = %.2f V\n",
			       vdev_t, vdev_omega * rpm_per_rad_s, vdev_v);
		} else {
			printf("c) |v_dc - v_ref| 全程未超过 2 V\n");
		}
	}

	if (cfg.sensorless) {
		const double rpm_per_rad_s = 9.549296585513721;	// rad/s -> rpm（30/π）
		const double iq_2x_steady = 112.6;			// 2×稳态 iq（20000 rpm 空载 ≈ 56.3 A）
		printf("\n=== 无感起动验收实测（I/F + STSMO + PLL -> 速度闭环） ===\n");
		if (sl_switch_t >= 0.0) {
			printf("切换（IF->过渡）：t = %.4f s，|theta_err| = %.3f deg（判据 <10 deg）\n",
			       sl_switch_t, sl_switch_theta_err_deg);
		} else {
			printf("切换（IF->过渡）：全程未触发（PLL 未锁定，检查判据/整定）\n");
		}
		if (sl_closed_t >= 0.0) {
			printf("进入闭环：t = %.4f s\n", sl_closed_t);
		}
		printf("过渡段 max|theta_err|：%.3f deg\n", sl_peak_err_blend_deg);
		printf("切换过程 max|iq|：%.2f A（判据 <2x稳态 %.1f A）\n", sl_peak_iq, iq_2x_steady);
		printf("稳态段（末 1 s）max|theta_err|：%.3f deg（判据 <5 deg）\n", sl_steady_err_deg);
		printf("稳态段 max|w_est-w_true|/w_true：%.3f %%（判据 <1 %%）\n", sl_steady_werr_pct);
		printf("末端转速：%.1f rpm（目标 %.1f rpm）\n",
		       shaft_view.omega_mech * rpm_per_rad_s, omega_ref_mech * rpm_per_rad_s);
		printf("nan/inf：%s\n", sl_bad ? "检测到（不合格）" : "无");

		// 任务剖面摘要：仅在切入发电（阶段 3）时打印。sensorless_start 不进阶段 3，
		// 故其 stdout 与改动前逐字节一致（回归）。各判据阈值见验收标准。
		if (sl_gen_t >= 0.0) {
			printf("\n=== 任务剖面摘要（无感全包线：起动 -> 速度闭环 -> 发电稳压）===\n");
			printf("阶段 0->1（I/F->过渡）   ：t = %.4f s，转速 %.1f rpm\n",
			       sl_switch_t, sl_switch_omega * rpm_per_rad_s);
			printf("阶段 1->2（进入速度闭环） ：t = %.4f s，转速 %.1f rpm\n",
			       sl_closed_t, sl_closed_omega * rpm_per_rad_s);
			printf("阶段 2->3（切入发电稳压） ：t = %.4f s，转速 %.1f rpm\n",
			       sl_gen_t, sl_gen_omega * rpm_per_rad_s);
			printf("交接瞬间 |Δiq_ref|：%.3f A（bumpless 判据 <10 A）\n", sl_gen_iq_jump);
			printf("发电段 v_dc 瞬态偏差（切入 0.5 s 内）：%.3f V（判据 <3 V）\n",
			       sl_gen_vdc_transient);
			printf("发电段 v_dc 稳态偏差（末 1 s）：%.3f V（判据 <1 V，v_ref = %.1f V）\n",
			       sl_gen_vdc_steady, cfg.v_ref);
			printf("发电段 max|theta_err|：%.3f deg（判据 <5 deg）；稳态末 1 s：%.3f deg（判据 <2 deg）\n",
			       sl_gen_theta_err_max, sl_gen_theta_err_steady);
			printf("发电段 min mod_scale：%.4f（=1 全程未触电压边界 25019 rpm）\n", sl_gen_mod_min);
			printf("平衡转速（末端）：%.1f rpm（判据 <24500 rpm）\n",
			       shaft_view.omega_mech * rpm_per_rad_s);
			printf("末端 p_inv：%.1f W（<0 = 回馈充电）\n", inverter_view.p_inv);
			printf("SOC：起点 %.5f -> 发电切入 %.5f -> 终态 %.5f（发电段 %s）\n",
			       soc_init_snapshot, sl_gen_soc, battery_view.soc,
			       battery_view.soc > sl_gen_soc ? "上升=充电 ✓" : "未上升 ✗");
		}
	}
	return 0;
}
