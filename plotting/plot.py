#!/usr/bin/env python3
"""仿真结果后处理绘图脚本。

用法示例:
    python plotting/plot.py                          # 默认读 results/simulation.csv
    python plotting/plot.py --csv path/to/sim.csv --outdir out --tmax 2.0 --show
"""
import argparse
import math
import os

import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)		# 仓库根目录（plotting/ 的上一级）
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')


def setup_style():
	"""注册中英文字体并设置全局绘图风格。"""
	for font in ('TimesNewRoman.ttf', 'SimSun.ttf'):
		path = os.path.join(SCRIPT_DIR, 'fonts', font)
		if os.path.exists(path):
			fm.fontManager.addfont(path)
	plt.rcParams['font.family'] = ['Times New Roman', 'SimSun', 'DejaVu Serif']
	plt.rcParams['axes.unicode_minus'] = False
	plt.rcParams['font.size'] = 12


def load_data(csv_path):
	"""加载仿真 CSV，返回按列名访问的结构化数组。"""
	data = np.genfromtxt(csv_path, delimiter=',', names=True)
	if data.ndim == 0:  # 只有一行数据时 genfromtxt 返回标量记录
		data = data.reshape(1)
	return data


def plot_series(x, y_list, labels, title, xlabel='t (s)', ylabel='',
				figsize=(6, 4), filename=None, outdir=None, dpi=300, show=False):
	"""统一绘图函数：多条曲线、网格、图例，可保存图片。"""
	plt.figure(figsize=figsize)
	for y, label in zip(y_list, labels):
		plt.plot(x, y, linewidth=2, label=label)
	plt.title(title, fontsize=15, fontweight='bold')
	plt.xlabel(xlabel, fontsize=13)
	plt.ylabel(ylabel, fontsize=13)
	plt.grid(True, alpha=0.2, linestyle='--')
	plt.legend(frameon=False, fontsize=12)
	plt.tight_layout()
	if filename and outdir:
		os.makedirs(outdir, exist_ok=True)
		plt.savefig(os.path.join(outdir, filename), dpi=dpi, bbox_inches='tight')
	if show:
		plt.show()
	else:
		plt.close()


# 声明式图表清单：每项描述一张图，新增图表只需在此追加。
# columns: CSV 列名（'omega_rpm' 为派生列）; tmax: 只画前 tmax 秒（None 为全程）
PLOTS = [
	dict(columns=['omega_rpm'], labels=['机械转速 omega_mech (rpm)'],
		 title='机械转速响应', ylabel='omega (rpm)', filename='omega_mech.png'),
	dict(columns=['theta_mech'], labels=['机械角度 theta_mech (rad)'],
		 title='机械角度响应', ylabel='theta (rad)', filename='theta_mech.png'),
	dict(columns=['Tem'], labels=['Tem (Nm)'],
		 title='电磁转矩响应', ylabel='Torque (Nm)', filename='Tem.png'),
	dict(columns=['id', 'iq'], labels=['id (A)', 'iq (A)'],
		 title='dq 电流响应', ylabel='Current (A)', filename='dq_current.png'),
	dict(columns=['Tem', 'T_turb', 'T_pump'],
		 labels=['Tem (Nm)', 'T_turb (Nm)', 'T_pump (Nm)'],
		 title='转矩响应', ylabel='Torque (Nm)', filename='T.png'),
	dict(columns=['soc'], labels=['State of Charge (SOC)'],
		 title='电池荷电状态响应', ylabel='SOC (%)', filename='soc.png'),
	dict(columns=['v_batt'], labels=['Battery Voltage (V)'],
		 title='电池电压响应', ylabel='Voltage (V)', filename='v_batt.png'),
	dict(columns=['ia', 'ib', 'ic'], labels=['ia (A)', 'ib (A)', 'ic (A)'],
		 title='三相电流响应', ylabel='Current (A)', filename='abc_current_full.png'),
	dict(columns=['ia', 'ib', 'ic'], labels=['ia (A)', 'ib (A)', 'ic (A)'],
		 title='三相电流响应', ylabel='Current (A)', filename='abc_current.png',
		 tmax=2.0),
	dict(columns=['p_inv'], labels=['Inverter Power (W)'],
		 title='逆变器功率响应', ylabel='Power (W)', filename='p_inv.png'),
]


def get_column(data, name):
	"""按名取列，支持派生列；列不存在返回 None。"""
	if name == 'omega_rpm':
		omega = get_column(data, 'omega_mech')
		return None if omega is None else omega * 60.0 / (2.0 * math.pi)
	return data[name] if name in (data.dtype.names or ()) else None


def render_plots(data, outdir, tmax=None, dpi=300, show=False):
	"""遍历图表清单绘图，缺失列的图跳过并提示。"""
	t = data['t']
	for spec in PLOTS:
		series = [get_column(data, name) for name in spec['columns']]
		missing = [name for name, y in zip(spec['columns'], series) if y is None]
		if missing:
			print(f"skip '{spec['filename']}': missing columns {missing}")
			continue
		# 时间窗：命令行 --tmax 全局生效，图表自带的 tmax 只约束该图
		t_lim = spec.get('tmax')
		if tmax is not None:
			t_lim = tmax if t_lim is None else min(t_lim, tmax)
		mask = slice(None) if t_lim is None else t <= t_lim
		plot_series(t[mask], [y[mask] for y in series], spec['labels'],
					spec['title'], ylabel=spec['ylabel'],
					filename=spec['filename'], outdir=outdir, dpi=dpi, show=show)
		print(f"saved {os.path.join(outdir, spec['filename'])}")


def _save_fig(filename, outdir, dpi, show):
	"""保存/显示当前图片并关闭（与 plot_series 的落盘行为一致）。"""
	plt.tight_layout()
	if filename and outdir:
		os.makedirs(outdir, exist_ok=True)
		plt.savefig(os.path.join(outdir, filename), dpi=dpi, bbox_inches='tight')
	if show:
		plt.show()
	else:
		plt.close()


def plot_gen_voltage_reg(data, outdir, vref=540.0, tmax=None, dpi=300, show=False):
	"""发电稳压能力边界图集（需控制链观测列，默认速度模式 CSV 自动跳过）。

	理论临界转速：反电动势 E = p·ω_m·ψ_f 达到 v_ref/√3（SVPWM 线性区上限），
	psi_f = 0.0595 Wb、p = 2、v_ref = 540 V 时为 25019 rpm。
	"""
	need = ('v_dc', 'iq', 'iq_ref_cmd', 'mod_scale', 'omega_mech')
	missing = [c for c in need if c not in (data.dtype.names or ())]
	if missing:
		print(f"skip gen-voltage-reg plots: missing columns {missing}")
		return
	t = data['t']
	mask = slice(None) if tmax is None else t <= tmax
	ts, v_dc = t[mask], data['v_dc'][mask]
	iq, iq_ref = data['iq'][mask], data['iq_ref_cmd'][mask]
	rpm = data['omega_mech'][mask] * 60.0 / (2.0 * math.pi)
	crit_rpm = (vref / math.sqrt(3.0)) / (2 * 0.0595) * 30.0 / math.pi

	plt.figure(figsize=(6, 4))
	plt.plot(ts, v_dc, linewidth=2, label='v_dc (V)')
	plt.axhline(vref, color='gray', linestyle='--', linewidth=1.5,
				label=f'v_ref = {vref:g} V')
	plt.title('发电稳压：母线电压跟踪', fontsize=15, fontweight='bold')
	plt.xlabel('t (s)', fontsize=13)
	plt.ylabel('Voltage (V)', fontsize=13)
	plt.grid(True, alpha=0.2, linestyle='--')
	plt.legend(frameon=False, fontsize=12)
	_save_fig('gen_v_dc.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'gen_v_dc.png')}")

	plt.figure(figsize=(6, 4))
	plt.plot(ts, iq, linewidth=2, label='iq (A)')
	plt.plot(ts, iq_ref, linewidth=2, linestyle='--', label='iq_ref (A)')
	plt.title('发电稳压：q 轴电流跟踪', fontsize=15, fontweight='bold')
	plt.xlabel('t (s)', fontsize=13)
	plt.ylabel('Current (A)', fontsize=13)
	plt.grid(True, alpha=0.2, linestyle='--')
	plt.legend(frameon=False, fontsize=12)
	_save_fig('gen_iq_tracking.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'gen_iq_tracking.png')}")

	# 能力边界主图（时间轴）：转速在过界后非单调（冲到峰值再折返振荡），按转速做横轴会
	# 出现多支重叠难以读；改按时间展开，左列全程、右列边界穿越放大，事件时刻从数据自动检测。
	sat_mask = data['mod_scale'][mask] < 1.0
	sat_i = int(np.argmax(sat_mask)) if sat_mask.any() else None
	vdev_i = None
	if sat_i is not None:
		vdev_mask = np.abs(v_dc[sat_i:] - vref) > 2.0
		if vdev_mask.any():
			vdev_i = sat_i + int(np.argmax(vdev_mask))
	t_zoom0 = ts[sat_i] - 0.05 if sat_i is not None else max(0.0, ts[-1] - 0.3)
	t_zoom1 = ts[vdev_i] + 0.06 if vdev_i is not None else t_zoom0 + 0.3
	fig, axes = plt.subplots(3, 2, figsize=(12, 9), sharex='col')
	for j, (ta, tb, tag) in enumerate([(ts[0], ts[-1], '全程'), (t_zoom0, t_zoom1, '边界穿越放大')]):
		win = (ts >= ta) & (ts <= tb)
		ax_ms, ax_iq, ax_v = axes[:, j]
		ax_ms.plot(ts[win], data['mod_scale'][mask][win], color='tab:blue', linewidth=1.2)
		ax_ms.axhline(1.0, color='k', linestyle='--', linewidth=0.8)
		ax_ms.set_ylabel('mod_scale', color='tab:blue', fontsize=12)
		ax_ms.set_ylim(0.95 if j else -0.05, 1.05)
		ax_iq.plot(ts[win], iq[win], label='iq', linewidth=1.2)
		ax_iq.plot(ts[win], iq_ref[win], label='iq_ref', linestyle='--', linewidth=1.2)
		ax_iq.set_ylabel('iq / iq_ref (A)', fontsize=12)
		ax_iq.legend(frameon=False, fontsize=10, loc='upper left')
		ax_v.plot(ts[win], v_dc[win], label='v_dc', linewidth=1.2)
		ax_v.axhline(vref, color='k', linestyle='--', linewidth=1.0, label='v_ref')
		ax_v.axhline(vref + 2.0, color='red', linestyle=':', linewidth=1.0)
		ax_v.axhline(vref - 2.0, color='red', linestyle=':', linewidth=1.0, label='v_ref ±2 V')
		ax_v.set_ylabel('v_dc (V)', fontsize=12)
		ax_v.set_xlabel('t (s)', fontsize=12)
		ax_v.legend(frameon=False, fontsize=10, loc='upper left')
		for ax in (ax_ms, ax_iq, ax_v):
			ax.grid(True, alpha=0.2, linestyle='--')
			if sat_i is not None and t_zoom0 <= ts[sat_i] <= tb:
				ax.axvline(ts[sat_i], color='tab:orange', linestyle='--', linewidth=1.2)
			if vdev_i is not None and t_zoom0 <= ts[vdev_i] <= tb:
				ax.axvline(ts[vdev_i], color='red', linestyle='--', linewidth=1.2)
		axes[0, j].set_title(tag, fontsize=13, fontweight='bold')
	if sat_i is not None:
		axes[0, 0].annotate(f"调制饱和起点\nn={rpm[sat_i]:.0f} rpm", (ts[sat_i], 1.0),
							xytext=(ts[sat_i] + 0.02 * ts[-1], 0.62),
							fontsize=10, color='tab:orange')
	if vdev_i is not None:
		axes[2, 1].annotate(f"|v_dc-v_ref|>2 V\nn={rpm[vdev_i]:.0f} rpm",
							(ts[vdev_i], v_dc[vdev_i]),
							xytext=(ts[vdev_i] - 0.12 * (t_zoom1 - t_zoom0), vref - 1.2),
							ha='right', fontsize=10, color='red')
	fig.suptitle(f'发电稳压能力边界（无弱磁，理论临界 {crit_rpm:.0f} rpm）',
				 fontsize=15, fontweight='bold')
	fig.tight_layout()
	_save_fig('gen_boundary.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'gen_boundary.png')}")

	# 电压上限图（升速段）：电压圆上限 U_max=v_dc/√3 基本不变，反电动势
	# E=p·ω·ψ_f 随转速线性升高；进入非线性（饱和）区后，电流控制可用裕度
	# U_max−E 随转速升高掉落并转负——mod_scale 同步跌落印证。只取升速段避免折返重叠。
	i_peak = int(np.argmax(rpm))
	rp_all = rpm[:i_peak + 1]
	Umax_all = v_dc[:i_peak + 1] / math.sqrt(3.0)
	E_all = (2 * 0.0595) * (rp_all * 2 * math.pi / 60.0)
	ms_all = data['mod_scale'][mask][:i_peak + 1]
	# 过界后转速存在极限环振荡（同一转速多次经过、mod_scale≈0），
	# 只保留“首次升速包络”（运行最大值滤波）保证曲线单值。
	env = rp_all >= np.maximum.accumulate(rp_all) - 1e-9
	rp, Umax, E, msv = rp_all[env], Umax_all[env], E_all[env], ms_all[env]
	room = Umax - E
	fig, (axA, axB) = plt.subplots(2, 1, figsize=(8, 7.5), sharex=True)
	for ax in (axA, axB):
		ax.axvspan(crit_rpm, rp[-1], color='red', alpha=0.06)
		ax.axvline(crit_rpm, color='red', linestyle=':', linewidth=1.2,
				   label=f'理论临界 {crit_rpm:.0f} rpm')
		ax.grid(True, alpha=0.2, linestyle='--')
	axA.plot(rp, Umax, color='k', linestyle='--', linewidth=1.5, label='电压上限 v_dc/√3')
	axA.plot(rp, E, color='red', linewidth=1.5, label='反电动势 E=p·ω·ψ_f')
	axA.set_ylabel('电压幅值 (V)', fontsize=12)
	axA.legend(frameon=False, fontsize=10)
	axB.plot(rp, room, color='tab:blue', linewidth=1.5, label='可用裕度 U_max−E')
	axB.axhline(0.0, color='k', linestyle='--', linewidth=0.8)
	axB.set_ylabel('可用电压裕度 (V)', fontsize=12, color='tab:blue')
	axB2 = axB.twinx()
	axB2.plot(rp, msv, color='tab:orange', linewidth=1.2, label='mod_scale')
	axB2.set_ylabel('mod_scale', fontsize=12, color='tab:orange')
	axB2.tick_params(axis='y', labelcolor='tab:orange')
	axB.set_xlabel('转速 (rpm)', fontsize=12)
	if sat_i is not None and sat_i <= i_peak:
		axB.annotate(f"进入非线性区\nn={rpm[sat_i]:.0f} rpm", (rpm[sat_i], room[sat_i]),
					 xytext=(rpm[sat_i] + 0.04 * rp[-1], room[sat_i] - 20),
					 fontsize=10, color='tab:orange')
	fig.suptitle('非线性区电压上限随转速升高掉落（无弱磁）', fontsize=15, fontweight='bold')
	fig.tight_layout()
	_save_fig('gen_v_limit.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'gen_v_limit.png')}")

	# 理论可稳压上限包络：每个转速下求电压环 PI 稳态还能维持的最高母线电压 v。
	# 稳态时电池充电电流 i_b=(v−v_oc)/R_tot（发电 iq=−i_b），所需电压
	# ud=ω_e·Lq·i_b、uq=ω_e·ψ_f−Rs·i_b 须落在电压圆 (v/√3)² 内；
	# 另 |iq|≤540 给出平顶 v_oc+540·R_tot。过临界后包络随转速掉落，
	# 无解即不存在稳压稳态（失控）。参数同 Battery.c / Motor.c。
	v_oc, r_tot = 540.0, 0.01 + 0.006 + 0.004
	rs, lq, psi_f, pp, iq_max = 0.0051, 34.6e-6, 0.0595, 2, 540.0
	v_top = v_oc + iq_max * r_tot

	def feas(v, we):
		i_b = (v - v_oc) / r_tot
		ud = we * lq * i_b
		uq = we * psi_f - rs * i_b
		return ud * ud + uq * uq <= v * v / 3.0

	rpm_grid = np.linspace(0.0, 30000.0, 601)
	env_v = np.full_like(rpm_grid, np.nan)
	for k, n in enumerate(rpm_grid):
		we = pp * n * 2 * math.pi / 60.0
		if not feas(v_oc, we):
			break  # 反电动势已超 v_oc/√3：任何 v 都无稳态解
		if feas(v_top, we):
			env_v[k] = v_top
			continue
		lo, hi = v_oc, v_top
		for _ in range(60):  # 二分求最大可行 v
			mid = 0.5 * (lo + hi)
			if feas(mid, we):
				lo = mid
			else:
				hi = mid
		env_v[k] = lo
	valid = ~np.isnan(env_v)
	n_end = rpm_grid[valid][-1] if valid.any() else crit_rpm
	fig, ax = plt.subplots(figsize=(8, 5))
	ax.plot(rpm_grid[valid], env_v[valid], color='tab:blue', linewidth=2,
			label='理论可稳压上限 max v_dc(n)')
	ax.axhline(v_top, color='gray', linestyle=':', linewidth=1.2,
			   label=f'电流限幅平顶 {v_top:.1f} V')
	ax.axhline(vref, color='k', linestyle='--', linewidth=1.0, label=f'v_ref={vref:.0f} V')
	ax.scatter(rpm, v_dc, s=2, alpha=0.25, color='tab:gray', label='实测 (n, v_dc)')
	ax.axvline(crit_rpm, color='red', linestyle=':', linewidth=1.2,
			   label=f'理论临界 {crit_rpm:.0f} rpm')
	ax.axvspan(n_end, rpm_grid[-1], color='red', alpha=0.06)
	ax.text(0.5 * (n_end + rpm_grid[-1]), v_oc + 2, '无稳压稳态\n（失控区）',
			ha='center', fontsize=11, color='red')
	cross = valid & (env_v <= vref)
	n_reg = rpm_grid[cross][0] if cross.any() else n_end
	ax.annotate(f'可稳压极限\n{n_reg:.0f} rpm', (n_reg, env_v[valid][-1] if valid.any() else vref),
				xytext=(n_reg - 0.16 * rpm_grid[-1], v_top + 2.5),
				fontsize=10, color='tab:blue')
	ax.set_xlim(0, rpm_grid[-1])
	ax.set_ylim(v_oc - 1, v_top + 6)
	ax.set_xlabel('转速 (rpm)', fontsize=12)
	ax.set_ylabel('稳态母线电压 (V)', fontsize=12)
	ax.legend(frameon=False, fontsize=9, loc='upper left')
	ax.grid(True, alpha=0.2, linestyle='--')
	fig.suptitle('转速升高 → 可稳压稳态电压掉落（无弱磁）', fontsize=15, fontweight='bold')
	fig.tight_layout()
	_save_fig('gen_v_envelope.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'gen_v_envelope.png')}")


def plot_sensorless_startup(data, outdir, tmax=None, dpi=300, show=False):
	"""无感起动图集（需无感诊断列 theta_err_rad/omega_est/startup_stage）。

	默认有感 CSV（25 列）缺这三列时整组自动跳过，符合“列存在才画”的约定。
	三张图对应验收判据的可视化：
	  1) sensorless_omega.png     —— 起动全程真值转速与估计转速对比（左全程、右切换放大）；
	  2) sensorless_theta_err.png —— 位置估计误差 theta_err 时程，标注切换/进入闭环时刻与判据线；
	  3) sensorless_stage.png     —— 起动状态机阶段阶梯叠加转速曲线。
	"""
	need = ('omega_mech', 'omega_est', 'theta_err_rad', 'startup_stage')
	missing = [c for c in need if c not in (data.dtype.names or ())]
	if missing:
		print(f"skip sensorless-startup plots: missing columns {missing}")
		return
	t = data['t']
	mask = slice(None) if tmax is None else t <= tmax
	ts = t[mask]
	rpm = data['omega_mech'][mask] * 60.0 / (2.0 * math.pi)
	rpm_est = data['omega_est'][mask] * 60.0 / (2.0 * math.pi)
	theta_err_deg = data['theta_err_rad'][mask] * 180.0 / math.pi
	stage = data['startup_stage'][mask]

	# 阶段跳变时刻：0->1 切换起点、1->2 进入闭环；用于竖线标注与放大窗右界
	def _first_at_or_above(level):
		idx = np.where(stage >= level - 0.5)[0]
		return float(ts[idx[0]]) if idx.size else None
	t_switch = _first_at_or_above(1)
	t_closed = _first_at_or_above(2)
	anchor = t_closed if t_closed is not None else (t_switch if t_switch is not None else 0.0)
	t_zoom1 = min(anchor + 0.6, float(ts[-1]))

	def _mark(ax):
		"""在给定坐标轴标注切换/进入闭环竖线（时刻落在窗内才画）。"""
		if t_switch is not None and ts[0] <= t_switch <= ts[-1]:
			ax.axvline(t_switch, color='tab:orange', linestyle='--', linewidth=1.2,
					   label=f'切换 {t_switch:.4f} s')
		if t_closed is not None and ts[0] <= t_closed <= ts[-1]:
			ax.axvline(t_closed, color='tab:green', linestyle=':', linewidth=1.2,
					   label=f'进入闭环 {t_closed:.4f} s')

	# --- 图 1：转速真值 vs 估计（左全程、右切换放大）---
	fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4))
	for ax, win, tag in ((axA, slice(None), '全程'), (axB, ts <= t_zoom1, '切换放大')):
		ax.plot(ts[win], rpm[win], linewidth=2, label='omega_mech 真值')
		ax.plot(ts[win], rpm_est[win], linewidth=1.6, linestyle='--', label='omega_est 估计')
		ax.set_title(tag, fontsize=13, fontweight='bold')
		ax.set_xlabel('t (s)', fontsize=12)
		ax.grid(True, alpha=0.2, linestyle='--')
		_mark(ax)
		ax.legend(frameon=False, fontsize=10, loc='upper left')
	axA.set_ylabel('omega (rpm)', fontsize=12)
	fig.suptitle('无感起动：转速真值与估计对比', fontsize=15, fontweight='bold')
	fig.tight_layout()
	_save_fig('sensorless_omega.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'sensorless_omega.png')}")

	# --- 图 2：theta_err 时程（左全程、右切换放大；±10° 切换判据、±5° 稳态判据）---
	fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4))
	for ax, win, tag in ((axA, slice(None), '全程'), (axB, ts <= t_zoom1, '切换放大')):
		ax.plot(ts[win], theta_err_deg[win], linewidth=1.4, color='tab:blue', label='theta_err')
		labeled = (ax is axA)  # 判据线仅在左图进图例，右图只画线避免重复条目
		ax.axhline(10.0, color='tab:orange', linestyle=':', linewidth=1.0,
				   label='±10° 切换判据' if labeled else None)
		ax.axhline(-10.0, color='tab:orange', linestyle=':', linewidth=1.0)
		ax.axhline(5.0, color='tab:green', linestyle=':', linewidth=1.0,
				   label='±5° 稳态判据' if labeled else None)
		ax.axhline(-5.0, color='tab:green', linestyle=':', linewidth=1.0)
		ax.set_title(tag, fontsize=13, fontweight='bold')
		ax.set_xlabel('t (s)', fontsize=12)
		ax.grid(True, alpha=0.2, linestyle='--')
		_mark(ax)
	axA.set_ylabel('theta_err (deg)', fontsize=12)
	if axA.get_legend_handles_labels()[1]:
		axA.legend(frameon=False, fontsize=10, loc='upper right')
	fig.suptitle('无感起动：位置估计误差 theta_err = wrap(θ̂ − p·theta_mech)',
				 fontsize=15, fontweight='bold')
	fig.tight_layout()
	_save_fig('sensorless_theta_err.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'sensorless_theta_err.png')}")

	# --- 图 3：状态机阶段阶梯叠加转速曲线（左轴阶段、右轴 rpm）---
	fig, axS = plt.subplots(figsize=(8, 4.5))
	axS.step(ts, stage, where='post', color='tab:red', linewidth=1.8, label='startup_stage')
	# 阶段刻度自适应：任务剖面工况（到达阶段 3 发电）扩展到 0/1/2/3，否则保持 0/1/2
	# （sensorless_start 不进阶段 3，此处行为与其原图逐像素一致）
	has_gen = float(np.max(stage)) >= 2.5
	axS.set_yticks([0, 1, 2, 3] if has_gen else [0, 1, 2])
	axS.set_yticklabels(['0 I/F', '1 过渡', '2 闭环', '3 发电'] if has_gen
						else ['0 I/F', '1 过渡', '2 闭环'])
	axS.set_ylim(-0.2, 3.3 if has_gen else 2.3)
	axS.set_xlabel('t (s)', fontsize=12)
	axS.set_ylabel('起动阶段', fontsize=12, color='tab:red')
	axS.tick_params(axis='y', labelcolor='tab:red')
	axS.grid(True, alpha=0.2, linestyle='--')
	_mark(axS)
	axR = axS.twinx()
	axR.plot(ts, rpm, color='tab:blue', linewidth=1.6, label='omega_mech (rpm)')
	axR.set_ylabel('omega (rpm)', fontsize=12, color='tab:blue')
	axR.tick_params(axis='y', labelcolor='tab:blue')
	lines = axS.get_lines() + axR.get_lines()
	axS.legend(lines, [ln.get_label() for ln in lines],
			   frameon=False, fontsize=10, loc='center right')
	fig.suptitle('无感起动：状态机阶段与转速', fontsize=15, fontweight='bold')
	fig.tight_layout()
	_save_fig('sensorless_stage.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'sensorless_stage.png')}")


def plot_mission_profile(data, outdir, vref=540.0, tmax=None, dpi=300, show=False):
	"""无感全包线任务剖面图集（需无感诊断列，且状态机到达阶段 3 发电）。

	仅当 CSV 含 theta_err_rad/startup_stage 等列且 startup_stage 出现过 3（发电）时才绘制；
	默认有感 CSV 或无感起动（不到阶段 3）自动跳过，符合“列存在才画”的约定。
	四张图对应任务剖面验收判据：
	  a) mission_speed_stage.png —— 全剖面转速 + 状态机阶段阶梯（0/1/2/3 标注）；
	  b) mission_iq.png          —— iq/iq_ref 全程（正转矩→零→负发电符号穿越）；
	  c) mission_v_dc.png        —— v_dc 与 v_ref（标注阶段 3 起点，验证稳压收敛）；
	  d) mission_theta_err.png   —— theta_err 全程（验证模式切换不扰观测器）。
	"""
	need = ('omega_mech', 'iq', 'iq_ref_cmd', 'v_dc', 'theta_err_rad', 'startup_stage')
	missing = [c for c in need if c not in (data.dtype.names or ())]
	if missing:
		print(f"skip mission-profile plots: missing columns {missing}")
		return
	if float(np.max(data['startup_stage'])) < 2.5:
		print("skip mission-profile plots: startup_stage 未到达阶段 3（非任务剖面工况）")
		return
	t = data['t']
	mask = slice(None) if tmax is None else t <= tmax
	ts = t[mask]
	stage = data['startup_stage'][mask]
	rpm = data['omega_mech'][mask] * 60.0 / (2.0 * math.pi)
	iq, iq_ref = data['iq'][mask], data['iq_ref_cmd'][mask]
	v_dc = data['v_dc'][mask]
	theta_err_deg = data['theta_err_rad'][mask] * 180.0 / math.pi

	# 阶段跳变时刻（0->1 过渡、1->2 闭环、2->3 发电）用于竖线标注
	def _first_at_or_above(level):
		idx = np.where(stage >= level - 0.5)[0]
		return float(ts[idx[0]]) if idx.size else None
	t_switch, t_closed, t_gen = _first_at_or_above(1), _first_at_or_above(2), _first_at_or_above(3)
	# 切换放大窗：发电切入点前 0.5 s 至后 2.0 s（覆盖交接瞬态与稳压收敛）
	z0 = (t_gen - 0.5) if t_gen is not None else float(ts[0])
	z1 = (t_gen + 2.0) if t_gen is not None else float(ts[-1])
	zoom = (ts >= z0) & (ts <= z1)

	def _mark(ax, gen_only=False):
		"""标注阶段切换竖线；gen_only 时只标发电切入点（放大图避免图例拥挤）。"""
		if not gen_only:
			if t_switch is not None and ts[0] <= t_switch <= ts[-1]:
				ax.axvline(t_switch, color='tab:orange', linestyle='--', linewidth=1.0,
						   label=f'0→1 过渡 {t_switch:.3f} s')
			if t_closed is not None and ts[0] <= t_closed <= ts[-1]:
				ax.axvline(t_closed, color='tab:green', linestyle=':', linewidth=1.0,
						   label=f'1→2 闭环 {t_closed:.3f} s')
		if t_gen is not None and ts[0] <= t_gen <= ts[-1]:
			ax.axvline(t_gen, color='tab:red', linestyle='--', linewidth=1.3,
					   label=f'2→3 切入发电 {t_gen:.3f} s')

	# --- 图 a：全剖面转速 + 状态机阶段阶梯（左轴阶段 0/1/2/3、右轴 rpm）---
	fig, axS = plt.subplots(figsize=(9, 4.8))
	axS.step(ts, stage, where='post', color='tab:red', linewidth=1.8, label='startup_stage')
	axS.set_yticks([0, 1, 2, 3])
	axS.set_yticklabels(['0 I/F', '1 过渡', '2 闭环', '3 发电'])
	axS.set_ylim(-0.2, 3.3)
	axS.set_xlabel('t (s)', fontsize=12)
	axS.set_ylabel('任务阶段', fontsize=12, color='tab:red')
	axS.tick_params(axis='y', labelcolor='tab:red')
	axS.grid(True, alpha=0.2, linestyle='--')
	_mark(axS)
	axR = axS.twinx()
	axR.plot(ts, rpm, color='tab:blue', linewidth=1.6, label='omega_mech (rpm)')
	axR.axhline(24500.0, color='gray', linestyle=':', linewidth=1.0, label='平衡上限 24500 rpm')
	axR.set_ylabel('omega (rpm)', fontsize=12, color='tab:blue')
	axR.tick_params(axis='y', labelcolor='tab:blue')
	lines = axS.get_lines() + axR.get_lines()
	axS.legend(lines, [ln.get_label() for ln in lines],
			   frameon=False, fontsize=9, loc='center right')
	fig.suptitle('任务剖面：全阶段转速与状态机（0 I/F→1 过渡→2 闭环→3 发电）',
				 fontsize=14, fontweight='bold')
	fig.tight_layout()
	_save_fig('mission_speed_stage.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'mission_speed_stage.png')}")

	# --- 图 b：iq/iq_ref 全程（左全程含正转矩->负发电符号穿越、右切换放大）---
	fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4))
	for ax, win, tag in ((axA, slice(None), '全程'), (axB, zoom, '切换放大')):
		ax.plot(ts[win], iq[win], linewidth=1.8, label='iq (A)')
		ax.plot(ts[win], iq_ref[win], linewidth=1.4, linestyle='--', label='iq_ref (A)')
		ax.axhline(0.0, color='k', linestyle='-', linewidth=0.8)
		ax.set_title(tag, fontsize=13, fontweight='bold')
		ax.set_xlabel('t (s)', fontsize=12)
		ax.grid(True, alpha=0.2, linestyle='--')
		_mark(ax, gen_only=(ax is axB))
		ax.legend(frameon=False, fontsize=10, loc='upper right')
	axA.set_ylabel('Current (A)', fontsize=12)
	fig.suptitle('任务剖面：q 轴电流 iq/iq_ref（正转矩 → 零 → 负发电符号穿越）',
				 fontsize=14, fontweight='bold')
	fig.tight_layout()
	_save_fig('mission_iq.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'mission_iq.png')}")

	# --- 图 c：v_dc 与 v_ref（左全程、右切入放大，标注阶段 3 起点与 ±1 V 稳态判据）---
	fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4))
	for ax, win, tag in ((axA, slice(None), '全程'), (axB, zoom, '切入发电放大')):
		ax.plot(ts[win], v_dc[win], linewidth=1.8, color='tab:blue', label='v_dc (V)')
		ax.axhline(vref, color='gray', linestyle='--', linewidth=1.4, label=f'v_ref = {vref:g} V')
		ax.axhline(vref + 1.0, color='tab:green', linestyle=':', linewidth=1.0)
		ax.axhline(vref - 1.0, color='tab:green', linestyle=':', linewidth=1.0, label='v_ref ±1 V')
		ax.set_title(tag, fontsize=13, fontweight='bold')
		ax.set_xlabel('t (s)', fontsize=12)
		ax.grid(True, alpha=0.2, linestyle='--')
		_mark(ax, gen_only=(ax is axB))
		ax.legend(frameon=False, fontsize=10, loc='lower right')
	axA.set_ylabel('Voltage (V)', fontsize=12)
	fig.suptitle('任务剖面：母线电压 v_dc 稳压（发电段收敛至 v_ref ±1 V）',
				 fontsize=14, fontweight='bold')
	fig.tight_layout()
	_save_fig('mission_v_dc.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'mission_v_dc.png')}")

	# --- 图 d：theta_err 全程（±5° 全程判据、±2° 发电稳态判据；验证模式切换不失锁）---
	fig, (axA, axB) = plt.subplots(1, 2, figsize=(12, 4))
	for ax, win, tag in ((axA, slice(None), '全程'), (axB, zoom, '切换放大')):
		ax.plot(ts[win], theta_err_deg[win], linewidth=1.3, color='tab:blue', label='theta_err')
		labeled = (ax is axA)
		ax.axhline(5.0, color='tab:green', linestyle=':', linewidth=1.0,
				   label='±5° 全程判据' if labeled else None)
		ax.axhline(-5.0, color='tab:green', linestyle=':', linewidth=1.0)
		ax.axhline(2.0, color='tab:orange', linestyle=':', linewidth=1.0,
				   label='±2° 发电稳态判据' if labeled else None)
		ax.axhline(-2.0, color='tab:orange', linestyle=':', linewidth=1.0)
		ax.set_title(tag, fontsize=13, fontweight='bold')
		ax.set_xlabel('t (s)', fontsize=12)
		ax.grid(True, alpha=0.2, linestyle='--')
		_mark(ax, gen_only=(ax is axB))
	axA.set_ylabel('theta_err (deg)', fontsize=12)
	if axA.get_legend_handles_labels()[1]:
		axA.legend(frameon=False, fontsize=10, loc='upper right')
	fig.suptitle('任务剖面：位置估计误差 theta_err（模式切换不失锁）',
				 fontsize=14, fontweight='bold')
	fig.tight_layout()
	_save_fig('mission_theta_err.png', outdir, dpi, show)
	print(f"saved {os.path.join(outdir, 'mission_theta_err.png')}")


def parse_args():
	parser = argparse.ArgumentParser(description='ATR_C 仿真结果绘图')
	parser.add_argument('--csv', default=os.path.join(RESULTS_DIR, 'simulation.csv'),
						help='仿真 CSV 路径（默认 results/simulation.csv）')
	parser.add_argument('--outdir', default=os.path.join(RESULTS_DIR, 'figures'),
						help='图片输出目录（默认 results/figures）')
	parser.add_argument('--tmax', type=float, default=None,
						help='只绘制前 tmax 秒（默认全程）')
	parser.add_argument('--vref', type=float, default=540.0,
						help='发电稳压图集的母线电压参考 V（默认 540，用于参考线与理论临界转速）')
	parser.add_argument('--dpi', type=int, default=300, help='保存图片 DPI（默认 300）')
	parser.add_argument('--show', action='store_true', help='交互式显示图片（默认只保存）')
	return parser.parse_args()


def main():
	args = parse_args()
	if not args.show:
		matplotlib.use('Agg')  # 批量模式无需 GUI 后端
	setup_style()
	data = load_data(args.csv)
	print(f"loaded {args.csv}: {len(data)} rows, columns: {list(data.dtype.names)}")
	render_plots(data, args.outdir, tmax=args.tmax, dpi=args.dpi, show=args.show)
	plot_gen_voltage_reg(data, args.outdir, vref=args.vref, tmax=args.tmax,
						 dpi=args.dpi, show=args.show)
	plot_sensorless_startup(data, args.outdir, tmax=args.tmax, dpi=args.dpi, show=args.show)
	plot_mission_profile(data, args.outdir, vref=args.vref, tmax=args.tmax,
						 dpi=args.dpi, show=args.show)


if __name__ == '__main__':
	main()
