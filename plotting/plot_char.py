#!/usr/bin/env python3
"""部件特性曲线绘图：读取 results/char_<部件>.csv，输出特性曲线图。

用法（通常由 tools/char.py --plot 自动调用）:
    python plotting/plot_char.py --csv results/char_turbine.csv [--outdir ...]
    python plotting/plot_char.py --csv a.csv --csv b.csv --show
"""
import argparse
import os

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 注册 3d 投影

from plot import setup_style, load_data

MAX_FAMILIES = 8	# 参数族曲线最多取的转速档数，避免图面拥挤


def family_values(data, key):
	"""取参数族坐标的档位（均匀抽样至最多 MAX_FAMILIES 档）。"""
	values = sorted(set(data[key]))
	if len(values) <= MAX_FAMILIES:
		return values
	idx = np.linspace(0, len(values) - 1, MAX_FAMILIES).astype(int)
	return [values[i] for i in idx]


def plot_family(data, x_col, y_col, fam_col, fam_fmt, fam_unit,
				title, xlabel, ylabel, filename, outdir, dpi=300, show=False):
	"""按 fam_col 分档画参数族曲线并保存。"""
	plt.figure(figsize=(6, 4))
	for fam in family_values(data, fam_col):
		mask = data[fam_col] == fam
		plt.plot(data[x_col][mask], data[y_col][mask], linewidth=2,
				 label=f'{fam_fmt.format(fam)} {fam_unit}')
	plt.title(title, fontsize=15, fontweight='bold')
	plt.xlabel(xlabel, fontsize=13)
	plt.ylabel(ylabel, fontsize=13)
	plt.grid(True, alpha=0.2, linestyle='--')
	plt.legend(frameon=False, fontsize=10)
	plt.tight_layout()
	os.makedirs(outdir, exist_ok=True)
	plt.savefig(os.path.join(outdir, filename), dpi=dpi, bbox_inches='tight')
	print(f"saved {os.path.join(outdir, filename)}")
	if show:
		plt.show()
	else:
		plt.close()


def plot_surface(data, x_col, y_col, z_col, title, xlabel, ylabel, zlabel,
				 filename, outdir, dpi=300, show=False, z_scale=1.0,
				 view=(28.0, -58.0)):
	"""二维扫描数据的三维曲面图（x×y → z，z 可选缩放换单位）。"""
	xs = sorted(set(data[x_col]))
	ys = sorted(set(data[y_col]))
	X, Y = np.meshgrid(xs, ys)
	Z = np.full(X.shape, np.nan)
	for row in data:
		i = ys.index(row[y_col])
		j = xs.index(row[x_col])
		Z[i, j] = row[z_col] / z_scale
	fig = plt.figure(figsize=(7.5, 5.5))
	ax = fig.add_subplot(projection='3d')
	if np.isnan(Z).any():
		# x 在不同 y 下不共用网格（如泵 Q 随转速变化）时矩形网格大量空洞；
		# 扫描行内点序对齐时按行对应构建结构化四边形曲面，否则退回三角网
		groups = [data[data[y_col] == yv] for yv in ys]
		if all(len(g) == len(groups[0]) for g in groups):
			Xs = np.array([[row[x_col] for row in g] for g in groups])
			Ys = np.array([[yv] * len(g) for yv, g in zip(ys, groups)])
			Zs = np.array([[row[z_col] / z_scale for row in g] for g in groups])
			ax.plot_surface(Xs, Ys, Zs, cmap='viridis',
							linewidth=0, antialiased=True)
		else:
			mask = ~np.isnan(Z)
			ax.plot_trisurf(X[mask], Y[mask], Z[mask], cmap='viridis',
							linewidth=0, antialiased=True)
	else:
		ax.plot_surface(X, Y, Z, cmap='viridis', linewidth=0, antialiased=True)
	ax.set_title(title, fontsize=14, fontweight='bold', pad=12)
	ax.set_xlabel(xlabel, fontsize=11, labelpad=8)
	ax.set_ylabel(ylabel, fontsize=11, labelpad=8)
	# 3D z 轴标签的 extent 不进入 tight 裁剪计算，savefig 时会被裁掉；
	# 改用 text2D 锚定在轴右侧（普通 Text 艺术家，extent 可被识别）
	ax.text2D(1.1, 0.5, zlabel, transform=ax.transAxes,
			  ha='left', va='center', fontsize=11)
	ax.tick_params(labelsize=9)
	ax.view_init(elev=view[0], azim=view[1])
	plt.tight_layout()
	os.makedirs(outdir, exist_ok=True)
	plt.savefig(os.path.join(outdir, filename), dpi=dpi, bbox_inches='tight')
	print(f"saved {os.path.join(outdir, filename)}")
	if show:
		plt.show()
	else:
		plt.close()


def plot_curve(x, y, title, xlabel, ylabel, filename, outdir, dpi=300, show=False):
	"""单曲线图（一维特性，如电池端口特性）。"""
	plt.figure(figsize=(6, 4))
	plt.plot(x, y, linewidth=2)
	plt.title(title, fontsize=15, fontweight='bold')
	plt.xlabel(xlabel, fontsize=13)
	plt.ylabel(ylabel, fontsize=13)
	plt.grid(True, alpha=0.2, linestyle='--')
	plt.tight_layout()
	os.makedirs(outdir, exist_ok=True)
	plt.savefig(os.path.join(outdir, filename), dpi=dpi, bbox_inches='tight')
	print(f"saved {os.path.join(outdir, filename)}")
	if show:
		plt.show()
	else:
		plt.close()


def render_turbine(data, outdir, dpi, show):
	xlabel = '膨胀比 π = p_plenum/p_down'
	rpm_label = '转速 (rpm)'
	plot_family(data, 'pi', 'm_dot_turb', 'omega_rpm', '{:.0f}', 'rpm',
				'涡轮流量特性', xlabel,
				'$\\dot{m}_{turb}$ (kg/s)', 'char_turbine_flow.png',
				outdir, dpi, show)
	plot_family(data, 'pi', 'power_W', 'omega_rpm', '{:.0f}', 'rpm',
				'涡轮功率特性', xlabel,
				'Power (W)', 'char_turbine_power.png',
				outdir, dpi, show)
	plot_family(data, 'pi', 'efficiency', 'omega_rpm', '{:.0f}', 'rpm',
				'涡轮等熵效率', xlabel,
				'Efficiency (-)', 'char_turbine_eff.png',
				outdir, dpi, show)
	plot_surface(data, 'pi', 'omega_rpm', 'm_dot_turb',
				 '涡轮流量特性（三维）', '$\\pi$', rpm_label,
				 '$\\dot{m}_{turb}$ (kg/s)', 'char_turbine_flow_3d.png',
				 outdir, dpi, show)
	plot_surface(data, 'pi', 'omega_rpm', 'power_W',
				 '涡轮功率特性（三维）', '$\\pi$', rpm_label,
				 'Power (kW)', 'char_turbine_power_3d.png',
				 outdir, dpi, show, z_scale=1e3)


def render_pump(data, outdir, dpi, show):
	if np.max(data['Q_pump']) <= 0.0:
		print("警告：Q_pump 全为 0 —— 扫描范围可能超出泵试验包线"
			  "（效率拟合为负的区域流量记 0），请检查参数或网格")
	rpm_label = '转速 (rpm)'
	plot_family(data, 'Q_pump', 'H_m', 'omega_rpm', '{:.0f}', 'rpm',
				'泵扬程-流量特性（H-Q）', 'Q (m$^3$/s)', 'H (m)',
				'char_pump_head.png', outdir, dpi, show)
	plot_family(data, 'Q_pump', 'T_pump', 'omega_rpm', '{:.0f}', 'rpm',
				'泵转矩-流量特性（T-Q）', 'Q (m$^3$/s)', 'T (N·m)',
				'char_pump_torque.png', outdir, dpi, show)
	plot_surface(data, 'Q_pump', 'omega_rpm', 'H_m',
				 '泵扬程特性（三维）', 'Q (m$^3$/s)', rpm_label,
				 'H (m)', 'char_pump_head_3d.png',
				 outdir, dpi, show)
	plot_surface(data, 'Q_pump', 'omega_rpm', 'efficiency',
				 '泵效率特性（三维）', 'Q (m$^3$/s)', rpm_label,
				 '$\\eta$', 'char_pump_eff_3d.png',
				 outdir, dpi, show)


def render_battery(data, outdir, dpi, show):
	p_kw = data['p_load_W'] / 1e3
	plot_curve(p_kw, data['v_batt_V'], '电池端电压-负载功率特性',
			   '$P_{load}$ (kW)', 'V (V)', 'char_battery_voltage.png',
			   outdir, dpi, show)
	plot_curve(p_kw, data['i_batt_A'], '电池电流-负载功率特性',
			   '$P_{load}$ (kW)', 'I (A)', 'char_battery_current.png',
			   outdir, dpi, show)
	plot_curve(p_kw, data['efficiency'], '电池效率-负载功率特性',
			   '$P_{load}$ (kW)', 'Efficiency (-)', 'char_battery_eff.png',
			   outdir, dpi, show)


def render_battery_dynamic(data, outdir, dpi, show):
	# 族分组用 kW 列、极化电压用两支路之和，追加派生列后统一出图
	combined = np.empty(len(data),
						dtype=np.dtype(data.dtype.descr +
									   [('p_load_kW', float), ('v_rc_V', float)]))
	for name in data.dtype.names:
		combined[name] = data[name]
	combined['p_load_kW'] = data['p_load_W'] / 1e3
	combined['v_rc_V'] = data['v_rc1_V'] + data['v_rc2_V']
	plot_family(combined, 't_s', 'v_batt_V', 'p_load_kW', '{:.0f}', 'kW',
				'电池阶跃响应-端电压', 't (s)', 'V (V)',
				'char_battery_dyn_voltage.png', outdir, dpi, show)
	plot_family(combined, 't_s', 'i_batt_A', 'p_load_kW', '{:.0f}', 'kW',
				'电池阶跃响应-电流', 't (s)', 'I (A)',
				'char_battery_dyn_current.png', outdir, dpi, show)
	plot_family(combined, 't_s', 'v_rc_V', 'p_load_kW', '{:.0f}', 'kW',
				'电池阶跃响应-极化电压', 't (s)', 'V (V)',
				'char_battery_dyn_polar.png', outdir, dpi, show)
	# 三维曲面：功率×时间满网格无空洞；视角选在扇面张开一侧，
	# 避免放电/充电两翼互相遮挡
	plot_surface(combined, 't_s', 'p_load_kW', 'v_batt_V',
				 '电池阶跃响应-端电压（三维）', 't (s)', '$P_{load}$ (kW)',
				 'V (V)', 'char_battery_dyn_voltage_3d.png',
				 outdir, dpi, show, view=(28.0, -145.0))
	plot_surface(combined, 't_s', 'p_load_kW', 'i_batt_A',
				 '电池阶跃响应-电流（三维）', 't (s)', '$P_{load}$ (kW)',
				 'I (A)', 'char_battery_dyn_current_3d.png',
				 outdir, dpi, show, view=(22.0, -115.0))
	plot_surface(combined, 't_s', 'p_load_kW', 'v_rc_V',
				 '电池阶跃响应-极化电压（三维）', 't (s)', '$P_{load}$ (kW)',
				 'V (V)', 'char_battery_dyn_polar_3d.png',
				 outdir, dpi, show, view=(22.0, -115.0))


# 按 CSV 特征列自动识别部件类型
DETECTORS = [
	(lambda cols: 'pi' in cols, render_turbine),
	(lambda cols: 'H_m' in cols, render_pump),
	(lambda cols: 't_s' in cols and 'v_rc1_V' in cols, render_battery_dynamic),
	(lambda cols: 'p_load_W' in cols, render_battery),
]


def parse_args():
	parser = argparse.ArgumentParser(description='ATR_C 部件特性曲线绘图')
	parser.add_argument('--csv', action='append', required=True,
						help='特性 CSV 路径，可重复指定多个')
	parser.add_argument('--outdir', required=True, help='图片输出目录')
	parser.add_argument('--dpi', type=int, default=300, help='保存图片 DPI（默认 300）')
	parser.add_argument('--show', action='store_true', help='交互式显示（默认只保存）')
	return parser.parse_args()


def main():
	args = parse_args()
	if not args.show:
		matplotlib.use('Agg')
	setup_style()
	for csv in args.csv:
		data = load_data(csv)
		cols = data.dtype.names or ()
		render = next((r for detect, r in DETECTORS if detect(cols)), None)
		if render is None:
			print(f"skip '{csv}': 无法识别部件类型（列：{list(cols)}）")
			continue
		print(f"plotting {csv} ({len(data)} rows)")
		render(data, args.outdir, args.dpi, args.show)


if __name__ == '__main__':
	main()
