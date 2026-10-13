#!/usr/bin/env python3
"""部件特性数据生成统一入口。

封装 build/char 二进制的编译、调用与可选绘图：
    python tools/char.py turbine                # 单部件（→ results/char_turbine.csv）
    python tools/char.py all                    # 全部部件
    python tools/char.py all --plot             # 生成 + 绘图（→ results/figures/char/）
    python tools/char.py pump -n 61 -c params.txt --plot
"""
import argparse
import os
import subprocess
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)		# 仓库根目录（tools/ 的上一级）
CHAR_BIN = os.path.join(PROJECT_ROOT, 'build', 'char')
RESULTS_DIR = os.path.join(PROJECT_ROOT, 'results')

COMPONENTS = ('turbine', 'pump', 'battery', 'battery-dynamic')
# 校验类子命令：不参与 all，输出 CSV 也不参与 --plot 识别
CHECKS = ('turbine-check', 'turbine-check-dynamic')


def ensure_binary():
	"""特性工具二进制不存在时先编译。"""
	if os.path.exists(CHAR_BIN):
		return
	print(f'[char] 未找到 {os.path.relpath(CHAR_BIN, PROJECT_ROOT)}，先执行 make char ...')
	subprocess.run(['make', 'char'], cwd=PROJECT_ROOT, check=True)


def run_component(comp, args):
	"""运行单部件扫描，返回输出 CSV 路径。"""
	out_csv = os.path.join(RESULTS_DIR, f'char_{comp}.csv')
	cmd = [CHAR_BIN, comp, '--grid', str(args.grid), '--output', out_csv]
	if args.config:
		cmd += ['--config', os.path.abspath(args.config)]
	subprocess.run(cmd, cwd=PROJECT_ROOT, check=True)
	return out_csv


def parse_args():
	parser = argparse.ArgumentParser(description='ATR_C 部件特性数据生成统一入口')
	parser.add_argument('component', choices=COMPONENTS + CHECKS + ('all',),
						help='目标部件、校验（turbine-check / turbine-check-dynamic），或 all 生成全部特性')
	parser.add_argument('-n', '--grid', type=int, default=41,
						help='每维网格点数（默认 41）')
	parser.add_argument('-c', '--config', default=None,
						help='组件参数文件（与主仿真 --config 同格式）')
	parser.add_argument('--plot', action='store_true',
						help='生成后调用 plotting/plot_char.py 出图')
	parser.add_argument('--outdir', default=os.path.join(RESULTS_DIR, 'figures', 'char'),
						help='图片输出目录（默认 results/figures/char）')
	return parser.parse_args()


def main():
	args = parse_args()
	ensure_binary()

	targets = COMPONENTS if args.component == 'all' else (args.component,)
	csvs = [run_component(comp, args) for comp in targets]

	if args.plot:
		plot_csvs = [csv for comp, csv in zip(targets, csvs) if comp not in CHECKS]
		if not plot_csvs:
			return
		plot_script = os.path.join(PROJECT_ROOT, 'plotting', 'plot_char.py')
		cmd = [sys.executable, plot_script, '--outdir', args.outdir]
		for csv in plot_csvs:
			cmd += ['--csv', csv]
		subprocess.run(cmd, cwd=PROJECT_ROOT, check=True)


if __name__ == '__main__':
	main()
