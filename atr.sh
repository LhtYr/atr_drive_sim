#!/usr/bin/env bash
# ATR_C 新手入口：对常用 Makefile 目标的薄封装。
# 除 char 外所有命令都原样转发给 make，不含任何编译逻辑；
# 熟悉后可直接使用对应的 make 目标（帮助中标注了等价命令）。
#
# 用法：./atr.sh [命令] [参数]

set -e
cd "$(dirname "$0")"

usage() {
	cat <<'EOF'
ATR_C 快捷命令（熟悉后可直接使用等价的 make 目标）

用法：./atr.sh <命令> [参数]

  build      编译 release 版（默认命令）        （= make all）
  run        编译并运行仿真                     （= make run）
  plot       绘制仿真结果图（需先 run）          （= make plot）
  run-plot   一键仿真 + 出图                    （= make run-plot）
  char ...   部件特性数据工具（转发 tools/char.py）
             例：./atr.sh char all --plot
             例：./atr.sh char turbine -n 61 --config params.txt
  debug      编译 debug 版                      （= make debug）
  gdb        进入 GDB 调试                      （= make rungdb）
  clean      清理编译产物                        （= make clean）
  help       显示本帮助
EOF
}

cmd="${1:-build}"
case "$cmd" in
	help|-h|--help)
		usage
		;;
	build)
		exec make all
		;;
	run|run-plot|plot|debug|clean)
		exec make "$cmd"
		;;
	gdb)
		exec make rungdb
		;;
	char)
		shift
		PY=python3
		[ -x .venv/bin/python ] && PY=.venv/bin/python
		exec "$PY" tools/char.py "$@"
		;;
	*)
		echo "未知命令：$cmd" >&2
		echo >&2
		usage >&2
		exit 1
		;;
esac
