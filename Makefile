# ============================================================
# ATR_C 顶层 Makefile
# 模块源文件清单拆分在 mk/*.mk 中，本文件只负责：
#   ① 聚合模块  ② 编译选项  ③ 构建规则  ④ 任务目标
# ============================================================

.DEFAULT_GOAL := all

# ---------- ① 模块源文件聚合（新增模块：新建 mk/模块名.mk 并在此 include） ----------
include mk/app.mk
include mk/component.mk
include mk/control.mk
include mk/solver.mk
include mk/tools.mk

SRCS := $(app_SRCS) $(component_SRCS) $(control_SRCS) $(solver_SRCS)

BIN       := build/main
BIN_DEBUG := build/main_debug
BIN_CHAR  := build/char
BIN_SELFCHECK := build/sl_selfcheck
# 特性工具链接对象：自带源码 + 复用 release 的组件/控制器/参数文件/边界默认值目标文件，不链求解器（参数绑定引用控制器全局实例，故控制链 .o 一并链入）
CHAR_OBJS := $(tools_SRCS:%.c=build/tools/%.o) $(component_SRCS:%.c=build/%.o) \
	$(control_SRCS:%.c=build/%.o) \
	build/src/param_file.o build/model/solver/BoundaryDefaults.o

# ---------- ② 编译选项（公共部分只定义一次） ----------
CC     := clang
INCLUDES := -Iinclude -Iinclude/component -Iinclude/control -Iinclude/solver
CFLAGS_COMMON := -std=c11 -Wall -Wextra -Wshadow -MMD -MP $(INCLUDES)
CFLAGS       := $(CFLAGS_COMMON) -O3 -march=native	# release
CFLAGS_DEBUG := $(CFLAGS_COMMON) -O0 -g3			# debug
LDFLAGS :=
LDLIBS  := -lm

# ---------- ③ 构建规则 ----------
# 编译器生成的头文件依赖（保证修改头文件后能正确重编译）
-include $(SRCS:%.c=build/%.d) $(SRCS:%.c=build/debug/%.d) $(tools_SRCS:%.c=build/tools/%.d)

# 构建规则模板：为指定目标文件目录生成一套「编译 .o + 链接」规则
# $(1)=目标文件目录  $(2)=可执行文件路径  $(3)=编译选项
define BUILD_TEMPLATE
$(1)/%.o: %.c
	@mkdir -p $$(dir $$@)
	$$(CC) $(3) -c $$< -o $$@

$(2): $$(SRCS:%.c=$(1)/%.o)
	@mkdir -p $$(dir $$@)
	$$(CC) $$(LDFLAGS) $$^ -o $$@ $$(LDLIBS)
endef

$(eval $(call BUILD_TEMPLATE,build,$(BIN),$(CFLAGS)))
$(eval $(call BUILD_TEMPLATE,build/debug,$(BIN_DEBUG),$(CFLAGS_DEBUG)))

# 特性工具：只编译 tools/ 自身源码，组件目标文件复用 release 版（见 CHAR_OBJS）
build/tools/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -c $< -o $@

$(BIN_CHAR): $(CHAR_OBJS)
	@mkdir -p $(dir $@)
	$(CC) $(LDFLAGS) $^ -o $@ $(LDLIBS)

# 运行仿真/绘图时静默 make 的命令回显，避免干扰进度条显示
ifneq ($(filter run run-plot,$(MAKECMDGOALS)),)
.SILENT:
endif

PYTHON ?= $(if $(wildcard .venv/bin/python),.venv/bin/python,python3)

# ---------- ④ 任务目标 ----------
BEAR ?= bear

.PHONY: all run debug compile-commands rungdb plot run-plot char selfcheck rebuild clean help

all: $(BIN)

# 运行仿真；ARGS="..." 可透传命令行参数给可执行文件（无感起动等）
# 例：make run ARGS="--sensorless --config config/sensorless_start.txt --duration 10"
run: $(BIN)
	./$(BIN) $(ARGS)

debug: $(BIN_DEBUG)

compile-commands:
	$(BEAR) -- $(MAKE) clean all debug

rungdb: debug
	lldb ./$(BIN_DEBUG)

# 绘制仿真结果图片（需先运行仿真生成 results/simulation.csv）
plot:
	@test -f results/simulation.csv || { echo "错误：未找到 results/simulation.csv，请先执行 make run"; exit 1; }
	$(PYTHON) plotting/plot.py

# 一键仿真 + 出图
run-plot: run plot

# 编译部件特性工具；运行入口为 tools/char.py（见 make help）
char: $(BIN_CHAR)

# 无感观测器离线自检：在合成 PMSM αβ 数据上验证 STSMO+PLL 收敛性（秒级，不依赖全系统仿真）
selfcheck:
	@mkdir -p build
	$(CC) -std=c11 -Wall -Wextra -Wshadow -O2 $(INCLUDES) \
		tools/sensorless_selfcheck.c model/control/STSMO.c model/control/PLL.c \
		-o $(BIN_SELFCHECK) $(LDLIBS)
	./$(BIN_SELFCHECK)

rebuild: clean all

clean:
	rm -rf build

help:
	@echo "可用目标："
	@echo "  all       编译 release 版（$(BIN)），默认目标"
	@echo "  run       编译并运行仿真（ARGS=\"...\" 透传参数，如 ARGS=\"--sensorless --config config/sensorless_start.txt\"）"
	@echo "  debug     编译 debug 版（$(BIN_DEBUG)）"
	@echo "  compile-commands  用 Bear 生成 clangd 所需 compile_commands.json"
	@echo "  rungdb    用 LLDB 调试 debug 版（CodeLLDB 外的命令行入口）"
	@echo "  plot      用 $(PYTHON) 绘制仿真结果图（需先 make run）"
	@echo "  run-plot  一键仿真 + 出图（同样支持 ARGS=\"...\" 透传）"
	@echo "  char      编译部件特性工具（$(BIN_CHAR)），运行入口："
	@echo "            $(PYTHON) tools/char.py <部件|all> [-n 网格数] [--config 参数文件] [--plot]"
	@echo "  selfcheck 编译并运行无感观测器离线自检（STSMO+PLL 收敛性，秒级）"
	@echo "  rebuild   clean 后重新编译"
	@echo "  clean     删除 build 目录"
	@echo "  help      显示本帮助"
