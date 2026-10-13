// 部件稳态特性扫描工具：各部件特性数据生成的统一 C 入口。
//
// 用法：
//   build/char <turbine|pump|battery> [选项]
//   -n, --grid N        每维网格点数（默认 41）
//   -o, --output PATH   输出 CSV（默认 results/char_<部件>.csv）
//   -c, --config FILE   组件参数文件（与主仿真 --config 同格式）
//
// 各部件的扫描逻辑见 tools/char_<部件>.c，统一经调度表分发；
// 新增部件：实现 char_sweep_xxx() 并在下方 COMPONENTS 追加一行。

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

#include "char_sweep.h"
#include "utils/param_file.h"

typedef struct CharComponent {
	const char* name;
	CharSweepFn sweep;
} CharComponent;

static const CharComponent COMPONENTS[] = {
	{ "turbine", char_sweep_turbine },
	{ "pump",    char_sweep_pump },
	{ "battery", char_sweep_battery },
	{ "battery-dynamic", char_sweep_battery_dynamic },	// 阶跃响应动态特性
	{ "turbine-check", char_check_turbine },	// 涡轮 7 点试验数据校验
	{ "turbine-check-dynamic", char_check_turbine_dynamic },	// 动态稳态校验
};
#define N_COMPONENTS ((int)(sizeof(COMPONENTS) / sizeof(COMPONENTS[0])))

void char_linspace(double* xs, double x0, double x1, int n) {
	for (int i = 0; i < n; i++) {
		xs[i] = x0 + (x1 - x0) * (double)i / (double)(n - 1);
	}
}

static void print_usage(const char* prog) {
	fprintf(stderr,
		"用法：%s <部件> [选项]\n"
		"部件：",
		prog);
	for (int i = 0; i < N_COMPONENTS; i++) {
		fprintf(stderr, "%s%s", COMPONENTS[i].name, i + 1 < N_COMPONENTS ? " | " : "\n");
	}
	fprintf(stderr,
		"选项：\n"
		"  -n, --grid N        每维网格点数（默认 41）\n"
		"  -o, --output PATH   输出 CSV（默认 results/char_<部件>.csv）\n"
		"  -c, --config FILE   组件参数文件（与主仿真 --config 同格式）\n"
		"  -h, --help          显示本帮助\n");
}

// 解析 "--flag value" 或 "--flag=value"（argv[*i] 为当前选项本体），
// 返回取到的值；缺失时报错并返回 NULL
static const char* option_value(int argc, char** argv, int* i) {
	const char* eq = strchr(argv[*i], '=');
	if (eq != NULL) {
		return eq + 1;
	}
	if (*i + 1 >= argc) {
		fprintf(stderr, "错误：选项 %s 缺少取值\n", argv[*i]);
		return NULL;
	}
	(*i)++;
	return argv[*i];
}

int main(int argc, char** argv) {
	if (argc < 2) {
		print_usage(argv[0]);
		return 1;
	}
	if (strcmp(argv[1], "-h") == 0 || strcmp(argv[1], "--help") == 0) {
		print_usage(argv[0]);
		return 0;
	}

	const CharComponent* comp = NULL;
	for (int i = 0; i < N_COMPONENTS; i++) {
		if (strcmp(argv[1], COMPONENTS[i].name) == 0) {
			comp = &COMPONENTS[i];
			break;
		}
	}
	if (comp == NULL) {
		fprintf(stderr, "错误：未知部件 '%s'\n\n", argv[1]);
		print_usage(argv[0]);
		return 1;
	}

	int grid = 41;
	const char* output_path = NULL;
	const char* config_path = NULL;
	for (int i = 2; i < argc; i++) {
		const char* value = NULL;
		if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0) {
			print_usage(argv[0]);
			return 0;
		} else if (strcmp(argv[i], "-n") == 0 || strncmp(argv[i], "--grid", 6) == 0) {
			if ((value = option_value(argc, argv, &i)) == NULL) {
				return 1;
			}
			grid = atoi(value);
			if (grid < 2) {
				fprintf(stderr, "错误：网格点数须 >= 2（得到 %d）\n", grid);
				return 1;
			}
		} else if (strcmp(argv[i], "-o") == 0 || strncmp(argv[i], "--output", 8) == 0) {
			if ((output_path = option_value(argc, argv, &i)) == NULL) {
				return 1;
			}
		} else if (strcmp(argv[i], "-c") == 0 || strncmp(argv[i], "--config", 8) == 0) {
			if ((config_path = option_value(argc, argv, &i)) == NULL) {
				return 1;
			}
		} else {
			fprintf(stderr, "错误：未知选项 '%s'\n\n", argv[i]);
			print_usage(argv[0]);
			return 1;
		}
	}

	if (config_path != NULL && param_file_apply(config_path) != 0) {
		return 1;
	}

	char default_path[256];
	if (output_path == NULL) {
		snprintf(default_path, sizeof(default_path), "results/char_%s.csv", comp->name);
		output_path = default_path;
	}

	// 确保输出目录存在（仅最后一级，与主仿真同策略）
	const char* slash = strrchr(output_path, '/');
	if (slash != NULL && slash != output_path) {
		char dir_buf[256];
		const size_t dir_len = (size_t)(slash - output_path);
		if (dir_len < sizeof(dir_buf)) {
			memcpy(dir_buf, output_path, dir_len);
			dir_buf[dir_len] = '\0';
			ATR_MKDIR(dir_buf, 0755);
		}
	}

	FILE* f = fopen(output_path, "w");
	if (f == NULL) {
		fprintf(stderr, "错误：无法打开输出文件 %s\n", output_path);
		return 1;
	}

	const int status = comp->sweep(f, grid);
	fclose(f);
	if (status != 0) {
		fprintf(stderr, "错误：%s 特性扫描失败\n", comp->name);
		return 1;
	}
	printf("Characteristic data written to %s\n", output_path);
	return 0;
}
