#include "Logger.h"

#include "Motor.h"
#include "Inverter.h"
#include "Shaft.h"
#include "Pump.h"
#include "Turbine.h"
#include "Battery.h"

// ---------------------------------------------------------------------------
// 列定义表：增删/调整日志列只需修改本表，表头与数据行自动保持一致。
// 列顺序与历史 CSV 格式保持一致（首列为仿真时间 t）。
// ---------------------------------------------------------------------------

static double col_t(double t, const LogViews* views) {
    (void)views;
    return t;
}

static double col_id(double t, const LogViews* views) {
    (void)t;
    return views->motor->id;
}

static double col_iq(double t, const LogViews* views) {
    (void)t;
    return views->motor->iq;
}

static double col_omega_mech(double t, const LogViews* views) {
    (void)t;
    return views->shaft->omega_mech;
}

static double col_theta_mech(double t, const LogViews* views) {
    (void)t;
    return views->shaft->theta_mech;
}

static double col_tem(double t, const LogViews* views) {
    (void)t;
    return views->motor->T_em;
}

static double col_t_turb(double t, const LogViews* views) {
    (void)t;
    return views->turbine->T_turb;
}

static double col_t_pump(double t, const LogViews* views) {
    (void)t;
    return views->pump->T_pump;
}

static double col_q_pump(double t, const LogViews* views) {
    (void)t;
    return views->pump->Q_pump;
}

static double col_p_discharge(double t, const LogViews* views) {
    (void)t;
    return views->pump->p_discharge;
}

static double col_v_dc(double t, const LogViews* views) {
    (void)t;
    return views->inverter->v_dc;
}

static double col_ua(double t, const LogViews* views) {
    (void)t;
    return views->inverter->ua;
}

static double col_ub(double t, const LogViews* views) {
    (void)t;
    return views->inverter->ub;
}

static double col_uc(double t, const LogViews* views) {
    (void)t;
    return views->inverter->uc;
}

static double col_ia(double t, const LogViews* views) {
    (void)t;
    return views->motor->ia;
}

static double col_ib(double t, const LogViews* views) {
    (void)t;
    return views->motor->ib;
}

static double col_ic(double t, const LogViews* views) {
    (void)t;
    return views->motor->ic;
}

static double col_i_inv_dc(double t, const LogViews* views) {
    (void)t;
    return views->inverter->i_inv_dc;
}

static double col_p_inv(double t, const LogViews* views) {
    (void)t;
    return views->inverter->p_inv;
}

static double col_v_batt(double t, const LogViews* views) {
    (void)t;
    return views->battery->v_batt;
}

static double col_soc(double t, const LogViews* views) {
    (void)t;
    return views->battery->soc;
}

static double col_mod_scale(double t, const LogViews* views) {
    (void)t;
    return views->control->mod_scale;
}

static double col_iq_ref_cmd(double t, const LogViews* views) {
    (void)t;
    return views->control->iq_ref_cmd;
}

static double col_ud_ref(double t, const LogViews* views) {
    (void)t;
    return views->control->ud_ref;
}

static double col_uq_ref(double t, const LogViews* views) {
    (void)t;
    return views->control->uq_ref;
}

// ---- 无感起动诊断列（仅 sensorless_columns 选用；默认 25 列不含，有感回归不变）----
static double col_theta_err_rad(double t, const LogViews* views) {
    (void)t;
    return views->control->theta_err_rad;
}

static double col_omega_est(double t, const LogViews* views) {
    (void)t;
    return views->control->omega_est;
}

static double col_startup_stage(double t, const LogViews* views) {
    (void)t;
    return views->control->startup_stage;
}

// 默认 25 列条目（ macro 复用，保证 default_columns 与 sensorless_columns 前 25 列完全一致）
#define DEFAULT_COLUMN_ENTRIES                  \
    {"t", col_t},                               \
    {"id", col_id},                             \
    {"iq", col_iq},                             \
    {"omega_mech", col_omega_mech},             \
    {"theta_mech", col_theta_mech},             \
    {"Tem", col_tem},                           \
    {"T_turb", col_t_turb},                     \
    {"T_pump", col_t_pump},                     \
    {"Q_pump", col_q_pump},                     \
    {"p_discharge", col_p_discharge},           \
    {"v_dc", col_v_dc},                         \
    {"ua", col_ua},                             \
    {"ub", col_ub},                             \
    {"uc", col_uc},                             \
    {"ia", col_ia},                             \
    {"ib", col_ib},                             \
    {"ic", col_ic},                             \
    {"i_inv_dc", col_i_inv_dc},                 \
    {"p_inv", col_p_inv},                       \
    {"v_batt", col_v_batt},                     \
    {"soc", col_soc},                           \
    {"mod_scale", col_mod_scale},               \
    {"iq_ref_cmd", col_iq_ref_cmd},             \
    {"ud_ref", col_ud_ref},                     \
    {"uq_ref", col_uq_ref}

static const LoggerColumn default_columns[] = {
    DEFAULT_COLUMN_ENTRIES,
};

// 无感起动扩展列（28）：默认 25 列 + 3 诊断列。仅 --sensorless 时经
// SimulationLogger_use_sensorless_columns 选用，默认路径仍为 25 列（回归逐字节一致）。
static const LoggerColumn sensorless_columns[] = {
    DEFAULT_COLUMN_ENTRIES,
    {"theta_err_rad", col_theta_err_rad},
    {"omega_est", col_omega_est},
    {"startup_stage", col_startup_stage},
};

// ---------------------------------------------------------------------------
// 通用写入原语
// ---------------------------------------------------------------------------

static int SimulationLogger_write_separator(SimulationLogger* logger) {
    if (logger == NULL || logger->stream == NULL) {
        return -1;
    }

    if (!logger->first_column) {
        if (fputc(',', logger->stream) == EOF) {
            return -1;
        }
    } else {
        logger->first_column = false;
    }

    return 0;
}

static void SimulationLogger_row_begin(SimulationLogger* logger) {
    if (logger == NULL || logger->stream == NULL) {
        return;
    }

    logger->first_column = true;
}

static int SimulationLogger_add_text(SimulationLogger* logger, const char* value) {
    if (logger == NULL || logger->stream == NULL || value == NULL) {
        return -1;
    }

    if (SimulationLogger_write_separator(logger) != 0) {
        return -1;
    }

    return fputs(value, logger->stream) == EOF ? -1 : 0;
}

static int SimulationLogger_add_double(SimulationLogger* logger, double value) {
    if (logger == NULL || logger->stream == NULL) {
        return -1;
    }

    if (SimulationLogger_write_separator(logger) != 0) {
        return -1;
    }

    return fprintf(logger->stream, "%.6f", value) < 0 ? -1 : 0;
}

static void SimulationLogger_row_end(SimulationLogger* logger) {
    if (logger == NULL || logger->stream == NULL) {
        return;
    }

    fputc('\n', logger->stream);
    logger->first_column = true;
}

// ---------------------------------------------------------------------------
// 公共接口
// ---------------------------------------------------------------------------

bool SimulationLogger_init(SimulationLogger* logger, const char* path) {
    if (logger == NULL || path == NULL) {
        return false;
    }

    logger->stream = fopen(path, "w");
    logger->first_column = true;
    logger->columns = default_columns;
    logger->column_count = (int)(sizeof(default_columns) / sizeof(default_columns[0]));
    return logger->stream != NULL;
}

// 覆盖默认列定义表（须在 write_header 之前调用）
void SimulationLogger_set_columns(SimulationLogger* logger,
                                  const LoggerColumn* columns,
                                  int column_count) {
    if (logger == NULL || columns == NULL || column_count < 1) {
        return;
    }

    logger->columns = columns;
    logger->column_count = column_count;
}

// 选用无感起动扩展列（28 列）。仅 --sensorless 时调用，默认路径仍为 25 列。
void SimulationLogger_use_sensorless_columns(SimulationLogger* logger) {
    SimulationLogger_set_columns(logger, sensorless_columns,
                                 (int)(sizeof(sensorless_columns) / sizeof(sensorless_columns[0])));
}

void SimulationLogger_write_header(SimulationLogger* logger) {
    if (logger == NULL || logger->stream == NULL) {
        return;
    }

    SimulationLogger_row_begin(logger);
    for (int i = 0; i < logger->column_count; ++i) {
        SimulationLogger_add_text(logger, logger->columns[i].name);
    }
    SimulationLogger_row_end(logger);
}

void SimulationLogger_log_snapshot(SimulationLogger* logger,
                                   double t,
                                   const LogViews* views) {
    if (logger == NULL || logger->stream == NULL || views == NULL ||
        views->motor == NULL || views->shaft == NULL || views->pump == NULL ||
        views->turbine == NULL || views->inverter == NULL || views->battery == NULL ||
        views->control == NULL) {
        return;
    }

    SimulationLogger_row_begin(logger);
    for (int i = 0; i < logger->column_count; ++i) {
        SimulationLogger_add_double(logger, logger->columns[i].extract(t, views));
    }
    SimulationLogger_row_end(logger);
}

void SimulationLogger_close(SimulationLogger* logger) {
    if (logger == NULL || logger->stream == NULL) {
        return;
    }

    fclose(logger->stream);
    logger->stream = NULL;
}
