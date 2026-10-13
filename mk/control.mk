# 控制器：电流环、速度环、电压环、SVPWM + 无感起动（STSMO/PLL/IF/管理器）
control_SRCS := \
	model/control/Current.c \
	model/control/Speed.c \
	model/control/SVPWM.c \
	model/control/Voltage.c \
	model/control/STSMO.c \
	model/control/PLL.c \
	model/control/IFStart.c \
	model/control/sensorless.c
