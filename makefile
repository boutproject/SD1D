
BOUT_TOP	= ../..

TARGET = sd1d
SOURCEC		= sd1d.cxx div_ops.cxx loadmetric.cxx radiation.cxx cubic_spline_local.cxx non-local_parallel.cxx non-local_parallel_integration.cxx non-local_parallel_toroidal_solver.cxx non-local_parallel_serial_matrixsolve.cxx non-local_parallel_sheath_solver.cxx

include $(BOUT_TOP)/make.config
