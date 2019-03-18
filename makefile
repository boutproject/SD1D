
BOUT_TOP	= ../..

TARGET = sd1d

DIRS = atomicpp 

SOURCEC		= sd1d.cxx div_ops.cxx loadmetric.cxx radiation.cxx \
                  reaction_impurity.cxx reaction_elastic.cxx reaction_recombination.cxx \
                  reaction_ionisation.cxx reaction_excitation.cxx \
                  reaction_cx.cxx reaction_neutralpressure.cxx \
                  species.cxx reaction_atomicpp.cxx reaction_neutral_diffusion.cxx \
                  reaction_loss_rate.cxx

include $(BOUT_TOP)/make.config
