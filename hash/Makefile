OBJS=	fmamp_subs.o \
        pol_subs.o \
        station_subs.o \
        uncert_subs.o \
    	util_subs.o \
        vel_subs.o

hash_hashpy1D: Makefile hash_hashpy1D.f $(OBJS)
	gfortran  hash_hashpy1D.f $(OBJS) -O -o hash_hashpy1D

pol_subs.o : Makefile pol_subs.f
	gfortran pol_subs.f -O -c

fmamp_subs.o : Makefile fmamp_subs.f
	gfortran fmamp_subs.f -O -c

station_subs.o : Makefile station_subs.f
	gfortran station_subs.f -O -c

uncert_subs.o : Makefile uncert_subs.f
	gfortran uncert_subs.f -O -c

util_subs.o : Makefile util_subs.f
	gfortran util_subs.f -O -c

vel_subs.o : Makefile vel_subs.f
	gfortran vel_subs.f -O -c



