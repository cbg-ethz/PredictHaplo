CFLAGS = -O3 -DSCYTHE_LAPACK -DSCYTHE_DEBUG=0 -DHAVE_TRUNC -fPIC -std=c++11 -Werror # -Wall -Wextra
SCYTHE = SCYTHE/include/
EXOBJS = PredictHaplo_externAlign.o


UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    CFLAGS += -stdlib=libc++
endif


all: PredictHaplo-Paired

clean:
	rm -f $(EXOBJS) PredictHaplo-Paired

PredictHaplo-Paired: $(EXOBJS)
	g++ $(CFLAGS) -o $@ $(EXOBJS) -lblas -llapack
PredictHaplo_externAlign.o: PredictHaplo_externAlign.cpp
	g++ $(CFLAGS) -I$(SCYTHE) -c -o $@ $<
