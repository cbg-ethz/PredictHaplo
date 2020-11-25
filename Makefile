SCYTHE = SCYTHE/include/
CPPFLAGS += -DSCYTHE_LAPACK -DSCYTHE_DEBUG=0 -DHAVE_TRUNC -I$(SCYTHE)
CXXFLAGS += -O3 -fPIC -std=c++11 -Werror # -Wall -Wextra
LDLIBS += -lblas -llapack
EXOBJS = PredictHaplo_externAlign.o


UNAME := $(shell uname)
ifeq ($(UNAME), Darwin)
    CXXFLAGS += -stdlib=libc++
endif


all: predicthaplo

clean:
	rm -f $(EXOBJS) predicthaplo

predicthaplo: $(EXOBJS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $^ $(LDLIBS)
