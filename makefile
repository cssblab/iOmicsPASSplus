mkdir /bin
all: bin/iOmicsPASS

objs = main globals subjectClass PAMClass geneClass
objs := $(addsuffix .o, $(objs))
objs := $(addprefix iOmicsPASS/src/, $(objs))
CPPFLAGS = -Wall -O3 -std=c++0x -I./TMPLIB
CPP = g++

bin/iOmicsPASS: $(objs)
	$(CPP) $(CPPFLAGS) -o $@ $(objs) 
	rm -f $(objs)

.PHONY : clean
clean :
	$(RM) bin/iOmicsPASS $(objs)
