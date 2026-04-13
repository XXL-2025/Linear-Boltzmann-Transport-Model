EXE=main01
# 合并编译选项，统一使用 -fPIE
GCC=g++ -std=c++17 -fPIE -O3 -w -pthread

#PYTHIALIBS="-lpythia8"
PYTHIAINC=$(shell pythia8-config --cxxflags)
PYTHIALIBS=$(shell pythia8-config --libs)

ROOTINC=$(shell root-config --cflags)
#ROOTLIBS=$(shell root-config --libs) -lMathMore
ROOTLIBS=$(shell root-config --libs) -ltbb -lvdt

FASTJETINC=$(shell fastjet-config --cxxflags)
FASTJETLIBS=$(shell fastjet-config --libs)

# 定义通用的编译和链接选项
INC_FLAGS = $(PYTHIAINC) $(ROOTINC) $(FASTJETINC)
LIB_FLAGS = $(PYTHIALIBS) $(ROOTLIBS) $(FASTJETLIBS)

all: $(EXE)

main01.o: main01.cc
	$(GCC) $(INC_FLAGS) -c main01.cc

$(EXE): main01.o
	$(GCC) $(INC_FLAGS) -pie -o $(EXE) main01.o $(LIB_FLAGS)

.PHONY: clean
clean:
	-rm *.o $(EXE)
