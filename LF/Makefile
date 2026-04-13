EXE=lbt
# 编译参数：注意这里不要加 -I，后续通过 $(PYTHIAINC) 等自动引入
GCC=g++ -std=c++17 -fPIE -O3 -w -pthread

# 获取 Pythia8 的编译参数（含 -I 头文件路径）和链接参数（含 -L 库路径和 -l 库名）
PYTHIAINC=$(shell pythia8-config --cxxflags)  # 例如：-I/path/to/pythia/include
PYTHIALIBS=$(shell pythia8-config --libs)      # 例如：-L/path/to/pythia/lib -lpythia8 -ldl

# 获取 FastJet 的编译和链接参数（同上）
FASTJETINC=$(shell fastjet-config --cxxflags)  # 例如：-I/path/to/fastjet/include
FASTJETLIBS=$(shell fastjet-config --libs)      # 例如：-L/path/to/fastjet/lib -lfastjet...

all: $(EXE)

# 编译 LBT.o：直接使用 PYTHIAINC 和 FASTJETINC（已包含 -I），无需额外加 -I
LBT.o: LBT.cpp
	$(GCC) $(PYTHIAINC) $(FASTJETINC) -c LBT.cpp -o $@

# 编译 main.o：同上
main.o: maintest.cpp
	$(GCC) $(PYTHIAINC) $(FASTJETINC) -c maintest.cpp -o $@

# 链接生成可执行文件：使用链接参数（PYTHIALIBS 和 FASTJETLIBS），无需 -I 选项
$(EXE): LBT.o maintest.o
	$(GCC) -o $(EXE) LBT.o maintest.o $(PYTHIALIBS) $(FASTJETLIBS)

.PHONY: clean
clean:
	-rm -f *.o $(EXE)  # 加 -f 避免删除不存在文件时的报错
