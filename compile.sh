cd Complete-Striped-Smith-Waterman-Library/src;
make default;
cd ../..;
cd spoa;
mkdir build;
cd build;
cmake -DCMAKE_BUILD_TYPE=Release ..;
make;
cd ../../;
cp Complete-Striped-Smith-Waterman-Library/src/libssw.so libssw.so
mpic++ --std=c++17 -Xpreprocessor -fopenmp -I$(brew --prefix libomp)/include -L$(brew --prefix libomp)/lib -lomp -o main main.cpp -lz -lzstd -O2 Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.cpp -Ispoa/include -IComplete-Striped-Smith-Waterman-Library/src -Lspoa/build/lib -LComplete-Striped-Smith-Waterman-Library/src -lspoa -lssw -pthread