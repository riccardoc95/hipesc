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
mpic++ --std=c++17 -fopenmp -o main main.cpp -lz -lzstd -O2 Complete-Striped-Smith-Waterman-Library/src/ssw_cpp.cpp -Ispoa/include -IComplete-Striped-Smith-Waterman-Library/src -Lspoa/build/lib -LComplete-Striped-Smith-Waterman-Library/src -lspoa -lssw -pthread

srun --time=15:00 --nodes=3 --ntasks=3 --cpus-per-task=16 --mem=64G --nodelist=cn1a,cn1b,cn1c ./main ../CONSENT/data/chr1/chr1.fastq ../CONSENT/data/chr1/alignment_chr1.paf