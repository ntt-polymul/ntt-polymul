This directory contains our multiplication implementations for the Saber, NTRU and LAC KEMs, which are using 16-bit NTT strategies that are optimized for amd64 CPUs supporting the AVX2 instruction set. Also included are the upstream multipliers to compare performance. More precisely, the code in the subdirectories `sabermul`, `ntrumul` and `lacmul` are copied from the official public-domain AVX2-optimized implementations of the Saber, NTRU and LAC KEM, respectively. On the other hand the code in the base directory constitutes our new multiplication code.

## Build Instructions

For running the range analysis, you will need to install the Berkeley Database Libraries (libdb-dev).
On Ubuntu/Debian, you can run
```
sudo apt install libdb-dev
```

The remainder of the source code in this directory is self-contained and does not require any external libraries. Just 
run `make` to build various test and benchmarking programs. The programs
``` 
./test_lightsabermul
./test_sabermul
./test_firesabermul
```
test the correctness of our multipliers for the Saber schemes by comparing random matrix-vector products over the polynomial ring in Saber against a naive multiplication implementation and the AVX2-optimized upstream multipliers. Similarly,
```
./test_ntruhps509mul
./test_ntruhps677mul
./test_ntruhps821mul
./test_ntruhrss701mul
```
test our NTRU multipliers by performing random polynomial products in the different NTRU rings, and
```
./test_lac128mul
./test_lac192mul
./test_lac256mul
```
do the same for the LAC multipliers.

## Benchmarking
We provide benchmarking programs for all of the multiplication implementations. These programs report the median and average cycle counts of 100 executions of the upstream multipliers and our NTT-based multipliers, and additionally of several inner functions. By default the time step counter (TSC) in the CPU is used. Hence the reported numbers are only actual cycle counts if the CPU runs at the constant frequency of the TSC. Alternatively, if the CPU is set to run at a constant frequency that is slightly different to the TSC frequency, the numbers can be converted to cycle counts. Note that the TSC frequency is often different to the rated nominal frequency of the CPU. In any case frequency scaling must be disabled in order for the reported numbers to be interpretable as cycle counts. The benchmarking programs also support the Performance Measurement Counters (PMC) that allow to measure actual cycle counts even when frequency scaling is used. For this the compiler flag `-DUSE_RDPMC` can be used when compiling the source.

For SABER there are the benchmarking programs
```
./test_speed256nx2
./test_speed256nx3
./test_speed256nx4
```
for LightSaber, Saber and FireSaber, respectively. For NTRU, the programs are
```
./test_speed1024
./test_speed1536
./test_speed1728
```
for ntruhps509, ntruhps677, ntruhrss701 and ntruhps821. Finally, for LAC, the programs are
```
./test_speed512n
./test_speed1024n
```
for LAC128, LAC192, and LAC256.

## Integration into the upstream KEM implementations
Our multiplication implementations can be compiled into static libraries by running `make lib`. The global symbols of the libraries are namespaced and the upstream KEM sources can be linked against them. The sub-directory `patches/` contains patches for the upstream sources that integrate our multipliers.

As an example, here is how to integrate our multiplier into Saber. The patch `patches/saber.patch` can be applied to the Saber implementation from the GitHub repository <https://github.com/KULeuven-COSIC/SABER>. It adds calls to our multiplication functions to the Saber code. They are used instead of the original multiplication functions when the macro `NTTMUL` is defined. Also the patch adds souces for two simple testing and benchmarking programs and corresponding make targets for every parameter set: `test_lightsaber`, `test_lightsaberspeed`, `test_saber`, `test_saberspeed`, `test_firesaber`, and `test_firesaberspeed`. To compile these programs the static libraries `liblightsabermul.a`, `libsabermul.a`, and `libfiresabermul.a` from this repository must be copied to the Saber repository, or their location must be added to the search path for the compiler by setting the environment variable `LDFLAGS` accordingly. Then the test and benchmarking programs can be compiled such that our NTT-based multiplier is used by adding `-DNTTMUL` to the compiler flags in the `EXTRAFLAGS` variable before running make.

