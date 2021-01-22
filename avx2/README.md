
Run `make` in this directory to create a number of binaries. 
The ones used for benchmarking and reproducing the paper results are listed below, the others are testing the correctness of the multiplications. 

# Benchmarking Saber Multiplication
```
# Lightsaber
./test_speed256nx2

# Saber
./test_speed256nx3

# Firesaber
./test_speed256nx4
```

# Benchmarking NTRU Multiplication
```
# ntruhps2048509
./test_speed1024

# ntruhps2048677
./test_speed1536

# ntruhps4096821
./test_speed1728

# ntruhrss701
./test_speed1536
```

# Benchmarking LAC Multiplication
```
# LAC128
./test_speed512n
# LAC192
./test_speed1024n
# LAC256
./test_speed1024n
```
