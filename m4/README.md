Our benchmarking code is mostly borrowed from [pqm4](https://github.com/mupq/pqm4) and the implementations in `crypto_kem` can simply be plugged into pqm4. For setting things up, please refer to [pqm4](https://github.com/mupq/pqm4).

We include implementations for the schemes
 - `lightsaber`
 - `saber`
 - `firesaber`
 - `ntruhps2048509`
 - `ntruhps2048677`
 - `ntruhrss701`
 - `ntruhps4096821`
 - `lac-128-v3a`
 - `lac-192-v3a`
 - `lac-256-v3a`


Our implementations are contained in the `m4` sub-directories under `crypto_kem/{scheme}/`. For convenience, the `toom` sub-directories contain the previously the fastest implementations which use Toom-Cook multiplications.

Similar to pqm4, this build system allows to build three different binaries for each implementation:
- `test`: for basic functional testing
- `testvectors`: for verifying that the provided implementations produce the same testvectors as ones submitted to NIST.
- `speed`: for benchmarking

You can build binaries by running
`make IMPLEMENTATION_PATH=crypto_kem/{scheme}/{implementation} bin/crypto_kem_{scheme}_{implementation}_{test}.bin`

Hence, the testing binaries can be built with
```
make IMPLEMENTATION_PATH=crypto_kem/lightsaber/m4 bin/crypto_kem_lightsaber_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/saber/m4 bin/crypto_kem_saber_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/firesaber/m4 bin/crypto_kem_firesaber_m4_test.bin

make IMPLEMENTATION_PATH=crypto_kem/ntruhps2048509/m4 bin/crypto_kem_ntruhps2048509_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhps2048677/m4 bin/crypto_kem_ntruhps2048677_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhrss701/m4 bin/crypto_kem_ntruhrss701_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhps4096821/m4 bin/crypto_kem_ntruhps4096821_m4_test.bin

make IMPLEMENTATION_PATH=crypto_kem/lac-128-v3a/m4 bin/crypto_kem_lac-128-v3a_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/lac-192-v3a/m4 bin/crypto_kem_lac-192-v3a_m4_test.bin
make IMPLEMENTATION_PATH=crypto_kem/lac-256-v3a/m4 bin/crypto_kem_lac-256-v3a_m4_test.bin

```

Similar for the testvector binaries
```
make IMPLEMENTATION_PATH=crypto_kem/lightsaber/m4 bin/crypto_kem_lightsaber_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/saber/m4 bin/crypto_kem_saber_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/firesaber/m4 bin/crypto_kem_firesaber_m4_testvectors.bin

make IMPLEMENTATION_PATH=crypto_kem/ntruhps2048509/m4 bin/crypto_kem_ntruhps2048509_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhps2048677/m4 bin/crypto_kem_ntruhps2048677_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhrss701/m4 bin/crypto_kem_ntruhrss701_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhps4096821/m4 bin/crypto_kem_ntruhps4096821_m4_testvectors.bin

make IMPLEMENTATION_PATH=crypto_kem/lac-128-v3a/m4 bin/crypto_kem_lac-128-v3a_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/lac-192-v3a/m4 bin/crypto_kem_lac-192-v3a_m4_testvectors.bin
make IMPLEMENTATION_PATH=crypto_kem/lac-256-v3a/m4 bin/crypto_kem_lac-256-v3a_m4_testvectors.bin

```

and benchmarks
```
make IMPLEMENTATION_PATH=crypto_kem/lightsaber/m4 bin/crypto_kem_lightsaber_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/saber/m4 bin/crypto_kem_saber_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/firesaber/m4 bin/crypto_kem_firesaber_m4_speed.bin

make IMPLEMENTATION_PATH=crypto_kem/ntruhps2048509/m4 bin/crypto_kem_ntruhps2048509_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhps2048677/m4 bin/crypto_kem_ntruhps2048677_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhrss701/m4 bin/crypto_kem_ntruhrss701_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/ntruhps4096821/m4 bin/crypto_kem_ntruhps4096821_m4_speed.bin

make IMPLEMENTATION_PATH=crypto_kem/lac-128-v3a/m4 bin/crypto_kem_lac-128-v3a_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/lac-192-v3a/m4 bin/crypto_kem_lac-192-v3a_m4_speed.bin
make IMPLEMENTATION_PATH=crypto_kem/lac-256-v3a/m4 bin/crypto_kem_lac-256-v3a_m4_speed.bin
```

For convenience, we also provide scripts that run tests (`tests.py`), check testvectors (`testvectors.py`), and run benchmarks (`benchmarks.py`) for all the implementations.