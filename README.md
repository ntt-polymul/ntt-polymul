This code accompanies the paper "NTT Multiplication for NTT-unfriendly Rings".

Authors:
    - Chi-Ming Marvin Chung `<marvin852316497@gmail.com>`
    - Vincent Hwang `<vincentvbh7@gmail.com>`
    - [Matthias J. Kannwischer](https://kannwischer.eu/) `<matthias@kannwischer.eu>`
    - Gregor Seiler `<gseiler@inf.ethz.ch>`
    - Cheng-Jhih Shi `<cs861324@gmail.com>`
    - [Bo-Yin Yang](https://www.iis.sinica.edu.tw/pages/byyang/) `<by@crypto.t>`


This repository contains our NTT-based implementations for Saber and NTRU for Cortex-M4 and AVX2.

Please clone this repository recursively to include [libopencm3](http://libopencm3.org/).
```
    git clone --recursive https://github.com/ntt-polymul/ntt-polymul
```

For details on how to build and use our Cortex-M4 see [m4/README.md](m4/README.md).
For our AVX2 implementations, see [avx2/README.md](avx2/README.md).

