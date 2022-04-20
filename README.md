<img src="img/wfa-logo.png" alt="WFA logo" width="120px" align="right"/>
<h1 align="center">Aligner testbed</h1>
<p align="center">Testbed for different aligners on CLR data</p>

***
# Goal

Goal is to benchmark KSW2 vs WFA2 on realistic PacBio CLR data. We mimic one
alignment step in the CCS algorithm, here the [draft
step](https://ccs.how/how-does-ccs-work.html#2-draft-generation). For this, all
subreads from one ZMW are aligned against a representative subread called
backbone. We perform mapping / seeding and extract alignment regions. The
histogram of those alignment regions has a sharp peak at ~250 bp, as we require
at least 200 bp long regions, but some are much longer due to the nature of
single-molecule sequencing.

# Data

Provided are two data files from a human sample from a few ZMWs, containing
alignment regions that go into the base-level alignment step:

 * `data/short.txt`, sequences with a max length of 500bp
 * `data/long.txt`, sequences with a lengths greater 500bp

# Algorithms

Three algorithms are currently available:

 * KSW2 with `ksw_extd2_sse41`
 * WFA2-lib with `WFAlignerGapAffine2Pieces`
 * miniwfa

# Compile

You need `meson`, `ninja`, `boost`, and a modern C++ compiler that supports C++20.

```
mkdir build && cd build
meson --prefix ~/mytools
ninja
ninja install
```

For debug builds
```
meson --buildtype debug
ninja
./at
```

# How to run

```
$ at ../data/short.txt
| 20220420 14:07:59.594 | INFO | Number of sequence pairs : xxx
| 20220420 14:07:59.658 | INFO | MWFA time xxx
| 20220420 14:07:59.732 | INFO | WFA2 time xxx
| 20220420 14:07:59.736 | INFO | KSW2 time xxx
```

Optionally, adjust the number of rounds to average timings and deactivate algos:
```
$ ./at ../data/short.txt --miniwfa=true --wfa2=false --ksw2=true --rounds 100
```

# Results

Results are per sequence pair on average running the data 10-times:

|  Compiler   |    Options     | Dataset | WFA2 | MWFA | KSW2 |
| ----------- | -------------- | ------- | ---- | ---- | ---- |
| Apple Clang | -              | Short   |      |      |      |
|             |                | Long    |      |      |      |
| GCC11       | -              | Short   |      |      |      |
|             |                | Long    |      |      |      |
| GCC11       | `march=native` | Short   |      |      |      |
|             |                | Long    |      |      |      |
