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

# Compile

You need `meson`, `ninja`, and a modern C++ compiler that supports C++20.

```
mkdir build && cd build
meson --prefix ~/mytools
ninja
ninja install
```

# How to run

```
$ cas ../data/short.txt
| 20220420 10:47:13.291 | INFO | Number of sequence pairs : 24670
| 20220420 10:47:15.943 | INFO | WFA2 time xxx
| 20220420 10:47:18.152 | INFO | KSW2 time xxx

 ./cas ../data/long.txt
| 20220420 10:47:22.864 | INFO | Number of sequence pairs : 2301
| 20220420 10:47:25.763 | INFO | WFA2 time xxx
| 20220420 10:47:34.653 | INFO | KSW2 time xxx
```

# Results

Results are per sequence pair on average running the data 10-times:

|  Compiler   |    Options     | Dataset | WFA2 | KSW2 |
| ----------- | -------------- | ------- | ---- | ---- |
| Apple Clang | -              | Short   |      |      |
|             |                | Long    |      |      |
| GCC11       | -              | Short   |      |      |
|             |                | Long    |      |      |
| GCC11       | `march=native` | Short   |      |      |
|             |                | Long    |      |      |
