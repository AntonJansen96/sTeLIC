#!/usr/bin/env python3

from science.cphmd import plotdVdl

#! LYST BPSI VS ANTON'S REWEIGHTS
# plotdVdl(
#     [
#         "-10957.1 44076.7 -67242.6 48014.3 -16026. 2075.64 6.98266 -707.24 526.209",
#         "62240.648 -234764.266 368714.986 -313597.193 156104.864 -45638.868 7523.61 -1331.275 547.595"
#     ],
#     [
#         'orig LYST from BP SI',
#         'reweighted LYST'
#     ],
#     (-0.1, 1.1)
# )

#! LYST BPSI VS PAVEL'S REWEIGHTS
plotdVdl(
    [
        "-10957.1 44076.7 -67242.6 48014.3 -16026. 2075.64 6.98266 -707.24 526.209",
        "-10985.52 44105.334 -67223.723 48005.063 -16043.649 2075.13 18.132 -710.34 526.185"
    ],
    [
        'orig LYST from BP SI',
        'reweighted Pavel'
    ],
    (-0.1, 1.1)
)
