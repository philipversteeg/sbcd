#!/bin/bash
# Copyright (c) 2021-2022, Philip Versteeg (p.j.j.p.versteeg@gmail.com). All rights reserved.
# Use of this source code is governed by a BSD-style license that can be found in the LICENSE file.

# logs
logdir=logs
mkdir -p "$logdir"

# setup
# Fig. 5b: fixed graph, stats
e=fixed nseed=20 N=10000 ct=stats
Rscript simulate.R "$e" "$nseed" "$N" "$ct" > "$logdir"/sim_"$e"_"$nseed"_"$N"_"$ct".log

# Fig. 6a: small graphs, ancestors
e=small nseed=20 N=10000 ct=ancestors
Rscript simulate.R "$e" "$nseed" "$N" "$ct" > "$logdir"/sim_"$e"_"$nseed"_"$N"_"$ct".log

# Fig. 6b: large graphs, ancestors
e=large nseed=20 N=10000 ct=ancestors
Rscript simulate.R "$e" "$nseed" "$N" "$ct" > "$logdir"/sim_"$e"_"$nseed"_"$N"_"$ct".log

# Fig. 8a: small graphs, pattern
e=small nseed=20 N=10000 ct=pattern
Rscript simulate.R "$e" "$nseed" "$N" "$ct" > "$logdir"/sim_"$e"_"$nseed"_"$N"_"$ct".log

# Fig. 8b: large graphs, pattern
e=large nseed=20 N=10000 ct=pattern
Rscript simulate.R "$e" "$nseed" "$N" "$ct" > "$logdir"/sim_"$e"_"$nseed"_"$N"_"$ct".log

# Fig. 9: small and large graphs, ancestors, varied n
nseed=20 ct=ancestors
for e in small large; do
    for N in 1000 5000 10000 20000; do
        Rscript simulate.R "$e" "$nseed" "$N" "$ct" > "$logdir"/sim_"$e"_"$nseed"_"$N"_"$ct".log
    done
done