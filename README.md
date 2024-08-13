# tempograph

Various algorithms in rust for computing optimum-cost temporal walks in temporal graphs, and related measures.

## Install

Intall rust and then:
```
cargo build --release 
```

Example:
```
./target/release/tgraph 4_03_bordeaux.patg --command src-opt-cost --criterion S --beta 3000 --source 1 
```


## Usage

```
USAGE:
    tgraph [FLAGS] [OPTIONS] <input>

FLAGS:
    -e, --exact      Exact computation for betweenness with big rationals (otherwise approximate with floats)
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    Log information about progress

OPTIONS:
    -b, --beta <beta>               Maximum waiting time: a temporal walk is restless if the maximum time between the
                                    arrival time of an edge and the departure time of the next edge in the walk is
                                    bounded by a value `beta`, non-restless is indicated by value `-1` [default: -1]
    -c, --command <command>         Computation to perform:
                                    `soc` or `src-opt-cost` for single-source optimal temporal walks:
                                       compute optimal walks from a source `s` (set with `-source`)
                                       and output for each node `t` the cost (see `-criterion`) of
                                       an optimal `st`-walk.
                                    `sc` or `shortest-closeness` for harmonic closeness of all nodes:
                                       the harmonic closeness of `s` is `\sum_{t != s} 1 / dist(s,t)`
                                       where `dist(s,t)` is the length in number of temporal edges
                                       of a shortest walk from `s` to `t` (shortest criterion).
                                    `tksc` or `top-k-shortest-closeness` for top-k harmonic closseness:
                                       same as `shortest-closeness`, but compute only the top `kappa` values
                                       (see `--kappa`).
                                    `c` or `closeness` for generalized harmonic closeness of all nodes:
                                       the harmonic closeness of `s` is `\sum_{t != s} 1 / dist_C(s,t)`
                                       where `dist_C(s,t)` is the length of an optimal walk from `s` to `t` 
                                       for criterion `C` used for optimality (see `-C`).
                                    `b` or `betweenness` for betweenness of all nodes:
                                       the betweenness of `v` is `\sum_{s,t != v} nwalks(s,v,t) / nwalks(s,t)`
                                       where `nwalks(s,t)` is the number of optimal walks from `s` to `t`
                                       and `nwalks(s,v,t)` is the number of optimal walks from `s` to `t` passing
                                       through `v` (the multiplicity of `v` in an optimal walk is counted).
                                    `rc` | `reachability` for reachability:
                                       the reachability of `s` is the number of nodes `t` such that there exists 
                                       an `st`-walk.
                                    `p` or `print` for temporal edges (sorted by arrival time). [default: print]
    -C, --criterion <criterion>     Criterion considered for optimal temporal walks: 
                                    `S` for shortest (in number of temporal edges),
                                    `Fo` for foremost (with earliest arrival time),
                                    `L` for latest departure time,
                                    `Fa` for fastest (minimum time difference between arrival and departure),
                                    `W` for minimum (overall) waiting time,
                                    `DS` for minimum sum of delays,  
                                    `SFo`, `SL`, `SFa`, `SW`, `SDS` for shortest foremost (shortest walks among 
                                    foremost ones), shortest latest, shortest fastest, shortest with minimum 
                                    waiting time, shortest with minimum sum of delays. [default: S]
    -k, --kappa <kappa>             Number of top nodes for top-k computation (see `top-k-shortest-closeness` command)
                                    [default: 100]
    -n, --num-threads <nthreads>    Number of threads to use (closeness from different source nodes are computed in
                                    parallel) [default: 3]
    -s, --source <source>           Source node (with commands considering a source node such as `min-cost-from`):
                                    consider only walks from that node [default: 1]

ARGS:
    <input>    Input file: a temporal graph in the following format:
               an optional first line with node maximum number `n` (nodes are numbered from `0`
               to `n`), followed by a quadruple `u v t delay` of unsigned integers per line for
               each temporal edge from node `u` to node `v` at time `t` with travel time `delay`.
```

