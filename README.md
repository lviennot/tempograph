# tempograph

Various algorithms in rust for computing optimum-cost temporal walks in temporal graphs, and related measures.

## Install

Intall rust and then:
```
cargo build --release
```

Example:
```
./target/release/tgraph 4_03_bordeaux.patg -c src-opt-cost --criterion S --beta 3000 --source 1 
```


## Usage

```
Compute optimum-cost temporal walks in a temporal graph, and related subjects.

USAGE:
    tgraph [FLAGS] [OPTIONS] <input>

FLAGS:
    -h, --help       Prints help information
    -V, --version    Prints version information
    -v, --verbose    Log information about progress

OPTIONS:
    -b, --beta <beta>              Maximum waiting time: a temporal walk is restless if the maximum time between the
                                   arrival time of an edge and the departure time of the next is bounded by a value
                                   `beta`, non-restless is indicated by value `-1` [default: -1]
    -c, --command <command>        Computation to perform:
                                   `soc` or `src-opt-cost` for single source optimum cost temporal walks:
                                      compute optimum-cost walks from a source `s` (set with `-source`)
                                      and output for each node `t` the cost (see `-criterion`) of
                                      an optimum-cost st-walk.
                                    `c` or `closeness` for harmonic closeness of all nodes:
                                      the harmonic closeness of `s` is `\sum_{t != s} 1 / dist(s,t)`
                                      where `dist(s,t)` is the length in number of temporal edges
                                      of a shortest walk from `s` to `t` (shortest criterion).
                                   `p` or `print` for temporal edges (sorted by arrival time). [default: print]
    -C, --criterion <criterion>    Criterion considered for optimal temporal walks: 
                                   `S` for shortest (in number of temporal edges),
                                   `Fo` for foremost (with earliest arrival time),
                                   `L` for latest departure time,
                                   `Fa` for fastest (minimum difference between arrival and departure times),
                                   `SFo`, `SL`, `SFa` for shortest foremost (shortest walks among foremost 
                                   ones), shortest latest, shortest fastest,
                                   `W` for minimum waiting time.   [default: SFa]
    -s, --source <source>          Source node (with commands considering a source node such as `min-cost-from`):
                                   consider only walks from that node [default: 1]

ARGS:
    <input>    Input file: a temporal graph in the following format: an optional first line with node maximum number
               `n` (nodes are numbered from `0` to `n`), followed by a quadruple `u v t delay` of unsigned integers
               per line for each temporal edge from node `u` to node `v` at time `t` with travel time `delay`
```

