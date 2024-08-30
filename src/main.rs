extern crate log;

pub mod graph;
pub mod tgraph;
pub mod cost;
pub mod tsweep;
pub mod tbetweenness;
pub mod tbfs;
pub mod tcloseness;
pub mod tbetw;

use tcloseness::*;
use tgraph::*;
use tsweep::TSweep;
use tbfs::TBFS;
use tbetweenness::*;

use std::fs::File;
use std::io::{BufRead, BufReader}; // use std::io::{self, prelude::*, BufReader};

use stderrlog;

use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
/// Compute optimum-cost temporal walks in a temporal graph, and related subjects.
/// 
#[structopt(name = "tempograph")]
struct Opt {
    /// Log information about progress.
    #[structopt(short, long)]
    verbose: bool,

    /// Computation to perform:
    /// `soc` or `src-opt-cost` for single-source optimal temporal walks:
    ///    compute optimal walks from a source `s` (set with `-source`)
    ///    and output for each node `t` the cost (see `-criterion`) of
    ///    an optimal `st`-walk.
    /// `sc` or `shortest-closeness` for harmonic closeness of all nodes:
    ///    the harmonic closeness of `s` is `\sum_{t != s} 1 / dist(s,t)`
    ///    where `dist(s,t)` is the length in number of temporal edges
    ///    of a shortest walk from `s` to `t` (shortest criterion).
    /// `tksc` or `top-k-shortest-closeness` for top-k harmonic closseness:
    ///    same as `shortest-closeness`, but compute only the top `kappa` values
    ///    (see `--kappa`).
    /// `c` or `closeness` for generalized harmonic closeness of all nodes:
    ///    the harmonic closeness of `s` is `\sum_{t != s} 1 / dist_C(s,t)`
    ///    where `dist_C(s,t)` is the length of an optimal walk from `s` to `t` 
    ///    for criterion `C` used for optimality (see `-C`).
    /// `b` or `betweenness` for betweenness of all nodes:
    ///    the betweenness of `v` is `\sum_{s,t != v} nwalks(s,v,t) / nwalks(s,t)`
    ///    where `nwalks(s,t)` is the number of optimal walks from `s` to `t`
    ///    and `nwalks(s,v,t)` is the number of optimal walks from `s` to `t` passing
    ///    through `v` (the multiplicity of `v` in an optimal walk is counted).
    /// `rc` | `reachability` for reachability:
    ///    the reachability of `s` is the number of nodes `t` such that there exists 
    ///    an `st`-walk.
    /// `p` or `print` for temporal edges (sorted by arrival time).
    #[structopt(short, long, default_value = "print", verbatim_doc_comment)]
    command: String,

    /// Criterion considered for optimal temporal walks: 
    /// `S` for shortest (in number of temporal edges),
    /// `Fo` for foremost (with earliest arrival time),
    /// `L` for latest departure time,
    /// `Fa` for fastest (minimum time difference between arrival and departure),
    /// `W` for minimum (overall) waiting time,
    /// `DS` for minimum sum of delays,  
    /// `SFo`, `SL`, `SFa`, `SW`, `SDS` for shortest foremost (shortest walks among 
    /// foremost ones), shortest latest, shortest fastest, shortest with minimum 
    /// waiting time, shortest with minimum sum of delays.
    #[structopt(short = "C", long, default_value = "S", verbatim_doc_comment)]
    criterion: String,

    /// Maximum waiting time: a temporal walk is restless if the maximum time between 
    /// the arrival time of an edge and the departure time of the next edge in the walk
    /// is bounded by a value `beta`, non-restless is indicated by value `-1`.
    #[structopt(short, long, default_value = "-1")]
    beta: i64,

    /// Source node (with commands considering a source node such as `min-cost-from`):
    /// consider only walks from that node.
    #[structopt(short, long, default_value = "1")]
    source: Node,

    /// Number of threads to use (closeness from different source nodes
    /// are computed in parallel).
    #[structopt(short, long = "num-threads", default_value = "3")]
    nthreads: u32,

    /// Input file: a temporal graph in the following format: 
    /// an optional first line with node maximum number `n` (nodes are numbered from `0`
    /// to `n`), followed by a quadruple `u v t delay` of unsigned integers per line for
    /// each temporal edge from node `u` to node `v` at time `t` with travel time `delay`.
    #[structopt(parse(from_os_str), verbatim_doc_comment)]
    input: PathBuf,

    /// Number of top nodes for top-k computation (see `top-k-shortest-closeness` command).
    #[structopt(short, long, default_value = "100")]
    kappa: usize,

    /// Exact computation for betweenness with big rationals (otherwise approximate with floats).
    #[structopt(short, long)]
    exact: bool, 
}

use std::io::{BufWriter, Write};

fn main() {
    let opt = Opt::from_args();

    if opt.verbose {
        stderrlog::new()
            .verbosity(2)
            //.timestamp(stderrlog::Timestamp::Second)
            .init().expect("log problem");
            log::info!("Options: {:?}", opt);
    }

    log::info!("Opening {:?} ...", opt.input); //.into_os_string().into_string()
    let file = File::open(&opt.input)
        .expect("Should have been able to read the file");
    let reader = BufReader::new(file);

    let tedges: Vec<tgraph::TEdge> = reader.lines()
        .enumerate()
        .filter(|(i, l)| -> bool {
            // skip first line if it is a single number:
            *i > 0 || l.as_ref().unwrap().parse::<Node>().is_err()
        }).map(|(i, l)| -> TEdge {
            let l = l.unwrap();
            match l.parse::<tgraph::TEdge>() {
                Err(err) => panic!("Error: {} at line {} : `{}`", err, i+1, l),
                Ok(e) => e                
            }
        }).collect();

    let tg = TGraph::new(tedges);
    log::info!("n={} M={}", tg.n, tg.m);

    // println! is not buffered and flushes after each line, use a BufWriter:
    let stdout = std::io::stdout();
    let lock = stdout.lock();
    let mut out = BufWriter::new(lock);
    // use `writeln!(out, "{}", smthg).expect("IO goes fine")` instead of `println!("{}", smthg)`

    match opt.command.as_str() {

        "size" | "sz" => {
            writeln!(out, "{} {}", tg.n, tg.m).expect("IO goes fine");
        }

        "earr" | "p" | "print" => {
            for &i in tg.earr.iter() {
                writeln!(out, "{}", &tg.edep[i]).expect("IO goes fine");
            }
        },

        "edep" => {
            for e in tg.edep.iter() {
                writeln!(out, "{}", e).expect("IO goes fine");
            }
        },

        "rc" | "reachability" => {
            let beta = if opt.beta == -1 { Time::MAX } else { Time::try_from(opt.beta).expect("unexpected beta value") };
            let rc = reachability_par(&tg, beta, opt.nthreads);
            for r in rc { writeln!(out, "{}", r).expect("IO goes fine"); }
        },

        "sc" | "shortest-closeness" => {
            let beta = if opt.beta == -1 { Time::MAX } else { Time::try_from(opt.beta).expect("unexpected beta value") };
            let hc = shortest_closeness(&tg, beta);
            for c in hc { writeln!(out, "{}", c).expect("IO goes fine"); }
        },

        "tsc" | "top-shortest-closeness" => {
            let beta = if opt.beta == -1 { Time::MAX } else { Time::try_from(opt.beta).expect("unexpected beta value") };
            let top_hc = if opt.nthreads == 1 { top_shortest_closeness(&tg, beta) } else { top_shortest_closeness_par(&tg, beta, opt.nthreads) };
            writeln!(out, "{}", top_hc).expect("IO goes fine");
        },

        "tksc" | "top-k-shortest-closeness" => {
            let beta = if opt.beta == -1 { Time::MAX } else { Time::try_from(opt.beta).expect("unexpected beta value") };
            let top_hc = top_k_shortest_closeness(&tg, beta, opt.kappa);
            writeln!(out, "{:?}", top_hc).expect("IO goes fine");
        },

        "test-clear" => {
            let mut tbfs = TBFS::new(&tg);
            for _ in 0..tg.n { tbfs.clear(); }
        }

        _ => {
            drop(out);
            match opt.criterion.as_str() {
                "foremost" | "Fo" => command::<cost::Foremost>(&tg, &opt),
                "latest" | "L" => command::<cost::Latest>(&tg, &opt),
                "fastest" | "Fa" => command::<cost::Fastest>(&tg, &opt),
                "delay-sum" | "DS" => command::<cost::DelaySum>(&tg, &opt),
                "waiting" | "W" => command::<cost::Waiting>(&tg, &opt),
                "shortest" | "S" => command::<cost::Shortest>(&tg, &opt),
                "shortest-foremost" | "SFo" => command::<cost::ShortestForemost>(&tg, &opt),
                "shortest-latest" | "SL" => command::<cost::ShortestLatest>(&tg, &opt),
                "shortest-fastest" | "SFa" => command::<cost::ShortestFastest>(&tg, &opt),
                "shortest-delay-sum" | "SDS" => command::<cost::ShortestDelaySum>(&tg, &opt),
                "shortest-waiting" | "SW" => command::<cost::ShortestWaiting>(&tg, &opt),
                _ => panic!("Cost '{}' not suppoted.", opt.criterion)
            }
        },

    };

}

fn command<C: cost::Cost>(tg: &TGraph, opt: &Opt) {

    let stdout = std::io::stdout();
    let lock = stdout.lock();
    let mut out = BufWriter::new(lock);

    let beta = if opt.beta == -1 { Time::MAX } else { Time::try_from(opt.beta).expect("unexpected beta value") };
    log::info!("use of beta={beta}");
 
    let mut tsweep: TSweep<C> = TSweep::new(&tg);

    match opt.command.as_str() {

        "soc" | "src-opt-cost" => {
            tsweep.scan(opt.source, beta);
            let opt_costs = tsweep.opt_costs(opt.source);
            for c in opt_costs {
                writeln!(out, "{:?}", c).expect("IO goes fine");
            }
        },

        "c" | "closeness" => {
            let hc = if opt.nthreads == 1 { closeness::<C>(&tg, beta) } else { closeness_par::<C>(&tg, beta, opt.nthreads) };
            for c in hc { writeln!(out, "{}", c).expect("IO goes fine"); }
        },

        "b" | "betweenness" => {
            if opt.exact {
                let bc = if opt.nthreads == 1 { betweenness_seq::<C, NumExact>(&tg, beta) } else { betweenness_par::<C, NumExact>(&tg, beta, opt.nthreads) };
                for c in bc { writeln!(out, "{}", c).expect("IO goes fine"); }
            } else {
                let bc = if opt.nthreads == 1 { betweenness_seq::<C, NumApprox>(&tg, beta) } else { betweenness_par::<C, NumApprox>(&tg, beta, opt.nthreads) };
                for c in bc { writeln!(out, "{}", c).expect("IO goes fine"); }
            }
        }

        "bd" | "betweenness-diff" => {
            let bc_ext = if opt.nthreads == 1 { betweenness_seq::<C, NumExact>(&tg, beta) } else { betweenness_par::<C, NumExact>(&tg, beta, opt.nthreads) };
            let bc_apx = if opt.nthreads == 1 { betweenness_seq::<C, NumApprox>(&tg, beta) } else { betweenness_par::<C, NumApprox>(&tg, beta, opt.nthreads) };
            assert_eq!(bc_ext.len(), bc_apx.len());
            for i in 0..bc_ext.len() {
                let ext = &bc_ext[i];
                let apx = bc_apx[i];
                let err = ext.to_f64() - apx;
                let err = err.abs();
                writeln!(out, "{}", err).expect("IO goes fine"); 
            }
        }

        _ => panic!("Unkown command '{}'.", opt.command)
    }

}
