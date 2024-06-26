extern crate log;

pub mod graph;
pub mod tgraph;
pub mod cost;
pub mod tsweep;
pub mod tcloseness;

use tcloseness::*;
use tgraph::*;
use tsweep::TSweep;

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
    /// `soc` or `src-opt-cost` for single source optimum cost temporal walks:
    ///    compute optimum-cost walks from a source `s` (set with `-source`)
    ///    and output for each node `t` the cost (see `-criterion`) of
    ///    an optimum-cost st-walk.
    ///  `c` or `closeness` for harmonic closeness of all nodes:
    ///    the harmonic closeness of `s` is `\sum_{t != s} 1 / dist(s,t)`
    ///    where `dist(s,t)` is the length in number of temporal edges
    ///    of a shortest walk from `s` to `t` (shortest criterion).
    /// `p` or `print` for temporal edges (sorted by arrival time).
    #[structopt(short, long, default_value = "print", verbatim_doc_comment)]
    command: String,

    /// Criterion considered for optimal temporal walks: 
    /// `S` for shortest (in number of temporal edges),
    /// `Fo` for foremost (with earliest arrival time),
    /// `L` for latest departure time,
    /// `Fa` for fastest (minimum difference between arrival and departure times),
    /// `SFo`, `SL`, `SFa` for shortest foremost (shortest walks among foremost 
    /// ones), shortest latest, shortest fastest,
    /// `W` for minimum waiting time.  
    #[structopt(short = "C", long, default_value = "S", verbatim_doc_comment)]
    criterion: String,

    /// Maximum waiting time: a temporal walk is restless if the maximum time between 
    /// the arrival time of an edge and the departure time of the next is bounded by
    /// a value `beta`, non-restless is indicated by value `-1`.
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
    /// an optional first line with node maximum number `n` 
    /// (nodes are numbered from `0` to `n`), followed by a quadruple
    /// `u v t delay` of unsigned integers per line for each temporal edge from node `u` 
    /// to node `v` at time `t` with travel time `delay`.
    #[structopt(parse(from_os_str))]
    input: PathBuf,

}

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

    match opt.command.as_str() {

        "size" | "sz" => {
            println!("{} {}", tg.n, tg.m);
        }

        "earr" | "p" | "print" => {
            for &i in tg.earr.iter() {
                println!("{}", &tg.edep[i]);
            }
        },

        "edep" => {
            for e in tg.edep.iter() {
                println!("{}", e);
            }
        },

        _ => {
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

    let beta = if opt.beta == -1 { Time::MAX } else { Time::try_from(opt.beta).expect("unexpected beta value") };
    log::info!("use of beta={beta}");
 
    let mut tsweep: TSweep<C> = TSweep::new(&tg);

    match opt.command.as_str() {

        "soc" | "src-opt-cost" => {
            tsweep.scan(opt.source, beta);
            let opt_costs = tsweep.opt_costs(opt.source);
            for c in opt_costs {
                println!("{:?}", c);
            }
        },

        "c" | "closeness" => {
            let hc = if opt.nthreads == 1 { closeness::<C>(&tg, beta) } else { closeness_par::<C>(&tg, beta, opt.nthreads) };
            for c in hc { println!("{}", c); }
        },

        _ => panic!("Unkown command '{}'.", opt.command)
    }

}
