extern crate log;

pub mod graph;
pub mod tgraph;
pub mod cost;
pub mod tsweep;
pub mod tcloseness;

use tgraph::*;
use tsweep::TSweep;

use std::fs::File;
use std::io::{BufRead, BufReader}; // use std::io::{self, prelude::*, BufReader};

use stderrlog;

use std::path::PathBuf;
use structopt::StructOpt;

#[derive(Debug, StructOpt)]
/// Compute temporal betweenness of all nodes of a temporal graph.
/// 
/// The betweenness of a node `v` evaluates the fraction of optimal walks
/// going through node `v`. Optimality refers to any criterion among shortest
/// (in hops), fastest, foremost,... (see option `--criterion`).
/// As the number of such walks can be very large, the default int type (UInt128)
/// might be too small and lead to a `attempt to add with overflow` error, use 
/// `--precision exact` in that case.
#[structopt(name = "twbc-rc")]
struct Opt {
    /// Log information about progress.
    #[structopt(short, long)]
    verbose: bool,

    /// Computation to perform:
    /// `mc` or `min-cost-from` for single source minimum cost walks:
    ///    compute minimum-cost walks from a source `s` (set with `-source`)
    ///    and output for each node `t` the cost (see `-criterion`) of
    ///    a minimum-cost st-walk.
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
    #[structopt(short = "C", long, default_value = "SFa", verbatim_doc_comment)]
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

        "mc" | "min-cost-from" => {
            match opt.criterion.as_str() {
                "foremost" | "Fo" => command::<cost::Foremost>(&tg, &opt),
                "latest" | "L" => command::<cost::Latest>(&tg, &opt),
                "fastest" | "Fa" => command::<cost::Fastest>(&tg, &opt),
                "waiting" | "W" => command::<cost::Waiting>(&tg, &opt),
                "shortest" | "S" => command::<cost::Shortest>(&tg, &opt),
                "shortest-foremost" | "SFo" => command::<cost::ShortestForemost>(&tg, &opt),
                "shortest-latest" | "SL" => command::<cost::ShortestLatest>(&tg, &opt),
                "shortest-fastest" | "SFa" => command::<cost::ShortestFastest>(&tg, &opt),
                _ => panic!("Cost '{}' not suppoted.", opt.criterion)
            }
        },

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


        _ => panic!("Unkown command '{}'.", opt.command)
    };

}

fn command<C: cost::Cost>(tg: &TGraph, opt: &Opt) {

    match opt.command.as_str() {

        "mc" | "min-cost-from" => {
            let mut tsweep: TSweep<C> = TSweep::new(&tg);
            tsweep.scan(opt.source, opt.beta);
            let opt_costs = tsweep.opt_costs();
            for c in opt_costs {
                println!("{:?}", c);
            }
        },

        _ => panic!("Unkown command '{}'.", opt.command)
    }

}