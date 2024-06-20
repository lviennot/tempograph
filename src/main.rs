extern crate log;

pub mod graph;
pub mod tgraph;
pub mod cost;
pub mod tbetweenness;

use tgraph::*;
use tbetweenness::*;

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

    /// Number of threads to use (betweenness from different source nodes
    /// are computed in parallel).
    #[structopt(short, long = "num-threads", default_value = "3")]
    nthreads: u32,

    /// Computation to perform:
    /// `b` or `betweenness` for betweenness values of all of nodes,
    /// `sb` or `src-betweenness` for betweenness values of all nodes and
    /// temporal edges restricted to walks from a given source (see `--source`),
    /// `p` or `print` for temporal edges (sorted by arrival time).
    #[structopt(short, long, default_value = "betweenness", verbatim_doc_comment)]
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

    /// Source node (with command `src-betweenness`):
    /// consider only walks from that node.
    #[structopt(short, long, default_value = "1")]
    source: Node,

    /// Precision: the computation of betweenness is done by computing the number
    /// optimal walks (for the choson criterion) between any pair of nodes.
    /// As such numbers can get very high, it may lead to overflow depending on
    /// the types chosen for storing integers and rationals. Smaller types lead
    /// to faster computation, but may fail due to overflow. Also, accumulation of
    /// imprecision with floats can lead to negative values. Possible choices are
    /// `approx` for f64 and f64, `low` for u64 and f64, `medium` for u128 and f64, 
    /// `high` for BigUInt and f256:f256, `exact` for BigUint and Ratio<BigUint> 
    /// (this latter choice provides exact computation). 
    #[structopt(short, long, default_value = "medium")]
    precision: String,

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
            for i in 0..tg.m as usize {
                println!("{}", &tg.earr[i]);
            }
        },

        "edep" => {
            for j in 0..tg.m as usize {
                let i = tg.edep[j];
                println!("{j} ({i}): {}", &tg.earr[i as usize]);
            }
        },

        "b" | "betweenness" | "sb" | "src-betweenness" | "berr" => {
            match opt.criterion.as_str() {
                "foremost" | "Fo" => precision::<cost::Foremost>(&tg, &opt),
                "latest" | "L" => precision::<cost::Latest>(&tg, &opt),
                "fastest" | "Fa" => precision::<cost::Fastest>(&tg, &opt),
                "waiting" | "W" => precision::<cost::Waiting>(&tg, &opt),
                "shortest" | "S" => precision::<cost::Shortest>(&tg, &opt),
                "shortest-foremost" | "SFo" => precision::<cost::ShortestForemost>(&tg, &opt),
                "shortest-latest" | "SL" => precision::<cost::ShortestLatest>(&tg, &opt),
                "shortest-fastest" | "SFa" => precision::<cost::ShortestFastest>(&tg, &opt),
                _ => panic!("Cost '{}' not suppoted.", opt.criterion)
            }
        },

        _ => panic!("Unkown command '{}'.", opt.command)
    };

}

/* rug::Rational appears to be faster than Ratio<BigUInt>,
 * but rug::Integer is slower than BigUint...    
 */
use num_bigint::{BigUint, BigInt};
use num_rational::Ratio;
use num_traits::ToPrimitive;
use rug::{Integer, Rational};

impl RatFrom<BigUint> for Ratio<BigUint> {
    fn rat_from(i: BigUint) -> Self { i.into() }
}

impl RatFrom<BigInt> for Ratio<BigInt> {
    fn rat_from(i: BigInt) -> Self { i.into() }
}

impl RatFrom<BigUint> for f64 {
    fn rat_from(i: BigUint) -> Self { i.to_f64().unwrap() }
}

impl RatFrom<f64> for f64 {
    fn rat_from(i: f64) -> Self { i }
}

fn precision<C : cost::Cost>(tg: &TGraph, opt: &Opt) {
    match opt.precision.as_str() {
        "approx" => command::<C, NumT<f64, f64, False, False>>(tg, opt),
        "low" => command::<C, NumT<u64, f64, True, False>>(tg, opt),
        "medium" => command::<C, NumT<u128, f64, True, False>>(tg, opt),
        "high-rug" => command::<C, NumT<Integer, f64, True, False>>(tg, opt),
        "high" => command::<C, NumT<BigUint, f64, True, False>>(tg, opt),
        "exact-num" => command::<C, NumT<BigUint, Ratio<BigUint>, True, True>>(tg, opt),
        "exact" => command::<C, NumT<Integer, Rational, True, True>>(tg, opt),
        _ => panic!("Unkown precision '{}'.", opt.precision)
    }
}

fn command<C : cost::Cost, N : Num>(tg: &TGraph, opt: &Opt) {

    let beta = if opt.beta >= 0 { Time::try_from(opt.beta).expect("unexpected beta value") } else { Time::MAX };

    match opt.command.as_str() {

        "b" | "betweenness" if opt.nthreads == 1 => {
            let betw = betweenness_seq::<C, N>(tg, beta);
            for v in 0..=tg.n as usize {
                //println!("{:.4} {}", betw[v].to_f64().unwrap(), betw[v])
                println!("{}", betw[v])
            }
        },

        "b" | "betweenness" => {
            let betw = betweenness_par::<C, N>(tg, beta, opt.nthreads);
            for v in 0..=tg.n as usize {
                //println!("{:.4} {}", betw[v].to_f64().unwrap(), betw[v])
                println!("{}", betw[v])
            }
        },

        "berr" => {
            let betw_aprx = betweenness_par::<C, N>(tg, beta, opt.nthreads);
            let betw_exct = betweenness_par::<C, NumT<Integer, Rational, True, True>>(tg, beta, opt.nthreads);
            let mut nb_diff = 0 as Node;
            let mut max_diff: Rational = Integer::from(0_u32).into();
            let mut max_diff_exact: Rational = Integer::from(0_u32).into();
            let mut max_rel: Rational = Integer::from(0_u32).into();
            let one: Rational = Integer::from(1_u32).into();
            let epsilon : Rational = Rational::from_f64(1e-6_f64).unwrap();
            for v in 0..=tg.n as usize {
                let apx: f64 = betw_aprx[v].to_string().as_str().parse::<f64>().unwrap();
                let apx = Rational::from_f64(apx).unwrap();
                let diff: Rational = betw_exct[v].clone() - apx;
                println!("{v} {}", diff.to_f64());
                let diff = diff.abs();
                if diff > epsilon {
                    nb_diff += 1;
                    if betw_exct[v] > one {
                        let rel = &diff / betw_exct[v].clone();
                        if rel > max_rel { max_rel = rel }
                    }
                }
                if diff > max_diff { 
                    max_diff = diff.clone();
                    max_diff_exact = betw_exct[v].clone();
                }
            }   
            log::info!("nb_diff > {}: {}, max_diff: {:.4e} (for value {:.4e}), max_rel(>1): {:+e}", 
                epsilon.to_f64(), nb_diff, 
                max_diff.to_f64(), max_diff_exact.to_f64(), max_rel.to_f64())
        },

        "sb" | "src-betweenness" => src_betweenness::<C, N>(tg, beta, opt.source),

        _ => panic!("Unkown command '{}'.", opt.command)
    }
}

fn src_betweenness<C: cost::Cost, N: Num>(tg: &TGraph, beta: Time, src: Node) {

    let mut sweep: TGraphSweep<C, N> = TGraphSweep::new(&tg);
    
    for s in src..=src as Node {
        if s > tg.n { break }
        //eprintln!(" --- {s} ---------------- ");
        sweep.clear();
        sweep.forward(s, beta);
        sweep.backward(s);
        /*
        let costs = sweep.opt_costs();
        let walks = sweep.opt_walk_counts();    
        let s_betw = sweep.nodes_betweenness();
        for v in 1..=tg.n as usize {
            eprintln!("{s} {v} {} {} {}", costs[v], walks[v], s_betw[v])
        }
        */
    
    }
    let costs = sweep.opt_costs();
    let walks = sweep.opt_walk_counts();
    let e_betw = sweep.edges_betweenness();
    //let betw: Vec<num::Rat> = sweep.betweenness();

    for v in 0..=tg.n as usize {
        //let b: f64 = betw[v].into();
        println!("{v} {:?} {}", costs[v], walks[v])
    }

    println!(" --- "); 
    for i in 0..tg.m as usize {
        println!("{} [{}] {}", i+1, tg.earr[i], e_betw[i]);
    }

}
