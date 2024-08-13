use crate::tgraph::*;
use crate::tsweep::*;
use crate::cost::*;
use crate::tbfs::*;

use dsi_progress_logger::*;
use std::time::Instant;

pub fn closeness<C: Cost>(tg: &TGraph, beta: Time) -> Vec<f64> {
    let start = Instant::now();
    let mut tsweep: TSweep<C> = TSweep::new(&tg);
    let mut pl = ProgressLogger::default();
    pl.expected_updates(Some(tg.n as usize))
        .display_memory(true)
        .item_name("node");
    pl.start("Computing closeness...");
    (0..tg.n).map(|s| -> f64 {
        tsweep.clear();
        tsweep.scan(s, beta);
        let hc = harmonic_centrality::<C>(s, &tsweep.opt_costs(s));
        pl.update();
        if s == tg.n - 1 { pl.stop(); let duration = start.elapsed(); println!("Time elapsed: {:?}", duration); }
        hc
    }).collect()
}

fn harmonic_centrality<C: Cost>(s: Node, dist: &Vec<C::TargetCost>) -> f64 {
    let mut hc = 0.;
    let infty = C::infinite_target_cost();
    for t in 0..dist.len() {
        if t != s && dist[t] < infty { hc += 1. / C::target_cost_to_distance(dist[t].clone()) as f64 }
    }
    hc
}

pub fn shortest_closeness(tg: &TGraph, beta: Time) -> Vec<f64> {
    let start = Instant::now();
    let succ = tg.extend_indexes(beta);
    let mut hc = vec![0.; tg.n];
    let mut tbfs = TBFS::new(tg);
    for s in 0..tg.n {
        if beta == Time::MAX { tbfs.tbfs_inf(tg, &succ, s) } else { tbfs.tbfs(tg, &succ, s) };
        hc[s] = harmonic_centrality::<Shortest>(s, &tbfs.u_hop);
    }
    let duration = start.elapsed(); println!("Time elapsed: {:?}", duration);
    hc
}

pub fn top_shortest_closeness(tg: &TGraph, beta: Time) -> f64 {
    let mut deg_ord: Vec<Node> = (0..tg.n).collect();
    deg_ord.sort_by(|&a, &b| (tg.u_fst[b+1]-tg.u_fst[b]).cmp(&(tg.u_fst[a+1]-tg.u_fst[a])));
    let start = Instant::now();
    let succ = tg.extend_indexes(beta);
    let mut hc_max = 0.;
    let mut tbfs = TBFS::new(tg);
    for &s in &deg_ord {
        let hc =
            if beta == Time::MAX { tbfs.tbfs_inf_prune(tg, &succ, s, hc_max) } 
            else { tbfs.tbfs_prune(tg, &succ, s, hc_max) };
        if hc > hc_max { hc_max = hc }
    }
    let duration = start.elapsed(); println!("Time elapsed: {:?}", duration);
    hc_max
}

pub fn top_k_shortest_closeness(tg: &TGraph, beta: Time, mut k: Node) -> Vec<Node> {
    let mut deg_ord: Vec<Node> = (0..tg.n).collect();
    deg_ord.sort_by(|a, b| (tg.u_fst[b+1]-tg.u_fst[*b]).cmp(&(tg.u_fst[a+1]-tg.u_fst[*a])));
    let start = Instant::now();
    let succ = tg.extend_indexes(beta);
    if k > tg.n { k = tg.n-1; }
    let mut tbfs = TBFS::new(tg);
    let mut hc = vec![0.; k+1];
    let mut top = vec![0; k+1];
    for s in 1..k+1 {
        let ss = deg_ord[s];
        if beta == Time::MAX { tbfs.tbfs_inf(tg, &succ, ss) } else { tbfs.tbfs(tg, &succ, ss) };
        let dist = & tbfs.u_hop;
        hc[s-1] = harmonic_centrality::<Shortest>(ss, dist);
        top[s-1] = ss;
        let mut cur = s-1;
        while cur > 0 && hc[cur-1] < hc[cur] {
            let a = hc[cur];
            let b = top[cur];
            hc[cur] = hc[cur-1];
            top[cur] = top[cur-1];
            hc[cur-1] = a;
            top[cur-1] = b;
            cur -= 1;
        } 
    }
    for s in k+1..tg.n {
        let ss = deg_ord[s];
        hc[k] = 
            if beta == Time::MAX { tbfs.tbfs_inf_prune(tg, &succ, ss, hc[k-1]) } 
            else { tbfs.tbfs_prune(tg, &succ, ss, hc[k-1]) };
        top[k] = ss;
        let mut cur = k;
        while cur > 0 && hc[cur-1] < hc[cur] {
            let a = hc[cur];
            let b = top[cur];
            hc[cur] = hc[cur-1];
            top[cur] = top[cur-1];
            hc[cur-1] = a;
            top[cur-1] = b;
            cur -= 1;
        } 
    }
    let duration = start.elapsed(); println!("Time elapsed: {:?}", duration);
    top
}

use std::thread;
use std::sync::{Arc, Mutex, mpsc};

pub fn top_shortest_closeness_par(tg: &TGraph, beta: Time, nthread: u32) -> f64 {
    let succ = tg.extend_indexes(beta);
    let succ_ref = & succ;
    let n_th_dft = std::thread::available_parallelism()
        .unwrap_or(std::num::NonZeroUsize::new(2).unwrap()).get() as u32;
    let nthread = if nthread > 0 { nthread } else { n_th_dft };
    log::info!("Using {nthread} threads.");
    let hc_max_all = Arc::new(Mutex::new(0.));
    thread::scope(|s| {
        for i_th in 0..nthread as Node {
            let hc_max_all = Arc::clone(&hc_max_all);
            s.spawn(move || {
                let mut tbfs = TBFS::new(tg);
                let mut hc_max = 0.;
                for s in 0..tg.n  {
                    if s % nthread as Node == i_th {
                        let hc_s = 
                            if beta == Time::MAX { tbfs.tbfs_inf_prune(tg, succ_ref, s, hc_max) } 
                            else { tbfs.tbfs_prune(tg, succ_ref, s, hc_max) };
                        if hc_s > hc_max {
                            hc_max = hc_s;
                            if hc_max > *hc_max_all.lock().unwrap() {
                                let mut hc_max_all = hc_max_all.lock().unwrap();
                                *hc_max_all = hc_max;
                            }
                        }
                    }
                }
            });
        }
    });
    let hc = *hc_max_all.lock().unwrap();
    hc
}

pub fn closeness_par<C: Cost>(tg: &TGraph, beta: Time, nthread: u32) -> Vec<f64> {
    let n_th_dft = std::thread::available_parallelism()
        .unwrap_or(std::num::NonZeroUsize::new(2).unwrap()).get() as u32;
    let nthread = if nthread > 0 { nthread } else { n_th_dft };
    log::info!("Using {nthread} threads.");
    let (tx, rx) = mpsc::channel();
    thread::scope(|s| {
        for i_th in 0..nthread as Node {
            let tx = tx.clone();
            //let tg = tg.clone();
            s.spawn(move || {
                let mut plopt = None;
                if i_th == 0 {
                    let mut pl = ProgressLogger::default();
                    pl.expected_updates(Some(tg.n as usize));
                    pl.display_memory(true);
                    pl.item_name("node");
                    pl.start("Computing closeness...");
                    plopt = Some(pl);
                }
                let mut tsweep: TSweep<C> = TSweep::new(&tg);
                let mut hc = vec![0.; tg.n];
                for s in 0..tg.n  {
                    if s % nthread as Node == i_th {
                        tsweep.clear();
                        tsweep.scan(s, beta);
                        hc[s] = harmonic_centrality::<C>(s, &tsweep.opt_costs(s))
                    }
                    if i_th == 0 { plopt.as_mut().expect("no logger").update() }
                }
                if i_th == 0 { plopt.unwrap().stop() }
                tx.send(hc).unwrap();
            });
        }
    });
    let mut hc = vec![0.; tg.n];
    for _ in 0..nthread {
        let h_th = rx.recv().unwrap();
        for v in 0..tg.n as usize {
            hc[v] += h_th[v];
        }
    }
    hc
}

pub fn reachability_par(tg: &TGraph, beta: Time, nthread: u32) -> Vec<Node> {
    let n_th_dft = std::thread::available_parallelism()
        .unwrap_or(std::num::NonZeroUsize::new(2).unwrap()).get() as u32;
    let nthread = if nthread > 0 { nthread } else { n_th_dft };
    log::info!("Using {nthread} threads.");
    let succ = tg.extend_indexes(beta);
    let rsucc = &succ;
    let (tx, rx) = mpsc::channel();
    thread::scope(|s| {
        for i_th in 0..nthread as Node {
            let tx = tx.clone();
            //let tg = tg.clone();
            s.spawn(move || {
                let mut plopt = None;
                if i_th == 0 {
                    let mut pl = ProgressLogger::default();
                    pl.expected_updates(Some(tg.n as usize));
                    pl.display_memory(true);
                    pl.item_name("node");
                    pl.start("Computing reachability...");
                    plopt = Some(pl);
                }
                let mut rch = vec![0 as Node; tg.n];
                let mut tbfs = TBFS::new(tg);
                for s in 0..tg.n {
                    if beta == Time::MAX { tbfs.tbfs_inf(tg, rsucc, s) } else { tbfs.tbfs(tg, rsucc, s) };
                    rch[s] = tbfs.u_hop.iter().filter(|&&h| h < Hop::MAX).count();
                    if i_th == 0 { plopt.as_mut().expect("no logger").update() }
                }
                if i_th == 0 { plopt.unwrap().stop() }
                tx.send(rch).unwrap();
            });
        }
    });
    let mut rch = vec![0 as Node; tg.n];
    for _ in 0..nthread {
        let r_th = rx.recv().unwrap();
        for v in 0..tg.n as usize {
            if r_th[v] > 0 { rch[v] = r_th[v]; }
        }
    }
    rch
}
