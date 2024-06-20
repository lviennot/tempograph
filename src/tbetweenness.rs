use crate::tgraph::*;
use crate::cost::Cost;

use std::thread;
use std::sync::mpsc;
use std::ops::{AddAssign, SubAssign, Mul, Div};

use log::info;
use dsi_progress_logger::*;

pub trait Num {
    /// Type for large integers
    type Int : PartialOrd + Clone + AddAssign + SubAssign + Div<Output = Self::Int> + From<u32> + std::fmt::Display;
    /// Type for rationals
    type Rat : PartialOrd + Clone + AddAssign + SubAssign 
                + Mul<Output = Self::Rat> + Div<Output = Self::Rat> + Send + std::fmt::Display;
    fn i_zero() -> Self::Int;
    fn i_one() -> Self::Int;
    fn r_zero() -> Self::Rat;
    fn r_one() -> Self::Rat;
    fn real_of_int(i: Self::Int) -> Self::Rat; 
    fn i_is_exact() -> bool; // are integers exact?
    fn r_is_exact() -> bool; // are rationals exact? (vs approximated)
}

pub trait RatFrom<I> {
    fn rat_from(i: I) -> Self;
}

pub trait Bool {
    fn val () -> bool;
}
pub struct True;
impl Bool for True { fn val () -> bool { true } }
pub struct False;
impl Bool for False { fn val () -> bool { false } }


pub struct NumT<I, R, IE, RE> { _i: I, _r: R, _ie: IE, _re: RE }

impl<I, R, IE: Bool, RE: Bool> Num for NumT<I, R, IE, RE>
where 
I: PartialOrd + Clone + AddAssign + SubAssign + Div<Output = I> + From<u32> + std::fmt::Display,
R: PartialOrd + Clone + AddAssign + SubAssign + RatFrom<I>
+ Mul<Output = R> + Div<Output = R> + Send + std::fmt::Display
{
    type Int = I;
    fn i_zero() -> Self::Int { 0u32.into() } 
    fn i_one() -> Self::Int { 1u32.into() }
    type Rat = R; // type for large unsigned rational
    fn r_zero() -> Self::Rat { Self::Rat::rat_from(Self::i_zero()) } 
    fn r_one() -> Self::Rat { Self::Rat::rat_from(Self::i_one()) }
    fn real_of_int(i: Self::Int) -> Self::Rat { Self::Rat::rat_from(i) }
    fn i_is_exact() -> bool { IE::val() }
    fn r_is_exact() -> bool { RE::val() }
}


impl RatFrom<u64> for f64 {
    fn rat_from(i: u64) -> f64 { i as f64 }
}

impl RatFrom<u128> for f64 {
    fn rat_from(i: u128) -> f64 { i as f64 }
}


use rug::{Integer, Rational, Float};

impl RatFrom<Integer> for f64 {
    fn rat_from(i: Integer) -> f64 { i.to_f64() }
}

impl RatFrom<Integer> for Float {
    fn rat_from(i: Integer) -> Self { 
        let f = Float::with_val(96, 0.);
        f + i 
    }
}

impl RatFrom<Integer> for Rational {
    fn rat_from(i: Integer) -> Self { 
        i.into() 
    }
}


//type NumExact = NumT<Integer, Rational, True>;
//type NumU128F64 = NumT<u128, f64, False>;


/// All information required to perform a sweep (i.e. a traversal) of a temporal graph.
pub struct TGraphSweep<'tg, C : Cost, N : Num> {
    tg: &'tg TGraph,
    e_inf: Vec<TEdgeInfo<C, N>>, // aligned with earr
    e_succ: Vec<TEdgeSucc>, // aligned with earr
    e_betw: Vec<N::Rat>, // aligned with earr
    u_inf: Vec<NodeInfo<C, N>>,
}



#[derive(Clone, Debug)]
struct TEdgeInfo<C : Cost, N : Num> { // info about a tedge e={u,v,t,d}
    idep: Eind, // index in edep
    min_cost: C, // minimum cost of a walk ending with e
    walks: N::Int, // number of minimum cost walks ending with e
}

/// index range (in edep) of successors of a tedge,
/// f is a successor of e if a minimum-cost walk ending with f has e right before f.
struct TEdgeSucc { 
    left: Eind,
    right: Eind,
}


use std::collections::VecDeque;

struct NodeInfo<C : Cost, N : Num> { // info about tedges from a node u (forward phase)
    intervs: VecDeque<Interv<C, N>>, // information about consecutive intervals of edges from u (in edgep)
    preds: VecDeque<Eind>, // predecessor tedges (their index in earr)
    l: Eind, // left bound in edep, intervals in intervs span l..r
    r: Eind, // right bound in edep (excluded)
    refresh: Eind, // index where we can recompute betw from suming in l..r-1 for better numerical stability
    opt_cost: C::TargetCost, // minimum target cost of a walk from source to u
    opt_walks: N::Int, // number of such walks
    betw: N::Rat, // sum of betweenness share for tedges in edep[l..r-1], and then betweenness 
}

#[derive(Debug)]
struct Interv<C : Cost, N : Num> { // an interval of edges from a given node v in edep
    l: Eind, // left bound in edep
    r: Eind, // right bound in edep (excluded)
    cost: C, // minimum cost of walks reaching this interval (and that edges in the interval extend)
    walks: N::Int, // count of walks
    preds: Eind, // number of pred. in pred fifo (a predecessor is an edge that terminates a counted walk)
}

use core::cmp::{min, max};

impl<'tg, C : Cost, N : Num> TGraphSweep<'tg, C, N>{

    pub fn new(tg: &'tg TGraph) -> TGraphSweep<'tg, C, N> {
        assert!(u128::try_from(tg.m).unwrap() <= <u32 as Into<u128>>::into(Eind::MAX));

        let mut e_inf: Vec<TEdgeInfo<C, N>> = Vec::with_capacity(tg.m as usize);
        for _ in 0..tg.m {
            e_inf.push(TEdgeInfo {
                idep: 0,
                min_cost: Cost::infinite_cost(),
                walks: N::i_zero(), 
            });
        }
        for (i, &j) in tg.edep.iter().enumerate() {
            e_inf[j as usize].idep = i as Eind;
        }

        let mut e_succ: Vec<TEdgeSucc> = Vec::with_capacity(tg.m as usize);
        for _ in 0..tg.m {
            e_succ.push(TEdgeSucc {
                left: 0,
                right: 0,
            });
        }

        let e_betw = vec![N::r_zero(); tg.m as usize];

        let mut u_inf: Vec<NodeInfo<C, N>> = Vec::with_capacity(tg.n as usize + 1);

        // in-degrees:
        let mut in_deg: Vec<usize> = vec![0; (tg.n+1) as usize];
        for e in &tg.earr {
            in_deg[e.v as usize] += 1;
        }
        for u in 0..=tg.n as usize {
            u_inf.push(NodeInfo { 
                intervs: VecDeque::with_capacity(in_deg[u]/64), // max size is in_deg[v], but try to be smaller than tg
                preds: VecDeque::with_capacity(in_deg[u]/64),
                l: tg.u_fst[u] as Eind, 
                r: tg.u_fst[u] as Eind,
                refresh: tg.u_fst[u+1] as Eind,
                opt_cost: C::infinite_target_cost(),
                opt_walks: N::i_zero(), 
                betw: N::r_zero(),
            })
        }

        TGraphSweep::<'tg, C, N> {
            tg, e_inf, e_succ, e_betw, u_inf
        }
    }

    pub fn clear(&mut self) {
        for ei in self.e_inf.iter_mut() {
            ei.min_cost = Cost::infinite_cost();
            ei.walks = N::i_zero();
        }
        for es in self.e_succ.iter_mut() {
            es.left = 0;
            es.right = 0; // empty interval
        }
        for eb in self.e_betw.iter_mut() {
            *eb = N::r_zero();
        }
        for (u, ui) in self.u_inf.iter_mut().enumerate() {
            ui.intervs.clear();
            ui.preds.clear();
            ui.l = self.tg.u_fst[u] as Eind;
            ui.r = self.tg.u_fst[u] as Eind;
            ui.refresh = self.tg.u_fst[u+1] as Eind;
            ui.opt_cost = C::infinite_target_cost();
            ui.opt_walks = N::i_zero();
            ui.betw = N::r_zero();
        }
    }

   pub fn forward(&mut self, s: Node, beta: Time) {
        let infty = Cost::infinite_cost();
        for (i, e) in self.tg.earr.iter().enumerate() {
            self.finalize(e.u, self.e_inf[i].idep+1);
            let ei = &self.e_inf[i];
            if e.u == s || ei.min_cost < infty {

                // source:
                if e.u == s {
                    let ei = &mut self.e_inf[i];
                    let ce = Cost::edge_cost(e);
                    if ei.min_cost == infty || ce < ei.min_cost {
                        ei.min_cost = ce;
                        ei.walks = N::i_one();
                    } else if ce == ei.min_cost {
                        ei.walks += N::i_one();
                    }
                }

                let mut l_ei = self.u_inf[e.v as usize].l;
                let arr = e.arr();
                while l_ei < self.tg.u_fst[e.v as usize + 1] 
                    && self.tg.earr[self.tg.edep[l_ei as usize] as usize].t < arr { l_ei += 1; }

                self.finalize(e.v, l_ei);
                let ei =  &self.e_inf[i];

                //println!("{:?} i={i}  c={:?} w={} l_u={}", e, ei.min_cost, ei.walks, self.u_inf[e.u as usize].l);

                let vi = &mut self.u_inf[e.v as usize];
                let mut l = vi.r;
                let mut r_ei = vi.r;
                while r_ei < self.tg.u_fst[e.v as usize + 1] 
                    && self.tg.earr[self.tg.edep[r_ei as usize] as usize].t - arr <= beta { r_ei += 1; }

                while let Some(itv) = vi.intervs.back_mut() {
                    if itv.cost > ei.min_cost {
                        assert!(vi.preds.len() >= itv.preds as usize);
                        for _ in 0..itv.preds { 
                            let j = vi.preds.pop_back().unwrap();
                            self.e_succ[j as usize].right = l_ei;
                        }
                        l = itv.l;
                        vi.intervs.pop_back();
                    } else {
                        break
                    }
                }
                match vi.intervs.back_mut() {
                    Some(itv) if itv.cost == ei.min_cost => { // extend interval to ei.right
                        self.e_succ[i].left = itv.l;
                        self.e_succ[i].right = r_ei;
                        vi.preds.push_back(i as Eind);
                        itv.r = r_ei;
                        itv.walks += ei.walks.clone();
                        itv.preds += 1;
                    }
                    _  if l < r_ei => { // add l..r_ei interval
                            self.e_succ[i].left = l;
                            self.e_succ[i].right = r_ei;
                            vi.preds.push_back(i as Eind);
                            vi.intervs.push_back(Interv { 
                                l, r: r_ei, cost: ei.min_cost.clone(), walks: ei.walks.clone(), preds: 1 
                            });
                        }
                    _ => {}
                }

                vi.r = r_ei; 
            }
        }

        self.count_opt_cost_walks(s);

    }

    fn finalize(&mut self, v: Node, j_end: Eind) {
        let vi = &mut self.u_inf[v as usize];
        if j_end <= vi.l { return }
        while let Some(itv) = vi.intervs.front_mut() {
            if itv.l < j_end {
                let c = itv.cost.clone();
                /*
                let mut fin_edges = 
                |walks: N::Int, l: usize, r: usize, i: usize| -> N::Int {
                    for &j in &self.tg.edep[l..r] {
                        self.e_inf[j as usize].min_cost = c.clone() + C::edge_cost(&self.tg.earr[j as usize]);
                        self.e_inf[j as usize].walks = walks.clone();
                    }
                    self.e_inf[i].walks.clone() // to avoid double borrow
                };
                */
                let r_end = min(itv.r, j_end);
                while let Some(&i) = vi.preds.front() {
                    let ri = self.e_succ[i as usize].right;
                    if ri <= r_end {
                        //let iwalks = fin_edges(itv.walks.clone(), 
                        //                                itv.l as usize, ri as usize, i as usize);
                        for &j in &self.tg.edep[itv.l as usize .. ri as usize] {
                            self.e_inf[j as usize].min_cost = c.clone() + C::edge_cost(&self.tg.earr[j as usize]);
                            self.e_inf[j as usize].walks = itv.walks.clone();
                        }
                        //itv.walks -= iwalks;
                        vi.preds.pop_front();
                        let pred_walks = itv.walks.clone();
                        itv.walks -= self.e_inf[i as usize].walks.clone();
                        let precision: N::Int = 1000u32.into();
                        if ( ! N::i_is_exact()) 
                                && itv.walks < pred_walks / precision {
                            //log::info!("recompute nb walks for {v}");
                            itv.walks = N::i_zero();
                            for &i in &vi.preds {
                                itv.walks += self.e_inf[i as usize].walks.clone();
                            }
                        }
                        assert!(itv.preds > 0);
                        itv.preds -= 1;
                        itv.l = self.e_succ[i as usize].right;
                    } else {
                        break;
                    }
                }
                if j_end >= itv.r {
                    vi.intervs.pop_front();
                } else {
                    let mut fin_edges = 
                |walks: N::Int, l: usize, r: usize, i: usize| -> N::Int {
                    for &j in &self.tg.edep[l..r] {
                        self.e_inf[j as usize].min_cost = c.clone() + C::edge_cost(&self.tg.earr[j as usize]);
                        self.e_inf[j as usize].walks = walks.clone();
                    }
                    self.e_inf[i].walks.clone() // to avoid double borrow
                };
                    fin_edges(itv.walks.clone(), 
                          itv.l as usize, r_end as usize, 0);
                    itv.l = j_end;
                }
            } else {
                break
            }
        }
        vi.l = j_end;
        vi.r = max(vi.r, vi.l);
        //self.check_intervs(v);
    }

    /* fn check_intervs(&self, v: Node) {
        assert_eq!(self.u_inf[v as usize].preds.len() as Eind,
                   self.u_inf[v as usize].intervs.iter().map(|itv| itv.preds).sum())
    } */

    fn count_opt_cost_walks(&mut self, s: Node) -> () {
        self.u_inf[s as usize].opt_cost = C::empty_target_cost();
        for (i, e) in self.tg.earr.iter().enumerate() {
            if self.e_inf[i].min_cost < Cost::infinite_cost() {
                let c = C::target_cost(&self.e_inf[i].min_cost, e);
                if c < self.u_inf[e.v as usize].opt_cost {
                    self.u_inf[e.v as usize].opt_cost = c;
                }
            }
        }
        self.u_inf[s as usize].opt_walks = N::i_one(); // empty walk
        for (i, e) in self.tg.earr.iter().enumerate() {
            if self.e_inf[i].min_cost < Cost::infinite_cost() {
                let c = C::target_cost(&self.e_inf[i].min_cost, e);
                if c == self.u_inf[e.v as usize].opt_cost {
                    self.u_inf[e.v as usize].opt_walks += self.e_inf[i].walks.clone();
                }
            }
        }
    }

    pub fn backward(&mut self, s: Node) -> u64 {
        for v in 0..=self.tg.n as usize {
            assert_eq!(self.u_inf[v].l, self.tg.u_fst[v+1]);
            assert_eq!(self.u_inf[v].r, self.tg.u_fst[v+1]);
            assert_eq!(self.u_inf[v].refresh, self.tg.u_fst[v+1]);
        }

        let mut r_walks: Vec<N::Rat> = Vec::with_capacity(self.tg.m as usize);
        for i in 0..self.tg.m as usize {
            r_walks.push(N::real_of_int(self.e_inf[i].walks.clone())); // costly
        }
        let r_walks = r_walks;

        // compute edges_betweenness:
        let mut _outdeg_sum: u64 = 0;
        for (i, e) in self.tg.earr.iter().rev().enumerate() {
            let i = self.tg.m as usize - 1 - i; // rev
            let e_l = self.e_succ[i].left;
            let e_r = self.e_succ[i].right;
            if e_l < e_r {
                _outdeg_sum += (e_r - e_l) as u64;
                // update vi.betw so that it stores sum_{j in edep[e_l..e_r]} ej.betw / ej.walks
                { // borrow mutability through vi
                    let ei = &self.e_inf[i];
                    let vi = &mut self.u_inf[e.v as usize];
                    if e_r <= vi.l { // e_l..e_r and vi.l..vi.r do not overlap
                        vi.l = e_r; 
                        vi.r = e_r; // empty interval
                        vi.betw = N::r_zero();
                    }
                    for i_j in (max(e_r,vi.l) as usize .. vi.r as usize).rev() {
                        let j = self.tg.edep[i_j];
                        let f = &self.tg.earr[j as usize];
                        let ej = &self.e_inf[j as usize];
                        if e.v != s || ei.min_cost.clone() + C::edge_cost(f) == ej.min_cost {
                            vi.betw -= self.e_betw[j as usize].clone() / r_walks[j as usize].clone(); //N::real_of_int(ej.walks); 
                        }
                    }
                    vi.r = e_r;
                    if e_r <= vi.refresh && ! N::r_is_exact() { // reset to zero and recompute sum for better numerical stability
                        vi.l = e_r; // empty interval
                        vi.betw = N::r_zero(); 
                    }
                    for i_j in (e_l as usize .. min(e_r, vi.l) as usize).rev() {
                        let j = self.tg.edep[i_j];
                        let f = &self.tg.earr[j as usize];
                        let ej = &self.e_inf[j as usize];
                        if e.v != s || ei.min_cost.clone() + C::edge_cost(f) == ej.min_cost {
                            vi.betw += self.e_betw[j as usize].clone() / r_walks[j as usize].clone(); //N::real_of_int(ej.walks); 
                        }
                    }
                    vi.l = e_l;
                }
                /*
                let ei = &self.e_inf[i]; // now borrow mutability through ei
                let mut bsum = N::r_zero();
                for j in (e_l as usize .. e_r as usize).rev().map(|i_j| self.tg.edep[i_j]) {
                    let f = &self.tg.earr[j];
                    let ej = &self.e_inf[j];
                    if e.v != s || ei.min_cost.clone() + C::edge_cost(f) == ej.min_cost {
            //eprintln!("{i} {j} s={} bj={} wj={}", bsum, ej.betw, N::real_of_int(ej.walks));
                        bsum += ej.betw / N::real_of_int(ej.walks); 
                    }
                }
                */
                // set the betweenness of e to sum_{j in edep[e_l..e_r]} ej.betw * ei.walks / ej.walks
                let vi = &self.u_inf[e.v as usize];
                self.e_betw[i as usize] = /* bsum */ vi.betw.clone() * r_walks[i].clone(); //N::real_of_int(ei.walks);
                
            }
            let ei = &self.e_inf[i]; // now borrow mutability through ei
            let vi = &self.u_inf[e.v as usize];
            if ei.walks > N::i_zero() && C::target_cost(&ei.min_cost, e) == vi.opt_cost { // count the optimal walk through e
                self.e_betw[i as usize] += N::real_of_int(ei.walks.clone()) / N::real_of_int(vi.opt_walks.clone());
            }
        }

        /*
        for i in 0..self.tg.m {
            let e = &self.tg.earr[i];
            let ei = &self.e_inf[i];
            let si = &self.e_succ[i];
            eprintln!("{}: [{}] {} {} {:?} {} {:?}", i, self.tg.earr[i], si.left, si.right, ei.min_cost, ei.betw,
            C::target_cost(&ei.min_cost, e) == self.u_inf[e.v as usize].opt_cost);
        }
        */

        // compute nodes betweenness:
        for vi in &mut self.u_inf { 
            vi.betw = N::r_zero();
        }
        for (i, e) in self.tg.earr.iter().enumerate() {
            if e.v != s {  // do not count optimal walks to the source
                self.u_inf[e.v as usize].betw += self.e_betw[i].clone()
            }
        }
        for (v, vi) in self.u_inf.iter_mut().enumerate() {
            if v as Node != s && vi.opt_walks > N::i_zero() { 
                assert!(( ! N::r_is_exact() ) || vi.betw >= N::r_one());
                vi.betw -= N::r_one(); // do not count optimal walks to the node itself
            }
        }

        _outdeg_sum
    }

    pub fn betweenness(&mut self, beta: Time) -> Vec<N::Rat> {
        let mut pl = ProgressLogger::default();
        pl.expected_updates(Some(self.tg.n as usize));
        pl.display_memory(true);
        pl.item_name("node");
        pl.start("Computing betweenness...");
        let mut betw: Vec<N::Rat> = vec![N::r_zero(); self.tg.n as usize + 1];
        let mut _outdeg_sum: u128 = 0;
        for s in 0..=self.tg.n as Node {
            self.clear();
            self.forward(s, beta);
            _outdeg_sum += self.backward(s) as u128;
            for v in 0..=self.tg.n as usize {
                betw[v] += self.u_inf[v].betw.clone();
            }
            pl.update();
        }
        pl.stop();
        log::info!("successors: outdeg_sum_avg={:.4}  outdeg_avg={:.4}", 
                   _outdeg_sum as f64 / self.tg.m as f64, 
                   _outdeg_sum as f64 / (self.tg.n as u128 * self.tg.m as u128) as f64);
        betw
    }

    pub fn min_costs(&self) -> Vec<C> {
        (0..self.tg.m as usize).map(|i| self.e_inf[i].min_cost.clone()).collect()
    }
    
    pub fn walk_counts(&self) -> Vec<N::Int> {
        (0..self.tg.m as usize).map(|i| self.e_inf[i].walks.clone()).collect()
    }
    
    pub fn target_costs(&self) -> Vec<C::TargetCost> {
        (0..self.tg.m as usize).map(|i| 
            Cost::target_cost(&self.e_inf[i].min_cost, &self.tg.earr[i])
        ).collect()
    }
    
    pub fn opt_costs(&mut self) -> Vec<C::TargetCost> {
        (0..=self.tg.n as usize).map(|v| self.u_inf[v].opt_cost.clone()).collect()
    }

    pub fn opt_walk_counts(&mut self) -> Vec<N::Int> {
        (0..=self.tg.n as usize).map(|v| self.u_inf[v].opt_walks.clone()).collect()
    }

    pub fn nodes_betweenness(&mut self) -> Vec<N::Rat> {
        (0..=self.tg.n as usize).map(|v| self.u_inf[v].betw.clone()).collect()
    }

    pub fn edges_betweenness(&mut self) -> Vec<N::Rat> {
        (0..self.tg.m as usize).map(|i| self.e_betw[i].clone()).collect()
    }

}

pub fn betweenness_seq<C : Cost, N : Num>(tg: &TGraph, beta: Time) -> Vec<N::Rat> {
    let mut sweep: TGraphSweep<'_, C, N> = TGraphSweep::new(tg);
    sweep.betweenness(beta)
}

pub fn betweenness_par<C : Cost, N : Num>(tg: &TGraph, beta: Time, nthread: u32) -> Vec<N::Rat> {
    let n_th = std::thread::available_parallelism()
        .unwrap_or(std::num::NonZeroUsize::new(2).unwrap()).get() as Node;
    let n_th = if n_th > 2 { n_th - 1 } else { n_th };
    let n_th = if nthread > 0 { nthread } else { n_th };
    info!("Using {n_th} threads.");
    let (tx, rx) = mpsc::channel();
    thread::scope(|s| {
        for i_th in 0..n_th as Node {
            let tx = tx.clone();
            //let tg = tg.clone();
            s.spawn(move || {
                let mut plopt = None;
                if i_th == 0 {
                    let mut pl = ProgressLogger::default();
                    pl.expected_updates(Some(tg.n as usize));
                    pl.display_memory(true);
                    pl.item_name("node");
                    pl.start("Computing betweenness in parallel...");
                    plopt = Some(pl);
                }
                let mut sweep: TGraphSweep<'_, C, N> = TGraphSweep::new(tg);
                let mut betw: Vec<N::Rat> = vec![N::r_zero(); tg.n as usize + 1];
                for s in 0..=tg.n  {
                    if s % n_th as Node == i_th {
                        sweep.clear();
                        sweep.forward(s, beta);
                        sweep.backward(s);
                        for v in 0..=tg.n as usize {
                            betw[v] += sweep.u_inf[v].betw.clone();
                        }
                    }
                    if i_th == 0 { plopt.as_mut().expect("no logger").update() }
                }
                if i_th == 0 { plopt.unwrap().stop() }
                tx.send(betw).unwrap();
            });
        }
    });
    let mut betw: Vec<N::Rat> = vec![N::r_zero(); tg.n as usize + 1];
    for _ in 0..n_th {
        let b = rx.recv().unwrap();
        for v in 0..=tg.n as usize {
            betw[v] += b[v].clone();
        }
    }
    betw
}



#[cfg(test)]
pub mod tests_tbet {
    use num_bigint::BigUint;
    use num_rational::Ratio;

    type NumExact = NumT<BigUint, Ratio<BigUint>, True, True>;
    type NumU128F64 = NumT<u128, f64, True, False>;

    use super::*;
    use crate::cost;

    #[test]
    fn test_sweep() {
        for tg_str in tests::TEST_TGRAPHS_STR {
            let tg: TGraph = tg_str.parse().unwrap();
            for e in tg.earr.iter() { println!("{:?}", e) }
            print!("{}", tg_str);
            for e in tg.edep.iter().map(|&i| &tg.earr[i as usize]) { println!("{:?}", e) }
            for beta in [0, 1, 2, 100] {
                println!("--- beta = {beta} ---");
                let mut sweep: TGraphSweep<cost::Shortest, NumU128F64> = TGraphSweep::new(&tg);
                sweep.forward(1, beta);
                println!("costs from 1: {:?}", sweep.min_costs());
                println!("nwalks from 1: {:?}", sweep.walk_counts());
                println!("targets from 1: {:?}", sweep.target_costs());
                sweep.clear();
                sweep.forward(2, beta);
                println!("costs from 2: {:?}", sweep.min_costs());
                println!("targets from 2: {:?}", sweep.target_costs());
            }
        }
    }

    #[test]
    fn examples() {
        let tg = tests::LZO.parse().unwrap();
        let big= |i: &u32| -> BigUint { (*i).into() };
        let br = |v : Vec<u32>| -> Vec<Ratio<BigUint>> {
            v.iter().map(|i| big(i).into()).collect()
        };

        assert_eq!(betweenness_seq::<cost::Shortest, NumExact>(&tg, 0), br(Vec::from([0, 0, 1, 4, 0, 1, 0])));
        assert_eq!(betweenness_seq::<cost::Shortest, NumExact>(&tg, 1), br(Vec::from([0, 0, 1, 6, 0, 2, 0])));
        assert_eq!(betweenness_seq::<cost::Shortest, NumExact>(&tg, 2), br(Vec::from([0, 0, 1, 6, 0, 1, 0])));
        assert_eq!(betweenness_seq::<cost::Shortest, NumExact>(&tg, 100), br(Vec::from([0, 0, 1, 6, 0, 0, 0])));
        
        assert_eq!(betweenness_seq::<cost::Foremost, NumExact>(&tg, 0), br(Vec::from([0, 0, 2, 4, 0, 2, 0])));
        assert_eq!(betweenness_seq::<cost::Foremost, NumExact>(&tg, 1), br(Vec::from([0, 0, 3, 6, 0, 3, 0])));
        assert_eq!(betweenness_seq::<cost::Foremost, NumExact>(&tg, 2), br(Vec::from([0, 0, 3, 6, 0, 3, 0])));
        assert_eq!(betweenness_seq::<cost::Foremost, NumExact>(&tg, 100), br(Vec::from([0, 0, 3, 6, 0, 3, 0])));

        let tg = tests::ZERO.parse().unwrap();
        let br = |v : Vec<(u32, u32)>| -> Vec<Ratio<BigUint>> {
            v.iter().map(|(p,q)| Into::<Ratio<BigUint>>::into(big(p)) / big(q)).collect()
        };

        assert_eq!(betweenness_seq::<cost::Shortest, NumExact>(&tg, 0), br(Vec::from([(0,1), (0,1), (1,1), (1,1), (5,1), (3,1), (2,1)])));
        assert_eq!(betweenness_seq::<cost::Latest, NumExact>(&tg, 1), br(Vec::from([(0,1), (0,1), (0,1), (1,1), (20,3), (13,3), (11,3)])));
        assert_eq!(betweenness_seq::<cost::Foremost, NumExact>(&tg, 2), br(Vec::from([(0,1), (0,1), (4,1), (2,1), (19,3), (8,3), (10,3)])));

        for tg_str in tests::TEST_TGRAPHS_STR {
            let tg: TGraph = tg_str.parse().unwrap();
            for e in tg.earr.iter() { println!("{:?}", e) }
            print!("{}", tg_str);
            for e in tg.edep.iter().map(|&i| &tg.earr[i as usize]) { println!("{:?}", e) }
            for beta in [0, 1, 2, 100] {
                assert_eq!(betweenness_seq::<cost::Shortest, NumExact>(&tg, beta), 
                           betweenness_par::<cost::Shortest, NumExact>(&tg, beta, 3));
                assert_eq!(betweenness_seq::<cost::Fastest, NumExact>(&tg, beta), 
                           betweenness_par::<cost::Fastest, NumExact>(&tg, beta, 3));
                assert_eq!(betweenness_seq::<cost::Foremost, NumExact>(&tg, beta), 
                           betweenness_par::<cost::Foremost, NumExact>(&tg, beta, 3));
            }
        }
    }
}