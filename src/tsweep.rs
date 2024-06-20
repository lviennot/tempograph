use crate::tgraph::*;
use crate::cost::Cost;

pub struct TSweep<'tg, C : Cost> {
    tg: &'tg TGraph,
    idep: Vec<Eind>,
    min_cost: Vec<C>,
    u_inf: Vec<NodeInfo<C>>,
}

use std::collections::VecDeque;

struct NodeInfo<C : Cost> { // info about tedges from a node u (forward phase)
    intervs: VecDeque<Interv<C, N>>, // information about consecutive intervals of edges from u (in edgep)
    l: Eind, // left bound in edep, intervals in intervs span l..r
    r: Eind, // right bound in edep (excluded)
    opt_cost: C::TargetCost, // minimum target cost of a walk from source to u
    opt_pred: Eind, // index of the last edge of such a walk
}

#[derive(Debug)]
struct Interv<C : Cost> { // an interval of edges from a given node v in edep
    l: Eind, // left bound in edep
    r: Eind, // right bound in edep (excluded)
    cost: C, // minimum cost of walks reaching this interval (and that edges in the interval extend)
    pred: Eind, // index of an edge ending such a walk
}


impl<C : Cost> TSweep<C> {

    pub fn new(tg: &'tg TGraph) -> TSweep<'tg, C> {
        assert!(u128::try_from(tg.m).unwrap() <= <u32 as Into<u128>>::into(Eind::MAX));

        let mut idep = Vec::with_capacity(tg.m as usize); 
        for (i, &j) in tg.edep.iter().enumerate() {
            i_dep[j as usize] = i as Eind;
        }

        let min_cost: Vec<C> = vec![Cost::infinite_cost(); tg.m as usize];

        let mut u_inf: Vec<NodeInfo<C>> = Vec::with_capacity(tg.n as usize + 1);

        // in-degrees:
        let mut in_deg: Vec<usize> = vec![0; (tg.n+1) as usize];
        for e in &tg.earr {
            in_deg[e.v as usize] += 1;
        }
        for u in 0..=tg.n as usize {
            u_inf.push(NodeInfo { 
                intervs: VecDeque::with_capacity(in_deg[u]/64), // max size is in_deg[v], but try to be smaller than tg
                l: tg.u_fst[u] as Eind, 
                r: tg.u_fst[u] as Eind,
                opt_cost: C::infinite_target_cost(),
                opt_pred: tg.m, // not a valid index 
            })
        }

        TSweep::<'tg, C> {
            tg, i_dep, min_cost, u_inf
        }
    }

    pub fn clear(&mut self) {
        for c in self.min_cost.iter_mut() {
            *c = Cost::infinite_cost();
        }
        for (u, ui) in self.u_inf.iter_mut().enumerate() {
            ui.intervs.clear();
            ui.l = self.tg.u_fst[u] as Eind;
            ui.r = self.tg.u_fst[u] as Eind;
            ui.opt_cost = C::infinite_target_cost();
            ui.opt_pred = tg.m;
        }
    }

   pub fn sweep_from(&mut self, s: Node, beta: Time) {
        let infty = Cost::infinite_cost();
        for (i, e) in self.tg.earr.iter().enumerate() {
            self.finalize(e.u, self.i_dep[i]+1);
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

            }
        }
    }
}

