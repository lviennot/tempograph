use crate::tgraph::*;
use crate::cost::Cost;

pub struct TSweep<'tg, C : Cost> {
    tg: &'tg TGraph,
    min_cost: Vec<C>, // min-cost of a walk ending with that tedge (edep order)
    pred: Vec<Eind>, // previous tedge of such a walk
    u_inf: Vec<NodeInfo<C>>,
}

use std::collections::VecDeque;

struct NodeInfo<C : Cost> { // info about tedges from a node u (forward phase)
    intervs: VecDeque<Interv<C>>, // information about consecutive intervals of edges from u (in edgep)
    l: Eind, // left bound in edep, intervals in intervs span l..r
    r: Eind, // right bound in edep (excluded)
}

#[derive(Debug)]
struct Interv<C : Cost> { // an interval of edges from a given node v in edep
    l: Eind, // left bound in edep
    r: Eind, // right bound in edep (excluded)
    cost: C, // minimum cost of walks reaching this interval (and that edges in the interval extend)
    pred: Eind, // index of an edge ending such a walk
}

use std::cmp::{min, max};

impl<'tg, C : Cost> TSweep<'tg, C> {

    pub fn new(tg: &'tg TGraph) -> TSweep<'tg, C> {

        let min_cost: Vec<C> = vec![C::infinite_cost(); tg.m];

        let pred = (0..tg.m).into_iter().collect();

        let mut u_inf: Vec<NodeInfo<C>> = Vec::with_capacity(tg.n + 1);

        // in-degrees:
        let mut in_deg: Vec<usize> = vec![0; tg.n+1];
        for e in &tg.edep {
            in_deg[e.v] += 1;
        }
        for u in 0..tg.n {
            u_inf.push(NodeInfo { 
                intervs: VecDeque::with_capacity(in_deg[u]/64), // max size is in_deg[v], but try to be smaller than tg
                l: tg.u_fst[u] as Eind, 
                r: tg.u_fst[u] as Eind,
            })
        }

        TSweep::<'tg, C> {
            tg, min_cost, pred, u_inf
        }
    }

    pub fn clear(&mut self) {
        for c in self.min_cost.iter_mut() {
            *c = C::infinite_cost();
        }
        for (u, ui) in self.u_inf.iter_mut().enumerate() {
            ui.intervs.clear();
            ui.l = self.tg.u_fst[u] as Eind;
            ui.r = self.tg.u_fst[u] as Eind;
        }
    }

   pub fn scan(&mut self, s: Node, beta: Time) {
        //let succ = self.tg.extend_indexes(beta);
        let infty = Cost::infinite_cost();
        for &i in self.tg.earr.iter() {
            let e = & self.tg.edep[i];
            let e_arr = e.arr();

            // ------------ finalize_tail e

            // finalize edges with tail e.u from self.u_inf[e.u].l up to e
            let ui = &mut self.u_inf[e.u];
            while let Some(itv) = ui.intervs.front_mut() {
                if itv.l <= i {
                    let r = min(itv.r, i+1);
                    for j in itv.l..r {
                        self.min_cost[j] = itv.cost.clone() + C::edge_cost(e);
                        self.pred[j] = itv.pred;
                    }
                    if r >= itv.r { ui.intervs.pop_front(); } // remove interval
                    else { itv.l = r; break; } // truncate interval
                } else {
                    break;
                }
            }

            // edges from the source
            if e.u == s {
                for j in ui.l..=i {
                    let c = C::edge_cost(&self.tg.edep[j]);
                    if c <= self.min_cost[j] { 
                        self.min_cost[j] = c;
                        self.pred[j] = j;
                    }
                }
            }

            // update left of union of intervals
            ui.l = i+1;

            //  ----------- relax_head e

            let e_min_cost = self.min_cost[i].clone();
            //eprintln!("[{}] {:?}", e, e_min_cost);

            if e_min_cost < infty {

                // Compute the interval edep[l_e..r_e] of edges extending e:
                let mut l_e = self.u_inf[e.v].l;
                while l_e < self.tg.u_fst[e.v + 1] 
                    && self.tg.edep[l_e].t < e_arr { l_e += 1; }

                let mut r_e = max(l_e, self.u_inf[e.v].r);
                while r_e < self.tg.u_fst[e.v + 1] 
                    && self.tg.edep[r_e].t - e_arr <= beta { r_e += 1; }
                //assert_eq!((l_e,r_e), succ[i]);

                // Compute the interval edep[l_c..r_e] where e provides min-cost:
                let vi = &mut self.u_inf[e.v];
                let mut l_c = max(l_e, vi.r);

                // Remove intervals with larger cost
                while let Some(itv) = vi.intervs.back_mut() {
                    if itv.r > l_e && itv.cost > e_min_cost { // larger cost
                        l_c = max(l_e, itv.l);
                        if itv.l >= l_e { vi.intervs.pop_back(); } // remove interval
                        else { itv.r = l_e; break; } // truncate interval
                    } else {
                        break
                    }
                }

                // add interval:
                if l_c < r_e { // if l_c..r_e is not empty
                    vi.intervs.push_back(Interv { 
                        l: l_c, r: r_e, cost: e_min_cost, pred: i 
                    });
                }

                // update right of union of intervals
                vi.r = r_e;

            }

        }
    }

    pub fn opt_costs(&self, src: Node) -> Vec<C::TargetCost> {
        let mut opt = vec![C::infinite_target_cost(); self.tg.n];
        opt[src] = C::empty_target_cost();
        let infty = C::infinite_cost();
        for (i, e) in self.tg.edep.iter().enumerate() {
            if self.min_cost[i] < infty {
                let c = C::target_cost(&self.min_cost[i], e);
                //eprintln!("[{}] {:?}", e, c);
                if c < opt[e.v] { opt[e.v] = c; }
            }
        }
        opt
    }

}


#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::cost::*;
    use crate::tgraph::tests as tgt;

    pub const TESTS_SHORTEST: [(&str, Time, [Hop; 6]); 8] = 
        [
            (tgt::BLACKBOARD, 1, [0, 1, 5, 2, 3, 777]),
            (tgt::BLACKBOARD, 2, [0, 1, 2, 2, 3, 777]),
            (tgt::LZO, 2, [0, 1, 1, 3, 2, 2]),
            (tgt::LZO, 3, [0, 1, 1, 3, 2, 2]),
            (tgt::ROMBUS_4, 0, [0, 1, 2, 2, 5, 777]),
            (tgt::ROMBUS_4, 2, [0, 1, 2, 2, 3, 777]),
            (tgt::ROMBUS_4, 10, [0, 1, 2, 2, 3, 777]),
            (tgt::WHITEBOARD, 4, [0, 1, 1, 3, 777, 777]),
        ];
        //[BLACKBOARD, LZO, ROMBUS_4, WHITEBOARD];

    #[test]
    fn test_shortest() {
        let src = 1;
        for (tg_str, beta, sol) in TESTS_SHORTEST {
            let tg: TGraph = tg_str.parse().unwrap();
            let mut tsweep: TSweep<Shortest> = TSweep::new(&tg);
            tsweep.scan(src, beta);
            let opt_costs = tsweep.opt_costs(src);
            for v in 1..tg.n { assert_eq!(opt_costs[v], sol[v-1]) }
        }
    }

}