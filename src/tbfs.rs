use crate::tgraph::*;
use crate::cost::*;

use std::collections::VecDeque;

pub fn tbfs(tg: &TGraph, succ: &Vec<(usize, usize)>, beta: Time, s: Node) -> Vec<Hop> {
    let mut hop: Vec<Hop> = vec![Hop::MAX; tg.m];
    let mut u_eat: Vec<Time> = vec![Time::MAX; tg.n];
    let mut u_hop: Vec<Hop> = vec![Hop::MAX; tg.n];
    let mut u_left = vec![Eind::MAX; tg.n];
    let mut queue: VecDeque<Eind> = VecDeque::new();
    u_hop[s] = 0; 
    u_eat[s] = Time::MIN;
    for i in tg.u_fst[s] .. tg.u_fst[s+1] {
        hop[i] = 1;
        let e = &tg.edep[i];
        if u_hop[e.v] == Hop::MAX || u_eat[e.v] > e.arr() {
            u_hop[e.v] = 1;
            u_eat[e.v] = e.arr();
            queue.push_back(i);
        }
    }
    while let Some(i) = queue.pop_front() {
        let (l, r) = succ[i];
        let vi = tg.edep[i].v;
        let r = std::cmp::min(r, u_left[vi]);
        let h = hop[i] + 1;
        for j in l..r {
            if hop[j] > h {
                hop[j] = h;
                let e = &tg.edep[j];
                if u_hop[e.v] == Hop::MAX || u_eat[e.v] > e.arr() {
                    if u_hop[e.v] == Hop::MAX { u_hop[e.v] = h }
                    u_eat[e.v] = e.arr();
                    queue.push_back(j);
                }
            }
        }
        u_left[vi] = l;
    }
    u_hop
}

