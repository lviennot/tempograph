use crate::tgraph::*;
use crate::cost::*;

use std::collections::VecDeque;

pub fn tbfs(tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node) -> Vec<Hop> {
    let mut hop: Vec<Hop> = vec![Hop::MAX; tg.m];
    let mut u_eat: Vec<Time> = vec![Time::MAX; tg.n]; // departure of edep[u_left[v]]
    let mut u_left = vec![Eind::MAX; tg.n]; // all tedges from v in edep from u_left[v] have been already scanned
    let mut u_hop: Vec<Hop> = vec![Hop::MAX; tg.n];
    let mut right_index: Vec<Eind> = Vec::with_capacity(tg.n); // right index of a segment of already scanned tedges
    let mut right: Vec<Eind> = vec![Eind::MAX; tg.m]; // for each a pointer to a right_index
    let mut queue: VecDeque<Eind> = VecDeque::new();
    u_hop[s] = 0; 
    u_eat[s] = Time::MIN;
    let rgt = right_index.len(); right_index.push(tg.u_fst[s+1]);
    for i in tg.u_fst[s] .. tg.u_fst[s+1] {
        right[i] = rgt;
        hop[i] = 1;
        let e = &tg.edep[i];
        if u_hop[e.v] == Hop::MAX {
            u_hop[e.v] = 1;
        }
        queue.push_back(i);
    }
    u_left[s] = tg.u_fst[s];
    while let Some(i) = queue.pop_front() {
        let (l, r) = succ[i];
        let vi = tg.edep[i].v;
        if l < r {
            let h = hop[i] + 1;
            let mut rgt = right[l];
            if rgt == Eind::MAX {
                rgt = right_index.len(); right_index.push(l);
            }
            if right_index[rgt] >= r { continue } // l..r window is included in the segment
            let mut r_of_segment = r;
            assert!(r_of_segment >= right_index[rgt]);
            for j in right_index[rgt]..std::cmp::min(r, u_left[vi]) {
                if right[j] != Eind::MAX { 
                    r_of_segment = right_index[right[j]];
                    assert!(r_of_segment >= r);
                    break 
                }
                right[j] = rgt;
                if hop[j] > h {
                    hop[j] = h;
                    let e = &tg.edep[j];
                    if true || e.arr() < u_eat[e.v] { // u_hop[e.v] == Hop::MAX { //
                        if u_hop[e.v] > h { u_hop[e.v] = h }
                        queue.push_back(j);
                    }
                }
            }
            right_index[rgt] = r_of_segment;
            if r >= u_left[vi] { 
                u_left[vi] = l; 
                if l < r { u_eat[vi] = tg.edep[l].t }
            }
        }
    }
    u_hop
}

pub fn tbfs_non_linear(tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node) -> Vec<Hop> {
    let mut hop: Vec<Hop> = vec![Hop::MAX; tg.m];
    let mut u_eat: Vec<Time> = vec![Time::MAX; tg.n]; // departure if edep[u_left[v]]
    let mut u_hop: Vec<Hop> = vec![Hop::MAX; tg.n];
    let mut u_left = vec![Eind::MAX; tg.n]; // all edges from v in edep from u_left[v] have been already scanned
    for v in 0..tg.n { u_left[v] = tg.u_fst[v+1] }
    let mut queue: VecDeque<Eind> = VecDeque::new();
    u_hop[s] = 0; 
    u_eat[s] = Time::MIN;
    for i in tg.u_fst[s] .. tg.u_fst[s+1] {
        hop[i] = 1;
        let e = &tg.edep[i];
        if u_hop[e.v] == Hop::MAX {
            u_hop[e.v] = 1;
        }
        queue.push_back(i);
    }
    u_left[s] = tg.u_fst[s];
    while let Some(i) = queue.pop_front() {
        let (l, r) = succ[i];
        let vi = tg.edep[i].v;
        let h = hop[i] + 1;
        for j in l..std::cmp::min(r, u_left[vi]) {
            if hop[j] > h {
                hop[j] = h;
                let e = &tg.edep[j];
                if e.arr() < u_eat[e.v] { // u_hop[e.v] == Hop::MAX { //
                    if u_hop[e.v] > h { u_hop[e.v] = h }
                    queue.push_back(j);
                }
            }
        }
        if r >= u_left[vi] { 
            u_left[vi] = l; 
            if l < r { u_eat[vi] = tg.edep[l].t }
        }
    }
    u_hop
}

pub fn tbfs_inf(tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node) -> Vec<Hop> {
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
    u_left[s] = tg.u_fst[s];
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

