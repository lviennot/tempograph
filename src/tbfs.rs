use crate::tgraph::*;
use crate::cost::*;

use std::collections::VecDeque;

pub struct TBFS {
    // for tbfs_inf_prune:
    hop: Vec<Hop>,
    u_eat: Vec<Time>, // departure if edep[u_left[v]]
    pub u_hop: Vec<Hop>,
    u_left: Vec<Eind>,  // all edges from v in edep from u_left[v] have been already scanned
    queue: VecDeque<Eind>,
    // for tbfs and tbfs_prune:
    right_index: Vec<Eind>, // right index of a segment of already scanned tedges
    right: Vec<Eind>,// for each a pointer to a right_index
}


impl TBFS {

    pub fn new(tg: &TGraph) -> Self {
        let hop = vec![Hop::MAX; tg.m];
        let u_eat = vec![Time::MAX; tg.n];
        let u_hop = vec![Hop::MAX; tg.n];
        let u_left = vec![Eind::MAX; tg.n];
        let queue = VecDeque::new();
        let right_index = Vec::with_capacity(tg.n);
        let right = vec![Eind::MAX; tg.m];
        TBFS { hop, u_eat, u_hop, u_left, queue, right_index, right }
    }

    pub fn clear(&mut self) {
        self.hop.fill(Hop::MAX);
        self.u_eat.fill(Time::MAX);
        self.u_hop.fill(Hop::MAX);
        self.u_left.fill(Eind::MAX);
        self.queue.clear();
        self.right_index.clear();
        self.right.fill(Eind::MAX);
    }

    fn check_size(&self, tg: &TGraph) {
        if tg.n != self.u_eat.len() || tg.m != self.hop.len() { 
            panic!("TGraph of {} nodes and {} tedges is too large for this TBFS struct of size {},{}!", 
                   tg.n, tg.m, self.u_eat.len(), self.hop.len()) 
        }
    }

    pub fn tbfs(&mut self, tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node) {
        self.check_size(tg);
        self.clear();
        self.u_hop[s] = 0; 
        let rgt = self.right_index.len(); self.right_index.push(tg.u_fst[s+1]);
        for i in tg.u_fst[s] .. tg.u_fst[s+1] {
            self.right[i] = rgt;
            self.hop[i] = 1;
            let e = &tg.edep[i];
            if self.u_hop[e.v] == Hop::MAX {
                self.u_hop[e.v] = 1;
            }
            self.queue.push_back(i);
        }
        while let Some(i) = self.queue.pop_front() {
            let (l, r) = succ[i];
            if l < r {
                let h = self.hop[i] + 1;
                let mut rgt = self.right[l];
                if rgt == Eind::MAX {
                    rgt = self.right_index.len(); self.right_index.push(l);
                }
                if self.right_index[rgt] >= r { continue } // l..r window is included in the segment
                let mut r_of_segment = r;
                assert!(r_of_segment >= self.right_index[rgt]);
                for j in self.right_index[rgt]..r {
                    if self.right[j] != Eind::MAX { 
                        r_of_segment = self.right_index[self.right[j]];
                        assert!(r_of_segment >= r);
                        break 
                    }
                    self.right[j] = rgt;
                    if self.hop[j] > h {
                        self.hop[j] = h;
                        let e = &tg.edep[j];
                        if self.u_hop[e.v] > h { self.u_hop[e.v] = h }
                        self.queue.push_back(j);
                    }
                }
                self.right_index[rgt] = r_of_segment;
            }
        }
    }    

    pub fn tbfs_inf(&mut self, tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node) {
        self.check_size(tg);
        self.clear();
        self.u_hop[s] = 0; 
        self.u_eat[s] = Time::MIN;
        for i in tg.u_fst[s] .. tg.u_fst[s+1] {
            self.hop[i] = 1;
            let e = &tg.edep[i];
            if self.u_hop[e.v] == Hop::MAX || self.u_eat[e.v] > e.arr() {
                self.u_hop[e.v] = 1;
                self.u_eat[e.v] = e.arr();
                self.queue.push_back(i);
            }
        }
        self.u_left[s] = tg.u_fst[s];
        while let Some(i) = self.queue.pop_front() {
            let (l, r) = succ[i];
            let vi = tg.edep[i].v;
            let r = std::cmp::min(r, self.u_left[vi]);
            let h = self.hop[i] + 1;
            for j in l..r {
                if self.hop[j] > h {
                    self.hop[j] = h;
                    let e = &tg.edep[j];
                    if self.u_hop[e.v] == Hop::MAX || self.u_eat[e.v] > e.arr() {
                        if self.u_hop[e.v] == Hop::MAX { self.u_hop[e.v] = h }
                        self.u_eat[e.v] = e.arr();
                        self.queue.push_back(j);
                    }
                }
            }
            self.u_left[vi] = l;
        }
    }
        
    pub fn tbfs_non_linear(&mut self, tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node) {
        self.check_size(tg);
        self.clear();
        for v in 0..tg.n { self.u_left[v] = tg.u_fst[v+1] }
        self.u_hop[s] = 0; 
        self.u_eat[s] = Time::MIN;
        for i in tg.u_fst[s] .. tg.u_fst[s+1] {
            self.hop[i] = 1;
            let e = &tg.edep[i];
            if self.u_hop[e.v] == Hop::MAX {
                self.u_hop[e.v] = 1;
            }
            self.queue.push_back(i);
        }
        self.u_left[s] = tg.u_fst[s];
        while let Some(i) = self.queue.pop_front() {
            let (l, r) = succ[i];
            let vi = tg.edep[i].v;
            let h = self.hop[i] + 1;
            for j in l..std::cmp::min(r, self.u_left[vi]) {
                if self.hop[j] > h {
                    self.hop[j] = h;
                    let e = &tg.edep[j];
                    if e.arr() < self.u_eat[e.v] { // u_hop[e.v] == Hop::MAX { //
                        if self.u_hop[e.v] > h { self.u_hop[e.v] = h }
                        self.queue.push_back(j);
                    }
                }
            }
            if r >= self.u_left[vi] { 
                self.u_left[vi] = l; 
                if l < r { self.u_eat[vi] = tg.edep[l].t }
            }
        }
    }


    /// Returns exact harmonic closeness of s if > min_harm_clos, otherwise returns an upper bound of it.
    pub fn tbfs_inf_prune(&mut self, tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node, min_harm_clos: f64) -> f64 {
        self.check_size(tg);
        self.clear();
        let mut s_hc = 0.;
        let mut n_unreached: Node = tg.n - 1; // do not count `s`
        self.u_hop[s] = 0; 
        self.u_eat[s] = Time::MIN;
        for i in tg.u_fst[s] .. tg.u_fst[s+1] {
            self.hop[i] = 1;
            let e = &tg.edep[i];
            if self.u_hop[e.v] == Hop::MAX || self.u_eat[e.v] > e.arr() {
                if self.u_hop[e.v] == Hop::MAX {
                    self.u_hop[e.v] = 1;
                    s_hc += 1.0; // 1 / 1
                    n_unreached -= 1;
                }
                self.u_eat[e.v] = e.arr();
                self.queue.push_back(i);
            }
        }
        self.u_left[s] = tg.u_fst[s];
        while let Some(i) = self.queue.pop_front() {
            let (l, r) = succ[i];
            let vi = tg.edep[i].v;
            let r = std::cmp::min(r, self.u_left[vi]);
            let h = self.hop[i] + 1;
            let hc_ub = s_hc + (n_unreached as f64) / (h as f64);
            if hc_ub < min_harm_clos { return hc_ub }
            for j in l..r {
                if self.hop[j] > h {
                    self.hop[j] = h;
                    let e = &tg.edep[j];
                    if self.u_hop[e.v] == Hop::MAX || self.u_eat[e.v] > e.arr() {
                        if self.u_hop[e.v] == Hop::MAX { 
                            self.u_hop[e.v] = h;
                            s_hc += 1.0 / (h as f64);
                            n_unreached -= 1;
                        }
                        self.u_eat[e.v] = e.arr();
                        self.queue.push_back(j);
                    }
                }
            }
            self.u_left[vi] = l;
        }
        s_hc
    }

    /// Returns exact harmonic closeness of s if > min_harm_clos, otherwise returns an upper bound of it.
    pub fn tbfs_prune(&mut self, tg: &TGraph, succ: &Vec<(usize, usize)>, s: Node, min_harm_clos: f64) -> f64 {
        self.check_size(tg);
        self.clear();
        let mut s_hc = 0.;
        let mut n_unreached: Node = tg.n - 1; // do not count `s`
        self.u_hop[s] = 0; 
        let rgt = self.right_index.len(); self.right_index.push(tg.u_fst[s+1]);
        for i in tg.u_fst[s] .. tg.u_fst[s+1] {
            self.right[i] = rgt;
            self.hop[i] = 1;
            let e = &tg.edep[i];
            if self.u_hop[e.v] == Hop::MAX {
                self.u_hop[e.v] = 1;
                s_hc += 1.0; // 1 / 1
                n_unreached -= 1;
            }
            self.queue.push_back(i);
        }
        while let Some(i) = self.queue.pop_front() {
            let (l, r) = succ[i];
            if l < r {
                let h = self.hop[i] + 1;
                let hc_ub = s_hc + (n_unreached as f64) / (h as f64);
                if hc_ub < min_harm_clos { return hc_ub }
                let mut rgt = self.right[l];
                if rgt == Eind::MAX {
                    rgt = self.right_index.len(); self.right_index.push(l);
                }
                if self.right_index[rgt] >= r { continue } // l..r window is included in the segment
                let mut r_of_segment = r;
                assert!(r_of_segment >= self.right_index[rgt]);
                for j in self.right_index[rgt]..r {
                    if self.right[j] != Eind::MAX { 
                        r_of_segment = self.right_index[self.right[j]];
                        assert!(r_of_segment >= r);
                        break 
                    }
                    self.right[j] = rgt;
                    if self.hop[j] > h {
                        self.hop[j] = h;
                        let e = &tg.edep[j];
                        if self.u_hop[e.v] > h { 
                            self.u_hop[e.v] = h;
                            s_hc += 1.0 / (h as f64);
                            n_unreached -= 1;
                        }
                        self.queue.push_back(j);
                    }
                }
                self.right_index[rgt] = r_of_segment;
            }
        }
        s_hc
    }


}





