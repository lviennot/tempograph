use crate::graph;
use std::collections::HashMap;

/// node index
pub type Node = usize; 

/// Time:MIN is reserved as infinite for fastest cost
pub type Time = i64; 

/// tedge index
pub type Eind = usize; 


/// A temporal graph as a doubly sorted list of temporal edges.
#[derive(Clone)]
pub struct TGraph {
    /// number of nodes (nodes are numbered from 0 to n-1)
    pub n: Node, 

    /// number of temporal edges
    pub m: Eind, 

    /// tedges sorted by tail and departure time
    pub edep: Vec<TEdge>, 

    /// indexes of tedges from u are in edep[u_fst[u]..u_fst[u+1]]
    pub u_fst: Vec<Eind>, 

    /// indexex of tedges (in edep) sorted by arrival time
    pub earr: Vec<Eind>, 
}


/// A temporal edge
#[derive(PartialEq, Eq, PartialOrd, Ord, Clone, Debug)]
pub struct TEdge {
    /// tail
    pub u: Node,
 
    /// head
    pub v: Node,
 
    /// departure time
    pub t: Time, 
 
    /// delay (or travel time)
    pub d: Time, 
}


use core::cmp::max;
use std::str::FromStr;
use std::cmp::Ordering;
use std::fmt;

impl fmt::Display for TEdge {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        //write!(f, "{} --{}+{}={}--> {}", self.u, self.t, self.d, self.arr(), self.v)
        write!(f, "{} {} {} {}", self.u, self.v, self.t, self.d)
    }
}
impl TEdge {
    /// Returns the arrival time of this [`TEdge`].
    pub fn arr(&self) -> Time {
        self.t + self.d
    }

    /// Returns the departure time of this [`TEdge`].
    pub fn dep(&self) -> Time {
        self.t
    }

    pub fn extends(&self, other: &TEdge, beta: Time) -> bool {
        self.u == other.v && self.t >= other.arr() && self.t - other.arr() <= beta 
    }

    fn parse<I : FromStr>(opt: Option<&str>) -> Result<I, String> {
        match opt {
            Some(s) => return s.parse::<I>().map_err(|_| format!("bad number format: `{}`", s)),
            None => return Err("number expected".to_string()),
        }
    }

    fn cmp_u_t(&self, other: &Self) -> Ordering {
        match self.u.cmp(&other.u) {
            Ordering::Equal => return self.t.cmp(&other.t),
            lg => return lg,
        } 
    }

    fn cmp_arr(&self, other: &Self) -> Ordering {
        match self.arr().cmp(&other.arr()) {
            Ordering::Equal =>
                match self.d.cmp(&other.d) { // zero delay edges after
                    Ordering::Equal => return self.cmp(&other),
                    Ordering::Less => return Ordering::Greater,
                    Ordering::Greater => return Ordering::Less,
                },
            lg => return lg,
        } 
    }
}

impl FromStr for TEdge {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let mut words = s.split_whitespace();
        let u: Node = TEdge::parse(words.next())?;
        let v: Node = TEdge::parse(words.next())?;
        let t: Time = TEdge::parse(words.next())?;
        if t <= Time::MIN { 
            log::info!("Time t={} is expected to be greater than Time::MIN={} \
                        (which is considered as infinite Fastest cost).", 
                       t, Time::MIN); 
        }
        let d: Time = 
            match words.next() {
                Some(s) => 
                    s.parse().map_err(|_| format!("bad number format in `{}`", s))?,
                None => 1,
            };
        if d < 0 {
            Err(format!("negative delay: {d}"))
        } else {
            Ok(TEdge { u, v, t, d })
        }
    }
}



impl TGraph {

    pub fn new(mut earr: Vec<TEdge>) -> Self {
        let m: usize = earr.len();
        if m > Eind::max as usize { panic!("graph too large, need larger Eind type!") } // if using Eind = u32
        let m = m as Eind;
        let mut n:  Node = 0; // number of nodes
        for e in &earr {
            if e.u + 1 > n { n = e.u + 1; }
            if e.v + 1 > n { n = e.v + 1; }
        }

        earr.sort_by(|e,f| e.cmp_arr(&f));

        // topological sort of zero delay edges:
        let mut index: HashMap<Node, graph::Node> = HashMap::new();
        let mut index_orig: Vec<Node> = vec![n as Node; n];
        let mut topord_rank: Vec<Node> = vec![n as Node; n];
        let mut i = 0;
        while i < earr.len() {
            let ei = &earr[i];
            if ei.d == 0 { // zero delay
                let mut j = i+1;
                while j < earr.len() && ei.arr() == earr[j].arr() {
                    assert_eq!(earr[j].d, 0);
                    j += 1;
                }
                if j > i+1 { // at least 2 zero delay edges
                    // Compute graph
                    let mut g = graph::LabDigraph::new(j-i);
                    index.clear();
                    let mut n: graph::Node = 0;
                    let mut index_nb = 
                        |v| -> graph::Node {
                            if ! index.contains_key(v) {
                                let i = n;
                                index.insert(*v, i);
                                index_orig[i] = *v; 
                                n += 1;
                                i
                            } else {
                                * index.get(v).unwrap()
                            }
                        };
                    for (k, e) in earr[i..j].iter().enumerate() {
                        g.add_arc(index_nb(&e.u), index_nb(&e.v), k);
                    }
                    // sort
                    let mut rank: Node = 0;
                    for u in graph::topological_sort(&g) {
                        let v = index_orig[u];
                        topord_rank[v] = rank;
                        rank += 1;
                    }
                    earr[i..j].sort_by_key(|e| topord_rank[e.u]);
                }
                i = j;
            } else {
                i += 1;
            }
        }
        // end topological sort

        let mut edep_ind: Vec<Eind> = (0..m).map(|i| i as Eind).collect();
        edep_ind.sort_by(|&i,&j| earr[i].cmp_u_t(&earr[j])); // stable sort
        
        let mut u_fst: Vec<Eind> = vec![0; n+1]; // store out-deg and then index of first out-edge
        for e in &earr {
            u_fst[e.u+1] += 1; // out-deg in next cell
        }
       for u in 0..n { // prefix sum
            u_fst[u+1] += u_fst[u];
        }
        assert_eq!(u_fst[n], m as Eind);

        // store edges in edep rather than earr
        let mut earr_ind = vec![0 as Eind; m];
        for (i, &j) in edep_ind.iter().enumerate() {
            earr_ind[j] = i as Eind;
        }
        let check = earr.clone();
        earr.sort_by(|e, f| e.cmp_u_t(f)); // stable sort
        for (i, &j) in edep_ind.iter().enumerate() {
            assert_eq!(&check[j], &earr[i]);
        }

        Self { n, m, edep: earr, u_fst, earr: earr_ind }
    }

    /// Returns for each tedge (in order of earr) a couple (l,r) where corresponding
    /// edges in edep (those in edep[l..r]) beta-extend it (can follow it in a in beta-restless walk).
    pub fn extend_indexes(&self, beta: Time) -> Vec<(Eind, Eind)> {
        let mut l_v = self.u_fst.clone(); // included
        let mut r_v = self.u_fst.clone(); // not included
        let mut ind: Vec<(Eind, Eind)> = Vec::with_capacity(self.m);
        for &i in self.earr.iter() {
            let e = & self.edep[i];
            let v = e.v;
            let arr = e.arr();
            let mut l = l_v[v];
            while l < self.u_fst[v+1] && self.edep[l].t < arr { l += 1; }
            l_v[v] = l;
            let mut r = max(r_v[v], l);
            while r < self.u_fst[v+1] && self.edep[r].t - arr <= beta { r += 1; }
            r_v[v] = r;
            ind.push((l,r)); 
        }
        ind
    }

}

impl FromStr for TGraph {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
            let tedges: Vec<TEdge> = s.lines()
                .map(|s| s.parse::<TEdge>())
                .collect::<Result<Vec<TEdge>, Self::Err>>()?;
            Ok(TGraph::new(tedges))   
    }
}


#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    fn check_usize() {
        //assert!(std::mem::size_of::<Node>() <= std::mem::size_of::<usize>());
        assert!(Node::BITS <= usize::BITS);
    }

    use super::TGraph;

    #[allow(dead_code)]
    pub const BLACKBOARD: &str = "\
1 2 1
1 2 2
2 3 5
2 4 2
4 5 3
5 2 4
";

    #[allow(dead_code)]
    pub const LZO: &str = "\
1 2 1
2 3 2
1 3 3
3 4 3
3 5 4
5 6 5
3 6 6
";

    #[allow(dead_code)]
    pub const ROMBUS_4: &str = "\
1 2 1
2 3 2
2 4 2
3 2 3
4 2 3
2 3 4
2 4 4
3 5 5
4 5 5
";

    #[allow(dead_code)]
    pub const WHITEBOARD: &str = "\
1 2 1 1
1 2 2 1
1 2 4 1
2 3 4 1
2 3 5 1
3 4 5 1
3 4 6 1
3 4 7 1
1 3 8 1
";

    #[allow(dead_code)]
    pub const ZERO: &str = "\
2 1 5 0
3 2 5 0
4 3 5 0
4 5 5 0
5 2 5 0
6 4 5 0
3 6 4 1
5 1 5 1
";

    #[allow(dead_code)]
    pub const TEST_TGRAPHS_STR: [&str; 4] = [BLACKBOARD, LZO, ROMBUS_4, WHITEBOARD];

    #[test]
    fn indexes () {
        for tg_str in TEST_TGRAPHS_STR {
            let tg: TGraph = tg_str.parse().unwrap();
            check_tgraph(&tg);
            for beta in [0, 1, 2, 100] {
                //eprintln!("graph={}, beta={}", tg_str, beta);
                let succ = tg.extend_indexes(beta);
                check_extend_indexes(&tg, beta, &succ);
            }
        }
    }

    #[allow(dead_code)]
    pub fn check_tgraph(tg: &TGraph) {
        let mut j: Eind = 0;
        for u in 0..tg.n {
            for i in tg.u_fst[u]..tg.u_fst[u+1] {
                assert_eq!(i, j);
                j += 1;
                let e = &tg.edep[i];
                assert_eq!(u, e.u);
            }
        }
    }

    #[allow(dead_code)]
    pub fn check_extend_indexes(tg: &TGraph, beta: Time, succ: &Vec<(Eind, Eind)>) -> () {
        eprintln!("{:?}", succ);
        // brute force check:
        for (i, &ei) in tg.earr.iter().enumerate() {
            let e = &tg.edep[ei];
            let &(l, r) = &succ[i];
            for j in 0..tg.m {
                let f = &tg.edep[j];
                assert_eq!(f.extends(e, beta), l <= j && j < r);
            }
        }
        // check non-decreasingness of l, and r:
        let mut l_v = tg.u_fst.clone(); 
        let mut r_v = tg.u_fst.clone();
        for (i, &ei) in tg.earr.iter().enumerate() {
            let &(l, r) = &succ[i];
            let v = tg.edep[ei].v;
            assert!(l >= l_v[v] && r >= r_v[v]);
            l_v[v] = l;
            r_v[v] = r;
        }
    }
}