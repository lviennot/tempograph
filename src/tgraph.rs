use crate::graph;
use std::collections::HashMap;

pub type Node = u32;
pub type Time = u32; // Time:MIN is reserved as infinite for fastest cost, use i32 if 0 is needed
pub type Eind = u32; // tedge index

/// a temporal edge
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

#[derive(Clone)]
pub struct TGraph {
    pub n: Node, // largest node number
    pub m: Eind, // number of temporal edges
    pub earr: Vec<TEdge>, // tedges sorted by arrival time
    pub edep: Vec<Eind>, // indexes (in earr) of tedges sorted by tail and departure time
    pub u_fst: Vec<Eind>, // indexes of tedges from u are in edep[u_fst[u]..u_fst[u+1]]
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
        let d: Time = 
            match words.next() {
                Some(s) => 
                    s.parse().map_err(|_| format!("bad number format in `{}`", s))?,
                None => 1,
            };
        Ok(TEdge { u, v, t, d })
    }
}



impl TGraph {

    pub fn new(mut earr: Vec<TEdge>) -> Self {
        let m: usize = earr.len();
        if m > Eind::max as usize { panic!("graph too large, need larger Eind type!") }
        let m = m as Eind;
        let mut n:  Node = 0; // max number of a node
        for e in &earr {
            if e.u > n { n = e.u; }
            if e.v > n { n = e.v; }
        }

        earr.sort_by(|e,f| e.cmp_arr(&f));

        // topological sort of zero delay edges:
        let mut index: HashMap<Node, graph::Node> = HashMap::new();
        let mut index_orig: Vec<Node> = vec![n as Node; n as usize + 1];
        let mut topord_rank: Vec<Node> = vec![n as Node; n as usize + 1];
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
                        topord_rank[v as usize] = rank;
                        rank += 1;
                    }
                    earr[i..j].sort_by_key(|e| topord_rank[e.u as usize]);
                }
                i = j;
            } else {
                i += 1;
            }
        }
        // end topological sort

        let mut edep: Vec<Eind> = (0..m).map(|i| i as Eind).collect();
        edep.sort_by(|&i,&j| earr[i as usize].cmp_u_t(&earr[j as usize]));
        
        let mut u_fst: Vec<Eind> = vec![0; (n+2) as usize]; // store out-deg and then index of first out-edge
        for e in &earr {
            u_fst[(e.u+1) as usize] += 1; // out-deg in next cell
        }
       for u in 0..=n as usize { // prefix sum
            u_fst[u+1] += u_fst[u];
        }
        assert_eq!(u_fst[(n+1) as usize], m as Eind);

        Self { n, m, earr, edep, u_fst }
    }

    fn edep_e(&self, i: Eind) -> &TEdge {
        &self.earr[self.edep[i as usize] as usize]
    } 

    /// Returns for each tedge a triplet (i,l,r) where i is its index in edep,
    /// and edges in edep (those in edep[l..r]) beta-extend it (can follow it in a in beta-restless walk).
    pub fn extend_indexes(&self, beta: Time) -> Vec<(Eind, Eind, Eind)> {
        let mut l_v = self.u_fst.clone(); // included
        let mut r_v = self.u_fst.clone(); // not included
        let mut ind: Vec<(Eind, Eind, Eind)> = Vec::with_capacity(self.m as usize);
        for e in &self.earr {
            let v = e.v as usize;
            let arr = e.arr();
            let mut l = l_v[v];
            while l < self.u_fst[v+1] && self.edep_e(l).t < arr { l += 1; }
            l_v[v] = l;
            let mut r = max(r_v[v], l);
            while r < self.u_fst[v+1] && self.edep_e(r).t - arr <= beta { r += 1; }
            r_v[v] = r;
            ind.push((0,l,r)); // idnex of e in edep not known yet
        }
        for (i, &j) in self.edep.iter().enumerate() {
            ind[j as usize].0 = i as Eind;
        }
        ind
    }

    // Counts for each tedge how many tedges it extends according to interval indexes given in `succ`.
    #[allow(dead_code)]
    fn extend_in_degrees(&self, succ: &Vec<(Eind, Eind, Eind)>) -> Vec<Eind> {
        let mut deg: Vec<Eind> = vec![0; self.m as usize];
        let mut minus: Vec<Eind> = vec![0; self.m as usize + 1];
        for &(_, l, r) in succ {
            if l < r {
                deg[l as usize] += 1; // one more interval over l..
                minus[r as usize] += 1; // one less interval from r..
            }
        }
        assert_eq!(minus[0], 0);
        // prefix sum:
        for i in 1..self.m as usize {
            deg[i] += deg[i-1];
            deg[i] -= minus[i];
        }
        assert_eq!(deg[self.m as usize - 1], minus[self.m as usize]);
        deg
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
                //println!("{:?}", tg.u_fst);
                let succ = tg.extend_indexes(beta);
                check_extend_indexes(&tg, beta, &succ);
                let deg = tg.extend_in_degrees(&succ);
                assert_eq!(deg.iter().sum::<Eind>(), succ.iter().map(|&(_,l,r)| r-l).sum::<Eind>());
            }
        }
    }

    #[allow(dead_code)]
    pub fn check_tgraph(tg: &TGraph) {
        let mut j: Eind = 0;
        for u in 0..tg.n as usize {
            for i in tg.u_fst[u]..tg.u_fst[u+1] {
                assert_eq!(i, j);
                j += 1;
                let e = tg.edep_e(i);
                assert_eq!(u, e.u as usize);
            }
        }
    }

    #[allow(dead_code)]
    pub fn check_extend_indexes(tg: &TGraph, beta: Time, succ: &Vec<(Eind, Eind, Eind)>) -> () {
        //for e in tg.edep.iter().map(|&i| &tg.earr[i]) { println!("{:?}", e) }
        // brute force check:
        for (i, e) in tg.earr.iter().enumerate() {
            let &(j, l, r) = &succ[i];
            assert_eq!(i, tg.edep[j as usize] as usize);
            for j in 0..tg.m as Eind {
                let f = tg.edep_e(j);
                assert_eq!(f.extends(e, beta), l <= j && j < r);
            }
        }
        // check non-decreasingness of l, and r:
        let mut l_v = tg.u_fst.clone(); 
        let mut r_v = tg.u_fst.clone();
        for (i, e) in tg.earr.iter().enumerate() {
            let &(_, l, r) = &succ[i];
            let v = e.v as usize;
           assert!(l >= l_v[v] && r >= r_v[v]);
            l_v[v] = l;
            r_v[v] = r;
        }
    }
}