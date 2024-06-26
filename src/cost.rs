use crate::tgraph::{Time, TEdge};

use std::ops::Add;
use std::cmp::Reverse;

pub trait Cost : Ord + Add<Output = Self> + Sized + Clone + std::fmt::Debug {
    type TargetCost : Ord + Clone + std::fmt::Debug;

    fn infinite_cost() -> Self;
    fn infinite_target_cost() -> Self::TargetCost;
    fn empty_target_cost() -> Self::TargetCost; // target cost of an empty walk
    fn edge_cost(e: &TEdge) -> Self;
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost;
    fn target_cost_to_distance(c : Self::TargetCost) -> f64;
}



/// Cost of a temporal walk Q as its arrival time arr(Q) (for minimization).
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct Foremost(pub bool); // true for infinite, false otherwise

impl Add for Foremost {
    type Output = Self; 
    fn add(self, _: Self) -> Self::Output { self }
}

impl Cost for Foremost {
    type TargetCost = Time; // arr(Q)
    fn infinite_cost() -> Self { Foremost(true) }
    fn infinite_target_cost() -> Self::TargetCost { Time::MAX }
    fn empty_target_cost() -> Self::TargetCost { Time::MIN }
    fn edge_cost(_: &TEdge) -> Self { Foremost(false) }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { e.arr() }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { c as f64 }
}



/// Cost of a temporal walk Q as its arrival time arr(Q) (for maximization).
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct Latest(pub bool); // true for not reachable, false otherwise

impl Add for Latest {
    type Output = Self; 
    fn add(self, _: Self) -> Self::Output { self }
}

impl Cost for Latest {
    type TargetCost = Reverse<Time>; // arr(Q)
    fn infinite_cost() -> Self { Latest(true) }
    fn infinite_target_cost() -> Self::TargetCost { Reverse(Time::MIN) }
    fn empty_target_cost() -> Self::TargetCost { Reverse(Time::MAX) }
    fn edge_cost(_: &TEdge) -> Self { Latest(false) }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { Reverse(e.arr()) }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { -(c.0 as f64) }
}



/// Cost of a temporal walk Q as its duration (for minimization).
///   (The duration of Q is arr(Q) - dep(Q)).
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct Fastest(pub Reverse<Time>); // Cost of Q is - dep(Q)

impl Add for Fastest {
    type Output = Self;
    fn add(self, _: Self) -> Self::Output { self }
}

impl Cost for Fastest {
    type TargetCost = Time; // TargetCost of Q is arr(Q) - dep(Q)
    fn infinite_cost() -> Self { Fastest(Reverse(Time::MIN)) }
    fn infinite_target_cost() -> Self::TargetCost { Time::MAX }
    fn empty_target_cost() -> Self::TargetCost { 0 as Time }
    fn edge_cost(e: &TEdge) -> Self { 
        if e.t <= Time::MIN { 
            panic!("Time {} is expected to be greater than Time::MIN={} \
                    (which is considered as infinite  Fastest cost).", 
                    e.t, Time::MIN); 
        }
        Fastest(Reverse(e.t))
    }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { 
        if self.0.0 == Time::MIN { Self::infinite_target_cost() } // Time::MIN is consdiered infinite
        else { e.arr() - self.0.0 }
    }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { c as f64 } // +1.0 if possible zero delay edges ?
}

/// Cost of a temporal walk Q as its sum of delays `sum_{e in Q} e.d`.
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct DelaySum(pub Time); // Cost of Q is sum_{e in Q} e.d

impl Add for DelaySum {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output { DelaySum(self.0 + rhs.0) }
}

impl Cost for DelaySum {
    type TargetCost = Time;
    fn infinite_cost() -> Self { DelaySum(Time::MAX) }
    fn infinite_target_cost() -> Self::TargetCost { Time::MAX }
    fn empty_target_cost() -> Self::TargetCost { 0 as Time }
    fn edge_cost(e: &TEdge) -> Self { DelaySum(e.d) }
    fn target_cost(&self, _: &TEdge) -> Self::TargetCost { self.0 }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { c as f64 } // +1.0 if possible zero delay edges ?
}

/// Cost of a temporal walk Q as its waiting time (for minimization). 
///   (The waiting time of Q is arr(Q) - dep(Q) - sum_{e \in Q} e.d).
#[derive(Copy, Clone, Debug)]
pub struct Waiting(pub Time, pub Time); // dep(Q), sum_{e \in Q} e.d

impl Add for Waiting {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output { Waiting(self.0 , self.1 + rhs.1) }
}

impl Cost for Waiting {
    type TargetCost = Time; // arr(Q) - dep(Q) - sum_{e \in Q}
    fn infinite_cost() -> Self { Waiting(Time::MIN, 0 as Time) }
    fn infinite_target_cost() -> Self::TargetCost { Time::MAX }
    fn empty_target_cost() -> Self::TargetCost { 0 as Time }
    fn edge_cost(e: &TEdge) -> Self { 
        if e.t <= Time::MIN { 
            panic!("Time {} is expected to be greater than Time::MIN={} \
                    (which is considered as infinite (negative) departure time).", 
                    e.t, Time::MIN); 
        }
        Waiting(e.t, e.d) 
    }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { 
        if self.0 == Time::MIN { Self::infinite_target_cost() } // Time::MIN is consdired infinite
        else { e.arr() - self.0 - self.1 }
    }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { 1.0 + c as f64 } // min. waiting time can be zero
}

impl PartialEq for Waiting {
    fn eq(&self, other: &Self) -> bool {
        self.0 + self.1 == other.0 + other.1
    }
}
impl Eq for Waiting {
}
impl PartialOrd for Waiting {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        let s: Time = self.0 + self.1;
        let so: Time = other.0 + other.1;
        Some(s.cmp(&so).reverse())
    }
}
impl Ord for Waiting {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        let s: Time = self.0 + self.1;
        let so: Time = other.0 + other.1;
        s.cmp(&so).reverse()
    }
}


pub type Hop = u32;

/// Cost of a temporal walk Q as its hop count (for minimization). 
///   (The hop count hop(Q) is its number of temporal edges).
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct Shortest(pub Hop); // hop(Q)

impl Add for Shortest {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output { Shortest(self.0 + rhs.0) }
}

impl Cost for Shortest {
    type TargetCost = Hop;
    fn infinite_cost() -> Self { Shortest(Hop::MAX) }
    fn infinite_target_cost() -> Self::TargetCost { Hop::MAX }
    fn empty_target_cost() -> Self::TargetCost { 0 }
    fn edge_cost(_: &TEdge) -> Self { Shortest(1) }
    fn target_cost(&self, _: &TEdge) -> Self::TargetCost { self.0 }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { c as f64 }
}


/// Cost of a temporal walk Q as (arr(Q), hop(Q)) (for minimization by lexicographic order). 
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct ShortestForemost(pub u32); // hop(Q)

impl Add for ShortestForemost {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output { ShortestForemost(self.0 + rhs.0) }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct ShortestForemostTarget(pub Time, pub Hop); // arr(Q), hop(Q)

impl Cost for ShortestForemost {
    type TargetCost = ShortestForemostTarget;
    fn infinite_cost() -> Self { ShortestForemost(Hop::MAX) }
    fn infinite_target_cost() -> Self::TargetCost { ShortestForemostTarget(Time::MAX, Hop::MAX) }
    fn empty_target_cost() -> Self::TargetCost { ShortestForemostTarget(Time::MIN, 0) }
    fn edge_cost(_: &TEdge) -> Self { ShortestForemost(1) }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { ShortestForemostTarget(e.arr(), self.0) }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { c.0 as f64 }
}



/// Cost of a temporal walk Q as (arr(Q), hop(Q)) (for minimization by lexicographic order). 
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct ShortestLatest(pub Hop); // hop(Q)

impl Add for ShortestLatest {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output { ShortestLatest(self.0 + rhs.0) }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct ShortestLatestTarget(pub Reverse<Time>, pub Hop); // arr(Q), hop(Q)

impl Cost for ShortestLatest {
    type TargetCost = ShortestLatestTarget;
    fn infinite_cost() -> Self { ShortestLatest(Hop::MAX) }
    fn infinite_target_cost() -> Self::TargetCost { ShortestLatestTarget(Reverse(Time::MIN), Hop::MAX) }
    fn empty_target_cost() -> Self::TargetCost { ShortestLatestTarget(Reverse(Time::MAX), 0) }
    fn edge_cost(_: &TEdge) -> Self { ShortestLatest(1) }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { ShortestLatestTarget(Reverse(e.arr()), self.0) }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { -(c.0.0 as f64) }
}



/// Generic shortest among minimum walks for another criterion
#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct ShortestOf<C: Cost>(pub C, pub Hop); // lexicographic min of (cost(Q), hop(Q))

impl<C: Cost> Add for ShortestOf<C> {
    type Output = Self;
    fn add(self, rhs: Self) -> Self::Output { ShortestOf(self.0 + rhs.0, self.1 + rhs.1) }
}

#[derive(PartialEq, Eq, PartialOrd, Ord, Copy, Clone, Debug)]
pub struct ShortestOfTarget<C: Cost>(pub C::TargetCost, pub Hop); // targetcost(Q), hop(Q)

impl<C: Cost> Cost for ShortestOf<C> {
    type TargetCost = ShortestOfTarget<C>;
    fn infinite_cost() -> Self { ShortestOf(C::infinite_cost(), Hop::MAX) }
    fn infinite_target_cost() -> Self::TargetCost { ShortestOfTarget(C::infinite_target_cost(), Hop::MAX) }
    fn empty_target_cost() -> Self::TargetCost { ShortestOfTarget(C::empty_target_cost(), 0) }
    fn edge_cost(e: &TEdge) -> Self { ShortestOf(C::edge_cost(e), 1) }
    fn target_cost(&self, e: &TEdge) -> Self::TargetCost { ShortestOfTarget(C::target_cost(&self.0, e), self.1) }
    fn target_cost_to_distance(c : Self::TargetCost) -> f64 { C::target_cost_to_distance(c.0) }
}

pub type ShortestFastest = ShortestOf<Fastest>;
pub type ShortestWaiting = ShortestOf<Waiting>;
pub type ShortestDelaySum = ShortestOf<DelaySum>;

#[cfg(test)]
mod cost_tests {
    use super::*;
    use crate::tgraph::*;

    #[test]
    fn cost_foremost() {
        let c = Foremost(false);
        let d = Foremost(false);
        assert!(c <= d);
        assert!( ! (c < d) );
        assert!(c == d);
        let e = TEdge{ u: 1, v: 2, t: 1, d: 0 };
        let c = Foremost::edge_cost(&e);
        let inf = Foremost::infinite_cost();
        assert!(c < inf);
    }

    #[test]
    #[should_panic] //(expected = "time of edge greater than Time::MIN")
    fn cost_shortest_fastest_edge_cost() {
        let e = TEdge{ u: 1, v: 2, t: Time::MIN, d: 0 };
        let _ = ShortestFastest::edge_cost(&e);
    }
    #[test]
    #[should_panic] //(expected = "time of edge greater than Time::MIN")
    fn cost_fastest_edge_cost() {
        let e = TEdge{ u: 1, v: 2, t: Time::MIN, d: 0 };
        let _ = Fastest::edge_cost(&e);
    }
}
