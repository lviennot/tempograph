use crate::tgraph::*;
use crate::tsweep::*;
use crate::cost::*;

pub fn closeness(tg: &TGraph, beta: Time) -> Vec<f64> {
    let mut tsweep: TSweep<Shortest> = TSweep::new(&tg);
    let infty = Shortest::infinite_target_cost();
    (0..tg.n).map(|s| -> f64 {
        tsweep.clear();
        tsweep.scan(s, beta);
        let opt_costs = tsweep.opt_costs();
        let mut hc = 0.;
        for t in 0..tg.n {
            if opt_costs[t] < infty { hc += 1. / opt_costs[t].0 as f64 }
        }
        hc
    }).collect()
}

