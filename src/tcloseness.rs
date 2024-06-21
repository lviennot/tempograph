use crate::tgraph::*;
use crate::tsweep::*;
use crate::cost::*;

use dsi_progress_logger::*;

pub fn closeness(tg: &TGraph, beta: Time) -> Vec<f64> {
    let mut tsweep: TSweep<Shortest> = TSweep::new(&tg);
    let infty = Shortest::infinite_target_cost();
    let mut pl = ProgressLogger::default();
    pl.expected_updates(Some(tg.n as usize))
        .display_memory(true)
        .item_name("node");
    pl.start("Computing closeness...");
    (0..tg.n).map(|s| -> f64 {
        tsweep.clear();
        tsweep.scan(s, beta);
        let opt_costs = tsweep.opt_costs(s);
        let mut hc = 0.;
        for t in 0..tg.n {
            if t != s && opt_costs[t] < infty { hc += 1. / opt_costs[t] as f64 }
        }
        pl.update();
        if s == tg.n - 1 { pl.stop() }
        hc
    }).collect()
}

