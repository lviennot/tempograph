pub type Node = usize;

pub trait Graph : std::fmt::Display {
    /// Returns an iterator over the neighbors of node `u`.
    fn neighbors<'a>(&'a self, u: Node) -> impl Iterator<Item = &Node> + 'a;

    /// Returns the number of nodes.
    fn n(&self) -> usize;

    /// Print the graph.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for u in 0..self.n() as Node {
            write!(f, "{u} :")?;
            for &v in self.neighbors(u) {
                write!(f, " {v}")?;
            }
            write!(f, "\n")?;
        }
        Ok(())
    }
}

pub struct LabDigraph<L> {
    adj: Vec<Vec<(Node, L)>>,
}

impl<L> Graph for LabDigraph<L> { 
    fn neighbors<'a>(&'a self, u: Node) -> impl Iterator<Item = &Node> + 'a {
        self.adj[u].iter().map(|(v, _)| v)
    }

    fn n(&self) -> usize {
        self.adj.len()
    }
}

impl<L> std::fmt::Display for  LabDigraph<L> { 
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Graph::fmt(self, f)
    }
}

impl<L> LabDigraph<L> {
    
    pub fn new(n: usize) -> Self {
        LabDigraph { adj: Vec::with_capacity(n) }
    }

    pub fn add_arc(&mut self, u: Node, v: Node, l: L) {
        let m = std::cmp::max(u, v) as usize;
        if m >= self.adj.len() { 
            self.adj.resize_with(m+1, || Vec::with_capacity(0));             
        }
        self.adj[u as usize].push((v, l));
    }

    pub fn neighbors_labels<'a>(&'a self, u: Node) -> impl Iterator<Item = &L> + 'a {
        self.adj[u].iter().map(|(_, l)| l)
    }

}


pub struct DFS {
    stack: Vec<Node>,
    in_stack: Vec<bool>,
    pub number: Vec<Node>,
    pub exit_order: Vec<Node>,
    num: usize,
}

impl DFS {

    const INFINITY: Node = Node::MAX;

    pub fn new(n: usize) -> Self {
        let stack = Vec::with_capacity(n);
        let in_stack = vec![false; n];
        let number = vec![Self::INFINITY; n];
        let exit_order = Vec::with_capacity(n);
        DFS { stack, in_stack, number, exit_order, num: 0 }
    }

    fn check_size(&self, n: usize) {
        if n > self.number.len() { 
            panic!("Graph of {} nodes is too large for this DFS struct of size {}!", 
                   n, self.number.len()) 
        }
    }

    pub fn dfs<G: Graph>(&mut self, g: &G, u: Node) {
        self.stack.push(u);
        while let Some(u) = self.stack.pop() {
            if self.number[u as usize] == Self::INFINITY {
                self.number[u as usize] = self.num as Node;
                self.num += 1;
                self.in_stack[u] = true;
                self.stack.push(u); // to later detect end of visit
                for &v in g.neighbors(u) {
                    if self.number[v as usize] == Self::INFINITY {
                        self.stack.push(v)
                    } 
                }
            } else if self.in_stack[u] { // end of visit
                self.in_stack[u] = false;
                self.exit_order.push(u)
            }
        }
    }

    pub fn dfs_all<G: Graph>(&mut self, g: &G) {
        self.check_size(g.n());
        for v in 0..g.n() as Node {
            if self.number[v as usize] == Self::INFINITY {
                self.dfs(g, v);
            }
        }
        assert_eq!(self.exit_order.len(), g.n());
    }

}

pub fn topological_sort<G: Graph>(g: &G) -> Vec<Node> {
    let mut dfs = DFS::new(g.n());
    dfs.dfs_all(g);
    dfs.exit_order.into_iter().rev().collect()
}

pub fn has_topological_order<G: Graph>(g: &G, ordering: &Vec<Node>) -> bool {
    let mut rank = vec![g.n(); g.n()];
    for (i, &v) in ordering.iter().enumerate() {
        rank[v] = i;
    } 
    for u in 0..g.n() as Node {
        for &v in g.neighbors(u) {
            if rank[u] >= rank[v] {
                return false
            }
        }
    }
    return true
}

pub fn acyclic<G: Graph>(g: &G) -> bool {
    let ordering = topological_sort(g);
    has_topological_order(g, &ordering)
}

// Tarjan, Robert E. (1972), "Depth-first search and linear graph algorithms", SIAM Journal on Computing, 1 (2): 146–160, CiteSeerX 10.1.1.327.8418, doi:10.1137/0201010
// TODO