//! API for working with POA graphs. Supports graph building  
//! from sequence iterators and consensus extraction from graphs.

// dependencies
use petgraph::graph::IndexType as PetgraphIndexType;
use petgraph::visit::EdgeRef;
use petgraph::prelude::NodeIndex;
use petgraph::Direction;
use serde::de::DeserializeOwned;
use crate::errors::PoastaError;
use crate::graphs::poa::POAGraph;
use crate::graphs::AlignableRefGraph;
use crate::aligner::config::{AffineMinGapCost};
use crate::aligner::scoring::{AlignmentType, GapAffine};
use crate::aligner::PoastaAligner;

impl<Ix> POAGraph<Ix>
where
    Ix: PetgraphIndexType + DeserializeOwned,
{

    /// Extract the consensus sequence from a POA graph. 
    pub fn traverse_heaviest_bundle(&self) -> Vec<u8>{

        // Rust adaptation of https://github.com/rvaser/spoa/blob/master/src/graph.cpp#L466
        let node_count = self.node_count_with_start_and_end();
        let mut predecessor_nodes: Vec<Option<NodeIndex<Ix>>> = vec!(None; node_count);
        let mut scores: Vec<usize> = vec!(0; node_count);
        let mut max_node: NodeIndex<Ix> = self.start_node();

        // forward pass to calculate cumulative weights and find max scoring node
        self.get_topological_sorted().iter().for_each(|node|{
            let node_i = node.index();
            self.graph.edges_directed(*node, Direction::Incoming).for_each(|edge_in|{
                let node_in = edge_in.source();
                let weight_in = edge_in.weight().weight; // always a positive non-zero integer
                if scores[node_i] < weight_in || (
                    scores[node_i] == weight_in && 
                    predecessor_nodes[node_i].is_some() &&
                    scores[predecessor_nodes[node_i].unwrap().index()] <= scores[node_in.index()]
                ){
                    scores[node_i] = weight_in;
                    predecessor_nodes[node_i] = Some(node_in);
                }
            });
            if let Some(predecessor_node) = predecessor_nodes[node_i]  {
                scores[node_i] += scores[predecessor_node.index()];
            }
            if scores[max_node.index()] < scores[node_i] {
                max_node = *node;
            }
        });

        // fill out from a max node that is not a terminal node
        if !self.is_terminal_node(max_node){
            let node_ranks = self.get_node_ranks();
            while !self.is_terminal_node(max_node) {
                max_node = self.branch_completion(
                    max_node,
                    node_ranks[max_node.index()], 
                    &mut scores,
                    &mut predecessor_nodes,
                );
            }
        }

        // traceback to finalize the best-path consensus
        let mut consensus: Vec<u8> = Vec::with_capacity(node_count);
        consensus.push(self.get_symbol(max_node)); // always a terminal node with no outgoing edges
        while let Some(predecessor_node) = predecessor_nodes[max_node.index()] {
            consensus.push(self.get_symbol(predecessor_node));
            max_node = predecessor_node;
        }
        consensus.reverse();
        consensus
    }

    /// Determine whether any given node is a right-sided terminal node.
    pub fn is_terminal_node(&self, node: NodeIndex<Ix>) -> bool {
        let edges_out: Vec<_> = self.graph.edges(node).collect();
        !edges_out.iter().any(|edge| edge.weight().weight > 0)
    }

    /// Fill out from max_node to a node with no outgoing edges.
    fn branch_completion(
        &self, 
        mut max_node: NodeIndex<Ix>,
        max_node_rank: usize,
        scores: &mut Vec<usize>,
        predecessor_nodes: &mut Vec<Option<NodeIndex<Ix>>>,
    ) -> NodeIndex<Ix> {

        // reset scores distal relevant to all paths from max_node
        self.graph.edges(max_node).for_each(|max_edge_out|{
            self.graph.edges_directed(max_edge_out.target(), Direction::Incoming).for_each(|edge_in|{ 
                let node_in = edge_in.source();
                if node_in != max_node {
                    scores[node_in.index()] = 0;
                }
            });
        });

        // re-run the logic of traverse_heaviest_bundle with appropriate resets and restrictions
        self.get_topological_sorted().iter().skip(max_node_rank + 1).for_each(|node|{
            let node_i = node.index();
            scores[node_i] = 0;
            predecessor_nodes[node_i] = None;
            self.graph.edges_directed(*node, Direction::Incoming).for_each(|edge_in|{
                let node_in = edge_in.source();
                if scores[node_in.index()] > 0 {
                    let weight_in = edge_in.weight().weight; // always a positive non-zero integer
                    if scores[node_i] < weight_in || (
                        scores[node_i] == weight_in && 
                        predecessor_nodes[node_i].is_some() &&
                        scores[predecessor_nodes[node_i].unwrap().index()] <= scores[node_in.index()]
                    ){
                        scores[node_i] = weight_in;
                        predecessor_nodes[node_i] = Some(node_in);
                    }                    
                }
            });
            if let Some(predecessor_node) = predecessor_nodes[node_i]  {
                scores[node_i] += scores[predecessor_node.index()];
            }
            if scores[max_node.index()] < scores[node_i] {
                max_node = *node;
            }
        });
        max_node
    }
}

/// Methods for building and extracting consensus sequences.
pub trait PoastaConsensus<'a> {

    /// Build a consensus sequence from a POA graph.
    fn build_consensus(
        &mut self,
        cost_mismatch:   u8,
        cost_gap_extend: u8,
        cost_gap_open:   u8,
        aln_type: AlignmentType,
    ) -> Result<Vec<u8>, PoastaError>;

    /// Build a consensus sequence from a POA graph using default parameters.
    fn build_consensus_with_defaults(&mut self) -> Result<Vec<u8>, PoastaError>{
        self.build_consensus(
            4, 
            2, 
            6, 
            AlignmentType::Global
        )
    }
}

impl<'a, I> PoastaConsensus<'a> for I 
where 
    I: Iterator<Item = (&'a [u8], Option<&'a [usize]>)> 
{
    fn build_consensus(
        &mut self,
        cost_mismatch:   u8,
        cost_gap_extend: u8,
        cost_gap_open:   u8,
        aln_type: AlignmentType,
    ) -> Result<Vec<u8>, PoastaError> {

        // check for something to do
        if let Some((seq, weights)) = self.next() {
            if seq.is_empty() {
                return Err(PoastaError::Other);
            }

            // initialize the aligner and graph
            let costs = AffineMinGapCost(GapAffine::new(
                cost_mismatch, 
                cost_gap_extend,
                cost_gap_open, 
            ));
            let aligner = PoastaAligner::new(
                costs, 
                aln_type
            );
            let mut graph: POAGraph<u32> = POAGraph::new();

            // add the first sequence to the graph
            let mut i = 1;
            if let Some(weights) = weights{
                graph.add_alignment_with_weights(&i.to_string(), seq, None, weights)?;
            } else {
                let weights = vec![1_usize; seq.len()];
                graph.add_alignment_with_weights(&i.to_string(), seq, None, &weights)?;
            }

            // add additional sequences to the graph
            while let Some((seq, weights)) = self.next() {
                i = i + 1;
                let result = aligner.align::<u32, _>(&graph, seq);
                if let Some(weights) = weights{
                    graph.add_alignment_with_weights(&i.to_string(), seq, Some(&result.alignment), weights)?;
                } else {
                    let weights = vec![1_usize; seq.len()];
                    graph.add_alignment_with_weights(&i.to_string(), seq, Some(&result.alignment), &weights)?;
                }
            }

            // extract and return the consensus
            Ok(graph.traverse_heaviest_bundle())

        // empty iterator throws an error, no sequences to process
        } else {
            return Err(PoastaError::Other);
        }
    }
}
