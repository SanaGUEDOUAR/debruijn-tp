#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""Perform assembly based on debruijn graph."""

import argparse
import os
import sys
import random
from random import randint
import statistics
import textwrap
from typing import Iterator, Dict, List
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
import networkx as nx
from networkx import (
    DiGraph,
    all_simple_paths,
    lowest_common_ancestor,
    has_path,
)
random.seed(9001)

matplotlib.use("Agg")

__author__ = "Sana GUEDOUAR"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Sana GUEDOUAR"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Sana GUEDOUAR"
__email__ = "sana.guedouar@etu.u-paris.fr"
__status__ = "Developpement"


def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments():  # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(
        description=__doc__, usage=f"{sys.argv[0]} -h"
    )
    parser.add_argument(
        "-i", dest="fastq_file", type=isfile, required=True, help="Fastq file"
    )
    parser.add_argument(
        "-k", dest="kmer_size", type=int, default=22, help="k-mer size (default 22)"
    )
    parser.add_argument(
        "-o",
        dest="output_file",
        type=Path,
        default=Path(os.curdir + os.sep + "contigs.fasta"),
        help="Output contigs in fasta file (default contigs.fasta)",
    )
    parser.add_argument(
        "-f", dest="graphimg_file", type=Path, help="Save graph as an image (png)"
    )
    return parser.parse_args()


def read_fastq(fastq_file: Path) -> Iterator[str]:
    """Extract reads from fastq files.

    :param fastq_file: (Path) Path to the fastq file.
    :return: A generator object that iterate the read sequences.
    """
    with fastq_file.open() as f:
        lines = iter(f)
        while True:
            try:
                next(lines)  # Skip the first line
                yield next(lines).strip()
                next(lines)  # Skip the third line
                next(lines)  # Skip the fourth line
            except StopIteration:
                break


def cut_kmer(read: str, kmer_size: int) -> Iterator[str]:
    """Cut read into kmers of size kmer_size.

    :param read: (str) Sequence of a read.
    :return: A generator object that provides the kmers (str) of size kmer_size.
    """
    for i in range(len(read) - kmer_size + 1):
        yield read[i : i + kmer_size]


def build_kmer_dict(fastq_file: Path, kmer_size: int) -> Dict[str, int]:
    """Build a dictionnary object of all kmer occurrences in the fastq file

    :param fastq_file: (str) Path to the fastq file.
    :return: A dictionnary object that identify all kmer occurrences.
    """
    kmer_dict = {}
    for read in read_fastq(fastq_file):
        for kmer in cut_kmer(read, kmer_size):
            if kmer in kmer_dict:
                kmer_dict[kmer] += 1
            else:
                kmer_dict[kmer] = 1
    return kmer_dict


def build_graph(kmer_dict: Dict[str, int]) -> DiGraph:
    """Build the de Bruijn graph from k-mers.

    :param kmer_dict: A dictionary object that identifies all k-mer occurrences.
    :return: A directed graph (DiGraph) with k-mer substrings 
            as nodes and occurrences as edge weights.
    """
    graph = nx.DiGraph()
    # Browse the dictionary of k-mers and their occurrences
    for kmer, occurrence in kmer_dict.items():
        # Calculate the prefix and suffix of the k-mer
        prefix = kmer[:-1]
        suffix = kmer[1:]
        # Add a link in the graph between the prefix
        # and the suffix with the weight of the occurrence
        graph.add_edge(prefix, suffix, weight=occurrence)
    return graph


def remove_paths(
    graph: DiGraph,
    path_list: List[List[str]],
    delete_entry_node: bool,
    delete_sink_node: bool,
) -> DiGraph:
    """Remove a list of path in a graph. A path is set of connected node in
    the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    for path in path_list:
        # Determine which nodes in the path to remove based on the flags
        nodes_to_remove = path[:]  # Start by assuming all nodes are to be removed
        if not delete_entry_node:
            # Keep the first node (entry node)
            nodes_to_remove = nodes_to_remove[1:]
        if not delete_sink_node:
            # Keep the last node (sink node)
            nodes_to_remove = nodes_to_remove[:-1]
        # Remove the specified nodes from the graph
        graph.remove_nodes_from(nodes_to_remove)
    return graph


def select_best_path(
    graph: DiGraph,
    path_list: List[List[str]],
    path_length: List[int],
    weight_avg_list: List[float],
    delete_entry_node: bool = False,
    delete_sink_node: bool = False,
) -> DiGraph:
    """Select the best path between different paths

    :param graph: (nx.DiGraph) A directed graph object
    :param path_list: (list) A list of path
    :param path_length_list: (list) A list of length of each path
    :param weight_avg_list: (list) A list of average weight of each path
    :param delete_entry_node: (boolean) True->We remove the first node of a path
    :param delete_sink_node: (boolean) True->We remove the last node of a path
    :return: (nx.DiGraph) A directed graph object
    """
    std_weight_avg = statistics.stdev(weight_avg_list) if len(weight_avg_list) > 1 else 0
    std_length = statistics.stdev(path_length) if len(path_length) > 1 else 0
    best_path_index = 0  # Initialize the variable
    if std_weight_avg > 0 :
        best_path_index = weight_avg_list.index(max(weight_avg_list))
    else:
        if std_length > 0:
            best_path_index = path_length.index(max(path_length))
        elif std_length == 0:
            best_path_index = randint(0, len(path_list) - 1)
    best_path = path_list[best_path_index]
    for path in path_list:
        if path != best_path:
            graph = remove_paths(graph, [path], delete_entry_node, delete_sink_node)
    return graph


def path_average_weight(graph: DiGraph, path: List[str]) -> float:
    """Compute the weight of a path

    :param graph: (nx.DiGraph) A directed graph object
    :param path: (list) A path consist of a list of nodes
    :return: (float) The average weight of a path
    """
    return statistics.mean(
        [d["weight"] for (u, v, d) in graph.subgraph(path).edges(data=True)]
    )


def solve_bubble(graph: DiGraph, ancestor_node: str, descendant_node: str) -> DiGraph:
    """Explore and solve bubble issue

    :param graph: (nx.DiGraph) A directed graph object
    :param ancestor_node: (str) An upstream node in the graph
    :param descendant_node: (str) A downstream node in the graph
    :return: (nx.DiGraph) A directed graph object
    """
    # Find all simple paths between ancestor and descendant
    paths = list(nx.all_simple_paths(graph, source=ancestor_node, target=descendant_node))
    all_weigths = [path_average_weight(graph, path) for path in paths]
    all_length = [len(path)-1 for path in paths]
    graph = select_best_path(graph, paths, all_length, all_weigths)
    return graph


def simplify_bubbles(graph: DiGraph) -> DiGraph:
    """Detect and explode bubbles

    :param graph: (nx.DiGraph) A directed graph object
    :return: (nx.DiGraph) A directed graph object
    """
    print("Detecting and exploding bubbles...")
    bubble_detected = False
    bubble_ancestor = None
    bubble_descendant = None
    cp = 0
    for node in graph:
        predecessors = list(graph.predecessors(node))
        if len(predecessors) > 1:
            for i, node_i in enumerate(predecessors):
                for node_j in predecessors[i + 1:]:
                    common_ancestor = lowest_common_ancestor(graph, node_i, node_j)
                    if common_ancestor is not None:
                        bubble_ancestor = common_ancestor
                        bubble_descendant = node
                        bubble_detected = True
                        break
                if bubble_detected:
                    break
    if bubble_detected:
        cp += 1
        graph = simplify_bubbles(solve_bubble(graph, bubble_ancestor, bubble_descendant))
    print(f"{cp} bubble(s) exploded sucessfully ! \n")
    return graph


def solve_entry_tips(graph: DiGraph, starting_nodes: List[str]) -> DiGraph:
    """Remove entry tips

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of starting nodes
    :return: (nx.DiGraph) A directed graph object
    """
    while True:
        entry_tips = [node for node in graph if len(list(graph.predecessors(node))) > 1]
        if not entry_tips:
            break
        for node in entry_tips:
            paths = [list(all_simple_paths(graph, node_start_i, node))
                     for node_start_i in starting_nodes]
            paths = [path[0] for path in paths if len(path) > 0]
            lengths = [len(path) - 1 for path in paths]
            weights = [path_average_weight(graph, path) if lengths[i] > 1
                       else graph.get_edge_data(*path)["weight"]
                        for i, path in enumerate(paths)]
            graph = select_best_path(graph, paths, lengths, weights,
                                        delete_entry_node=True, delete_sink_node=False)
    return graph


def solve_out_tips(graph: DiGraph, ending_nodes: List[str]) -> DiGraph:
    """Remove out tips

    :param graph: (nx.DiGraph) A directed graph object
    :param ending_nodes: (list) A list of ending nodes
    :return: (nx.DiGraph) A directed graph object
    """
    print("Removing out tips...")
    while True:
        found_tip = False
        for node in graph:
            node_success = list(graph.successors(node))
            if len(node_success) > 1:
                paths = [list(all_simple_paths(graph, node, node_end_i))\
                         for node_end_i in ending_nodes]
                paths = [path[0] for path in paths if len(path) > 0]
                lengths = [len(path) - 1 for path in paths]
                weights = [path_average_weight(graph, path) if lengths[i] > 1 else \
                           graph.get_edge_data(*path)["weight"]
                           for i, path in enumerate(paths)]
                graph = select_best_path(graph, paths, lengths, weights,
                                         delete_entry_node=False, delete_sink_node=True)
                found_tip = True
                break
        if not found_tip:
            break
    print("Out tips removed successfully! \n")
    return graph


def get_starting_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without predecessors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without predecessors
    """
    input_nodes = []
    for node in graph.nodes():
        predecessors = list(graph.predecessors(node))
        if not predecessors:
            input_nodes.append(node)
    return input_nodes


def get_sink_nodes(graph: DiGraph) -> List[str]:
    """Get nodes without successors

    :param graph: (nx.DiGraph) A directed graph object
    :return: (list) A list of all nodes without successors
    """
    output_nodes = []
    for node in graph.nodes():
        successors = list(graph.successors(node))
        if not successors:
            output_nodes.append(node)
    return output_nodes


def get_contigs(
    graph: DiGraph, starting_nodes: List[str], ending_nodes: List[str]
) -> List:
    """Extract the contigs from the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param starting_nodes: (list) A list of nodes without predecessors
    :param ending_nodes: (list) A list of nodes without successors
    :return: (list) List of [contiguous sequence and their length]
    """
    contigs = []
    for start_node in starting_nodes:
        for end_node in ending_nodes:
            if start_node in graph and end_node in graph:
                if has_path(graph, start_node, end_node):
                    for path in all_simple_paths(graph, start_node, end_node):
                        contig = path[0]  # Start with the first kmer
                        for node in path[1:]:
                            contig += node[-1]  # Add the last character of the kmer
                        contigs.append((contig, len(contig)))
                else:
                    print(f"No path between {start_node} and {end_node}")
            else:
                print(f"Node not found - Start: {start_node}, End: {end_node}")
    return contigs


def save_contigs(contigs_list: List[str], output_file: Path) -> None:
    """Write all contigs in fasta format

    :param contigs_list: (list) List of [contiguous sequence and their length]
    :param output_file: (Path) Path to the output file
    """
    with output_file.open("w") as f:
        for i, (contig, length) in enumerate(contigs_list):
            f.write(f">contig_{i} len={length}\n")
            f.write(textwrap.fill(contig, width=80)) #  Fasta format
            f.write("\n")


def draw_graph(graph: DiGraph, graphimg_file: Path) -> None:  # pragma: no cover
    """Draw the graph

    :param graph: (nx.DiGraph) A directed graph object
    :param graphimg_file: (Path) Path to the output file
    """
    elarge = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] > 3]
    esmall = [(u, v) for (u, v, d) in graph.edges(data=True) if d["weight"] <= 3]
    # Draw the graph with networkx
    pos = nx.random_layout(graph)
    nx.draw_networkx_nodes(graph, pos, node_size=6)
    nx.draw_networkx_edges(graph, pos, edgelist=elarge, width=6)
    nx.draw_networkx_edges(
        graph, pos, edgelist=esmall, width=6, alpha=0.5, edge_color="b", style="dashed"
    )
    # save image
    plt.savefig(graphimg_file.resolve())


# ==============================================================
# Main program
# ==============================================================
def main() -> None:  # pragma: no cover
    """
    Main program function
    """
    # Step 1: Parse arguments
    args = get_arguments()
    # Retrieve parameters
    fastq_file = args.fastq_file
    kmer_size = args.kmer_size
    output_file = Path(args.output_file)
    # Step 2: Read the FASTQ file and extract k-mers
    print("Reading the FASTQ file...")
    # Build the k-mer dictionary
    print(f"Building the k-mer dictionary with k={kmer_size}...")
    kmer_dict = build_kmer_dict(fastq_file, kmer_size)
    # Display a few k-mers to verify
    print("Displaying the first 10 k-mers from the dictionary:")
    for kmer, occurrence in list(kmer_dict.items())[:10]:
        print(f"{kmer}: {occurrence}")
    # Step 3: Build the graph from k-mers
    print("Building the De Bruijn graph from the k-mers...")
    graph = build_graph(kmer_dict)
    # Display a few edges to verify
    print("Displaying a few edges from the graph:")
    for u, v, d in list(graph.edges(data=True))[:10]:
        print(f"{u} -> {v} (weight: {d['weight']})")
    # Step 4: Simplify bubbles in the graph
    print("Simplifying bubbles in the graph...")
    graph = simplify_bubbles(graph)
    # Step 5: Extract starting and ending nodes
    print("Extracting starting and ending nodes...")
    starting_nodes = get_starting_nodes(graph)
    ending_nodes = get_sink_nodes(graph)
    # Step 6: Simplify entry and exit tips in the graph
    print("Removing unwanted entry tips...")
    graph = solve_entry_tips(graph, starting_nodes)
    print("Removing unwanted exit tips...")
    graph = solve_out_tips(graph, ending_nodes)
    # Step 7: Extract contigs from starting and ending nodes
    print("Extracting contigs...")
    contigs_list = get_contigs(graph, starting_nodes, ending_nodes)
    # Display a few contigs to verify
    print("Displaying the first 5 extracted contigs:")
    for contig, length in contigs_list[:5]:
        print(f"Contig of length {length}: {contig[:50]}...")  # Displays first 50 characters
    # Step 8: Save contigs to a FASTA file
    print(f"Saving contigs to the file: {output_file}")
    save_contigs(contigs_list, output_file)
    print("Process completed successfully.")
    # Step 9: Save the graph as an image
    if args.graphimg_file:
        print(f"Saving the graph as an image: {args.graphimg_file}")
        draw_graph(graph, args.graphimg_file)


if __name__ == "__main__":  # pragma: no cover
    main()
