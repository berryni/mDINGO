\name{gen_graphs_sf}
\alias{gen_graphs_sf}

\title{
Generate Covariate Dependent Scale-Free Graphs
}
\description{
Use igraph to generate synthetic underlying scale-free graph structure for simulated OTU data. A list of graphs is outputted, the disease, change, control, and global graphs. See ``values" for full info. The disease graph is a scale free network generated with Barab{\'a}si-Albert algorithm. A hub node is then selected and expanded until it is at least as large as `hubsize` of the full graph. This hub is then permuted to form the control graph.
}

\arguments{
  \item{numNodes}{
    The number of nodes to include in created graphs.
}
  \item{edgePower}{
    Default value of 1. The power of the preferential attachment. Determines strength of scale-free structure.
}
  \item{hubgraphchangeprop}{
    Defaults to .6. Probability of altering each edge in the graph subset that is determined to be the hub (disease generating phenomena).
}
\item{nonhubgraphchangeprop}{
    Defaults to .05. Probability of altering each edge not in the graph subset that is determined to be the hub.
}
\item{hubsize}{
    Defaults to 1/3. Must be between 0 and 1. Dictates size of disease generating hub node. 1/3 indicates that at least 1/3 of the entire graph must be detached from the disease generation node. Larger number -> smaller hub node.
}
}
\value{
    List of 4 igraph objects.
	\item{disease}{
		Dependence graph for disease having patient subset.
	}
	\item{change}{
		Graph outlining what changed between disease and control graphs.
	}
	\item{control}{
		Dependence graph for control patients. Result of applying change graph to disease graph.
	}
	\item{global}{
		Intersection of disease and control graph edges.
	}
}
\author{
    Nick Berry berryni@iastate.edu
}
