## Using SpiecEasi data generation, but with underlying graph being created as of the DINGO supplementary materials.

#### Create scale-free graph set for use in mDINGO.

### Remove and add nodes to get the control graph. Will probably no longer be scale-free, that is okay.
### There is also a bunch of weird changing of values, and upper.tris. This is to make sure that the network
### is symmetric and to get the change network to appear correctly. Mostly bookkeeping. Hopefully I never need
### to change it.

gen_graphs_sf = function(numNodes, edgePower = 1, hubgraphchangeprop = .6, nonhubgraphchangeprop = .05, hubsize = 1/3)
{
  graph_dis = sample_pa(n=numNodes, power=edgePower, m=1,  directed=F)
  hub = which.max(degree(graph_dis))
  detached = rep(NA, degree(graph_dis, v = hub))
  cutlist = rep(NA, degree(graph_dis, v = hub)*2)
  cutcount = 1
  while(sum(distances(graph_dis, v=hub) == Inf) < hubsize*length(V(graph_dis)))
  {
    ne = neighborhood(graph_dis, 1, nodes = hub)[[1]][-1]
    for(nd in 1:length(ne))
    {
      tempgraph = delete_edges(graph_dis, paste0(hub, "|", ne[nd]))
      detached[nd] = sum(distances(tempgraph, v=hub) == Inf)
    } 
    graph_dis = delete_edges(graph_dis, paste0(hub, "|", ne[which.max(detached)]))
    cutlist[cutcount:(cutcount+1)] = c(hub, ne[which.max(detached)])
    cutcount = cutcount+2
  }

  cutlist = cutlist[!is.na(cutlist)]

  reachfromhub = which(!(distances(graph_dis, v=hub) == Inf))
  nonreachfromhub = which(distances(graph_dis, v=hub) == Inf)
 
  graph_dis = add_edges(graph_dis, cutlist)
  
  hubgraph = induced_subgraph(graph_dis, reachfromhub)
  nonhubgraph = induced_subgraph(graph_dis, nonreachfromhub)
  
  hubgraph = rewire(hubgraph, each_edge(prob = hubgraphchangeprop))
  nonhubgraph = rewire(nonhubgraph, each_edge(prob = nonhubgraphchangeprop))
  
  graph_con = matrix(NA, numNodes, numNodes)
  graph_con[reachfromhub,][,reachfromhub] = as.matrix(hubgraph[])
  graph_con[nonreachfromhub,][,nonreachfromhub] = as.matrix(nonhubgraph[])
  graph_con[is.na(graph_con)] = 0
  graph_con = graph_from_adjacency_matrix(adjmatrix = graph_con, weighted = NULL, mode = "undirected")
  graph_con = add_edges(graph_con, cutlist)
  
  #V(graph_dis)$name = V(graph_dis)
  change_graph1 = (graph_dis - graph_con)
  change_graph2 = (graph_con - graph_dis)
  
  change_graph1 <- set_edge_attr(change_graph1, name = "color", value = rep("red", length(E(change_graph1))))
  change_graph2 <- set_edge_attr(change_graph2, name = "color", value = rep("green", length(E(change_graph2))))
  
  change_graph = change_graph1 %u% change_graph2
  if(length(E(change_graph)) > 0)
  {
    col = data.frame("color_1" = get.edge.attribute(change_graph)$color_1, "color_2" = get.edge.attribute(change_graph)$color_2)
    levels(col$color_1) = c("red", "green")
    col$color_1[is.na(col$color_1)] = "green"
    change_graph <- set_edge_attr(change_graph, name = "color", value = as.character(col$color_1))
  }
  
  global_graph = intersection(graph_dis, graph_con)
  
  retList = list(disease = graph_dis, change = change_graph, control = graph_con, global = global_graph)
  return(retList)
}

plotGraphList = function(graphList)
{
  l_fr = layout_with_fr(graphList$disease)
  plot(graphList$disease, vertex.size= 0, main = "Disease Graph", layout = l_fr, vertex.label = NA)
  plot(graphList$change, vertex.size= 0, main = "Change Graph", layout = l_fr, vertex.label = NA)
  plot(graphList$control, vertex.size= 0, main = "Control Graph", layout = l_fr, vertex.label = NA)
  plot(graphList$global, vertex.size= 0, main = "Global Graph", layout = l_fr, vertex.label = NA)
}