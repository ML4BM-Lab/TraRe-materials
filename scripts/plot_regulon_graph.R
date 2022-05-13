# SCRIPT TO PLOT SUB-FIGURES 3D (REGULON GRAPH)

library(visNetwork)
library(igraph)

# Load data
regulons <- readRDS('./output/sm_3500_50b_regulons.rds')

plot_regulon_graph <-function(data,
                              regulon_name = NULL,
                              filter = FALSE,
                              weighted = TRUE,
                              tf_in_color = "orange",
                              tf_out_color = "darkblue",
                              tf_shape = "box",
                              target_in_color = "lightblue",
                              target_shape = "ellipse",
                              legend = T,
                              nodeDistance = 130 ,
                              springLength = 50,
                              springConstant = 0.05) {
  if (!regulon_name %in% names(data$graph)) {
    stop("Error: Regulon does not exist")
  }
  
  g <- data$graph[[regulon_name]]
  
  if (filter > max(edge.attributes(g)$weight)) {
    stop("Error: filter > max weight")
  }
  
  if (filter) {
    g <-igraph::subgraph.edges(graph = g, eids = E(g)[edge.attributes(g)$weight >=filter])
    weight_dim <- dim(table(edge.attributes(g)$weight))
    myfactor <- seq(from = 1,to = 2 * weight_dim,by = 2)
  } else {
    myfactor <- unique(edge.attributes(g)$weight)
  }
  
  g_vis <- visNetwork::toVisNetworkData(g)
  
  if (weighted) {
    g_vis$edges$width <- factor(E(g)$weight, labels = myfactor)
    message("Weights: ", paste(names(table(edge.attributes(g)$weight)), collapse = ", "))
  } else {
    g_vis$edges$width <- 1
  }
  
  g_vis$nodes$group <- c("TF", rep("target", nrow(g_vis$nodes) - 1))
  message("Number of targets: ", sum(g_vis$nodes$group == "target"))
  
  p <- visNetwork::visNetwork(g_vis$nodes, g_vis$edges) %>%
    visGroups(groupname = "TF",
              color = list(background = tf_in_color, border = tf_out_color),
              shape = tf_shape,
              shadow = T) %>%
    visGroups( groupname = "target",
               color = list(background = target_in_color),
               shape = target_shape) %>%
    visLegend(enabled = legend) %>% 
    visPhysics(solver = "repulsion",
               repulsion = list(springConstant = springConstant,
                                nodeDistance = nodeDistance,
                                springLength = springLength))
  
  return(p)
}

# PLOT the interactive graphs

plot_regulon_graph(regulons, 
                   regulon_name = "ELK3",
                   filter = 7,
                   weighted = T,
                   springConstant = 0.10, 
                   nodeDistance = 200,
                   springLength = 25)
plot_regulon_graph(regulons, 
                   regulon_name = "MYB",
                   filter = F,
                   weighted = T,
                   springConstant = 0.10, 
                   nodeDistance = 200,
                   springLength = 25)
plot_regulon_graph(regulons, 
                   regulon_name = "MXD1",
                   filter = F,
                   weighted = T,
                   springConstant = 0.10, 
                   nodeDistance = 200,
                   springLength = 25)

# For manuscript sub-figures 3D Adjust the vertices in teh interactive graph to
# avoid overlap and export the image as png file or screenshot.
################################ EOF #############################################