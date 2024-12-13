## -----------------------------------------------------------------------------
##
## [ PROJ ] Replication and Extension for PLSC 508 Final Project
## [ FILE ] data_analysis.R
## [ AUTH ] Maya Dalton; mad6821@psu.edu
## [ INIT ] 22 October 2024
##
## -----------------------------------------------------------------------------

## libraries
libs <- c("dplyr", "tidyverse", "RColorBrewer", "network", "igraph", "sna", "ggraph", 
          "blockmodeling", "RSiena", "texreg", "intergraph", "PAFit", "ergm", "ContagionTest")
sapply(libs, require, character.only = TRUE)

## paths - set wd as scripts folder
args <- commandArgs(trailingOnly = TRUE)
root <- ifelse(length(args) == 0, file.path(".."), args)
dat_dir <- file.path(root, "data")
fig_dir <- file.path(root, "figures")
scr_dir <- file.path(root, "scripts")
tab_dir <- file.path(root, "tables")

color <- brewer.pal(4, "RdBu") # choosing a specific palette

## -----------------------------------------------------------------------------
## Analysis data frame
## -----------------------------------------------------------------------------

## Load in dataset ------------------------------
drugnet <- read_csv(file.path(dat_dir, "drugnet.csv"))

## Initializing network objects ------------------------------
edge_df <- drugnet %>% # Create edge list
  select(sideA, sideB, weight) 

combined_graph <- graph_from_data_frame(edge_df) # turn into igraph

combined_matrix <- as.matrix(as_adjacency_matrix(combined_graph)) # adjacency matrix

main_graph <- readRDS(file.path(dat_dir, "cartel_network.rds"), refhook = NULL) # load in full network object

## Split into two graphs for pre and post arrest ------------------------------
preNet <- subgraph.edges(main_graph, E(main_graph)[E(main_graph)$period < 2])

postNet <- subgraph.edges(main_graph, E(main_graph)[E(main_graph)$period > 1])

## -----------------------------------------------------------------------------
## Plot of the Networks (for 10/22 Update)
## -----------------------------------------------------------------------------

## Before El Chapo arrest ------------------------------
png(file.path(fig_dir, "pre_arrest_net.png"), res=100)
set.seed(5)
xy_cool <- layout_with_fr(preNet)
xy <- plot(preNet, # network
           vertex.label="", # don't display the labels
           coord = xy_cool,  # set the coordinates
           vertex.color=color,
           edge.arrow.size = 0.5) # color the nodes
dev.off()

## After El Chapo arrest ------------------------------
png(file.path(fig_dir, "post_arrest_net.png"), res=100)
set.seed(5)
xy_cool <- layout_with_fr(postNet)
xy <- plot(postNet, # network
           vertex.label="", # don't display the labels
           coord = xy_cool,  # set the coordinates
           vertex.color=color,
           edge.arrow.size = 0.5) # color the nodes
dev.off()

## -----------------------------------------------------------------------------
## Plot of the Vertex Attributes (for 10/22 Update)
## -----------------------------------------------------------------------------

## Aggression Attribute ----------------------------------------
png(file.path(fig_dir, "aggression.png"), res=100)
hist(V(main_graph)$aggression, col="palegreen3",
     main="",
     xlab = "Aggression Attribute")
dev.off()

## Militia Attribute ----------------------------------------
png(file.path(fig_dir, "militia.png"), res=100)
hist(V(main_graph)$militia, col="palegreen3",
     main="",
     xlab = "Militia Status")
dev.off()

## Subfaction Attribute ----------------------------------------
png(file.path(fig_dir, "subfaction.png"), res=100)
hist(V(main_graph)$subfaction, col="palegreen3",
     main="",
     xlab = "Subfaction Status")
dev.off()

## Role Attribute ----------------------------------------
roles <- V(main_graph)$role
counts <- table(roles)

png(file.path(fig_dir, "roles.png"), res=100)
barplot(counts, col="palegreen3", width=0.5, cex.names = 0.6,
        names.arg = c("Jalisco NG", "Gulf Cartel", "Los Zetas", "Rising Challengers", "Sinaloa Cartel", 
                      "Small Cartels", "White Dwarfs"),
        main="")
dev.off()

## Time Period Attribute ----------------------------------------
png(file.path(fig_dir, "time.png"), res=100)
hist(E(main_graph)$period, col="palegreen3",
     main="",
     xlab = "Time Period (1=Pre-Arrest, 2=Post-Arrest)")
dev.off()

## -----------------------------------------------------------------------------
## Reproduction of Results (for 10/22 Update)
## -----------------------------------------------------------------------------
## -------------------------------------
## In- and Out-Degree Plots
## -------------------------------------

## Pre-Arrest In Degree ----------------------------------------
pre_in_degree <- igraph::degree(preNet, mode="in")
png(file.path(fig_dir, "pre_in_degree.png"), res=100)
hist(pre_in_degree, breaks = 30, 
     main="",
     xlab = "In-Degree", col = "lightblue")
dev.off()

## Post-Arrest In Degree ----------------------------------------
post_in_degree <- igraph::degree(postNet, mode="in")
png(file.path(fig_dir, "post_in_degree.png"), res=100)
hist(post_in_degree, breaks = 30, 
     main="",
     xlab = "In-Degree", col = "lightblue")
dev.off()

## Pre-Arrest Out Degree ----------------------------------------
pre_out_degree <- igraph::degree(preNet, mode="out")
png(file.path(fig_dir, "pre_out_degree.png"), res=100)
hist(pre_out_degree, breaks = 30, 
     main="",
     xlab = "Out-Degree", col = "lightblue")
dev.off()

## Post-Arrest Out Degree ----------------------------------------
post_out_degree <- igraph::degree(postNet, mode="out")
png(file.path(fig_dir, "post_out_degree.png"), res=100)
hist(post_out_degree, breaks = 30, 
     main="",
     xlab = "Out-Degree", col = "lightblue")
dev.off()

## -------------------------------------
## SAOM Estimation
## -------------------------------------
## Estimate SAOM Function ----------------------------------------

# Function to estimate a stochastic actor oriented model
# @Params: effect, the effect to estimate; seed, the seed to set the machine to
# @Returns: A sienaFit object
estimate_saom <- function(effect, seed){
  
  effect <- deparse(substitute(effect))
  
  # Create the basic data
  data <- sienaDataCreate(dv, aggression_var, subfaction_var, militia_var, 
                          role_var, joiners)
  
  # Effects for rate, reciprocity, and out-degree
  effects <- getEffects(data)
  
  # Effects for in-degree popularity, in-degree activity, out-degree popularity,
  # and balance. These reflect preferential attachment and similarity for outgoing 
  # ties between actors
  if (effect == "Jout") effects <- includeEffects(effects, Jout)
  if (effect == "inPop") effects <- includeEffects(effects, inPop)
  if (effect == "outInAss") effects <- includeEffects(effects, outInAss)
  if (effect == "transTrip1") effects <- includeEffects(effects, transTrip1)
  
  # Effects for homophily
  effects <- includeEffects(effects, sameX, interaction1 = "aggression_var")
  effects <- includeEffects(effects, sameX, interaction1 = "subfaction_var")
  effects <- includeEffects(effects, sameX, interaction1 = "militia_var")
  effects <- includeEffects(effects, sameX, interaction1 = "role_var")
  
  # Define algorithm and get results
  algorithm <- sienaAlgorithmCreate(projname = "output", seed = seed, 
                                    n3 = 3000, nsub = 7, firstg = 0.01)
  results <- siena07(algorithm, data = data, effects = effects, 
                     returnDeps = TRUE)
  
  # Check for convergence
  if (max(results$tstat > 1) | abs(results$tconv.max) > 0.25){
    
    results <- siena07(algorithm, data = data, effects = effects, 
                       returnDeps = TRUE, prevAns = results)
    
  }
  
  return(results)
  
}

## Prepare the data for SAOM estimation ------------------------------------

# Edges in the first period
t0_network <- main_graph %>% # Edges in the first period
  delete_edges(as.vector(E(postNet)))

# Edges in the second period
t1_network <- subgraph.edges(main_graph, as.vector(E(postNet)), 
                             delete.vertices = FALSE)

# Make a list for joiner nodes
joiner_nodes <- as_tibble(as.vector(V(t0_network)$name)) %>% 
  mutate(node = row_number()) %>% 
  ungroup() %>% 
  relocate(node) %>% 
  filter(!(value %in% as.vector(V(preNet)$name))) %>% 
  pull(node)

# Make a list for all nodes and code them as being present in both periods
joiners <- rep(list(c(1, 2)), 151)

# Update periods for nodes that were only present in the second period
for (i in seq_along(joiners)){
  if (i %in% joiner_nodes){ # Check to see if a node was a joiner in the second period
    joiners[[i]][1] <- 2     # Recode the list to show joiners only in period 2
  }
}

## Creates an RSiena object for joiners
joiners <- sienaCompositionChange(joiners)

## Create a dependent variable
dv <- sienaDependent(array(c(as.matrix(as_adjacency_matrix(t0_network)), 
                             as.matrix(as_adjacency_matrix(main_graph))), 
                           dim = c(151, 151, 2)))

## Create constant covariates objects
blocks <- readRDS(file.path(dat_dir, "blocks.rds"))

role_var <- coCovar(as.vector(blocks$best$best1$clu))
militia_var <- coCovar(as.vector(V(main_graph)$militia))
aggression_var <- coCovar(as.vector(V(main_graph)$aggression))
subfaction_var <- coCovar(as.vector(V(main_graph)$subfaction))


## Run models with different effects ---------------------------------------
# Effects for Jaccard out-degree, in degree popularity, out-in-degree 
# assortativity, and transitive triples type 1
m1_results <- estimate_saom(Jout, 100)
m2_results <- estimate_saom(inPop, 200)
m3_results <- estimate_saom(outInAss, 300)
m4_results <- estimate_saom(transTrip1, 400)

## Create a regression table -----------------------------------------------
texreg(l = list(m1_results, m2_results, m3_results, m4_results),
          omit.coef = "basic rate parameter dv",
          custom.coef.map = list("out-Jaccard similarity" = "Jaccard Similarity",
                                 "indegree - popularity" = "In-degree Popularity",
                                 "out-in degree^(1/2) assortativity" 
                                 = "Out-in-degree Assortativity",
                                 "transitive triplets (1)" = "Transitive Closure",
                                 "same aggression_var" = "Aggression Homophily",
                                 "same subfaction_var" = "Subfaction Homophily",
                                 "same militia_var" = "Militia Homophily",
                                 "same role_var" = "Role Homophily",
                                 "reciprocity" = "Reciprocity"),
          custom.model.names = c("Alliances", "Reputation", "Strong-vs-weak",
                                 "Clustering"),
          main = "Table 1. SAOM Model Estimates - Replicated from Colby (2021)",
          file = file.path(tab_dir, "saom.tex"))


## GOF Diagnostics --------------------------------------------------
png(file.path(fig_dir, "gof_m1.png"), res=100)
plot(btergm::gof(m1_results), main=NULL)
dev.off()

png(file.path(fig_dir, "gof_m2.png"), res=100)
plot(btergm::gof(m2_results), main=NULL)
dev.off()

png(file.path(fig_dir, "gof_m3.png"), res=100)
plot(btergm::gof(m3_results), main=NULL)
dev.off()

png(file.path(fig_dir, "gof_m4.png"), res=100)
plot(btergm::gof(m4_results), main=NULL)
dev.off()

## -----------------------------------------------------------------------------
## Descriptive Extension (for 11/12 Update)
## -----------------------------------------------------------------------------

## -------------------------------------
## Reciprocity and Transitivity
## -------------------------------------

## Reciprocity -----------------------------------------------
pre_recip <- igraph::reciprocity(preNet) # 0.0
post_recip <- igraph::reciprocity(postNet) # 0.243
main_recip <- igraph::reciprocity(main_graph) # 0.288

## Transitivity -----------------------------------------------
pre_trans <- igraph::transitivity(preNet) # 0.0
post_trans <- igraph::transitivity(postNet) # 0.039
main_trans <- igraph::transitivity(main_graph) # 0.033

# bind all scores into a matrix
scores <- rbind(pre_recip, post_recip, main_recip, 
                pre_trans, post_trans, main_trans)
scores <- as.data.frame(scores)

## Plot of Scores -----------------------------------------------
png(file.path(fig_dir, "recip_trans_scores.png"), res=100)
ggplot(scores) + 
  geom_bar(aes(x=rownames(scores), y=V1), stat='identity', fill="plum4") +
  labs(x="", y="Reciprocity and Transitivity Scores") + 
  theme_minimal()
dev.off()

# rename columns
rownames(scores) <- c("Pre-Arrest Recip", "Post-Arrest Recip", "Full Graph Recip", 
                      "Pre-Arrest Trans", "Post-Arrest Trans", "Full Graph Trans")

## Create table of reciprocity and transitivity scores -----------------------------------------------
recip_trans_tab <- knitr::kable(scores, booktabs = TRUE, digits=3, format = "latex", 
                                table.envir = "table", position = "!h",
                                col.names = NULL, 
                                caption="Reciprocity and Transitivity Scores") %>%
  kableExtra::pack_rows(index = c("Reciprocity" = 3, "Transitivity" = 3))

writeLines(recip_trans_tab, file.path(tab_dir, 'recip_trans_tab.tex')) ## save

## -------------------------------------
## Preferential Attachment
## -------------------------------------
drug_EL <- as_data_frame(main_graph, what = "edges")
drug_EL <- drug_EL %>% 
  mutate(from = as.numeric(as.factor(from)),  # converting side names to unique IDs
         to = as.numeric(as.factor(to)))

drug_PAFit <- as.PAFit_net(as.matrix(drug_EL))

## Plotting degree distribution -----------------------------------------------

# Distribution of degree centrality at t=0 (pre-2017)
png(file.path(fig_dir, "pre_PA_deg.png"), res=100)
plot(drug_PAFit, slice = 0, plot = "degree", cex = 2, cex.axis = 1, cex.lab = 1) 
dev.off()

# Distribution of degree centrality at t=9 (2020)
png(file.path(fig_dir, "post_PA_deg.png"), res=100)
plot(drug_PAFit, slice = 1, plot = "degree", cex = 2, cex.axis = 1, cex.lab = 1) 
dev.off()

## Plotting Preferrential Attachment -----------------------------------------------
# Network Plot at t=0 (pre-2017)
png(file.path(fig_dir, "pre_PA.png"), res=100)
plot(drug_PAFit, slice = 0, 
     arrowhead.cex = 2, vertex.cex = 3,
     vertex.col=color) 
dev.off()

# Network Plot at t=1 (post-2017)
png(file.path(fig_dir, "post_PA.png"), res=100)
plot(drug_PAFit, slice = 1, 
     arrowhead.cex = 2, vertex.cex = 3,
     vertex.col=color) 
dev.off()

## Preferrential Attachment Table -----------------------------------------------
result_OS <- PAFit_oneshot(drug_PAFit)
summary(result_OS)

PA_tab <- knitr::kable(result_OS$alpha, booktabs = TRUE, digits=3, format = "latex", 
                         table.envir = "table", position = "!h",
                         col.names = NULL, 
                         caption="Preferential Attachment") 

writeLines(PA_tab, file.path(tab_dir, 'PA_tab.tex')) ## save

## -------------------------------------
## Community Detection
## -------------------------------------
## Comparing Modularities -----------------------------------------------
compare_modularity <- function(graph) {
    # Compute each clustering method
    infomap <- cluster_infomap(graph)            # Infomap
    lead_eig <- cluster_leading_eigen(graph)     # Leading Eigenvector
    label_prop <- cluster_label_prop(graph)       # Lable Prop
    walktrap <- cluster_walktrap(graph)          # Walktrap
    
    # Store the modularity for each method
    modularities <- c(
      modularity(infomap),
      modularity(lead_eig),
      modularity(label_prop),
      modularity(walktrap)
    )
    
  # Convert list to a data frame for easy comparison
  mod_df <- data.frame(modularities)
  rownames(mod_df) <- c("Infomap", "Leading Eigenvector", "Spinglass", "Walktrap")
  
  return(mod_df)
}

drugMods <- compare_modularity(main_graph) # Infomap performs the best!

## Create table of modularity scores -----------------------------------------------
mods_tab <- knitr::kable(drugMods, booktabs = TRUE, digits=3, format = "latex", 
                                table.envir = "table", position = "!h",
                                col.names = NULL, 
                                caption="Modularity Scores") 

writeLines(mods_tab, file.path(tab_dir, 'mods_tab.tex')) ## save

## Plotting Infomap for Pre-Arrest -----------------------------------------------
infomap <- cluster_infomap(preNet) #Perform info map community detection
membership(infomap)  # Shows the community each node belongs to
sizes(infomap)       # Shows the sizes of the detected communities

set.seed(5) #set seed for replication
png(file.path(fig_dir, "pre_community.png"), res=100)
plot(infomap, preNet, #visualize the communities identified by the infomap method
     layout = layout_with_fr, #set layout using FR algorithm
     vertex.label = "", #remove labels
     vertex.size = 5, #set size of nodes
     edge.arrow.size = 0.2 #set size of arrows (if directed network)
) 
dev.off()

## Plotting Infomap for Post-Arrest -----------------------------------------------
infomap <- cluster_infomap(postNet) #Perform info map community detection
membership(infomap)  # Shows the community each node belongs to
sizes(infomap)       # Shows the sizes of the detected communities

set.seed(5) #set seed for replication
png(file.path(fig_dir, "post_community.png"), res=100)
plot(infomap, postNet, #visualize the communities identified by the infomap method
     layout = layout_with_fr, #set layout using FR algorithm
     vertex.label = "", #remove labels
     vertex.size = 5, #set size of nodes
     edge.arrow.size = 0.2 #set size of arrows (if directed network)
) 
dev.off()

## -----------------------------------------------------------------------------
## Inferential Extension (for 12/17 Submission)
## -----------------------------------------------------------------------------
## -------------------------------------
## ERGM Specification
## -------------------------------------

## Switching to network objects ------------------------------

# Convert igraph to adjacency matrix
adj_matrix <- as_adjacency_matrix(main_graph, sparse = FALSE)

# Create a network object
drugNet <- network(adj_matrix, directed = is_directed(main_graph))

# Transfer vertex attributes
network::set.vertex.attribute(drugNet, "name", V(main_graph)$name)
network::set.vertex.attribute(drugNet, "aggression", V(main_graph)$aggression)
network::set.vertex.attribute(drugNet, "role", V(main_graph)$role)
network::set.vertex.attribute(drugNet, "militia", V(main_graph)$militia)
network::set.vertex.attribute(drugNet, "subfaction", V(main_graph)$subfaction)

# Transfer edge attributes
edge_list <- as_edgelist(main_graph)
edge_period <- E(main_graph)$period
network::set.edge.attribute(drugNet, "period", edge_period)

## Standardizing aggression attributes ------------------------------
# Aggression
aggression <- network::get.vertex.attribute(drugNet, "aggression")
scaled_aggression <- as.vector(scale(aggression))

# Assign the standardized attribute back to the network
network::set.vertex.attribute(drugNet, "aggression", scaled_aggression)


## Baseline ERGM ------------------------------
control <- control.ergm(MCMLE.maxit=50, seed=47)
set.seed(123)
ergm_mod <- ergm(drugNet ~ edges + 
                   nodeicov("aggression") + 
                   nodeocov("aggression") +
                   absdiff("aggression") + 
                   nodeicov("subfaction") + 
                   nodeocov("subfaction") +
                   nodematch("subfaction") + 
                   nodeicov("militia") + 
                   nodeocov("militia") +
                   nodematch("militia") + 
                   nodematch("role"),
                 control=control)

summary(ergm_mod)

## node fixed effects term - sender() and receiver()
## nodeifactor and nodeofactor 
## delete bottom rows from textreg 

## Dyad-Dependent ERGM ------------------------------
control <- control.ergm(MCMLE.maxit=50, seed=47)
set.seed(123)
ergm_mod_DD <- ergm::ergm(drugNet ~ edges + 
                            nodeicov("aggression") + 
                            nodeocov("aggression") +
                            absdiff("aggression") + 
                            nodeicov("subfaction") + 
                            nodeocov("subfaction") +
                            nodematch("subfaction") + 
                            nodeicov("militia") + 
                            nodeocov("militia") +
                            nodematch("militia") + 
                            nodematch("role") + 
                            mutual + 
                            gwesp(0.5, fixed=T), 
                          control = control)

screenreg(ergm_mod_DD)


## Create a regression table -----------------------------------------------
texreg(list(ergm_mod, ergm_mod_DD),
       omit.coef = "basic rate parameter dv",
       custom.coef.map = list("edges" = "Edges",
                              "nodeicov.aggression" = "Receiver Aggression",
                              "nodeocov.aggression" = "Sender Aggression",
                              "absdiff.aggression" = "Aggression Homophily",
                              "nodeicov.subfaction" = "Receiver Subfaction",
                              "nodeocov.subfaction" = "Sender Subfaction",
                              "nodematch.subfaction" = "Subfaction Homophily",
                              "nodeicov.militia" = "Receiver Militia",
                              "nodeocov.militia" = "Sender Militia",
                              "nodematch.militia" = "Militia Homophily",
                              "nodematch.role" = "Role Homophily",
                              "mutual" = "Reciprocity",
                              "gwesp.OTP.fixed.0.5" = "Transitivity (gwesp)"),
       label = "tab:ergms",
       file = file.path(tab_dir, "ergms.tex"))


## GOF Diagnostics --------------------------------------------------
png(file.path(fig_dir, "mcmc_dd.png"))
mcmc.diagnostics(ergm_mod_DD, which="plots")
dev.off()

png(file.path(fig_dir, "gof_ergm.png"))
plot(btergm::gof(ergm_mod), main=NULL)
dev.off()

png(file.path(fig_dir, "gof_dd.png"))
plot(btergm::gof(ergm_mod_DD))
dev.off()

## -----------------------------------------------------------------------------
################################################################################









