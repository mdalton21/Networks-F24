## -----------------------------------------------------------------------------
##
## [ PROJ ] Replication and Extension for PLSC 508 Final Project
## [ FILE ] network_init.R
## [ AUTH ] Maya Dalton; mad6821@psu.edu
## [ INIT ] 20 November 2024
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

main_graph <- edge_df %>% # distinct observations for igraph
  distinct(sideA, sideB) %>% 
  graph_from_data_frame()

## -------------------------------------
## Add structural equivalence group attribute
## -------------------------------------

## Create a 7 block solution ------------------------------
blocks <- optRandomParC(M = combined_matrix, 
                        k = 7,  
                        rep = 100,  
                        approaches = "ss", 
                        blocks = "com",  
                        seed = 42)

## Save blocks for easy loading ------------------------------
saveRDS(blocks, file.path(dat_dir, "blocks.rds"))

## Assign actors to blocks
V(main_graph)$role <- blocks$best$best1$clu

# Convert roles to names ------------------------------
V(main_graph)$role <- case_when(
  V(main_graph)$role == 1 ~ "Gulf Cartel",
  V(main_graph)$role == 2 ~ "Rising Challengers",
  V(main_graph)$role == 3 ~ "Sinaloa Cartel",
  V(main_graph)$role == 4 ~ "Small Cartels and Militias",
  V(main_graph)$role == 5 ~ "Los Zetas",
  V(main_graph)$role == 6 ~ "White Dwarfs",
  V(main_graph)$role == 7 ~ "Cartel Jalisco Nueva Generacion",
  TRUE ~ as.character(V(main_graph)$role)
)

## -------------------------------------
## Add attribute for militia status
## -------------------------------------

## List of militia groups from replication code
militias <- c("malinaltepec communal militia (mexico)", "chamula militia", 
              "san pablo cuatro venados communal militia (mexico)", 
              "dzan communal militia (mexico)", 
              "el pinar communal militia (mexico)", 
              "loma de la cruz communal militia (mexico)", 
              "molinos los arcos communal militia (mexico)", 
              "el carrizal de bravo communal militia (mexico)", 
              "huahua communal militia (mexico)", 
              "tlacotepec communal militia (mexico)", 
              "maya indigenous militia (mexico)", 
              "mixe indigenous militia (mexico)", 
              "san agustin oapan communal militia (mexico)", 
              "santa isabel de la reforma communal militia (mexico)", 
              "aldama indigenous militia (mexico)", 
              "san mateo yucutindoo communal militia (mexico)", 
              "apetlanca communal militia (mexico)", 
              "rio santiago communal militia (mexico)",
              "yautepec communal militia (mexico)", 
              "upoeg self defense group", 
              "faction of front for security and development vigilante group", 
              "totolapan self-defense force", "autodefensas unidas de michoacan", 
              "alacatlatzala communal militia (mexico)", 
              "chocomanatlan communal militia (mexico)", 
              "coahuayutla de guerrero communal militia (mexico)", 
              "cotija de la paz communal militia (mexico)", 
              "cuilapam de guerrero communal militia (mexico)", 
              "el cipresal communal militia (mexico)", 
              "el nith communal militia (mexico)", 
              "ezln: zapatista army of national liberation", 
              "leonardo bravo communal militia (mexico)", 
              "los reyes communal militia (mexico)", 
              "san cristobal de las casas communal militia (mexico)", 
              "san miguel tecuiciapan communal militia (mexico)", 
              "santa maria nativitas coatlan communal militia (mexico)", 
              "santa martha indigenous militia (mexico)", 
              "santiago amoltepec communal militia (mexico)", 
              "santiago miltepec communal militia (mexico)", 
              "santo reyes zochiquilazala communal militia (mexico)", 
              "tangancicuaro communal militia (mexico)", 
              "tenquizolco communal militia (mexico)", 
              "tepalcatepec communal militia (mexico)", 
              "tinguindin communal militia (mexico)", 
              "tocumbo communal militia (mexico)", 
              "villa morelos communal militia (mexico)", 
              "villa victoria communal militia (mexico)", 
              "xochiltepec communal militia (mexico)",
              "rival faction of front for security and development vigilante group")

## Assign militia status to nodes ------------------------------
V(main_graph)$militia <- ifelse(V(combined_graph)$name %in% militias, 1, 0)

## -------------------------------------
## Add aggressiveness attribute 
## -------------------------------------

## Convert the graph to a dataframe ------------------------------
temp_vector <- as_data_frame(main_graph) %>% 
  as_tibble() %>% 
  select(from) %>% 
  
  ## Merge the vertices with the the aggression data for sending actors
  left_join(drugnet %>% 
              select(sideA, agression) %>% 
              distinct(), 
            by = c("from" = "sideA")) %>% 
  distinct(from, agression) %>% 
  
  ## Add the sideB actors to the dataframe. 
  full_join(drugnet %>% 
              select(sideB, agression) %>% 
              distinct(),
            by = c("from" = "sideB")) %>% 
  select(-agression.y) %>% 
  rename("agression" = agression.x) %>% 
  distinct(from, agression) %>% 
  select(agression) %>% 
  
  ## Codes actors that never directed an attack as having an aggression of 0
  mutate(agression = ifelse(is.na(agression), 0, agression)) %>% 
  ungroup() %>% 
  pull()

## Assign aggression scores to vertices ------------------------------
V(main_graph)$aggression <- temp_vector

## -------------------------------------
## Add subfaction attribute 
## -------------------------------------

subfaction <- vector(length = 151, mode = "numeric")

## The indices of groups that are subfactions ------------------------------
subfaction[c(8, 72, 74, 75, 77:79, 113, 134, 136:138, 142:146, 148, 149)] <- 1
V(main_graph)$subfaction <- subfaction

## -------------------------------------
## Create time period before and after arrest
## -------------------------------------

## Create a vector of periods in which each tie existed ------------------------------
period <- drugnet %>% 
  mutate(period = ifelse(year < 2017, 1, 2)) %>% 
  group_by(period, sideA, sideB) %>% 
  ungroup() %>% 
  distinct(sideA, sideB, period) %>% 
  group_by(sideA, sideB) %>% 
  summarise(period = max(period)) %>% 
  ungroup() %>% 
  pull(period)

## Add a year attribute to the combined graph ------------------------------
E(main_graph)$period <- period

## Save network graph for easy loading ----------------------------------------
saveRDS(main_graph, file = file.path(dat_dir, "cartel_network.rds"), ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

## -----------------------------------------------------------------------------
################################################################################


