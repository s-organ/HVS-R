library(dplyr)
library(glmnet)
library(gor)
library(circlize)
library(igraph)
library(SetSel)

####### INPUTS ##############################

# x: a design matrix
# y: response vector
# q = 0.05 (or alternative desired FDR control level)

# 1) hypothesis testing

p = ncol(x)
n = nrow(x)

reg = lm(y ~ x)
pairs = t(combn(c(1:p),2))
p_val = c()

# pairs #

for(i in 1:nrow(pairs)){
  var1 = pairs[i,1]
  var2 = pairs[i,2]
  p1 = x[,var1]
  p2 = x[,var2]
  x2 = x[,-c(var1,var2)]
  nested = lm(y ~ x2)
  p_val[i] = anova(reg,nested)[2,6]
  print(i)
}

p_val_pairs = p_val

# singles 

lreg <- summary(lm(y ~ x))$coefficients
p_val_singles <- lreg[-1, 4]  # Remove intercept

# Combine edges
singles <- cbind(1:p, 1:p)
colnames(singles) <- c("V1", "V2")
results <- rbind(pairs, singles)
S <- apply(results, 1, function(row) as.integer(row))  # list of edges
S <- lapply(1:nrow(results), function(i) unname(unlist(results[i, ])))  # as list

# Combined p-values 
p_vals_combined <- c(p_val_pairs, p_val_singles)

# 2) alpha functions ########################

alpha_HVSPRDS <- function(q, p) {
  log_ratio <- log(4) - log(3)
  return((1/q) * p + p * (log_ratio / log(2)))
}
alpha <- alpha_HVSPRDS(q, p)


alpha_HVS <- function(q, p0, p) {
  term1 <- q * p0 * (1 + log(p))
  term2 <- (p0 / 2) * (log(4) - log(3)) / log(2)
  numerator <- 1 + log(p) - log(log(4)) - log(3)
  
  alpha <- numerator / (term1 + term2)
  return(alpha)
}

# 3) sample independent sets function #######

approximate_independent_sets <- function(graph, N = 1000) {
  n <- vcount(graph)
  total_weight <- 0
  for (i in 1:N) {
    remaining <- 1:n
    indep_set <- c()
    while (length(remaining) > 0) {
      v <- sample(remaining, 1)
      indep_set <- c(indep_set, v)
      neighbors_v <- neighbors(graph, v)
      remaining <- setdiff(remaining, c(v, neighbors_v))
    }
    w <- 2^(n - length(indep_set))
    total_weight <- total_weight + w
  }
  return(total_weight / N)
}

# 4) get graph ###################

build_full_graph <- function(S) {
  edge_list <- lapply(S, function(e) {
    if (length(e) == 2 && e[1] != e[2]) {
      return(c(e[1], e[2]))
    } else {
      return(NULL)
    }
  })
  edge_list <- do.call(rbind, Filter(Negate(is.null), edge_list))
  if (!is.null(edge_list) && nrow(edge_list) > 0) {
    g <- graph_from_edgelist(edge_list, directed = FALSE)
  } else {
    g <- make_empty_graph()
  }
  return(g)
}

# 5) sigma function #############

sigma <- function(selected_indices, S, p, full_graph) {
  selected_edges <- S[selected_indices]
  nodes <- sort(unique(unlist(selected_edges)))
  if (length(nodes) == 0) return(-Inf)
  
  if (length(nodes) == 1) {
    num_vertex_covers <- 1  # singleton graph: only one vertex
  } else {
    g_sub <- induced_subgraph(full_graph, vids = as.character(nodes))
    n <- vcount(g_sub)
    if (n == 0) return(-Inf)
    ind_sets <- approximate_independent_sets(g_sub, N = 1000)  # reduce N for speed
    num_vertex_covers <- 2^n - ind_sets
    if (num_vertex_covers <= 0) return(-Inf)
  }
  
  return(p - log2(num_vertex_covers))
}

# 6) selction ######################

full_graph <- build_full_graph(S)

# Step-up procedure with incremental selection
ord <- order(p_vals_combined)
sorted_pvals <- p_vals_combined[ord]
best_c <- 0
best_set <- NULL
current_set <- c()

for (i in seq_along(sorted_pvals)) {
  current_set <- c(current_set, ord[i])
  c <- sorted_pvals[i]
  
  sig_val <- sigma(current_set, S, p, full_graph)
  
  if (sig_val >= alpha * c && c > best_c) {
    best_c <- c
    best_set <- current_set
  }
  if (i %% 10 == 0) cat("Checked", i, "of", length(sorted_pvals), "\n")
}

rejected = best_set 

rejected_sets <- S[rejected]

# 7) get graph of final selection ############

build_graph_from_S <- function(S) {
  # Extract edges: keep only pairs with two distinct nodes
  edges <- do.call(rbind, lapply(S, function(e) {
    if (length(e) == 2 && e[1] != e[2]) {
      return(e)
    } else {
      return(NULL)
    }
  }))
  
  # Create graph from edges (will auto-create nodes)
  if (!is.null(edges) && nrow(edges) > 0) {
    g <- graph_from_edgelist(edges, directed = FALSE)
  } else {
    # No edges, just create empty graph
    g <- make_empty_graph()
  }
  
  # Add singleton nodes (those in S with length 1) if not present as vertices
  singletons <- unique(unlist(S[sapply(S, length) == 1]))
  missing_singletons <- setdiff(as.character(singletons), V(g)$name)
  if (length(missing_singletons) > 0) {
    g <- add_vertices(g, length(missing_singletons), name = missing_singletons)
  }
  
  return(g)
}

g_final <- build_graph_from_S(rejected_sets) # final selected graph

########### POWER AND FDR IF DESIRED ####################

# assuming truth is a vector true variables

false = which(!c(1:p) %in% truth)

false_sets = Filter(function(vars) all(vars %in% false), rejected_sets)

false_graph = build_graph_from_S(false_sets) 

FDR = (p-log2(approximate_independent_sets(false_graph,1000)))/(p-log2(approximate_independent_sets(g_final,1000)))

true_sets = Filter(function(vars) any(vars %in% truth), S)

true_graph = build_graph_from_S(true_sets)

power = (p-log2(approximate_independent_sets(true_graph,1000)))/(length(truth))

############## IF WANTING MSE ####################

# build cover greedy is from the gor package, builds the vertex cover from the greedy algorithm
# this is just an example, any algorithm can be used

cover = build_cover_greedy(g_final)

cover$set