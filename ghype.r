# 1. Carregar bibliotecas
library(ghypernet)
library(igraph)

# 2. Ler os arquivos de dados
cat("Lendo os arquivos 'nodes.txt' e 'edges.txt'...\n")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Por favor, forneça o parâmetro 'size'. Uso: Rscript ghype.r <size>")
}
size <- as.numeric(args[1])
num_networks <- 100

nodes_df <- read.table(paste0("./input/gHypEG/nodes_", size, ".txt"), header = FALSE, col.names = c("node", "category"), stringsAsFactors = FALSE)
edges_df <- read.table(paste0("./input/gHypEG/edges_", size, ".txt"), header = FALSE, col.names = c("from", "to"), stringsAsFactors = FALSE)
number_of_edges <- nrow(edges_df)
cat("Número de arestas lidas:", number_of_edges, "\n")
cat("Degree esperado: ", number_of_edges / nrow(nodes_df), "\n")
# Garantir que os labels sejam tratados como strings
nodes_df$category <- as.character(nodes_df$category)

# 3. Preparar os dados para o modelo
all_nodes <- sort(unique(nodes_df$node))
num_nodes <- length(all_nodes)
adj_matrix <- matrix(0, nrow = num_nodes, ncol = num_nodes,
                     dimnames = list(as.character(all_nodes), as.character(all_nodes)))

for (i in 1:nrow(edges_df)) {
  adj_matrix[as.character(edges_df$from[i]), as.character(edges_df$to[i])] <- 1
}
node_labels <- setNames(nodes_df$category, nodes_df$node)
node_labels <- node_labels[rownames(adj_matrix)]

# 4. Estimar os parâmetros com ghypernet
model_params <- bccm(
  adj = adj_matrix,      # Note que o argumento muda de 'graph' para 'adj' em algumas versões, ou mantenha graph se sua versão aceitar
  labels = node_labels,  # <--- O INGREDIENTE ESSENCIAL
  directed = TRUE,
  selfloops = FALSE
)
rm(adj_matrix)

# Função para calcular proporções 'de -> para'. O uso de unname() é mantido como boa prática.
calculate_dyadic_proportions <- function(g, all_labels) {
  results <- list(
    IN = sapply(all_labels, function(x) sapply(all_labels, function(y) list(), simplify = FALSE), simplify = FALSE),
    OUT = sapply(all_labels, function(x) sapply(all_labels, function(y) list(), simplify = FALSE), simplify = FALSE)
  )

  for (v_id in V(g)) {
    source_label_of_v <- V(g)$category[v_id]

    # --- ANÁLISE DE SAÍDA (OUT) ---
    total_out_degree <- degree(g, v_id, mode = "out")
    if (total_out_degree > 0) {
      out_neighbors <- neighbors(g, v_id, mode = "out")
      out_neighbor_labels <- V(g)$category[out_neighbors]
      counts <- table(factor(out_neighbor_labels, levels = all_labels))
      props <- counts / total_out_degree
      for (target_label in names(props)) {
        results$OUT[[source_label_of_v]][[target_label]] <- c(results$OUT[[source_label_of_v]][[target_label]], unname(props[target_label]))
      }
    }

    # --- ANÁLISE DE ENTRADA (IN) ---
    total_in_degree <- degree(g, v_id, mode = "in")
    if (total_in_degree > 0) {
      in_neighbors <- neighbors(g, v_id, mode = "in")
      in_neighbor_labels <- V(g)$category[in_neighbors]
      counts <- table(factor(in_neighbor_labels, levels = all_labels))
      props <- counts / total_in_degree
      for (source_label_of_neighbor in names(props)) {
        results$IN[[source_label_of_neighbor]][[source_label_of_v]] <- c(results$IN[[source_label_of_neighbor]][[source_label_of_v]], unname(props[source_label_of_neighbor]))
      }
    }
  }
  return(results)
}

# 5. Preparar para o Loop de Simulação
dir.create("output/gHypEG", recursive = TRUE, showWarnings = FALSE)
all_possible_labels <- sort(unique(nodes_df$category))

aggregated_proportions <- list(
  IN = sapply(all_possible_labels, function(x) sapply(all_possible_labels, function(y) list(), simplify = FALSE), simplify = FALSE),
  OUT = sapply(all_possible_labels, function(x) sapply(all_possible_labels, function(y) list(), simplify = FALSE), simplify = FALSE)
)

clust_coefs <- numeric(num_networks); avg_clusts <- numeric(num_networks); diams <- numeric(num_networks); mean_paths <- numeric(num_networks); avg_degrees <- numeric(num_networks)
corr_degree_clustering <- numeric(num_networks)
# 6. Loop de Simulação
for(i in 1:num_networks) {
  if(i %% 10 == 0) cat("Modelo", i, ": Processando...\n")
  graph_sim_matrix <- rghype(1, model_params, seed = i+42*size,multinomial = TRUE,m =number_of_edges)
  g_sim <- graph_from_adjacency_matrix(graph_sim_matrix, mode="directed")
  cat("Rede simulada", i, "- Nós:", vcount(g_sim), "- Arestas:", ecount(g_sim), "\n")
  V(g_sim)$category <- node_labels
  V(g_sim)$name <- as.character(1:vcount(g_sim))
  total_loops <- sum(which_loop(g_sim))

  total_multiples <- sum(which_multiple(g_sim))

  cat(" - Self-loops encontrados:", total_loops, "\n")
  cat(" - Arestas múltiplas encontradas:", total_multiples, "\n")
  g_sim <- simplify(g_sim)
  isolated_nodes <- which(degree(g_sim, mode = "all") == 0)
  g_sim <- delete_vertices(g_sim, isolated_nodes)
  # Ajustar node_labels para corresponder aos nós restantes no grafo
  remaining_node_names <- V(g_sim)$name
  new_node_labels <- node_labels[remaining_node_names]
  V(g_sim)$category <- new_node_labels
  rm(graph_sim_matrix)
  
  proportions_this_sim <- calculate_dyadic_proportions(g_sim, all_possible_labels)
  
  for (from_label in all_possible_labels) {
    for (to_label in all_possible_labels) {
      aggregated_proportions$IN[[from_label]][[to_label]] <- c(aggregated_proportions$IN[[from_label]][[to_label]], proportions_this_sim$IN[[from_label]][[to_label]])
      aggregated_proportions$OUT[[from_label]][[to_label]] <- c(aggregated_proportions$OUT[[from_label]][[to_label]], proportions_this_sim$OUT[[from_label]][[to_label]])
    }
  }
  
  clust_coefs[i] <- transitivity(g_sim, type="global"); avg_clusts[i] <- transitivity(g_sim, type="average"); diams[i] <- diameter(g_sim); mean_paths[i] <- mean_distance(g_sim); avg_degrees[i] <- mean(degree(g_sim))
  node_degrees <- degree(g_sim, mode = "all")
  local_clustering <- transitivity(g_sim, type = "local")
  
  # Calculate correlation, handling NaN values from local_clustering for nodes with degree < 2
  valid_indices <- !is.na(local_clustering)
  if(sum(valid_indices) > 2) { # Correlation requires at least 2 valid pairs
    corr_degree_clustering[i] <- cor(node_degrees[valid_indices], local_clustering[valid_indices])
  } else {
    corr_degree_clustering[i] <- NA # Not enough data to compute correlation
  }
  rm(g_sim)
}

# 7. Salvar resultados das métricas
write.table(data.frame(avg_degrees, clust_coefs, avg_clusts, diams, mean_paths, corr_degree_clustering), 
            file = paste0("output/gHypEG/metrics/gHypEG_", size, ".txt"), 
            row.names = FALSE, col.names = FALSE, sep = " ")

# ==============================================================================
# SEÇÃO 8 CORRIGIDA: Usamos unlist() antes de chamar a função hist().
# ==============================================================================

hist_breaks <- seq(0, 1.01, by = 0.01)

for (file_label in all_possible_labels) {
  out_row_list <- list()
  for (target_label in all_possible_labels) {
    # CORREÇÃO: "Achata" a lista de listas em um único vetor numérico
    data_vector <- unlist(aggregated_proportions$OUT[[file_label]][[target_label]])
    hist_data <- hist(data_vector, breaks = hist_breaks, plot = FALSE, right = FALSE)
    out_row_list[[target_label]] <- hist_data$counts
  }
  matrix_out <- do.call(rbind, out_row_list)
  filename_out <- paste0("output/gHypEG/matrix/Matrix_gHypEG_out_", file_label, "_", size, ".txt")
  write.table(matrix_out, file = filename_out, col.names = FALSE, row.names = FALSE, sep = " ")

  in_row_list <- list()
  for (source_label in all_possible_labels) {
    # CORREÇÃO: "Achata" a lista de listas em um único vetor numérico
    data_vector <- unlist(aggregated_proportions$IN[[source_label]][[file_label]])
    hist_data <- hist(data_vector, breaks = hist_breaks, plot = FALSE, right = FALSE)
    in_row_list[[source_label]] <- hist_data$counts
  }
  matrix_in <- do.call(rbind, in_row_list)
  filename_in <- paste0("output/gHypEG/matrix/Matrix_gHypEG_in_", file_label, "_", size, ".txt")
  write.table(matrix_in, file = filename_in, col.names = FALSE, row.names = FALSE, sep = " ")
}

