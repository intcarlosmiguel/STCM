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

nodes_df <- read.table(paste0("./input/gHypEG/nodes_", size, ".txt"), header = FALSE, col.names = c("node", "category"), stringsAsFactors = FALSE)
edges_df <- read.table(paste0("./input/gHypEG/edges_", size, ".txt"), header = FALSE, col.names = c("from", "to"), stringsAsFactors = FALSE)

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
model_params <- ghype(
  graph = adj_matrix,
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

clust_coefs <- numeric(100); avg_clusts <- numeric(100); diams <- numeric(100); mean_paths <- numeric(100); avg_degrees <- numeric(100)

# 6. Loop de Simulação
for(i in 1:100) {
  if(i %% 10 == 0) cat("Modelo", i, ": Processando...\n")
  graph_sim_matrix <- rghype(1, model_params, seed = i+42*size)
  g_sim <- graph_from_adjacency_matrix(graph_sim_matrix, mode="directed")
  V(g_sim)$category <- node_labels 
  rm(graph_sim_matrix)
  
  proportions_this_sim <- calculate_dyadic_proportions(g_sim, all_possible_labels)
  
  for (from_label in all_possible_labels) {
    for (to_label in all_possible_labels) {
      aggregated_proportions$IN[[from_label]][[to_label]] <- c(aggregated_proportions$IN[[from_label]][[to_label]], proportions_this_sim$IN[[from_label]][[to_label]])
      aggregated_proportions$OUT[[from_label]][[to_label]] <- c(aggregated_proportions$OUT[[from_label]][[to_label]], proportions_this_sim$OUT[[from_label]][[to_label]])
    }
  }
  
  clust_coefs[i] <- transitivity(g_sim, type="global"); avg_clusts[i] <- transitivity(g_sim, type="average"); diams[i] <- diameter(g_sim); mean_paths[i] <- mean_distance(g_sim); avg_degrees[i] <- mean(degree(g_sim))
  rm(g_sim)
}

# 7. Salvar resultados das métricas
write.table(data.frame(avg_degrees, clust_coefs, avg_clusts, diams, mean_paths), 
            file = paste0("output/gHypEG/gHypEG_results_", size, ".txt"), 
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
  filename_out <- paste0("output/gHypEG/Matrix_gHypEG_out_", file_label, "_", size, ".txt")
  write.table(matrix_out, file = filename_out, col.names = FALSE, row.names = FALSE, sep = " ")

  in_row_list <- list()
  for (source_label in all_possible_labels) {
    # CORREÇÃO: "Achata" a lista de listas em um único vetor numérico
    data_vector <- unlist(aggregated_proportions$IN[[source_label]][[file_label]])
    hist_data <- hist(data_vector, breaks = hist_breaks, plot = FALSE, right = FALSE)
    in_row_list[[source_label]] <- hist_data$counts
  }
  matrix_in <- do.call(rbind, in_row_list)
  filename_in <- paste0("output/gHypEG/Matrix_gHypEG_in_", file_label, "_", size, ".txt")
  write.table(matrix_in, file = filename_in, col.names = FALSE, row.names = FALSE, sep = " ")
}

