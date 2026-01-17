# 1. Carregar bibliotecas
library(ghypernet)
library(igraph)

# 2. Ler os arquivos de dados
cat("Lendo os arquivos 'nodes.txt' e 'edges.txt'...\n")

nodes_df <- read.table(paste0("./input/gHypEG/polymod_nodes.txt"), header = FALSE, col.names = c("node", "category"), stringsAsFactors = FALSE)
edges_df <- read.table(paste0("./input/gHypEG/polymod_edges.txt"), header = FALSE, col.names = c("from", "to"), stringsAsFactors = FALSE)
test_vector <- scan("test.txt", what = integer(), quiet = TRUE)

my_graph_igraph <- graph_from_data_frame(
  d = edges_df,        # Data frame com as arestas (colunas 'from' e 'to')
  vertices = nodes_df, # Data frame com os nós e seus atributos (colunas 'node', 'category')
  directed = FALSE      # Mude para FALSE se o seu grafo não for direcionado
)
# Garantir que os labels sejam tratados como strings
nodes_df$category <- as.character(nodes_df$category)



# Tente passar a matriz esparsa para ghype
# Não há garantia de que ghype aceitará uma matriz esparsa,
# mas vale a pena tentar, pois é uma matriz
model_params <- ghype(
  graph = my_graph_igraph,
  directed = FALSE,
  selfloops = FALSE,
)
cat("Parâmetros do modelo gHypEG calculados.\n")
graph_sim_matrix <- rghype(1, model_params, seed = 422,multinomial = TRUE,m =13783 )
cat("Grafo simulado com gHypEG carregado.\n")
g_sim <- graph_from_adjacency_matrix(graph_sim_matrix, mode="undirected")

cat("Mean Degree:", mean(degree(g_sim)), "\n")
cat("Number of Nodes:", vcount(g_sim), "\n")
cat("Mean Clustering:", transitivity(g_sim, type = "average"), "\n")
cat("Total Clustering:", transitivity(g_sim, type = "global"), "\n")

quit(save = "no", status = 0)
# 3. Extrair informações do grafo original
size <- vcount(my_graph_igraph)
node_labels <- V(my_graph_igraph)$category
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
corr_degree_clustering <- numeric(100)
# 6. Loop de Simulação
for(i in 1:100) {
  if(i %% 10 == 0) cat("Modelo", i, ": Processando...\n")
  graph_sim_matrix <- rghype(1, model_params, seed = i+42*size,multinomial = TRUE,m =test_vector[i] )
  g_sim <- graph_from_adjacency_matrix(graph_sim_matrix, mode="directed")

  V(g_sim)$category <- node_labels
  V(g_sim)$name <- as.character(1:vcount(g_sim))
  
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
            file = paste0("output/gHypEG/metrics/gHypEG_polymod.txt"), 
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
  filename_out <- paste0("output/gHypEG/matrix/Matrix_gHypEG_out_", file_label, "_polymod.txt")
  write.table(matrix_out, file = filename_out, col.names = FALSE, row.names = FALSE, sep = " ")

  in_row_list <- list()
  for (source_label in all_possible_labels) {
    # CORREÇÃO: "Achata" a lista de listas em um único vetor numérico
    data_vector <- unlist(aggregated_proportions$IN[[source_label]][[file_label]])
    hist_data <- hist(data_vector, breaks = hist_breaks, plot = FALSE, right = FALSE)
    in_row_list[[source_label]] <- hist_data$counts
  }
  matrix_in <- do.call(rbind, in_row_list)
  filename_in <- paste0("output/gHypEG/matrix/Matrix_gHypEG_in_", file_label, "_polymod.txt")
  write.table(matrix_in, file = filename_in, col.names = FALSE, row.names = FALSE, sep = " ")
}

