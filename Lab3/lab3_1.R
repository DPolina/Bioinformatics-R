library("graph")
library("stringdist")
library("eulerian")
library("ggplot2")

# Исходная последовательность
true_sequence <- "CGCTACCGTCTGAAGAA"

# Риды
reads <- c("CTACCGTCT", "ACCGTCTGA", "CGCTACCGT", "GCTACCGTC", "GTCTGAAGA", "TACCGTCTG", "TCTGAAGAA")

# Функция для извлечения k-меров
get_kmers <- function(reads, k) {
  kmers <- c()
  for (read in reads) {
    for (i in 1:(nchar(read) - k + 1)) {
      kmers <- c(kmers, substr(read, i, i + k - 1))
    }
  }
  return(unique(kmers))
}

# Функция для построения графа де Брюйна
build_de_bruijn_graph <- function(kmers) {
  prefix <- substr(kmers, 1, nchar(kmers) - 1)
  suffix <- substr(kmers, 2, nchar(kmers))
  nodes <- unique(c(prefix, suffix))
  gr <- graphNEL(nodes = nodes, edgemode = "directed")
  for (kmer in kmers) {
    from <- substr(kmer, 1, nchar(kmer) - 1)
    to <- substr(kmer, 2, nchar(kmer))
    gr <- addEdge(graph = gr, from = from, to = to)
  }
  return(gr)
}

# Функция для сборки генома
assemble_genome <- function(path) {
  genome <- path[1]
  for (i in 2:length(path)) {
    genome <- paste0(genome, substr(path[i], nchar(path[i]), nchar(path[i])))
  }
  return(genome)
}

# Основной цикл по k
results <- data.frame(k = integer(), deviation = numeric())
for (k in 2:9) {
  kmers <- get_kmers(reads, k)
  graph <- build_de_bruijn_graph(kmers)
  eulerian_path <- tryCatch(eulerian(graph), error = function(e) NULL)
  if (!is.null(eulerian_path)) {
    genome <- assemble_genome(eulerian_path)
    # Вычисляем расстояние Левенштейна как отклонение
    deviation <- stringdist(genome, true_sequence, method = "lv")
    results <- rbind(results, data.frame(k = k, deviation = deviation))
  } else {
    # Если граф не имеет эйлерова пути, задаем большое отклонение
    results <- rbind(results, data.frame(k = k, deviation = NA))
  }
}
results_clean <- results[!is.na(results$deviation), ]

# Замена NA на большие значения отклонений (чтобы отобразить на графике)
results$deviation[is.na(results$deviation)] <- Inf

# Построение графика
ggplot(results, aes(x = k, y = deviation)) +
  geom_line(data = results[is.finite(results$deviation), ], aes(x = k, y = deviation)) +
  geom_point() +
  scale_y_continuous(
    limits = c(0, max(results$deviation[is.finite(results$deviation)])),
    breaks = c(seq(0, max(results$deviation[is.finite(results$deviation)]), by = 2), Inf),
    labels = c(seq(0, max(results$deviation[is.finite(results$deviation)]), by = 2), "Нет пути")
  ) +
  labs(
    title = "Отклонения восстановленной последовательности от длины k-меров",
    x = "Длина k-меров (k)",
    y = "Отклонение (расстояние Левенштейна)"
  ) +
  theme_minimal()