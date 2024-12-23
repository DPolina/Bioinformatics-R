library(graph)
library(eulerian)
library(stringdist)


# paste: Объединяет строки в одну строку с заданным разделителем
# substr: Извлекает или заменяет подстроки в строках
# nchar: Возвращает количество символов в строке
# unique: Удаляет дубликаты из вектора или списка
# graphNEL: Создаёт пустой направленный или ненаправленный граф
# addEdge: Добавляет рёбра между вершинами в объекте graphNEL
# eulerian: Находит эйлеров путь (или цикл) в графе.

# Шаг 1. Разбиение ридов на k-меры

# Риды — это короткие фрагменты ДНК, которые получаются в результате секвенирования. 
# Чтобы собрать длинную последовательность ДНК, риды разбиваются на k-меры — 
# последовательности длиной k. Каждый k-мер состоит из префикса (первых k-1 символов) 
# и суффикса (последних k-1 символов). Эти префиксы и суффиксы используются как вершины графа.

reads <- c("CTACCGTCT", "ACCGTCTGA", "CGCTACCGT", "GCTACCGTC", "GTCTGAAGA", "TACCGTCTG", "TCTGAAGAA")
k <- 7 # длина k-меров (2-9)

# Для каждого рида извлекаются все возможные последовательности длиной k с помощью substr.
get_kmers <- function(reads, k) {
  kmers <- c()
  for (read in reads) {
    for (i in 1:(nchar(read) - k + 1)) {
      kmers <- c(kmers, substr(read, i, i + k - 1))
    }
  }
  return(kmers)
}

kmers <- get_kmers(reads, k)
print(kmers)
kmers <- unique(get_kmers(reads, k))
print(kmers)

# Шаг 2: Построение графа де Брюйна

# Граф де Брюйна — это ориентированный граф, где:
  # - Вершины — это уникальные префиксы и суффиксы длиной k-1.
  # - Ребра — это переходы от префикса k-мера к его суффиксу.

# Такой граф используется для сборки генома, так как он отражает перекрытие 
# последовательностей. Эйлеров путь в этом графе (путь, проходящий через каждое 
# ребро ровно один раз) представляет итоговую последовательность.

build_de_bruijn_graph <- function(kmers){
  
  # Создаем уникальные вершины (префиксы и суффиксы)
  prefix <- substr(kmers, 1, nchar(kmers) - 1)
  suffix <- substr(kmers, 2, nchar(kmers))
  nodes <- unique(c(prefix, suffix))
  
  # Создаем пустой граф
  gr <- graphNEL(nodes = nodes, edgemode = "directed")
  
  # Добавляем ребра
  for(kmer in kmers){
    from <- substr(kmer, 1, nchar(kmer) - 1)
    to <- substr(kmer, 2, nchar(kmer))
    gr <- addEdge(graph = gr, from = from, to = to)
  }
  return(gr)
}

# Построение графа
graph <- build_de_bruijn_graph(kmers)
print(graph)
#plot(grph, 'circo')
plot(graph)

# Шаг 3: Нахождение Эйлерова пути

# Эйлеров путь — это путь, который проходит через каждое ребро графа ровно один раз.

library(eulerian)

eulerian_path <- eulerian(graph)
print(eulerian_path)

# Шаг 4: Сборка последовательности

# После нахождения Эйлерова пути восстанавливаем длинную последовательность ДНК. Для этого:
  # - Берется первая вершина пути.
  # - К ней добавляется последний символ каждой последующей вершины.

assemble_genome <- function(path) {
  genome <- path[1]
  for (i in 2:length(path)) {
    genome <- paste0(genome, substr(path[i], nchar(path[i]), nchar(path[i])))
  }
  return(genome)
}

genome <- assemble_genome(eulerian_path)
print(genome)
