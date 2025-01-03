library(stringdist)
library(cluster)
library(dendextend)
library(factoextra)

# Шаг 1. Загрузка данных
file_path <- "data2.txt"
data <- as.matrix(read.table(file_path))

# Генерируем 10 новых последовательностей
sequence_length <- 90
num_sequences <- 10

generate_sequence <- function(length) {
  paste0(sample(c("A", "C", "G", "T"), length, replace = TRUE), collapse = "")
}

new_sequences <- replicate(num_sequences, generate_sequence(sequence_length))
cat("Сгенерированные последовательности:\n")
print(new_sequences)

# Приведение формата добавленных строк к виду существующих
new_sequences <- gsub("\"", "", new_sequences)  # Удаляем кавычки из новых строк

# Добавляем новые последовательности к data
data <- c(data, new_sequences)
cat("length(data):", length(data), "\n")
print(data)

# Шаг 2. Вычисление расстояний (различий между строками)

# Левенштейн измеряет минимальное количество операций (вставка, удаление, замена), 
# чтобы преобразовать одну строку в другую.
levenshtein_dist <- stringdistmatrix(data, data, method = "lv")

# Дамерау-Левенштейн учитывает также транспозицию двух соседних символов.
damerau_levenshtein_dist <- stringdistmatrix(data, data, method = "dl")

# Q-грамм измеряет схожесть строк на основе последовательностей символов длиной 
# q (обычно 2 или 3).
qgram_dist <- stringdistmatrix(data, data, method = "qgram")

# Шаг 3. Иерархическая кластеризация
# Преобразование матриц расстояний в объекты dist
levenshtein_dist_obj <- as.dist(levenshtein_dist)
damerau_levenshtein_dist_obj <- as.dist(damerau_levenshtein_dist)
qgram_dist_obj <- as.dist(qgram_dist)

# Левенштейн - методы кластеризации
hclust_lev_single <- hclust(levenshtein_dist_obj, method = "single")
hclust_lev_complete <- hclust(levenshtein_dist_obj, method = "complete")
hclust_lev_centroid <- hclust(levenshtein_dist_obj, method = "centroid")

# Дамерау-Левенштейн - методы кластеризации
hclust_dl_single <- hclust(damerau_levenshtein_dist_obj, method = "single")
hclust_dl_complete <- hclust(damerau_levenshtein_dist_obj, method = "complete")
hclust_dl_centroid <- hclust(damerau_levenshtein_dist_obj, method = "centroid")

# q-грамм - методы кластеризации
hclust_qgram_single <- hclust(qgram_dist_obj, method = "single")
hclust_qgram_complete <- hclust(qgram_dist_obj, method = "complete")
hclust_qgram_centroid <- hclust(qgram_dist_obj, method = "centroid")

# Шаг 4. Анализ качества кластеризации (кофенетический коэффициент)

# Кофенетический коэффициент оценивает, насколько хорошо кластеризация совпадает 
# с исходными расстояниями между объектами.

# Левенштейн
coph_lev_single <- cor(cophenetic(hclust_lev_single), levenshtein_dist_obj)
coph_lev_complete <- cor(cophenetic(hclust_lev_complete), levenshtein_dist_obj)
coph_lev_centroid <- cor(cophenetic(hclust_lev_centroid), levenshtein_dist_obj)

# Дамерау-Левенштейн
coph_dl_single <- cor(cophenetic(hclust_dl_single), damerau_levenshtein_dist_obj)
coph_dl_complete <- cor(cophenetic(hclust_dl_complete), damerau_levenshtein_dist_obj)
coph_dl_centroid <- cor(cophenetic(hclust_dl_centroid), damerau_levenshtein_dist_obj)

# q-грамм
coph_qgram_single <- cor(cophenetic(hclust_qgram_single), qgram_dist_obj)
coph_qgram_complete <- cor(cophenetic(hclust_qgram_complete), qgram_dist_obj)
coph_qgram_centroid <- cor(cophenetic(hclust_qgram_centroid), qgram_dist_obj)

# Таблица кофенетических коэффициентов
coph_table <- data.frame(
  Metric = c("Levenshtein", "Damerau-Levenshtein", "Q-gram"),
  Single = c(coph_lev_single, coph_dl_single, coph_qgram_single),
  Complete = c(coph_lev_complete, coph_dl_complete, coph_qgram_complete),
  Centroid = c(coph_lev_centroid, coph_dl_centroid, coph_qgram_centroid)
)
print("Кофенетические коэффициенты:")
print(coph_table)

# Шаг 5. Определение наиболее и наименее эффективных способов
# Преобразуем таблицу в числовую матрицу для вычислений
numeric_table <- as.matrix(coph_table[, 2:4])  # Извлекаем только числовые столбцы

# Определение наилучшего метода
max_index <- which(numeric_table == max(numeric_table), arr.ind = TRUE)
best_method <- coph_table[max_index[1], ]

# Определение наихудшего метода
min_index <- which(numeric_table == min(numeric_table), arr.ind = TRUE)
worst_method <- coph_table[min_index[1], ]

print("Наиболее эффективный способ кластеризации:")
print(best_method)

print("Наименее эффективный способ кластеризации:")
print(worst_method)

# Шаг 6. Построение дендрограммы для наиболее эффективного метода
best_hclust <- hclust_qgram_centroid  # Q-gram и Centroid
dend <- as.dendrogram(best_hclust)
plot(dend, main = "Дендрограмма (Q-gram, Centroid)", sub = "", xlab = "")
rect.hclust(best_hclust, k = k, border = "red")  # Рисуем разрез для 3 кластеров

# Шаг 7. Определение количества кластеров

# Метод: Использование фиксированного числа кластеров (например, k = 3)
k <- 4  # Количество кластеров
# Разделяем на кластеры
clusters <- cutree(best_hclust, k = k)
cat("Определено", length(unique(clusters)), "кластера.\n")
print(table(clusters))
#save(clusters, file = "clusters.RData")

# Пункт 8. Рассчет медоидов и отображение их на дендрограмме

# Медоид — это объект, который представляют собой центральные элементы кластеров 
# и имеет наименьшее среднее расстояние до всех остальных объектов в своем кластере. 
# В отличие от центроидов, которые являются средним арифметическим точек в 
# пространстве, медоиды всегда являются реальными объектами данных.

# Рассчёт медоидов вручную

# Функция для расчёта медоидов
calculate_medoids <- function(dist_matrix, clusters) {
  unique_clusters <- unique(clusters)  # Получаем уникальные кластеры
  medoids <- vector("list", length(unique_clusters))  # Для хранения медоидов
  
  for (cluster_id in unique_clusters) {
    # Выбираем объекты текущего кластера
    cluster_indices <- which(clusters == cluster_id)
    
    # Субматрица расстояний для объектов текущего кластера
    cluster_distances <- as.matrix(dist_matrix)[cluster_indices, cluster_indices]
    
    # Сумма расстояний для каждого объекта кластера
    total_distances <- rowSums(cluster_distances)
    
    # Индекс объекта с минимальной суммой расстояний
    medoid_index <- cluster_indices[which.min(total_distances)]
    
    # Сохраняем медоид
    medoids[[cluster_id]] <- medoid_index
  }
  
  return(unlist(medoids))
}

# Используем матрицу расстояний q-грамм и результаты кластеризации
medoid_indices <- calculate_medoids(qgram_dist_obj, clusters)

# Выводим индексы медоидов
cat("Индексы медоидов кластеров:\n")
print(medoid_indices)

# Создаём вектор цветов для визуализации медоидов
colors <- rep(2, length(data))  # Все объекты будут иметь красный цвет
colors[medoid_indices] <- 3  # Медоиды будут иметь зеленый цвет

# Рисуем дендрограмму с подсветкой медоидов
dend %>% plot(main = "Дендрограмма с медоидами (Q-gram, Centroid)")
dend %>% colored_bars(colors = colors, dend = dend)

# 
