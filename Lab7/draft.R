# Шаг 1: Загрузка данных
file_path <- "data2.txt"
data <- as.matrix(read.table(file_path, stringsAsFactors = FALSE))
n <- nrow(data)  # Количество объектов
head(data)

# Шаг 2: Вычисление матрицы G (двойное центрирование матрицы D2)
D2 <- as.matrix(dist(data, method = "euclidean"))  # Пример меры расстояния (оптимальное выравнивание можно подставить)
D2_j <- colMeans(D2)  # Средние по столбцам
D2_i <- rowMeans(D2)  # Средние по строкам
D2_all <- mean(D2)    # Общее среднее значение

G <- matrix(0, n, n)  # Инициализация матрицы G
for (i in 1:n) {
  for (j in 1:n) {
    G[i, j] <- 0.5 * (D2_all - D2_j[j] - D2_i[i] + D2[i, j])
  }
}

# Шаг 3: Степенной метод для нахождения собственных значений и векторов
power_method <- function(A, max_iter = 1000, tol = 1e-6) {
  b <- rep(1, nrow(A))  # Начальный вектор
  for (iter in 1:max_iter) {
    b_new <- A %*% b
    b_new <- b_new / sqrt(sum(b_new^2))  # Нормализация
    if (sqrt(sum((b_new - b)^2)) < tol) break
    b <- b_new
  }
  eigenvalue <- sum(b * (A %*% b)) / sum(b^2)
  return(list(eigenvalue = eigenvalue, eigenvector = b))
}

eigen_results <- power_method(G)
lambda <- eigen_results$eigenvalue
v <- eigen_results$eigenvector

# Шаг 4: Проекции объектов на главные координаты
M <- 10  # Количество главных координат для анализа
eigen_values <- c()
eigen_vectors <- matrix(0, n, M)

for (k in 1:M) {
  res <- power_method(G)
  eigen_values <- c(eigen_values, res$eigenvalue)
  eigen_vectors[, k] <- res$eigenvector
  G <- G - res$eigenvalue * (res$eigenvector %*% t(res$eigenvector))  # Дефляция
}

# Проекции на первые M главных координат
projections <- eigen_vectors[, 1:M]

# Шаг 5: Визуализация результатов
# Диаграмма рассеяния для первых двух главных координат
plot(projections[, 1], projections[, 2], 
     xlab = "Первая главная координата", 
     ylab = "Вторая главная координата", 
     main = "Диаграмма рассеяния для первых двух главных координат", 
     col = "blue", pch = 19)

# График зависимости собственных значений от порядкового номера
plot(1:M, eigen_values, 
     type = "b", 
     xlab = "Порядковый номер главной координаты", 
     ylab = "Собственное значение", 
     main = "График собственных значений")

# Определение оптимального числа главных координат по "локтю"
