# diag: создает или извлекает диагонали матриц.
# rownames / colnames: устанавливают или извлекают названия строк/столбцов матриц и таблиц.
# strsplit: разделяет строки на подстроки по заданному разделителю.
# unlist: преобразует списки в векторы.
# length: определяет длину объекта (вектора, строки или списка).
# paste: объединяет строки или символы.

# 1. Создание матрицы замен S
create_substitution_matrix <- function(match = 1, mismatch = -1) {
  nucleotides <- c("A", "C", "G", "T")
  S <- matrix(mismatch, nrow = 4, ncol = 4, dimnames = list(nucleotides, nucleotides))
  diag(S) <- match  # Совпадения имеют балл match
  return(S)
}

# 2. Инициализация точечной матрицы сходства F
initialize_matrices <- function(seq1, seq2, gap = -1) {
  n <- nchar(seq1) + 1
  m <- nchar(seq2) + 1
  F <- matrix(0, nrow = n, ncol = m)  # Точечная матрица сходства
  traceback <- matrix("", nrow = n, ncol = m)  
  # traceback - матрица направлений, которая указывает, откуда пришло значение в 
  # каждой ячейке точечной матрицы сходства.
  
  # Заполнение первой строки и столбца штрафами за удаление
  for (i in 2:n) {
    F[i, 1] <- F[i - 1, 1] + gap
    traceback[i, 1] <- "U"  # Up
  }
  for (j in 2:m) {
    F[1, j] <- F[1, j - 1] + gap
    traceback[1, j] <- "L"  # Left
  }
  
  # Установка подписей строк и столбцов
  rownames(F) <- c("-", strsplit(seq1, "")[[1]])
  colnames(F) <- c("-", strsplit(seq2, "")[[1]])
  
  list(F = F, traceback = traceback)
}

# 3. Заполнение точечной матрицы F с использованием матрицы замен S
fill_matrices <- function(seq1, seq2, F, traceback, S, gap = -1) {
  n <- nrow(F)
  m <- ncol(F)
  
  for (i in 2:n) {
    for (j in 2:m) {
      char1 <- rownames(F)[i]
      char2 <- colnames(F)[j]
      score <- S[char1, char2]  # Получение значения из матрицы замен
      
      # Рассчитываем варианты
      diagonal <- F[i - 1, j - 1] + score # Совпадение
      up <- F[i - 1, j] + gap # 
      left <- F[i, j - 1] + gap # 
      
      # Выбираем максимум
      max_score <- max(diagonal, up, left)
      F[i, j] <- max_score
      
      # Заполняем матрицу направлений
      if (max_score == diagonal) {
        traceback[i, j] <- "D"  # Diagonal
      } else if (max_score == up) {
        traceback[i, j] <- "U"  # Up
      } else {
        traceback[i, j] <- "L"  # Left
      }
    }
  }
  
  list(F = F, traceback = traceback)
}

# 4. Обратный проход для восстановления оптимального выравнивания
traceback <- function(seq1, seq2, traceback) {
  # Создаем пустые строки aligned_seq1 и aligned_seq2 для хранения итоговых 
  # последовательностей после выравнивания:
  aligned_seq1 <- ""
  aligned_seq2 <- ""
  # Начинаем с правого нижнего угла матрицы направлений traceback, который соответствует 
  # финальному состоянию выравнивания:
  i <- nrow(traceback)
  j <- ncol(traceback)
  
  # Двигаемся назад по матрице traceback до тех пор, пока не дойдем до верхнего 
  # левого угла (ячейка (1,1)):
  while (i > 1 || j > 1) {
    # На каждом шаге проверяем значение в текущей ячейке traceback[i, j]:
    if (traceback[i, j] == "D") {
      # Добавляем текущие символы из обеих последовательностей (из seq1[i-1] и 
      # seq2[j-1]) в выравненные строки:
      aligned_seq1 <- paste0(substr(seq1, i - 1, i - 1), aligned_seq1)
      aligned_seq2 <- paste0(substr(seq2, j - 1, j - 1), aligned_seq2)
      # Смещаемся по диагонали:
      i <- i - 1
      j <- j - 1
    } else if (traceback[i, j] == "U") {
      # Добавляем текущий символ из seq1[i-1] и символ пропуска "-" в aligned_seq2:
      aligned_seq1 <- paste0(substr(seq1, i - 1, i - 1), aligned_seq1)
      aligned_seq2 <- paste0("-", aligned_seq2)
      # Смещаемся вверх:
      i <- i - 1
    } else if (traceback[i, j] == "L") {
      # Добавляем символ пропуска "-" в aligned_seq1 и текущий символ из seq2[j-1]:
      aligned_seq1 <- paste0("-", aligned_seq1)
      aligned_seq2 <- paste0(substr(seq2, j - 1, j - 1), aligned_seq2)
      # Смещаемся влево:
      j <- j - 1
    }
  }
  # Возврат результата:
  list(aligned_seq1 = aligned_seq1, aligned_seq2 = aligned_seq2)
}

# 5. Основная функция
needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -1) {
  # Создаем матрицу замен S
  S <- create_substitution_matrix(match, mismatch)
  
  # Инициализируем матрицы F и traceback
  matrices <- initialize_matrices(seq1, seq2, gap)
  
  # Заполняем точечную матрицу сходства F
  filled_matrices <- fill_matrices(seq1, seq2, matrices$F, matrices$traceback, S, gap)
  
  # Обратный проход
  traceback_result <- traceback(seq1, seq2, filled_matrices$traceback)
  
  list(F = filled_matrices$F, 
       traceback = filled_matrices$traceback,
       aligned_seq1 = traceback_result$aligned_seq1, 
       aligned_seq2 = traceback_result$aligned_seq2)
}

# Вариант 2
seq1 <- "CCCTATATA"
seq2 <- "TACCCTATATACC"

result <- needleman_wunsch(seq1, seq2)

cat("Матрица замен S:\n")
print(create_substitution_matrix(1, -1))

cat("\nТочечная матрица сходства F:\n")
print(result$F)

# ------------------------------------------------ #

visualize_alignment <- function(aligned_seq1, aligned_seq2) {
  # Создаем строку для отображения совпадений
  matches <- ""
  for (i in seq_len(nchar(aligned_seq1))) {
    char1 <- substr(aligned_seq1, i, i)
    char2 <- substr(aligned_seq2, i, i)
    if (char1 == char2 && char1 != "-") {
      matches <- paste0(matches, "|")  # Совпадение
    } else {
      matches <- paste0(matches, " ")  # Несовпадение или пропуск
    }
  }
  
  cat("Выравнивание:\n")
  cat(aligned_seq1, "\n")
  cat(matches, "\n")
  cat(aligned_seq2, "\n")
}

# Визуальный анализ
result <- needleman_wunsch(seq1, seq2)
visualize_alignment(result$aligned_seq1, result$aligned_seq2)

