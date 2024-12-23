# Алгоритм Needleman-Wunsch для глобального выравнивания последовательностей

# 1. Создание матрицы замен S
# Матрица замен задает баллы за совпадения и несоответствия.
create_substitution_matrix <- function(match = 1, mismatch = -1) {
  nucleotides <- c("A", "C", "G", "T")  # Нуклеотиды
  S <- matrix(mismatch, nrow = 4, ncol = 4, dimnames = list(nucleotides, nucleotides))
  diag(S) <- match  # Устанавливаем баллы за совпадение
  return(S)
}

# 2. Инициализация матриц F (точечная матрица сходства) и traceback (матрица направлений)
initialize_matrices <- function(seq1, seq2, gap = -1) {
  n <- nchar(seq1) + 1  # Количество строк (длина первой последовательности + 1)
  m <- nchar(seq2) + 1  # Количество столбцов (длина второй последовательности + 1)
  
  F <- matrix(0, nrow = n, ncol = m)  # Матрица сходства
  traceback <- matrix("", nrow = n, ncol = m)  # Матрица направлений
  # Диагональ: Совпадение/замена
  # Вверх: удаление из первой последовательности
  # Влево: вставка во вторую последовательность
  
  # Инициализация первой строки и столбца штрафами за пропуски
  for (i in 2:n) {
    F[i, 1] <- F[i - 1, 1] + gap
    traceback[i, 1] <- "U"  # Направление вверх
  }
  for (j in 2:m) {
    F[1, j] <- F[1, j - 1] + gap
    traceback[1, j] <- "L"  # Направление влево
  }
  
  # Установка подписей строк и столбцов
  rownames(F) <- c("-", strsplit(seq1, "")[[1]])
  colnames(F) <- c("-", strsplit(seq2, "")[[1]])
  
  list(F = F, traceback = traceback)
}

# 3. Заполнение матриц F и traceback
fill_matrices <- function(seq1, seq2, F, traceback, S, gap = -1) {
  n <- nrow(F)  # Количество строк
  m <- ncol(F)  # Количество столбцов
  
  # Заполнение каждой ячейки
  for (i in 2:n) {
    for (j in 2:m) {
      char1 <- rownames(F)[i]  # Текущий символ из первой последовательности
      char2 <- colnames(F)[j]  # Текущий символ из второй последовательности
      score <- S[char1, char2]  # Балл из матрицы замен
      
      # Вычисляем три варианта:
      diagonal <- F[i - 1, j - 1] + score  # балл за совпадение или замену
      up <- F[i - 1, j] + gap              # штраф за пропуск сверху
      left <- F[i, j - 1] + gap            # штраф за пропуск слева
      
      # Выбираем максимальный балл
      max_score <- max(diagonal, up, left)
      F[i, j] <- max_score
      
      # Определяем направление
      if (max_score == diagonal) {
        traceback[i, j] <- "D"  # Диагональ
      } else if (max_score == up) {
        traceback[i, j] <- "U"  # Вверх
      } else {
        traceback[i, j] <- "L"  # Влево
      }
    }
  }
  
  list(F = F, traceback = traceback)
}

# 4. Обратный проход для восстановления оптимального выравнивания
traceback <- function(seq1, seq2, traceback) {
  aligned_seq1 <- ""  # Выравненная первая последовательность
  aligned_seq2 <- ""  # Выравненная вторая последовательность
  
  # Начинаем с правого нижнего угла
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
  
  list(aligned_seq1 = aligned_seq1, aligned_seq2 = aligned_seq2)
}

# 5. Основная функция
needleman_wunsch <- function(seq1, seq2, match = 1, mismatch = -1, gap = -1) {
  # Создаем матрицу замен
  S <- create_substitution_matrix(match, mismatch)
  
  # Инициализируем матрицы
  matrices <- initialize_matrices(seq1, seq2, gap)
  
  # Заполняем матрицы
  filled_matrices <- fill_matrices(seq1, seq2, matrices$F, matrices$traceback, S, gap)
  
  # Восстанавливаем выравнивание
  traceback_result <- traceback(seq1, seq2, filled_matrices$traceback)
  
  list(F = filled_matrices$F, 
       traceback = filled_matrices$traceback,
       aligned_seq1 = traceback_result$aligned_seq1, 
       aligned_seq2 = traceback_result$aligned_seq2)
}

# Визуализация выравнивания
visualize_alignment <- function(aligned_seq1, aligned_seq2) {
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

# Функция для удаления неопознанных символов из последовательности
remove_unidentified_characters <- function(sequence) {
  cleaned_sequence <- gsub("[^ACGT]", "", sequence)
  return(cleaned_sequence)
}

# Пример использования
#seq1 <- "CCCTATATA"
#seq2 <- "TACCCTATATACC"
seq1 <- "NNCCCTNNATATNA"
seq2 <- "TACCCMNNTATNATACCN"
cat("seq1: ", seq1, ", seq2: ", seq2, "\n")

# Удаляем неопознанные символы
seq1 <- remove_unidentified_characters(seq1)
seq2 <- remove_unidentified_characters(seq2)
cat("seq1: ", seq1, ", seq2: ", seq2, "\n")

result <- needleman_wunsch(seq1, seq2)

# Печать матрицы сходства
cat("\nМатрица сходства F:\n")
print(result$F)

# Печать выравнивания
visualize_alignment(result$aligned_seq1, result$aligned_seq2)
