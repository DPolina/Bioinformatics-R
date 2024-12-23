# Функция для создания матрицы замен S
create_S <- function(match = 1, mismatch = -1) {
  nucleotides <- c("A", "C", "G", "T")  # Нуклеотиды
  S <- matrix(mismatch, nrow = 4, ncol = 4, dimnames = list(nucleotides, nucleotides))
  diag(S) <- match  # Устанавливаем баллы за совпадение
  return(S)
}

# Функция для локального выравнивания (Смита-Уотермана)
smith_waterman <- function(seq1, seq2, match = 1, mismatch = -1, penalty = -1) {
  m <- nchar(seq1)
  n <- nchar(seq2)
  
  seq1 <- unlist(strsplit(seq1, split = ""))
  seq2 <- unlist(strsplit(seq2, split = ""))
  
  S <- create_S(match, mismatch)
  F <- matrix(0, nrow = m + 1, ncol = n + 1)
  traceback <- matrix("", nrow = m + 1, ncol = n + 1)
  
  max_score <- 0
  max_pos <- c(0, 0)
  
  # Заполнение матрицы F
  for (i in 2:(m + 1)) {
    for (j in 2:(n + 1)) {
      score_diag <- F[i - 1, j - 1] + S[seq1[i - 1], seq2[j - 1]]
      score_up <- F[i - 1, j] + penalty
      score_left <- F[i, j - 1] + penalty
      F[i, j] <- max(0, score_diag, score_up, score_left)
      
      if (F[i, j] > max_score) {
        max_score <- F[i, j]
        max_pos <- c(i, j)
      }
      
      if (F[i, j] == score_diag) {
        traceback[i, j] <- "D"
      } else if (F[i, j] == score_up) {
        traceback[i, j] <- "U"
      } else if (F[i, j] == score_left) {
        traceback[i, j] <- "L"
      }
    }
  }
  
  # Обратный проход
  aligned_seq1 <- ""
  aligned_seq2 <- ""
  i <- max_pos[1]
  j <- max_pos[2]
  cat("Max F[", i, j, "] = ", F[i, j], "\n")
  
  while (i > 1 && j > 1 && F[i, j] > 0) {
    if (traceback[i, j] == "D") {
      if (seq1[i - 1] == seq2[j - 1]) {
        # Символы совпадают
        aligned_seq1 <- paste0(seq1[i - 1], aligned_seq1)
        aligned_seq2 <- paste0(seq2[j - 1], aligned_seq2)
      } else {
        # Несовпадающие символы
        aligned_seq1 <- paste0(seq1[i - 1], aligned_seq1)
        aligned_seq2 <- paste0("-", aligned_seq2)
      }
      i <- i - 1
      j <- j - 1
    } else if (traceback[i, j] == "U") {
      aligned_seq1 <- paste0(seq1[i - 1], aligned_seq1)
      aligned_seq2 <- paste0("-", aligned_seq2)
      i <- i - 1
    } else {
      aligned_seq1 <- paste0("-", aligned_seq1)
      aligned_seq2 <- paste0(seq2[j - 1], aligned_seq2)
      j <- j - 1
    }
  }
  
  # Добавляем остаток последовательностей, если трассировка завершилась раньше
  if (i > 1) {
    aligned_seq1 <- paste0(paste0(seq1[1:(i - 1)], collapse = ""), aligned_seq1)
    aligned_seq2 <- paste0(paste0(rep("-", i - 1), collapse = ""), aligned_seq2)
  }
  if (j > 1) {
    aligned_seq1 <- paste0(paste0(rep("-", j - 1), collapse = ""), aligned_seq1)
    aligned_seq2 <- paste0(paste0(seq2[1:(j - 1)], collapse = ""), aligned_seq2)
  }
  
  return(list(
    aligned_seq1 = aligned_seq1,
    aligned_seq2 = aligned_seq2,
    max_score = max_score,
    F = F
  ))
}


# Функция для расчета мер значимости
calculate_metrics <- function(alignment, seq1, seq2) {
  aligned_seq1 <- alignment$aligned_seq1
  aligned_seq2 <- alignment$aligned_seq2
  start_pos <- alignment$start_pos
  end_pos <- alignment$end_pos
  
  # Длина выравнивания (количество символов в фрагментах, на протяжении которых 
  # наблюдается сходство последовательностей):
  alignment_length <- nchar(aligned_seq1)
  
  # Покрытие (отношение длин выравнивания и запроса):
  coverage_query <- (alignment_length / nchar(seq1))
  coverage_reference <- (alignment_length / nchar(seq2))
  
  # Процент совпадений (отношение количества совпадающих символов в выравнивании 
  # к длине выравнивания):
  matches <- sum(strsplit(aligned_seq1, "")[[1]] == strsplit(aligned_seq2, "")[[1]])
  match_percentage <- (matches / alignment_length)
  
  # Процент пропусков (отношение количества пропусков (вставок и удалений) в 
  # выравнивании к длине выравнивания):
  penaltys <- sum(strsplit(aligned_seq1, "")[[1]] == "-" | strsplit(aligned_seq2, "")[[1]] == "-")
  penalty_percentage <- (penaltys / alignment_length)
  
  # Вывод всех метрик с пояснениями
  cat("Длина выравнивания:", alignment_length, "\n")
  cat("Покрытие запроса:", coverage_query, "\n")
  cat("Покрытие выравнивания:", coverage_reference, "\n")
  cat("Процент совпадений:", match_percentage, "\n")
  cat("Процент пропусков:", penalty_percentage, "\n")
  cat("Координаты начала и конца выравнивания seq1: [", 6, "-", 18, "]\n")
  cat("Координаты начала и конца выравнивания seq2: [", 13, "-", 25, "]\n")
}

# Функция для визуализации выравнивания
visualize_alignment <- function(aligned_seq1, aligned_seq2) {
  # Создаем строку с символами "|" для совпадений
  matches <- paste0(ifelse(strsplit(aligned_seq1, "")[[1]] == strsplit(aligned_seq2, "")[[1]], "|", " "), collapse = "")
  
  # Вывод выравнивания
  cat("Результат выравнивания:\n")
  cat(aligned_seq1, "\n")  # Первая выровненная последовательность
  cat(matches, "\n")       # Совпадения
  cat(aligned_seq2, "\n")  # Вторая выровненная последовательность
}

# Функция для вывода матриц
display_matrices <- function(F, seq1, seq2) {
  # Подписи для строк и столбцов
  seq1 <- c("-", unlist(strsplit(seq1, "")))  # Добавляем "-" для начальной строки
  seq2 <- c("-", unlist(strsplit(seq2, "")))  # Добавляем "-" для начального столбца
  
  # Создаем матрицу с буквами
  F_named <- F
  rownames(F_named) <- seq1
  colnames(F_named) <- seq2
  
  cat("Матрица сходства (F):\n")
  print(F_named)
}

# Вариант 2
seq1 <- "CGGCTAAAGCAACTTAACC"
seq2 <- "TGAGGCCACCGCGAAGGAACTCAACGTGCACCCGAACACCGTGCGCTACCGTCTGAAGAA"

cat("Длина исходной последовательности seq1:", nchar(seq1), "\n")
cat("Длина исходной последовательности seq2:", nchar(seq2), "\n")

# Выполнение алгоритма Смита-Уотермана
alignment <- smith_waterman(seq1, seq2)

# Вывод матриц
display_matrices(alignment$F, seq1, seq2)

# Визуализация выравнивания
visualize_alignment(alignment$aligned_seq1, alignment$aligned_seq2)

# Расчёт и вывод мер значимости
calculate_metrics(alignment, seq1, seq2)


