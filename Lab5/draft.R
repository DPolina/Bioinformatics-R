# Функция для создания матрицы замен S
create_S <- function(match = 1, mismatch = -1) {
  nucleotides <- c("A", "C", "G", "T")  # Нуклеотиды
  S <- matrix(mismatch, nrow = 4, ncol = 4, dimnames = list(nucleotides, nucleotides))
  diag(S) <- match  # Устанавливаем баллы за совпадение
  return(S)
}

# Функция для локального выравнивания (Смита-Уотермана) с несколькими результатами
smith_waterman_multiple <- function(seq1, seq2, match = 1, mismatch = -1, penalty = -1, top_n = 5) {
  m <- nchar(seq1)
  n <- nchar(seq2)
  
  seq1 <- unlist(strsplit(seq1, split = ""))
  seq2 <- unlist(strsplit(seq2, split = ""))
  
  S <- create_S(match, mismatch)
  F <- matrix(0, nrow = m + 1, ncol = n + 1)
  traceback <- matrix("", nrow = m + 1, ncol = n + 1)
  
  results <- list()
  
  for (k in 1:top_n) {
    max_score <- 0
    max_pos <- c(0, 0)
    
    # Заполнение матрицы F
    for (i in 2:(m + 1)) {
      for (j in 2:(n + 1)) {
        if (F[i, j] < 0) next  # Пропускаем уже обработанные ячейки
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
    if (max_score == 0) break  # Если больше нет значимых выравниваний
    aligned_seq1 <- ""
    aligned_seq2 <- ""
    i <- max_pos[1]
    j <- max_pos[2]
    end_pos <- c(i - 1, j - 1)
    
    while (i > 1 && j > 1 && F[i, j] > 0) {
      if (traceback[i, j] == "D") {
        aligned_seq1 <- paste0(seq1[i - 1], aligned_seq1)
        aligned_seq2 <- paste0(seq2[j - 1], aligned_seq2)
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
    start_pos <- c(i - 1, j - 1)
    
    results[[k]] <- list(
      aligned_seq1 = aligned_seq1,
      aligned_seq2 = aligned_seq2,
      max_score = max_score,
      start_pos = start_pos,
      end_pos = end_pos
    )
    
    # Обнуление текущего выравнивания в матрице F
    i <- max_pos[1]
    j <- max_pos[2]
    while (i > 1 && j > 1 && F[i, j] > 0) {
      F[i, j] <- -1  # Помечаем клетки как обработанные
      if (traceback[i, j] == "D") {
        i <- i - 1
        j <- j - 1
      } else if (traceback[i, j] == "U") {
        i <- i - 1
      } else {
        j <- j - 1
      }
    }
  }
  
  return(results)
}

# Функция для расчета метрик
calculate_metrics <- function(alignment, seq1, seq2) {
  aligned_seq1 <- alignment$aligned_seq1
  aligned_seq2 <- alignment$aligned_seq2
  start_pos <- alignment$start_pos
  end_pos <- alignment$end_pos
  
  # Длина выравнивания
  alignment_length <- nchar(aligned_seq1)
  
  # Покрытие
  coverage_query <- (alignment_length / nchar(seq1))
  coverage_reference <- (alignment_length / nchar(seq2))
  
  # Процент совпадений
  matches <- sum(strsplit(aligned_seq1, "")[[1]] == strsplit(aligned_seq2, "")[[1]])
  match_percentage <- (matches / alignment_length)
  
  # Процент пропусков
  gaps <- sum(strsplit(aligned_seq1, "")[[1]] == "-" | strsplit(aligned_seq2, "")[[1]] == "-")
  gap_percentage <- (gaps / alignment_length)
  
  # Вывод метрик
  cat("Длина выравнивания:", alignment_length, "\n")
  cat("Покрытие запроса:", coverage_query, "\n")
  cat("Покрытие выравнивания:", coverage_reference, "\n")
  cat("Процент совпадений:", match_percentage, "\n")
  cat("Процент пропусков:", gap_percentage, "\n")
  cat("Координаты начала и конца выравнивания seq1: [", start_pos[1], "-", end_pos[1], "]\n")
  cat("Координаты начала и конца выравнивания seq2: [", start_pos[2], "-", end_pos[2], "]\n")
}

# Функция для визуализации выравнивания
visualize_alignment <- function(aligned_seq1, aligned_seq2) {
  # Создаем строку с символами "|" для совпадений, "." для несовпадений
  matches <- strsplit(aligned_seq1, "")[[1]] == strsplit(aligned_seq2, "")[[1]]
  match_line <- paste0(ifelse(matches, "|", " "), collapse = "")
  
  # Вывод выравнивания
  cat("Результат выравнивания:\n")
  cat(aligned_seq1, "\n")  # Первая выровненная последовательность
  cat(match_line, "\n")    # Совпадения / несовпадения
  cat(aligned_seq2, "\n")  # Вторая выровненная последовательность
}

# Пример использования
seq1 <- "CGGCTAAAGCAACTTAACC"
seq2 <- "TGAGGCCACCGCGAAGGAACTCAACGTGCACCCGAACACCGTGCGCTACCGTCTGAAGAA"

alignments <- smith_waterman_multiple(seq1, seq2, top_n = 5)

# Вывод результатов
for (i in seq_along(alignments)) {
  cat("\n=== Выравнивание", i, "===\n")
  visualize_alignment(alignments[[i]]$aligned_seq1, alignments[[i]]$aligned_seq2)
  calculate_metrics(alignments[[i]], seq1, seq2)
}
