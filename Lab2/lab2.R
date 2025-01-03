# Вариант 2

directory <- "C:\\Users\\dapol\\OneDrive\\Документы\\БГУ\\Биоинформатика\\Labs\\Data\\Lab2\\specdata"

# Функция для вычисления среднего значения загрязнителя сульфатами или нитратами 
# для заданного количества источников: 
pollutantmean <- function(directory, pollutant, id) {
  
  # directory – вектор типа character, указывающий путь к папке, содержащей  CSV файлы. 
  # pollutant – вектор типа character, указывающий имя загрязнителя ("sulfate" or "nitrate"), для которого необходимо вычислить среднее значение. 
  # id – вектор типа integer, идентификационные номера файлов данных.  
  
  values <- c() # Создаем пустой вектор для хранения значений загрязнителя
  
  for (i in id) {
    # Формируем имя файла
    filename <- paste(directory, "/", sprintf("%03d", i), ".csv", sep = "") 
    data <- read.csv(filename) # Читаем данные из CSV-файла
    values <- c(values, data[[pollutant]]) # Добавляем значения загрязнителя в вектор
  }
  mean_value <- mean(values, na.rm = TRUE) # Вычисляем среднее значение, игнорируя NA
  return(mean_value) # Возвращаем среднее значение
}

# Вычисляем средний уровень сульфата для файлов с ID от 50 до 99
mean_level_sulfate <- pollutantmean(directory, "sulfate", 50:99)
+
print(mean_level_sulfate) # Выводим средний уровень сульфата

# Функция для подсчета полных случаев в данных:
complete <- function(directory, id) {
  
  # id – вектор типа integer, идентификационные номера файлов данных. 
  
  res <- data.frame(id = integer(), nobs = integer()) # Создаем пустой data frame для результатов
  
  for (i in id) {
    # Формируем имя файла с нулями для индексов
    filename <- paste(directory, "/", sprintf("%03d", i), ".csv", sep = "")
    data <- read.csv(filename) # Читаем данные из CSV-файла
    nobs <- sum(complete.cases(data)) # Подсчитываем количество полных случаев
    res <- rbind(res, data.frame(id = i, nobs = nobs)) # Добавляем результаты в таблицу
  }
  return(res) # Возвращаем таблицу с результатами
}

# Подсчитываем полные случаи для файлов с определенными ID
result <- complete(directory, c(3, 5, 9, 11, 13))
print(result)

# Функция для вычисления коэффициента корреляции Пирсона между сульфатом и нитратом:
corr <- function(directory, threshold, method) {
  
  # threshold – вектор типа numeric, содержащий пороговое значение, при превышении 
  # которого вычисляется коэффициент корреляции Пирсона между переменными "sulfate" 
  # или "nitrate" для каждого из источников данных.
  
  corr_values <- c() # Создаем пустой вектор для хранения значений корреляции
  files <- list.files(directory, full.names = TRUE) # Получаем список файлов в директории
  
  for (file in files) {
    data <- read.csv(file) # Читаем данные из файла
    full_data <- data[complete.cases(data), ] # Оставляем только полные случаи
    
    if (nrow(full_data) > threshold) { # Проверяем, превышает ли количество полных случаев порог
      correlation <- cor(full_data$sulfate, full_data$nitrate, method = method) # Вычисляем корреляцию
      corr_values <- c(corr_values, correlation) # Добавляем значение корреляции в вектор
    }
  }
  return(corr_values) # Возвращаем вектор корреляций
}

# Вычисляем корреляции для файлов, где количество полных случаев превышает 950
corr_result <- corr(directory, threshold = 950, "pearson")
print(corr_result)

library(ggplot2)

plot_histograms <- function(directory, ids) {
  data_list <- list()
  
  for (i in ids) {
    filename <- paste(directory, "/", sprintf("%03d", i), ".csv", sep = "")
    data <- read.csv(filename)
    
    data$file_id <- as.factor(i)
    data_list[[as.character(i)]] <- data
  }
  
  all_data <- do.call(rbind, data_list)
  
  # Построение гистограммы
  ggplot(all_data, aes(x = sulfate, fill = file_id)) +
    geom_histogram(position = "identity", alpha = 0.4, bins = 50) +
    labs(title = "Гистограммы уровней сульфатов для каждого файла", x = "Уровень загрязнения сульфатами", y = "Частота") +
    theme_minimal() +
    scale_fill_manual(values = rainbow(length(ids)))
    #scale_fill_manual(values = c("orange", "yellow", "lightgreen", "lightblue", "lightpink"))
}

plot_histograms(directory, c(3, 5, 9, 11, 13))