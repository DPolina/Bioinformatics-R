my_gauss_gen <- function(N, m, D) {
  # Генерация случайных величин
  random_numbers <- rnorm(N, mean = m, sd = sqrt(D))
  
  # Возвращаем сгенерированные величины
  return(random_numbers)
}