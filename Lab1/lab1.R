# Вариант 2

# Задание 3
help(t)
help(asin)
help(plot)

# Задание 4
rm(list=ls())

# Задание 5
V <- c(1,2,3,4,5) # Создание вектора 
M <- matrix(1:15,nrow=3,ncol=5) # Создание матрицы 

V[3] # - третий элемент вектора V
M[2,3] # - элемент матрицы из второй строки и третьего столбца матрицы M 
M[,3] # - доступ ко всему третьему столбцу матрицы M
M[2,] # - доступ ко всей второй строке матрицы M

# Добавление строки к матрице: 
M <- matrix(1:15,nrow=3,ncol=5) 
S <- 1:5 
M <- rbind(M,S) 
M

# Добавление столбца к матрице: 
M <- matrix(1:15,nrow=3,ncol=5) 
S <- 1:3 
M <- cbind(M,S) 
M

M <- M[,-3] # - Удаление всего третьего столбца матрицы M
M <- M[-2,] # - Удаление всей второй строки матрицы M
M

C <- matrix(1:15,nrow=3,ncol=5); C
apply(C,2,min) # – вектор минимальных элементов, найденных по столбцам
apply(C,1,min) # – вектор минимальных элементов, найденных по строкам 

# Задание 6
C <- matrix(1:25,nrow=5,ncol=5); C
apply(C,2,min) # – вектор минимальных элементов, найденных по столбцам
apply(C,1,min) # – вектор минимальных элементов, найденных по строкам 
apply(C,2,max) # – вектор максимальных элементов, найденных по столбцам
apply(C,1,max) # – вектор максимальных элементов, найденных по строкам 
min(C) # - минимальный элемент по всей матрице
max(C) # - максимальный элемент по всей матрице

# Задание 7
A <- matrix(c(2, 1, 3, 4), nrow = 2); A
B <- c(4, 11); B
x <- solve(A,B); x
A %*% x == B

# Задание 8
x <- seq(0,6.28,0.01)   # Создание вектора аргумента 0<x<6.28 с шагом 0.1 
y <- sin(x)
y2 <- sqrt(abs(sin(x)))
plot(x,y,type = "p",xlab = 'Argument x',ylab = 'Function y')  # Построение графика функции 
lines(x,y2,col="red")
grid() # Отображение сетки на рисунке 
title(main='Function  sin(x)') 

tiff(filename = "plot.tiff",res = 100) 
plot(x,y) 
dev.off() 

# Задание 9
Lx <- -5 # Левая граница для x 
Rx <- 5  # Правая граница для x 
stepx <- 0.05 # Шаг по оси x 
Ly <- -5 # Левая граница для y 
Ry <- 5 # Правая граница для y 
stepy <- 0.05 # Шаг по оси y 

# Создание сетки координат 
xs <- seq(Lx,Rx,stepx) 
ys <- seq(Ly,Ry,stepy) 
z <- outer(xs, ys, function(x, y) 100*(y - x^2)^2 + (1-x)^2) 
persp(xs, ys, z, phi = 30, theta = -30,col = "lightblue") 

# Задание 10-11
# Задание диапазона изменения X 
X_left <- -2 
X_right <- 2 
# Задание диапазона изменения Y 
Y_left <- -3 
Y_right <- 3 
N <- 1000 # Задание количества сгенерированных точек 
source("my_func.R") # Вызов функции 
X <- my_func(X_left, X_right, N) 
Y <- my_func(Y_left, Y_right, N)
plot(X,Y) # Построение графика функции 

# Инициализация параметров для построения гистограммы 
BinNumber <- 20 
k <- 0:BinNumber 
# Вычисление границ карманов на оси X 
X_bins <- X_left + k*(X_right - X_left)/BinNumber 
# Вычисление границ карманов на оси Y 
Y_bins <- Y_left + k*(Y_right - Y_left)/BinNumber 
# Построение гистограмм для X и Y 
hist(X,X_bins) 
hist(Y,Y_bins)

# Задание 12
set.seed(42)
source("my_gauss_gen.R") # Вызов функции 
data <- my_gauss_gen(N = 10000, m = 0, D = 2)
data2 <- my_gauss_gen(N = 10000, m = 5, D = 1)

# Построение гистограммы с 40 каналами
hist(data, 
     breaks = 40, 
     main = "Гистограмма гауссового распределения (m = 0, D = 2)", 
     xlab = "Значения", 
     ylab = "Частота", 
     col = "lightblue", 
     border = "black")

combined_data <- c(data, data2)

hist(data_sum, 
     breaks = 40, 
     main = "Гистограмма суммы гауссовых распределений (m = 5, D = 3)", 
     xlab = "Значения", 
     ylab = "Частота", 
     col = rgb(0.1, 0.2, 0.5, 0.7), # Прозрачный цвет 
     border = "black")

# Задание 13
n <- 5  # Количество строк (можно изменить)
m <- 10 # Количество столбцов (можно изменить)
# Формируем нулевую матрицу Q размером n x m
Q <- matrix(0, nrow = n, ncol = m)

# Заполняем матрицу Q случайными целыми числами
for (k in 1:n) {
  for (j in 1:m) {
    Q[k, j] <-  round(10*runif(1)) # Используем sample для целых чисел от 1 до 30
  }
}
print(Q)

# Задание 14
sum_Q <- 0

# Используем цикл, чтобы вычислить сумму
for (k in 1:n) {
  for (j in 1:m) {
    sum_Q <- sum_Q + Q[k, j]
  }
}
print(paste("Сумма элементов матрицы Q:", sum_Q))

# Задание 15
x <- seq(from = 0, to = 5, by = 0.1)
y <- 2 * x^2 + x - 1
M <- cbind(x, y)

write.table(M, "task15.txt", row.names = FALSE, col.names = TRUE)

# Задание 16
A <- read.table("task15.txt", header = TRUE)
print(A)

plot(A[, 1], A[, 2], xlab = 'X', ylab = 'Y')
title(main = 'F(x)')