# Статические функции и переменные

#--- Auxiliary functions ---
#--- Оператор опережения/запаздывания ---
shift <- function(x, n_periods) {
  if (n_periods == 0) return(x);
  
  if (n_periods > 0) {
    out <- c(x[(n_periods + 1):length(x)], rep(NA, n_periods));
    return(out);
  }
  
  if (n_periods < 0) {
    out <- c(rep(NA, -n_periods), x[1:(length(x) + n_periods)]);
    return(out);
  }
}

#--- One-sided Hodrick-Prescott filter ---
hp1 <- function(x, lambda = 1600L) {
  dx.init <- x[2] - x[1];
  
  return (sapply(seq_along(x), 
                 function(t) {
                   return (hp2(c(x[1] - 2*dx.init, 
                                 x[1] - dx.init, 
                                 x[1:t]), lambda)[t]);
                 }));
}

#--- Two-sided Hodrick-Prescott filter ---
hp2 <- function(x, lambda = 1600L) {
  m <- length(x);
  E <- diag(nrow = m, ncol = m);
  
  D <- cbind(0L, diag(nrow = m - 2L, ncol = m - 2L)) - 
    cbind(diag(nrow = m - 2L, ncol = m - 2L), 0L);
  D <- cbind(0L, D) - cbind(D, 0L);
  
  return(solve(E + lambda * t(D) %*% D, x));
}

#--- Расчётные параметры ---
eps <- .Machine$double.eps^0.5; # допустимая погрешность
max.steps <- 200; # ограничение на число итераций
lambda <- eps; # регулирующий параметр Левенберга-Марквардта
