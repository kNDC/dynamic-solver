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

#--- Calculation parameters ---
eps <- .Machine$double.eps^0.5; # допустимая погрешность
max.steps <- 200; # ограничение на число итераций
lambda <- eps; # регулирующий параметр Левенберга-Марквардта

#--- Метаданные модели ---
write.metadata.eggs <- function() {
  if (!exists("metadata")) metadata <- list();
  if (is.null(metadata$Id)) metadata$Id <- jsonlite::unbox(ids::uuid());
  metadata$name <- jsonlite::unbox("Яйца");
  
  # Настройка переменных
  metadata$vars <- list();
  
  metadata$vars$Y <- list();
  metadata$vars$Y$name <- jsonlite::unbox("Выпуск");
  metadata$vars$Y$formula <- jsonlite::unbox("data()");
  metadata$vars$Y$nonneg <- jsonlite::unbox(T);
  metadata$vars$Y$dial <- jsonlite::unbox(F);
  metadata$vars$Y$unit <- jsonlite::unbox("млн ед.");
  
  metadata$vars$S.in <- list();
  metadata$vars$S.in$name <- jsonlite::unbox("Запасы, начальные");
  metadata$vars$S.in$formula <- jsonlite::unbox("dd(S.in) - S.in.base - shift((1 - own_use - loss) * Y, -1) + shift(C + Ex - Im, -1)");
  metadata$vars$S.in$nonneg <- jsonlite::unbox(T);
  metadata$vars$S.in$dial <- jsonlite::unbox(F);
  metadata$vars$S.in$unit <- jsonlite::unbox("млн ед.");
  
  metadata$vars$S.out <- list();
  metadata$vars$S.out$name <- jsonlite::unbox("Запасы, конечные");
  metadata$vars$S.out$formula <- jsonlite::unbox("S.out - S.in - (1 - own_use - loss) * Y + C + Ex - Im");
  metadata$vars$S.out$nonneg <- jsonlite::unbox(T);
  metadata$vars$S.out$dial <- jsonlite::unbox(F);
  metadata$vars$S.out$unit <- jsonlite::unbox("млн ед.");
  
  metadata$vars$C <- list();
  metadata$vars$C$name <- jsonlite::unbox("Потребление");
  metadata$vars$C$formula <- jsonlite::unbox("data()");
  metadata$vars$C$nonneg <- jsonlite::unbox(T);
  metadata$vars$C$dial <- jsonlite::unbox(F);
  metadata$vars$C$unit <- jsonlite::unbox("млн ед.");
  
  metadata$vars$Ex <- list();
  metadata$vars$Ex$name <- jsonlite::unbox("Вывоз");
  metadata$vars$Ex$formula <- jsonlite::unbox("data()");
  metadata$vars$Ex$nonneg <- jsonlite::unbox(T);
  metadata$vars$Ex$dial <- jsonlite::unbox(F);
  metadata$vars$Ex$unit <- jsonlite::unbox("млн ед.");
  
  metadata$vars$Im <- list();
  metadata$vars$Im$name <- jsonlite::unbox("Привоз");
  metadata$vars$Im$formula <- jsonlite::unbox("data()");
  metadata$vars$Im$nonneg <- jsonlite::unbox(T);
  metadata$vars$Im$dial <- jsonlite::unbox(F);
  metadata$vars$Im$unit <- jsonlite::unbox("млн ед.");
  
  # Настройка параметров
  metadata$pars <- list();
  
  jsonlite::write_json(x = metadata, 
                       path = "./data/metadata.json", 
                       pretty = T);
};

#--- Метаданные модели ---
write.metadata.milk <- function() {
  if (!exists("metadata")) metadata <- list();
  if (is.null(metadata$Id)) metadata$Id <- jsonlite::unbox(ids::uuid());
  metadata$name <- jsonlite::unbox("Молоко");
  
  # Настройка переменных
  metadata$vars <- list();
  
  metadata$vars$Ls1 <- list();
  metadata$vars$Ls1$name <- jsonlite::unbox("Поголовье, СХО");
  metadata$vars$Ls1$formula <- jsonlite::unbox("dd(Ls1 - Ls1.base) - I1 - I1.base");
  metadata$vars$Ls1$nonneg <- jsonlite::unbox(T);
  metadata$vars$Ls1$dial <- jsonlite::unbox(F);
  metadata$vars$Ls1$unit <- jsonlite::unbox("тыс. голов");
  
  metadata$vars$I1 <- list();
  metadata$vars$I1$name <- jsonlite::unbox("Новые скотоместа, СХО");
  metadata$vars$I1$formula <- jsonlite::unbox("data()");
  metadata$vars$I1$nonneg <- jsonlite::unbox(T);
  metadata$vars$I1$dial <- jsonlite::unbox(F);
  metadata$vars$I1$unit <- jsonlite::unbox("тыс.");
  
  metadata$vars$Ls2 <- list();
  metadata$vars$Ls2$name <- jsonlite::unbox("Поголовье, вне СХО");
  metadata$vars$Ls2$formula <- jsonlite::unbox("dd(Ls2 - Ls2.base) - I2 - I2.base");
  metadata$vars$Ls2$nonneg <- jsonlite::unbox(T);
  metadata$vars$Ls2$dial <- jsonlite::unbox(F);
  metadata$vars$Ls2$unit <- jsonlite::unbox("тыс. голов");
  
  metadata$vars$I2 <- list();
  metadata$vars$I2$name <- jsonlite::unbox("Новые скотоместа, вне СХО");
  metadata$vars$I2$formula <- jsonlite::unbox("data()");
  metadata$vars$I2$nonneg <- jsonlite::unbox(T);
  metadata$vars$I2$dial <- jsonlite::unbox(F);
  metadata$vars$I2$unit <- jsonlite::unbox("тыс.");
  
  metadata$vars$S.in <- list();
  metadata$vars$S.in$name <- jsonlite::unbox("Запасы, начальные");
  metadata$vars$S.in$formula <- jsonlite::unbox("dd(S.in) - S.in.base - shift(q1 * (1 - own_use1 - loss1) * Ls1, -1) - shift(q2 * (1 - own_use2 - loss2) * Ls2, -1) + shift(C + Ex - Im, -1)");
  metadata$vars$S.in$nonneg <- jsonlite::unbox(T);
  metadata$vars$S.in$dial <- jsonlite::unbox(F);
  metadata$vars$S.in$unit <- jsonlite::unbox("млн тонн");
  
  metadata$vars$S.out <- list();
  metadata$vars$S.out$name <- jsonlite::unbox("Запасы, конечные");
  metadata$vars$S.out$formula <- jsonlite::unbox("S.out - S.in - q1 * (1 - own_use1 - loss1) * Ls1 - q2 * (1 - own_use2 - loss2) * Ls2) + C + Ex - Im");
  metadata$vars$S.out$nonneg <- jsonlite::unbox(T);
  metadata$vars$S.out$dial <- jsonlite::unbox(F);
  metadata$vars$S.out$unit <- jsonlite::unbox("млн тонн");
  
  metadata$vars$C <- list();
  metadata$vars$C$name <- jsonlite::unbox("Потребление и переработка");
  metadata$vars$C$formula <- jsonlite::unbox("data()");
  metadata$vars$C$nonneg <- jsonlite::unbox(T);
  metadata$vars$C$dial <- jsonlite::unbox(F);
  metadata$vars$C$unit <- jsonlite::unbox("млн тонн");
  
  metadata$vars$Ex <- list();
  metadata$vars$Ex$name <- jsonlite::unbox("Вывоз");
  metadata$vars$Ex$formula <- jsonlite::unbox("data()");
  metadata$vars$Ex$nonneg <- jsonlite::unbox(T);
  metadata$vars$Ex$dial <- jsonlite::unbox(F);
  metadata$vars$Ex$unit <- jsonlite::unbox("млн тонн");
  
  metadata$vars$Im <- list();
  metadata$vars$Im$name <- jsonlite::unbox("Привоз");
  metadata$vars$Im$formula <- jsonlite::unbox("data()");
  metadata$vars$Im$nonneg <- jsonlite::unbox(T);
  metadata$vars$Im$dial <- jsonlite::unbox(F);
  metadata$vars$Im$unit <- jsonlite::unbox("тыс.");
  
  # Настройка параметров
  metadata$pars <- list();
  
  jsonlite::write_json(x = metadata, 
                       path = "./data/metadata.json", 
                       pretty = T);
};
