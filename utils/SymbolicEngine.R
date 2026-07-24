library(R6)

SymbolicEngine <- R6Class(
  classname = "SymbolicEngine", 
  cloneable = F, 
  
  public = list(
    zMx = as.symbol("zMx"), 
    uMx = as.symbol("uMx"), 
    sub_uMx = as.symbol("sub_uMx"),
    sup_uMx = as.symbol("sup_uMx"),
    
    add_sources = function(expr, sources) {
      if (is.expression(expr) || is.list(expr)) {
        return(as.expression(lapply(
          expr, 
          function(node) private$append_source(node, sources)
        )));
      }
      
      return(private$append_source(expr, sources));
    },
    
    parse = function(expr) {
      if (!is.recursive(expr)) return(expr);
      
      if (is.call(expr)) {
        if (length(expr) == 1) return(expr);
        
        arg <- self$parse(expr[[2]]);
        if (expr[[1]] == quote(dd)) {
          return(bquote((.(self$uMx) - .(self$sub_uMx)) %*% .(arg)));
        } else if (expr[[1]] == quote(shift)) {
          return(bquote(.(self$sub_uMx) %*% .(arg)));
        }
        
        return(as.call(lapply(expr, self$parse)));
      }
      
      return(as.expression(lapply(expr, self$parse)));
    },
    
    D = function(expr, var, flat_matching = F) {
      .D <- function(expr, var, flat_matching) {
        # Производная самого себя или постоянной величины
        if (!is.recursive(expr)) {
          if (is.symbol(expr) && expr == var) {
            return(self$uMx);
          }
          
          return(self$zMx);
        }
        
        if (!is.call(expr)) return(NA);
        
        op <- expr[[1]];
        op.name <- as.character(op);
        
        # Случай выражения в скобках
        if (op.name == "(") {
          return(bquote((.(.D(expr[[2]], var, flat_matching)))));
        }
        
        # Случай элементов data.frame, list и т.п.
        # src$x -> $, src, x
        if (op.name == "$") {
          if (identical(expr, var)) {
            return(self$uMx);
          } else {
            if (flat_matching && is.symbol(var)) {
              while (is.call(expr) && expr[[1]] == quote(`$`)) {
                expr <- expr[[3]];
              }
              
              return(if (expr == var) self$uMx else self$zMx);
            }
            return(self$zMx);
          }
        }
        
        # Матричные операции
        # A %*% x -> %*%, A, x
        if (op.name == "%*%") {
          symbs <- list(
            A = expr[[2]], 
            dx = .D(expr[[3]], var, flat_matching)
          );
          return(substitute(A %*% dx, symbs));
        }
        
        # Алгебраические операции
        # x op y -> op, x, y
        if (op.name == "+") {
          if (length(expr) == 2) {
            return(.D(expr[[2]], var, flat_matching));
          } else if (length(expr) == 3) {
            return(bquote(
              .(.D(expr[[2]], var, flat_matching)) + 
                .(.D(expr[[3]], var, flat_matching))
            ));
          }
        } else if (op.name == "-") {
          if (length(expr) == 2) {
            return(bquote(-.(.D(expr[[2]], var, flat_matching))));
          } else if (length(expr) == 3) {
            return(bquote(
              .(.D(expr[[2]], var, flat_matching)) - 
                .(.D(expr[[3]], var, flat_matching))
            ));
          }
        } else if (op.name == "*") {
          u = expr[[2]]; v = expr[[3]];
          du = .D(expr[[2]], var, flat_matching);
          dv = .D(expr[[3]], var, flat_matching);
          
          return(bquote(
            .(private$diag(v)) %*% (.(du)) + 
              .(private$diag(u)) %*% (.(dv))
          ));
        } else if (op.name == "/") {
          u = expr[[2]]; v = expr[[3]];
          du = .D(expr[[2]], var, flat_matching);
          dv = .D(expr[[3]], var, flat_matching);
        }
        
        if (op.name == "^") {
          u = expr[[2]]; v = expr[[3]];
          du = .D(expr[[2]], var, flat_matching);
          dv = .D(expr[[3]], var, flat_matching);
          
          return(bquote(
            (.(private$diag(bquote(.(v) * .(u)^(.(v) - 1)))) %*% (.(du)) + 
               .(private$diag(bquote(.(u)^.(v) * log(.(u))))) %*% (.(dv)))
          ));
        } else if (op.name == "sqrt") {
          return(.D(bquote(.(expr[[2]])^0.5), var, flat_matching));
        }
        
        # Введённые дополнительные операторы
        if (op.name == "dd") {
          return(bquote(
            (.(self$uMx) - .(self$sub_uMx)) %*% 
              (.(.D(expr[[2]], var, flat_matching)))
          ));
        } else if (op.name == "shift") {
          return(bquote(
            .(self$sub_uMx) %*% 
              (.(.D(expr[[2]], var, flat_matching)))
          ));
        }
        
        # Производные элементарных функций
        if (identical(expr[[2]], var)) {
          # Логарифмы
          if (op.name == "log") { # Логарифм с произвольным основанием
            x <- expr[[2]]; a <- expr[[3]];
            return(bquote(
              .(private$diag(1/x/log(a)))
            ));
          }
          
          if (op.name == "ln") { # Естественный логарифм
            return(.D(bquote(log(.(expr[[2]]), exp(1))), var, flat_matching));
          }
          
          if (op.name == "lg") { # Десятичный логарифм
            return(.D(bquote(log(.(expr[[2]]), 10)), var, flat_matching));
          }
          
          # Тригонометрические
          if (op.name == "sin") {
            x <- expr[[2]];
            return(bquote(
              .(private$diag(bquote(cos(.(x)))))
            ));
          }
          
          if (op.name == "cos") {
            x <- expr[[2]];
            return(bquote(
              -.(private$diag(bquote(sin(.(x)))))
            ));
          }
          
          if (op.name == "tg") {
            x <- expr[[2]];
            return(bquote(
              .(private$diag(1 + bquote(tan(.(x))^2)))
            ));
          }
          
          if (op.name == "ctg") {
            x <- expr[[2]];
            return(bquote(
              -.(private$diag(1 + 1/tan(x)^2))
            ));
          }
        }
        
        # Производная сложной функции 
        # Общий случай нескольких переменных
        summands <- lapply(
          seq_along(expr)[-1], 
          function(i) {
            df.dx <- .D(expr, expr[[i]], flat_matching);
            dx.dvar <- .D(expr[[i]], var, flat_matching);
            return(bquote(
              (.(df.dx)) %*% (.(dx.dvar))
            ));
          }
        );
        
        return(Reduce(f = function(x1, x2) call("+", x1, x2), summands));
      };
      
      # Если переменная - строка, перейти к символам
      if (is.character(var)) var <- parse(text = var)[[1]];
      
      if (is.recursive(expr) && !is.call(expr)) {
        return(private$simplify(as.expression(lapply(
          expr, function(expr) .D(expr, var, flat_matching))
        )));
      }
      
      if (is.call(expr)) {
        return(private$simplify(.D(expr, var, flat_matching)));
      }
      
      return(expr);
    }
  ), 
  
  private = list(
    diag_op = ".sparseDiagonal",
    
    # Диагонализация множителей для матрицы Якоби
    diag = function(expr) {
      
      return(bquote(.(as.symbol(private$diag_op))(x = .(expr))));
    },
    
    # Извлекает название источника для переменной
    find_source = function(node, sources) {
      name <- as.character(node);
      
      for (i in seq_along(sources)) {
        if (name %in% sources[[i]]) return(names(sources)[[i]]);
      }
      
      return(NULL);
    },
    
    append_source = function(node, sources) {
      # Дальнейшие уровни вложенности узлов выражения
      if (is.recursive(node)) {
        if (is.call(node)) {
          return(as.call(lapply(
            node, 
            function(sub_node) private$append_source(sub_node, sources)
          )));
        }
        
        return(node);
      }
      
      # Случай если узел - не переменная
      if (!is.symbol(node)) return(node);
      
      src.name <- private$find_source(node, sources);
      if (is.null(src.name)) return(node);
      src.symb <- parse(text = src.name)[[1]];
      
      symbs <- list(src = src.symb, var = node);
      return(substitute(src$var, symbs));
    },
    
    simplify = function(expr) {
      .simplify <- function(expr) {
        if (!is.recursive(expr)) return(expr);
        
        for (i in seq_along(expr)) {
          expr[[i]] <- .simplify(expr[[i]]);
        }
        
        op <- expr[[1]];
        op.name <- as.character(op)[[1]];
        
        # x + 0 -> x
        if (op.name == "+") {
          if ((is.integer(expr[[2]]) || is.numeric(expr[[2]])) && 
              (is.integer(expr[[3]]) || is.numeric(expr[[3]]))) {
            once_more <<- T;
            return(expr[[2]] + expr[[3]]);
          } else if (is.symbol(expr[[2]]) && expr[[2]] == self$zMx) {
            once_more <<- T;
            return(expr[[3]]);
          } else if (is.symbol(expr[[3]]) && expr[[3]] == self$zMx) {
            once_more <<- T;
            return(expr[[2]]);
          }
        }
        
        # x - 0 -> x
        if ((op.name == "-")) {
          # Упрощение арифметики
          if (length(expr) == 3) {
            if ((is.integer(expr[[2]]) || is.numeric(expr[[2]])) && 
                (is.integer(expr[[3]]) || is.numeric(expr[[3]]))) {
              once_more <<- T;
              return(expr[[2]] - expr[[3]]);
            } else if (is.symbol(expr[[2]]) && expr[[2]] == self$zMx) {
              once_more <<- T;
              return(bquote(-.(expr[[3]])));
            } else if (is.symbol(expr[[3]]) && expr[[3]] == self$zMx) {
              once_more <<- T;
              return(expr[[2]]);
            }
          } else if (length(expr) == 2) {
            if (is.symbol(expr[[2]]) && expr[[2]] == self$zMx) {
              return(expr[[2]]);
            }
          }
        }
        
        # 0 * x -> 0; 1 * x -> x
        if (op.name == "*") {
          # Упрощение арифметики
          if ((is.integer(expr[[2]]) || is.numeric(expr[[2]])) && 
              (is.integer(expr[[3]]) || is.numeric(expr[[3]]))) {
            once_more <<- T;
            return(expr[[2]] * expr[[3]]);
          } else if (((is.integer(expr[[2]]) || is.numeric(expr[[2]])) && 
                      (expr[[2]] == 0)) || 
                     ((is.integer(expr[[3]]) || is.numeric(expr[[3]])) && 
                      (expr[[3]] == 0))) {
            once_more <<- T;
            return(0);
          } else if (((is.integer(expr[[2]]) || is.numeric(expr[[2]])) && 
                      (expr[[2]] == 1))) {
            once_more <<- T;
            return(expr[[3]]);
          } else if (((is.integer(expr[[3]]) || is.numeric(expr[[3]])) && 
                      (expr[[3]] == 1))) {
            once_more <<- T;
            return(expr[[2]]);
          }
        }
        
        # x^1 -> x
        if (op.name == "^") {
          # Упрощение арифметики
          if ((is.integer(expr[[2]]) || is.numeric(expr[[2]])) && 
              (is.integer(expr[[3]]) || is.numeric(expr[[3]]))) {
            once_more <<- T;
            return(expr[[2]] ^ expr[[3]]);
          } else if ((is.integer(expr[[3]]) || is.numeric(expr[[3]])) && 
                     expr[[3]] == 1) {
            once_more <<- T;
            return(expr[[2]]);
          }
        }
        
        # zMx %*% A -> zMx; uMx %*% A -> A
        # diag(x) %*% diag(y) -> diag(x * y)
        if (op.name == "%*%") {
          if ((is.symbol(expr[[2]]) && expr[[2]] == self$zMx) || 
              (is.symbol(expr[[3]]) && expr[[3]] == self$zMx)) {
            once_more <<- T;
            return(self$zMx);
          } else if (is.symbol(expr[[2]]) && expr[[2]] == self$uMx) {
            once_more <<- T;
            return(expr[[3]]);
          } else if (is.symbol(expr[[3]]) && expr[[3]] == self$uMx) {
            once_more <<- T;
            return(expr[[2]]);
          }
          
          if ((is.call(expr[[2]]) && (expr[[2]][[1]] == private$diag_op)) && 
              (is.call(expr[[3]]) && (expr[[3]][[1]] == private$diag_op))) {
            once_more <<- T;
            return(bquote(
              .(private$diag(bquote(.(expr[[2]][[2]]) * .(expr[[3]][[2]]))))
            ));
          }
        }
        
        # (x) -> x; ((x)) -> (x)
        if (op.name == "(") {
          if (length(expr) == 2) {
            if (is.integer(expr[[2]]) || is.numeric(expr[[2]]) || 
                is.symbol(expr[[2]])) {
              once_more <<- T;
              return(expr[[2]]);
            }
          }
          
          if (is.call(expr[[2]]) && expr[[2]][[1]] == quote(`(`)) {
            once_more <<- T;
            return(expr[[2]]);
          }
        }
        
        # diag((x)) -> diag(x)
        if (is.call(op) && as.character(op[[1]]) == private$diag_op) {
          if (is.call(op[[2]]) && as.character(op[[2]][[1]]) == "(") {
            once_more <<- T;
            return(bquote(
              .(private$diag(op[[2]][[2]]))
            ));
          }
        }
        
        return(expr);
      };
      
      once_more <- F;
      
      if (is.recursive(expr) && !is.call(expr)) {
        repeat {
          expr <- as.expression(lapply(expr, .simplify));
          if (!once_more) break;
          once_more <- F;
        }
      } else if (is.call(expr)) {
        repeat {
          expr <- .simplify(expr);
          if (!once_more) break;
          once_more <- F;
        }
      }
      
      return(expr);
    }
  )
);

engine <- new.env(parent = emptyenv())
engine$get_instance <- function() {
  # Опасно в многопоточном окружении!!!
  if (is.null(engine$instance)) {
    engine$instance <- SymbolicEngine$new();
  }
  
  return(engine$instance);
};
