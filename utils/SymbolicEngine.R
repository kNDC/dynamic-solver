library(R6)

SymbolicEngine <- R6Class(
  classname = "SymbolicEngine", 
  cloneable = F, 
  
  public = list(
    zMx = as.symbol("zMx"), 
    uMx = as.symbol("uMx"), 
    
    ddMx = as.symbol("ddMx"), 
    adjMx = as.symbol("adjMx"), 
    
    initialize = function() {
      private$adjMx0 <- private$process_indices(
        call("[", self$uMx, call("-", quote(`T`)), quote(expr = ))
      );
      
      private$adjMx1 <- private$process_indices(
        call("[", self$uMx, call("-", 1), quote(expr = ))
      );
    },
    
    add_sources = function(expr, ...) {
      if (is.recursive(expr)) {
        return(as.expression(lapply(
          expr, 
          function(node) private$append_source(node, ...)
        )));
      }
      
      return(private$append_source(expr, ...));
    },
    
    parse = function(expr) {
      if (!is.recursive(expr)) return(expr);
      
      if (is.call(expr)) {
        if (length(expr) == 1) return(expr);
        
        arg <- self$parse(expr[[2]]);
        
        # Добавленные операторы
        if (expr[[1]] == quote(dd)) {
          return(call("%*%", self$ddMx, arg));
        } else if (expr[[1]] == quote(adj0) || 
                   expr[[1]] == quote(ip)) {
          return(call("%*%", private$adjMx0, arg));
        } else if (expr[[1]] == quote(adj1) || 
                   expr[[1]] == quote(fp)) {
          return(call("%*%", private$adjMx1, arg));
        } else if (expr[[1]] == quote(adj) || 
                   expr[[1]] == quote(mp)) {
          return(call("%*%", self$adjMx, arg));
        }
        
        if (expr[[1]] == quote(ln)) {
          return(call("log", arg))
        } else if (expr[[1]] == quote(lg)) {
          return(call("log", arg, 10))
        }
        
        if (expr[[1]] == quote(tg)) {
          return(call("tan", arg))
        } else if (expr[[1]] == quote(ctg)) {
          return(call("/", 1, call("tan", arg)))
        }
        
        if (expr[[1]] == quote(arctg)) {
          return(call("atan", arg))
        } else if (expr[[1]] == quote(arcctg)) {
          return(call("atan", call("/", 1, arg)))
        }
        
        # Нормализация индексов
        if (expr[[1]] == quote(`[`) || expr[[1]] == quote(`[[`)) {
          return(private$process_indices(expr));
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
        
        # Индексирующие записи
        # A[i.rows, i.cols] -> [, A, i.rows, i.cols
        if (op.name == "[" || op.name == "[[") {
          index.op <- call(op.name, self$uMx, expr[[3]], quote(expr = ));
          
          return(call("%*%", private$process_indices(index.op), 
                      .D(expr[[2]], var, flat_matching)
          ));
        }
        
        # Матричные операции
        # A %*% x -> %*%, A, x
        if (op.name == "%*%") {
          return(call("%*%", expr[[2]], 
                      .D(expr[[3]], var, flat_matching)));
        }
        
        # Алгебраические операции
        # x op y -> op, x, y
        if (op.name == "+") {
          if (length(expr) == 2) {
            return(.D(expr[[2]], var, flat_matching));
          } else if (length(expr) == 3) {
            return(call("+", 
              .D(expr[[2]], var, flat_matching), 
              .D(expr[[3]], var, flat_matching))
            );
          }
        } else if (op.name == "-") {
          if (length(expr) == 2) {
            return(call("-", .D(expr[[2]], var, flat_matching)));
          } else if (length(expr) == 3) {
            return(call("-", 
              .D(expr[[2]], var, flat_matching), 
              .D(expr[[3]], var, flat_matching))
            );
          }
        } else if (op.name == "*") {
          u = expr[[2]]; v = expr[[3]];
          du = .D(expr[[2]], var, flat_matching);
          dv = .D(expr[[3]], var, flat_matching);
          
          return(call(
            "+", 
            call("%*%", private$diag(v), du), 
            call("%*%", private$diag(u), dv)
          ));
        } else if (op.name == "/") {
          u = expr[[2]]; v = expr[[3]];
          du = .D(expr[[2]], var, flat_matching);
          dv = .D(expr[[3]], var, flat_matching);
          
          return(call(
            "-", 
            call("%*%", private$diag(call("/", 1, v)), du), 
            call("%*%", private$diag(call("/", u, call("^", v, 2))), dv)
          ));
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
          return(.D(call("^", expr[[2]], 0.5), var, flat_matching));
        } else if (op.name == "exp") {
          return(call("%*%", private$diag(expr), 
                      .D(expr[[2]], var, flat_matching)));
        }
        
        # Введённые дополнительные операторы
        if (op.name == "dd") {
          return(call("%*%", self$ddMx, 
                      .D(expr[[2]], var, flat_matching)));
        } else if (op.name == "adj" || op.name == "mp") {
          return(call("%*%", self$adjMx, 
                      .D(expr[[2]], var, flat_matching)));
        } else if (op.name == "adj0" || op.name == "ip") {
          return(call("%*%", private$adjMx0, 
                      .D(expr[[2]], var, flat_matching)));
        } else if (op.name == "adj1" || op.name == "fp") {
          return(call("%*%", private$adjMx1, 
                      .D(expr[[2]], var, flat_matching)));
        }
        
        # Производные элементарных функций
        if (identical(expr[[2]], var)) {
          # Логарифмы
          if (op.name == "log") { # Логарифм с произвольным основанием
            x <- expr[[2]]; a <- expr[[3]];
            return(private$diag(call("/", call("/", 1, x), 
                                     call("log", a))));
          }
          
          if (op.name == "ln") { # Естественный логарифм
            x <- expr[[2]];
            return(private$diag(call("/", 1, x)));
          }
          
          if (op.name == "lg") { # Десятичный логарифм
            x <- expr[[2]];
            return(private$diag(call("/", call("/", 1, x), 
                                     call("log", 10))));
          }
          
          # Тригонометрические
          if (op.name == "sin") {
            x <- expr[[2]];
            return(private$diag(call("cos", x)));
          }
          
          if (op.name == "cos") {
            x <- expr[[2]];
            return(call("-", private$diag(call("sin", x))));
          }
          
          if (op.name == "tg") {
            x <- expr[[2]];
            return(private$diag(call("+", 1, call("^", call("tan", x), 2))));
          }
          
          if (op.name == "ctg") {
            x <- expr[[2]];
            return(call(
              "-", 
              private$diag(call(
                "+", 1, call("/", 1, call("^", call("tan", x), 2)))
            )));
          }
          
          if (op.name == "arctg" || op.name == "atan") {
            x <- expr[[2]];
            return(private$diag(call("/", 1, call("+", 1, call("^", x, 2)))));
          }
          
          if (op.name == "arcctg") {
            x <- expr[[2]];
            return(private$diag(call("/", call("-", 1), 
                                     call("+", 1, call("^", x, 2)))));
          }
        }
        
        # Производная сложной функции 
        # Общий случай нескольких переменных
        summands <- lapply(
          seq_along(expr)[-1], 
          function(i) {
            return(call(
              "%*%", 
              .D(expr, expr[[i]], flat_matching), 
              .D(expr[[i]], var, flat_matching)
            ));
          }
        );
        
        return(Reduce(f = function(x1, x2) call("+", x1, x2), summands));
      };
      
      # Если переменная - строка, перейти к символам
      if (is.character(var)) var <- parse(text = var)[[1]];
      
      if (is.call(expr)) {
        out <- private$simplify(.D(expr, var, flat_matching));
      } else if (is.recursive(expr)) {
        out <- private$simplify(as.expression(lapply(
          expr, function(expr) .D(expr, var, flat_matching))
        ));
      } else return(expr);
      
      if (is.symbol(out) && out == self$zMx) {
        dimension_op <- self$dimension_op(expr);
        if (dimension_op != self$uMx) {
          out <- call("%*%", self$dimension_op(expr), out);
        }
      } else if (is.expression(out) && 
                 is.symbol(out[[1]]) && out[[1]] == self$zMx) {
        dimension_op <- self$dimension_op(expr);
        if (dimension_op != self$uMx) {
          out[[1]] <- call("%*%", self$dimension_op(expr), out[[1]]);
        }
      }
      
      return(out);
    },
    
    dimension_op = function(expr) {
      mult_chains <- list();
      
      .dimension_op <- function(node, chain) {
        if (is.symbol(node)) {
          mult_chains[[length(mult_chains) + 1]] <<- chain;
          return();
        }
        
        if (is.expression(node)) {
          for (i in seq_along(node)) {
            .dimension_op(node[[i]], chain);
          }
        }
        
        if (!is.call(node)) return();
        
        if (node[[1]] == quote(dd)) {
          return(.dimension_op(node[[2]], 
            call("%*%", chain, self$ddMx)
          ));
        } else if (node[[1]] == quote(adj) || 
                   node[[1]] == quote(mp)) {
          return(.dimension_op(node[[2]], 
            call("%*%", chain, self$adjMx)
          ));
        } else if (node[[1]] == quote(adj0) || 
                   node[[1]] == quote(ip)) {
          index.op <- call("[", self$uMx, 
                           call("-", quote(`T`)), quote(expr = ));
          return(.dimension_op(node[[2]], 
            call("%*%", chain, index)
          ));
        } else if (node[[1]] == quote(adj1) || 
                   node[[1]] == quote(fp)) {
          return(.dimension_op(node[[2]], 
                               call("%*%", chain, self$adjMx)
          ));
        }
        
        if (node[[1]] == quote(`[`)) {
          index.op <- call("[", self$uMx, node[[3]], quote(expr = ));
          mult_chains[[length(mult_chains) + 1]] <<- 
            call("%*%", chain, private$process_indices(index.op));
          return();
        }
        
        if (node[[1]] == quote(`%*%`)) {
          return(.dimension_op(
            node[[3]], 
            call("%*%", chain, node[[2]])
          ));
        }
        
        if (length(node) >= 2) {
          for (i in 2:length(node)) {
            .dimension_op(node[[i]], chain);
          }
        }
      }
      
      .dimension_op(expr, self$uMx);
      
      mult_chains <- lapply(
        mult_chains, private$simplify
        );
      
      # Так как все переменные одной длины, то 
      # выбирается просто самая короткая цепочка
      i_min <- 1;
      for (i in seq_along(mult_chains)) {
        if (is.symbol(mult_chains[[i]]) && 
            mult_chains[[i]] == self$uMx) {
          return(self$uMx);
        } else if (length(mult_chains[[i]]) <= 
                   length(mult_chains[[i_min]])) {
          i_min <- i;
        }
      }
      
      if (!length(mult_chains)) return(NULL);
      return(mult_chains[[i_min]]);
    }
  ), 
  
  private = list(
    diag_op = ".sparseDiagonal",
    adjMx0 = NULL,
    adjMx1 = NULL,
    
    # Диагонализация множителей для матрицы Якоби
    diag = function(expr) {
      return(bquote(.(as.symbol(private$diag_op))(x = .(expr))));
    },
    
    # Извлекает название источника для переменной
    find_source = function(node, ...) {
      name <- as.character(node);
      env.names <- as.character(substitute(...()));
      
      i <- 1;
      for (env in list(...)) {
        if (exists(name, envir = env, inherits = F)) return(env.names[i]);
        i <- i + 1;
      }
      
      return(NULL);
    },
    
    append_source = function(node, ...) {
      # Дальнейшие уровни вложенности узлов выражения
      if (is.recursive(node)) {
        if (is.call(node)) {
          return(as.call(lapply(
            node, 
            function(sub_node) private$append_source(sub_node, ...)
          )));
        }
        
        return(node);
      }
      
      # Случай если узел - не переменная
      if (!is.symbol(node)) return(node);
      
      src.name <- private$find_source(node, ...);
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
        
        # zMx %*% A -> zMx; uMx %*% A -> A; -uMx %*% A -> -A
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
          } else if (is.call(expr[[2]]) && length(expr[[2]]) == 2 && 
                     expr[[2]][[1]] == quote(`-`) && 
                     is.symbol(expr[[2]][[2]]) && 
                     expr[[2]][[2]] == self$uMx) {
            once_more <<- T;
            return(call("-", expr[[3]]));
          } else if (is.call(expr[[3]]) && length(expr[[3]]) == 2 && 
                     expr[[3]][[1]] == quote(`-`) && 
                     is.symbol(expr[[3]][[2]]) && 
                     expr[[3]][[2]] == self$uMx) {
            once_more <<- T;
            return(call("-", expr[[2]]));
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
    },
    
    # x[1, , drop = FALSE] -> `[`, x, 1, quote(expr = ), FALSE
    process_indices = function(expr) {
      var <- expr[[2]];
      if (length(expr) == 2) return(var);
      
      if(length(expr) > 3) {
        T.subst <- call("nrow", var);
      } else {
        T.subst <- call("length", var);
      }
      
      # Передаются только вызовы => нет выражений
      .replace_Ts <- function(node) {
        if (is.call(node)) {
          return(as.call(lapply(node, .replace_Ts)));
        }
        
        if (is.symbol(node) && node == quote(`T`)) return(T.subst);
        return(node);
      };
      
      .group_inds <- function(node) {
        .extract_set <- function(node) {
          if (is.call(node) && node[[1]] == quote(`{`)) {
            return(as.call(c(quote(c), as.list(node[-1]))));
          }
          
          return(node);
        }
        
        return(as.call(lapply(node, function(node) {
          if (is.call(node) && node[[1]] == quote(`-`) && 
              length(node) == 2) {
            return(call("-", .extract_set(node[[2]])));
          } else {
            return(.extract_set(node));
          }
        })));
      }
      
      out <- .replace_Ts(expr);
      out <- .group_inds(out);
      if (is.null(out$drop)) out$drop <- FALSE;
      return(out);
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
