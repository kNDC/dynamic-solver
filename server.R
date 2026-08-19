## Обрабатывающая составляющая
# https://shiny.posit.co/

library(shiny)
library(bslib)
library(yyjsonr)
library(Matrix)
# library(plotly)

source("utils/SymbolicEngine.R")

function(input, output, session) {
  # Предзаготовка решения системы
  # x_0 определяется на основе размерности данных
  # x_k пересчитывается итеративно Гауссом-Ньютоном
  x <- reactiveVal(numeric());
  x.df <- reactive({
    out <- list();
    
    from <- 1; to <- n.pers() * n.exprs$dials;
    
    out$main <- data.frame(t = core()$data$t.labels);
    out[[1]] <- cbind(out[[1]], matrix(x()[from:to], 
                                       ncol = n.exprs$dials));
    
    names(out[[1]]) <- c("t", unlist(lapply(
      core()$metadata$endogenous,
      function(v) return(v$name)
    )));
    
    from <- from + to; to <- to + n.pers() * n.exprs$nonneg;
    if (from <= to) {
      out$dual <- data.frame(t = core()$data$t.labels);
      out[[2]] <- cbind(out[[2]], matrix(x()[from:to], 
                                         ncol = n.exprs$nonneg));
      
      names(out[[2]]) <- c("t", unlist(lapply(
        core()$metadata$endogenous,
        function(v) {
          if (v$nonneg %||% FALSE) {
            return(paste(v$name, "(дв.)"));
          } else return(NULL);
        }
      )));
    }
    
    return(out);
  });
  
  # Хранилище опорного решения
  x.df.bl <- reactiveVal(); # bl = baseline
  
  # Неупорядоченное хранилище и перечень имён переменных
  vars.dict <- new.env(hash = TRUE, parent = emptyenv());
  vars.names <- new.env(hash = TRUE, parent = emptyenv());
  vars.names$primals <- c();
  vars.names$endogenous <- c();
  vars.names$exogenous <- c();
  
  # Стандартные подблоки (в ошибках и матрице Якоби)
  mx.blocks <- reactive({
    out <- new.env(parent = emptyenv());
    
    out$zMx <- Matrix(0, nrow = n.pers(), ncol = n.pers(), 
                            sparse = TRUE);
    out$uMx <- .sparseDiagonal(n.pers());
    
    out$ddMx <- out$uMx[-1,] - out$uMx[-n.pers(),];
    out$adjMx <- (out$uMx[-1,] + out$uMx[-n.pers(),])/2;
    
    return(out);
  });
  
  # Для отслеживания вызовов
  depth <- 0;
  
  # Болванка вектора ошибок
  # err_k пересчитывается итеративно Гауссом-Ньютоном
  err <- numeric();
  
  # Болванка матрицы Якоби
  # J_k пересчитывается итеративно Гауссом-Ньютоном
  J <- Matrix(0, 1, 1, sparse = TRUE);
  
  # Болванка вектора весов условий системы
  # Настраивается пользователем через ползунки в правом поле
  w <- numeric();
  
  # Вектор весов в условии выхода
  w.metric <- reactive({
    out <- Matrix(
      c(rep(1, n.pers() * n.exprs$dials), 
        rep(0, n.pers() * n.exprs$nonneg), 
        rep(0, n.pers() * n.exprs$neqs)), 
      sparse = TRUE
    );
    
    return(out);
  });
  
  # Болванка хранилища векторов целевых значений в ограничениях системы
  # Настраивается пользователем через указание значений в правом поле
  trg <- new.env(hash = TRUE, parent = emptyenv());
  
  # Динамические списки Id для 
  # управляющих элементов-источников данных
  ids <- reactiveValues(trg = list(), w.pwr = list());
  
  # Основа для последующих расчётов:
  # Подгрузка исходных данных и управляющих параметров
  core <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    metadata <- yyjsonr::read_json_file("./data/metadata.json");
    core.data <- read.table("./data/data.csv", sep = ",", 
                           fill = TRUE, header = TRUE);
    core.data[is.na(core.data)] <- 0;
    names(core.data)[1] <- "t.labels";
    
    return(list(
      metadata = metadata, 
      data = core.data,
      params = list(
        elasticity = 0.1,
        markup = 0.3, 
        storage_cost_factor = 0.3
      )
    ));
  });
  
  # Индекс размерностей
  n.exprs <- new.env(hash = TRUE, parent = emptyenv());
  n.pers <- reactive(nrow(core()$data));
  
  n.rows.err <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    zMx <- mx.blocks()$zMx;
    uMx <- mx.blocks()$uMx;
    
    ddMx <- mx.blocks()$ddMx;
    adjMx <- mx.blocks()$adjMx;
    
    env <- environment();
    
    out <- 0;
    
    for (expr in err.symb()) {
      out <- out + nrow(eval(
        engine$get_instance()$dimension_op(expr), env
        ));
    }
    
    return(out);
  });
  
  # Преднастройка по чтении исходных данных и метаданных
  observeEvent(
    core(), {
      message(paste0(strrep("  ", depth), "Called ", 
                     deparse(sys.call(0)[[1]])));
      depth <<- depth + 1;
      on.exit(depth <<- depth - 1);
      
      # Наполнение хранилищ текущих и целевых значений
      vars.names$exogenous <- names(core()$metadata$exogenous);
      for (name in vars.names$exogenous) {
        vars.dict[[name]] <- core()$data[[name]];
      }
      
      vars.dict$t <- seq_len(n.pers()) - 1;
      vars.names$exogenous <- c(vars.names$exogenous, "t");
      
      vars.names$primals <- names(core()$metadata$endogenous);
      
      duals.names <- vars.names$primals[unlist(lapply(
        core()$metadata$endogenous, 
        function(v) return(v$nonneg %||% FALSE)
      ))];
      
      if (length(duals.names)) {
        duals.names <- paste(duals.names, "dual", sep = ".");
      }
      
      # Наполнение реестра различных размерностей
      n.exprs$dials <- length(core()$metadata$endogenous);
      n.exprs$nonneg <- length(duals.names);
      n.exprs$neqs <- length(core()$metadata$neqs);
      
      for (i in seq_along(vars.names$primals)) {
        vals <- c();
        
        if (core()$metadata$endogenous[[i]]$initialiser %||% FALSE) {
          vals <- core()$data[[vars.names$primals[i]]];
        } else {
          # Мнимые значения для неинициализирующих переменных
          vals <- rep(1e0, n.pers());
        }
        
        vars.dict[[vars.names$primals[i]]] <- vals;
        trg[[vars.names$primals[i]]] <- vals;
      }
      
      if (!is.null(core()$metadata$neqs) && 
          length(core()$metadata$neqs)) {
        duals.names <- c(
          duals.names, 
          paste(names(core()$metadata$neqs), "dual", sep = ".")
          );
      }
      
      if (length(duals.names)) {
        for (name in duals.names) {
          vars.dict[[name]] <- rep(1e0, n.pers());
        }
      }
      
      vars.names$endogenous <- c(vars.names$primals, duals.names);
      
      # Определение начального вектора решений с учётом размерности
      x(do.call(
        c, lapply(vars.names$endogenous, 
                  function(name) return(vars.dict[[name]]))
        ));
      
      # Id динамических элементов страницы
      ids.indices <- unlist(lapply(
        seq_len(n.exprs$dials), 
        function(idx) paste0(idx, "_", seq_len(n.pers()))
        ));
      ids$trg <- paste0("trg_", ids.indices);
      ids$w.pwr <- paste0("pwr_", ids.indices);
      
      # Подготовка весов
      get.w.init();
      
      # Решение системы
      solve.system();
      
      # Сохранение опорного решения
      x.df.bl(x.df());
  });
  
  # Начальное сведение вектора весов
  get.w.init <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    # Показатели степеней весов ограничения
    # Сначала идут тривиальные ограничения вида x = x.trg
    # Пока не подгрузилась пользовательская сторона, источник - 
    # исходные данные для переменных, которые в них содержатся; 
    # для прочих переменных - произвольные значения с нулевыми весами.
    w <<- rep(0, n.pers() * n.exprs$dials);
    from <- 1; to <- n.pers();
    
    for (v in core()$metadata$endogenous) {
      # Ограничение через целевые значения - 
      # назначаются веса пониже
      if (v$initialiser %||% FALSE) w[from:to] <<- rep(1e6, n.pers());
      from <- from + n.pers(); to <- to + n.pers();
    }
    
    w <<- c(w, 1e1^rep(1e1, n.rows.err() - n.pers() * n.exprs$dials));
  };
  
  # Сведение вектора весов
  get.w <- reactive({
    req(length(ids$w.pwr));
    req(!is.null(n.exprs$dials));
    
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    no.nulls <- TRUE;
    
    w[seq_len(n.pers() * n.exprs$dials)] <<- 
      unlist(lapply(
        ids$w.pwr, 
        function(id) {
          if (is.null(id) || is.na(id) || is.null(input[[id]])) {
            no.nulls <<- FALSE;
            return(NA);
          }
          
          return(10^(input[[id]] - 4));
        }
      ));
    
    # NULLs => пользовательская составляющая 
    # ещё не прогрузилась
    req(no.nulls);
    
    names(w) <<- NULL;
    return(TRUE);
  });
  
  # Сведение вектора целевых значений
  get.trg <- reactive({
    req(length(ids$trg));
    
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    no.nulls <- TRUE;
    from <- 1; to <- n.pers();
    
    .ids <- ids$trg;
    
    for (name in vars.names$primals) {
      trg[[name]] <<- unlist(lapply(
        .ids[from:to], 
        function(id) {
          if (is.null(id) || is.na(id) || is.null(input[[id]])) {
            no.nulls <<- FALSE;
            return(NA);
          }
          
          return(input[[id]]);
        }
      ));
      
      from <- from + n.pers(); to <- to + n.pers();
    }
    
    # NULLs => пользовательская составляющая 
    # ещё не прогрузилась
    req(no.nulls);
    return(TRUE);
  });
  
  # Разбор выражений в векторе ошибок -> reactive!
  err.symb <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    out <- lapply(
      seq_along(core()$metadata$endogenous), 
      function(i) {
        expr <- as.expression(call(
          "-", 
          call("$", as.symbol("vars.dict"), 
               as.symbol(vars.names$primals[i])), 
          call("$", as.symbol("trg"), 
               as.symbol(vars.names$primals[i]))
          ));
        
        if (core()$metadata$endogenous[[i]]$nonneg %||% FALSE) {
          expr[[1]] <- call(
            "-", 
            expr[[1]], 
            call("$", as.symbol("vars.dict"), 
                 as.symbol(paste(vars.names$primals[i],
                                 "dual", sep = ".")))
            );
        }
        
        return(expr);
      }
    );
    
    out <- c(out, unlist(lapply(
      seq_along(core()$metadata$endogenous), 
      function(i) {
        if (!core()$metadata$endogenous[[i]]$nonneg %||% FALSE) return(NULL);
        
        v <- call("$", as.symbol("vars.dict"), 
               as.symbol(vars.names$endogenous[[i]]));
        
        v.dual <- call(
          "$", as.symbol("vars.dict"), 
          as.symbol(paste(vars.names$endogenous[[i]], "dual", sep = ".")));
        
        expr <- as.expression(call(
          "-", 
          call(
            "sqrt", 
            call(
              "+", 
              call("+", call("^", v, 2), call("^", v.dual, 2)), 
              eps^2
            )
          ), 
          call("+", v, v.dual)));
        
        return(expr);
      }
    )));
    
    out <- c(out, unlist(lapply(
      core()$metadata$eqs, 
      function(eq) {
        expr <- base::parse(text = eq$formula);
        
        expr <- engine$get_instance()$add_sources(
          expr = expr, vars.dict
        );
        expr <- engine$get_instance()$parse(expr);
        
        return(expr);
      }
    )));
    
    neqs.names <- names(core()$metadata$neqs);
    out <- c(out, unlist(lapply(
      seq_along(core()$metadata$neqs), 
      function(i) {
        neq <- core()$metadata$neqs[[i]]$formula;
        neq.dual <- paste(neqs.names[i], "dual", sep = ".");
        expr <- parse(text = sprintf(
          "%s + %s", neq, neq.dual
        ));
        
        expr <- engine$get_instance()$add_sources(
          expr = expr, vars.dict
        );
        expr <- engine$get_instance()$parse(expr);
        
        return(expr);
      }
    )));
    
    out <- c(out, unlist(lapply(
      seq_along(core()$metadata$neqs), 
      function(i) {
        neq <- core()$metadata$neqs[[i]]$formula;
        neq.dual <- paste(neqs.names[i], "dual", sep = ".");
        expr <- parse(text = sprintf(
          "sqrt((%s)^2 + %s^2) - (%s) - %s", 
          neq, neq.dual, neq, neq.dual
        ));
        
        expr <- engine$get_instance()$add_sources(
          expr = expr, vars.dict
        );
        
        return(expr);
      }
    )));
  });
  
  # Подсчёт вектора ошибок
  calc.err <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    zMx <- mx.blocks()$zMx;
    uMx <- mx.blocks()$uMx;
    
    ddMx <- mx.blocks()$ddMx;
    adjMx <- mx.blocks()$adjMx;
    
    env <- environment();
    err <<- do.call(rbind, lapply(
      err.symb(), function(expr) {
        return(Matrix(eval(expr, envir = env), sparse = TRUE));
      }
    ));
  };
  
  # Разбор выражений в матрице Якоби -> reactive!
  J.symb <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    vars <- lapply(
      vars.names$endogenous,
      function(name) {
        engine$get_instance()$add_sources(
          base::parse(text = name)[[1]], vars.dict
        )
      }
    );
    
    out <- lapply(
      err.symb(),
      function(expr) {
        lapply(
          vars,
          function(v) engine$get_instance()$D(expr, v)
        )
      }
    );
    
    return(out);
  });
  
  # Подсчёт и сборка матрицы Якоби
  calc.J <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    zMx <- mx.blocks()$zMx;
    uMx <- mx.blocks()$uMx;
    
    ddMx <- mx.blocks()$ddMx;
    adjMx <- mx.blocks()$adjMx;
    
    env <- environment();
    
    J <<- do.call(
      rbind, 
      lapply(J.symb(), function(expr_row) {
        do.call(cbind, lapply(
          expr_row,
          function(expr) return(eval(expr, envir = env))
        ))
      })
    );
  };
  
  solve.system <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    x.adj <- x(); # adj = adjusted
    
    from <- 1 + n.pers() * length(core()$metadata$endogenous);
    to <- length(x.adj);
    
    if (from <= to) x.adj[from:to] <- 1e-1;
    x(x.adj);
    
    # Для проверки условия остановки 
    metric <- 0;
    n.steps <- 0;
    
    withProgress(
      message = "Пересчёт решения",
      value = 0, {
        repeat {
          n.steps <- n.steps + 1;
          
          # Разбиение столбца-решения на переменные
          from <- 1; to <- n.pers();
          for (v.name in vars.names$endogenous) {
            vars.dict[[v.name]] <- x()[from:to];
            from <- from + n.pers(); to <- to + n.pers();
          }
          
          # Пересчёт вектора ошибок
          calc.err();
          
          # Пересчёт матрицы Якоби
          calc.J();
          
          # Пересчёт составляющих системы уравнений
          rhs <- t(w * J) %*% err;
          JJ <- t(J) %*% (w * J);
          
          # Предыдущее решение для сопоставления
          x.prev <- x();
          
          # Пересчёт решения
          # x.step <- NULL; lambda.adj <- lambda;
          # while (is.null(x.step)) {
          #   x.step <- tryCatch(
          #     solve(JJ + lambda.adj * .sparseDiagonal(ncol(JJ)), rhs), 
          #     error = function(e) NULL
          #   );
          #   
          #   lambda.adj <- 1.5 * lambda.adj;
          # }
          
          # x(x() - x.step);
          x(x() - solve(JJ, rhs));
          
          # Условие останова: решение перестало меняться
          metric <- max(abs(x() - x.prev));
          
          incProgress(amount = NULL, 
                      detail = paste0(
                        "[", n.steps, "] ",  
                        "Критерий: ", metric)
                      );
          
          # Останов
          if (!is.finite(metric) || 
              metric <= eps || 
              n.steps >= max.steps) break;
        }
    });
    
    # from <- 1 + n.pers() * n.exprs$dials;
    # to <- length(x());
    # 
    # if (from <= to) {
    #   x.adj <- x();
    #   x.adj[from:to] <- pmin(x.adj[from:to], 1e-1);
    #   x(x.adj);
    # }
    
    if (!is.finite(metric) || metric > eps) {
      showNotification(
        ui = "Ошибка поиска решения: алгоритм не сошёлся",
        closeButton = TRUE, 
        type = "error",
        duration = NULL
      );
    } else {
      showNotification(
        ui = "Решение найдено",
        type = "message",
        duration = 3
      );
    }
  };
  
  # Отрисовка заголовка
  output$title <- renderUI(
    h3("Макро-планировщик", 
       tags$small(class = "d-block text-muted", 
                  core()$metadata$name), 
       class = "text-center")
  );
  
  # Отрисовка таблиц
  lapply(
    list(
      list(frame = "main"),
      list(frame = "dual")
    ), 
    function(config) local({
      id <- paste("tables", config$frame, sep = "_");
      
      output[[id]] <- renderTable({
        x.df()[[config$frame]]
      });
    })
  );
  
  # Отрисовка графиков
  lapply(
    list(
      list(frame = "main", colour = "orange"),
      list(frame = "dual", colour = "steelblue")
    ),
    function(config) local({
      src <- config$frame;
      id <- paste("plots", src, sep = "_");
      colour <- config$colour;

      output[[id]] <- renderPlot({
        req(!is.null(x.df()[[src]]), length(x.df()[[src]]) > 1);

        par(mfrow = c(ceiling((ncol(x.df()[[src]]) - 1) / 2), 2),
            mar = c(3, 3, 2, 1));

        for (i in seq_along(x.df()[[src]])[-1]) {
          limits <- range(range(x.df.bl()[[src]][[i]]),
                          range(x.df()[[src]][[i]]));
          limits <- c(limits[1], max(limits[2], 1e2*eps));

          plot(x.df.bl()[[src]]$t, x.df.bl()[[src]][[i]],
               main = names(x.df()[[src]])[i], ylim = limits,
               type = "l", lty = "dashed",
               lwd = 1.5, col = colour);

          lines(x.df()[[src]]$t, x.df()[[src]][[i]],
                ylim = limits,
                type = "l", lty = "solid",
                lwd = 2.0, col = colour);
        }

        par(mfrow = c(1, 1));
      },

      res = 96,

      height = function() {
        if (is.null(x.df()[[src]]) ||
            length(x.df()[[src]]) <= 1) return(1);
        width <- session$clientData[[paste("output", id,
                                           "width", sep = "_")]] / 2;
        return(4 * width / 5 * ceiling((ncol(x.df()[[src]]) - 1) / 2));
      });
    })
  );
  
  # lapply(
  #   list(
  #     list(frame = "main", colour = "orange"),
  #     list(frame = "dual", colour = "steelblue")
  #   ),
  #   function(config) local({
  #     src <- config$frame;
  #     id <- paste("plots", src, sep = "_");
  #     colour <- config$colour;
  # 
  #     output[[id]] <- renderPlotly({
  #       req(!is.null(x.df()[[src]]), length(x.df()[[src]]) > 1);
  # 
  #       plots <- lapply(
  #         seq_along(x.df()[[src]])[-1],
  #         function(i) {
  #           plot_ly() |>
  #             add_lines(
  #               x = x.df.bl()[[src]]$t, y = x.df.bl()[[src]][[i]],
  #               line = list(dash = "dash", color = "orange", width = 1.5),
  #               name = "Опорное решение"#, showlegend = (i == 2)
  #             ) |>
  #             add_lines(
  #               x = x.df()[[src]]$t, y = x.df()[[src]][[i]],
  #               line = list(color = "orange", width = 2.0),
  #               name = "Текущее решение"#, showlegend = (i == 2)
  #             ) |>
  #             layout(
  #               annotations = list(
  #                 list(text = names(x.df()[[src]])[i],
  #                      xref = "paper", yref = "paper",
  #                      x = 0.5, y = 1.05, showarrow = FALSE)
  #               )
  #             );
  #         });
  # 
  #       subplot(plots, nrows = ceiling(length(plots) / 2),
  #               shareX = TRUE, titleY = TRUE) |>
  #         layout(height = 350 * ceiling(length(plots) / 2));
  #     });
  #   })
  # );
  
  output$dials <- renderUI({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    panels <- lapply(seq_len(n.exprs$dials), function(i) {
      accordion_panel(
        title = names(isolate(x.df())$main)[i + 1], 
        
        fluidRow(
          column(2, p(strong("Срок"))), 
          column(3, p(strong("Цель"))), 
          column(3, p(strong("Решение"))), 
          column(4, p(strong("Важность")))
        ),
        
        lapply(
          seq_len(nrow(isolate(x.df())$main)), 
          function(t) {
            fluidRow(
              column(2, {
                tags$span(isolate(x.df())$main[t, 1], 
                          id = paste("period", i, t, sep = "_"))
              }), 
              column(3, {
                item <- numericInput(
                  paste("trg", i, t, sep = "_"), 
                  label = NULL, 
                  value = round(isolate(x.df())$main[t, i + 1], 3)
                );
                
                item$children[[2]] <- tagAppendAttributes(
                  item$children[[2]], 
                  style = "text-align: right;"
                );
                
                item
              }), 
              column(3, {
                item <- numericInput(
                  paste("sln", i, t, sep = "_"), 
                  label = NULL, 
                  value = round(isolate(x.df())$main[t, i + 1], 3)
                );
                
                item$children[[2]] <- tagAppendAttributes(
                  item$children[[2]], 
                  readonly = "readonly", 
                  style = "text-align: right;"
                );
                
                item
              }), 
              column(4, sliderInput(paste("pwr", i, t, sep = "_"), 
                                    label = NULL, 0, 10, 5, 
                                    ticks = FALSE)), 
              class = "compact-row"
            )
          })
        )
      });
    
    return(do.call(accordion, c(list(open = FALSE), panels)));
  });
  
  # Отражение изменений решения на странице
  observe({
    req(x.df());
    for (i in seq_len(n.exprs$dials)) {
      for (t in seq_len(nrow(x.df()$main))) {
        updateNumericInput(session = session, 
                           inputId = paste0("sln_", i, "_", t), 
                           value = round(x.df()$main[t, i + 1], 3));
      }
    }
  });
  
  observeEvent(input$btnReset, {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    # Пересчёт решения происходит автоматически
    for (id in ids$trg) {
      id.parts <- strsplit(id, split = "_")[[1]];
      i <- as.integer(id.parts[2]);
      t <- as.integer(id.parts[3]);
      
      updateNumericInput(session = session, 
                         inputId = id, 
                         value = round(x.df.bl()$main[t, i + 1], 3));
      updateSliderInput(session = session, 
                        inputId = paste("pwr", id.parts[2], id.parts[3], sep = "_"), 
                        value = 5);
    }
  });
  
  observeEvent(input$btnShowAll, {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    message("Бабах!");
  });
  
  solver.launcher <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    req(get.w(), get.trg());
    return(TRUE);
  }) |> debounce(250);
  
  observeEvent(solver.launcher(), {
    message(paste0(strrep("  ", depth), "Called ",
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);

    solve.system();
  }, ignoreInit = TRUE);
};
