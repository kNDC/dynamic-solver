## Серверная сторона работы макропланировщика

library(shiny)
library(bslib)
library(Matrix)
# library(plotly)

source("utils/SymbolicEngine.R")

# Обрабатывающая составляющая
function(input, output, session) {
  # Предзаготовка решения системы
  # x_0 определяется на основе размерности данных
  # x_k пересчитывается итеративно Гауссом-Ньютоном
  x <- reactiveVal(numeric());
  x.df <- reactive({
    out <- list(
      main = lapply(seq_len(n.vars() + 1), function(v) numeric()), 
      dual = lapply(seq_len(n.neqs() + 1), function(v) numeric())
    );
    from <- 1; to <- n.pers();
    
    out[[1]][[1]] <- core()$data$t;
    main.pos <- 2;
    
    out[[2]][[1]] <- core()$data$t;
    dual.pos <- 2;
    
    for (v in core()$metadata$vars) {
      out[[1]][[main.pos]] <- x()[from:to];
      from <- from + n.pers(); to <- to + n.pers();
      main.pos <- main.pos + 1;
      
      if (v$nonneg) {
        out[[2]][[dual.pos]] <- x()[from:to];
        from <- from + n.pers(); to <- to + n.pers();
        dual.pos <- dual.pos + 1;
      }
    }
    
    out <- lapply(out, as.data.frame);
    
    names(out[[1]]) <- c("t", unlist(lapply(
      core()$metadata$vars,
      function(v) return(v$name)
    )));
    
    names(out[[2]]) <- c("t", unlist(lapply(
      core()$metadata$vars,
      function(v) return(if (v$nonneg) paste(v$name, "(дв.)") else NULL)
    )));
    
    return(out);
  });
  
  # Хранилище опорного решения
  x.df.bl <- reactiveVal(data.frame()); # bl = baseline
  
  # Стандартные подблоки (в ошибках и матрице Якоби)
  mx.blocks <- reactiveValues(
    zMx = Matrix(0, nrow = 1, ncol = 1, 
                 sparse = T), 
    uMx = Matrix(0, nrow = 1, ncol = 1, 
                 sparse = T), 
    sub_uMx = Matrix(0, nrow = 1, ncol = 1, 
                     sparse = T)
  );
  
  # Для отслеживания вызовов
  depth <- 0;
  
  # Болванка вектора ошибок
  # err_k пересчитывается итеративно Гауссом-Ньютоном
  err <- numeric();
  err.bp <- list(); # bp = blueprint
  
  # Болванка вектора весов условий системы
  # Настраивается пользователем через ползунки в правом поле
  w <- numeric();
  
  # Вектор весов в условии выхода
  w.metric <- reactive({
    out <- Matrix(
      rep(1, n.pers() * (n.vars() + n.neqs())), 
      sparse = T
    );
    
    from <- 1; to <- n.pers();
    
    for (v in core()$metadata$vars) {
      from <- from + n.pers(); to <- to + n.pers();
      
      if (v$nonneg) {
        out[from:to, 1] <- rep(0, n.pers());
        from <- from + n.pers(); to <- to + n.pers();
      }
    }
    
    return(out);
  });
  
  # Болванка вектора целевых значений в ограничениях системы
  # Настраивается пользователем через указание значений в правом поле
  trg <- list();
  
  # Болванка матрицы Якоби
  # J_k пересчитывается итеративно Гауссом-Ньютоном
  J <- Matrix(0, 1, 1, sparse = T);
  J.bp <- list();
  
  # Динамические списки Id для 
  # управляющих элементов-источников данных
  ids <- reactiveValues(
    trg = list(),
    w.pwr = list()
  );
  
  # Словарь переменных и параметров
  vars.index <- list();
  
  # Основа для последующих расчётов:
  # Подгрузка исходных данных и управляющих параметров
  core <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    metadata <- jsonlite::read_json("./data/metadata.json");
    core.data <- read.table("./data/data.csv", sep = ",", 
                           fill = T, header = T);
    core.data[is.na(core.data)] <- 0;
    names(core.data)[1] <- "t";
    
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
  
  # Размерности набора переменных
  n.vars <- reactive(length(core()$metadata$vars));
  n.pers <- reactive(nrow(core()$data));
  
  # Количество нетривиальных уравнений
  n.eqs <- reactive({
    out <- 0;
    for (v in core()$metadata$vars) {
      out <- out + (v$formula != "data()");
    }
    
    return(out);
  });
  
  # Количество условий неотрицательности
  n.neqs <- reactive({
    out <- 0;
    for (v in core()$metadata$vars) {
      out <- out + v$nonneg;
    }
    
    return(out);
  });
  
  # Начальные значения составляющих задачи
  observeEvent(
    core(), {
      message(paste0(strrep("  ", depth), "Called ", 
                     deparse(sys.call(0)[[1]])));
      depth <<- depth + 1;
      on.exit(depth <<- depth - 1);
      
      # Определение начального вектора решений с учётом размерности
      x(rep(rep(1, n.pers()), n.vars() + n.neqs()));
      
      # Определение строительных блоков с учётом размерности
      mx.blocks$zMx <- Matrix(0, nrow = n.pers(), ncol = n.pers(), 
                              sparse = T);
      mx.blocks$uMx <- .sparseDiagonal(n.pers());
      mx.blocks$sub_uMx <- bandSparse(n = n.pers(), k = -1, 
                                      diagonals = list(rep(1, n.pers() - 1)));
      
      # Определение Id динамических элементов страницы
      ids.indices <- unlist(lapply(seq_len(n.vars()), 
                               function(idx) {
                                 return(paste0(idx, "_", seq_len(n.pers())));
                               }));
      ids$trg <- paste0("trg_", ids.indices);
      ids$w.pwr <- paste0("pwr_", ids.indices);
      
      # Словарь переменных и параметров
      vars.names <- unlist(lapply(
        seq_along(core()$metadata$vars), 
        function(i) {
          out <- names(core()$metadata$vars)[[i]];
          if (core()$metadata$vars[[i]]$nonneg) {
            out[[2]] <- paste0(names(core()$metadata$vars)[[i]], ".dual")
          }
          
          return(out);
        }
      ));
      
      pars.names <- unlist(lapply(
        names(core()$data)[-1],
        function(name) {
          return(if (name %in% vars.names) NULL else name);
        }
      ));
      
      vars.index <<- list(
        "vars" = vars.names, 
        "core()$data" = pars.names
      );
      
      # Разбор выражений в метаданных
      parse.err.bp(UI = F);
      parse.J.bp();
      
      # Подготовка весов
      get.w.init();
      
      # Подготовка целевых значений
      get.trg.init();
      
      # Решение системы
      solve.system();
      
      # Разбор выражений в метаданных
      parse.err.bp();
      parse.J.bp();
      
      # Сохранение опорного решения
      x.df.bl(x.df());
  });
  
  # Начальное сведение вектора весов
  get.w.init <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    
    # Показатели степеней весов ограничения
    pwr <- rep(1e1, n.pers() * (n.vars() + n.neqs()));
    from <- 1; to <- n.pers();
    
    for (i in seq_along(core()$metadata$vars)) {
      # Ограничение через целевые значения - 
      # назначаются веса пониже
      if (core()$metadata$vars[[i]]$formula == "data()") {
        pwr[from:to] <- rep(1, n.pers());
      }
      
      from <- from + n.pers(); to <- to + n.pers();
      
      # Неотрицательность - остаются высокие веса
      if (core()$metadata$vars[[i]]$nonneg) {
        from <- from + n.pers(); to <- to + n.pers();
      }
    }
    
    names(pwr) <- NULL;
    
    w <<- 10^pwr;
  };
  
  # Сведение вектора весов
  get.w <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    
    no.nulls <- T;
    
    # Показатели степеней весов ограничения
    pwr <- rep(1e1, n.pers() * (n.eqs() + n.neqs() + n.vars()));
    from <- 1; to <- n.pers();
    
    .ids <- ids$w.pwr;
    ids.from <- 1; ids.to <- n.pers();
    
    for (i in seq_along(core()$metadata$vars)) {
      # Высокие веса для соответствующего уравнения системы
      if (core()$metadata$vars[[i]]$formula != "data()") {
        from <- from + n.pers(); to <- to + n.pers();
      }
      
      pwr[from:to] <- unlist(lapply(
        .ids[ids.from:ids.to], 
        function(id) {
          if (is.null(id) || is.na(id) || is.null(input[[id]])) {
            no.nulls <<- F;
            return(NA);
          }
          
          return(input[[id]] - 4);
        }
      ));
      
      from <- from + n.pers(); to <- to + n.pers();
      ids.from <- ids.from + n.pers();
      ids.to <- ids.to + n.pers();
      
      # Высокие веса для условия неотрицательности
      if (core()$metadata$vars[[i]]$nonneg) {
        from <- from + n.pers(); to <- to + n.pers();
      }
    }
    
    # NULLs => пользовательская составляющая 
    # ещё не прогрузилась
    shiny::req(no.nulls);
    
    names(pwr) <- NULL;
    
    w <<- 10^pwr;
  });
  
  # Начальное сведение вектора целевых значений
  get.trg.init <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    
    trg <<- lapply(seq_along(core()$metadata$vars), function(i) {});
    names(trg) <<- names(core()$metadata$vars);
    
    for (i in seq_along(core()$metadata$vars)) {
      if (core()$metadata$vars[[i]]$formula == "data()") {
        trg[[i]] <<- core()$data[[names(core()$metadata$vars)[i]]];
      }
    }
  };
  
  # Сведение вектора целевых значений
  get.trg <- reactive({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    
    no.nulls <- T;
    from <- 1; to <- n.pers();
    
    .ids <- ids$trg;
    
    for (i in seq_along(core()$metadata$vars)) {
      trg[[i]] <<- unlist(lapply(
        .ids[from:to], 
        function(id) {
          if (is.null(id) || is.na(id) || is.null(input[[id]])) {
            no.nulls <<- F;
            return(NA);
          }
          
          return(input[[id]]);
        }
      ));
      
      names(trg)[[i]] <<- names(core()$metadata$vars)[[i]];
      from <- from + n.pers(); to <- to + n.pers();
    }
    
    shiny::req(no.nulls);
  });
  
  # Разбор выражений в векторе ошибок
  parse.err.bp <- function(UI = T) {
    err.bp <<- unlist(lapply(
      seq_along(core()$metadata$vars), 
      function(i) {
        out <- list();
        # Условие, порождаемое уровнением
        if (core()$metadata$vars[[i]]$formula != "data()") {
          expr <- base::parse(text = core()$metadata$vars[[i]]$formula);
          
          expr <- engine$get_instance()$add_sources(
            expr = expr, 
            sources = vars.index
          );
          
          out[1] <- engine$get_instance()$parse(expr);
        }
        
        # Тривиальное управляющее условие
        if (UI || core()$metadata$vars[[i]]$formula == "data()") {
          out[length(out) + 1] <- base::parse(text = sprintf(
            "%s - %s",
            paste0("vars$", names(core()$metadata$vars)[i]), 
            paste0("trg$", names(core()$metadata$vars)[i])
          ))
        }
        
        # Условие, порождаемое требованием неотрицательности
        if (core()$metadata$vars[[i]]$nonneg) {
          v <- parse(text = names(core()$metadata$vars)[i])[[1]];
          v.sc <- parse(text = paste(names(core()$metadata$vars)[i][1], 
                                     "dual", sep = "."))[[1]];
          expr <- as.expression(bquote(
            sqrt(.(v)^2 + .(v.sc)^2) - .(v) - .(v.sc)
          ));
          
          expr <- engine$get_instance()$add_sources(
            expr = expr, 
            sources = vars.index
          );
          
          out[length(out) + 1] <- expr;
        }
        
        return(out);
      }
    ), recursive = F);
  };
  
  # Подсчёт вектора ошибок
  calc.err <- function(vars) {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    zMx <- mx.blocks$zMx;
    uMx <- mx.blocks$uMx;
    sub_uMx <- mx.blocks$sub_uMx;
    
    env <- environment();
    
    err <<- do.call(rbind, lapply(
      err.bp, function(expr) Matrix(eval(expr, envir = env), sparse = T)
    ));
  };
  
  # Разбор выражений в матрице Якоби
  parse.J.bp <- function() {
    vars <- lapply(
      vars.index$vars,
      function(name) {
        engine$get_instance()$add_sources(
          base::parse(text = name)[[1]], vars.index
        )
      }
    );
    
    J.bp <<- lapply(
      err.bp,
      function(expr) {
        lapply(
          vars,
          function(v) engine$get_instance()$D(expr, v)
        )
      }
    );
  };
  
  # Подсчёт и сборка матрицы Якоби
  calc.J <- function(vars) {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    zMx <- mx.blocks$zMx;
    uMx <- mx.blocks$uMx;
    sub_uMx <- mx.blocks$sub_uMx;
    
    env <- environment();
    
    J <<- do.call(
      rbind, 
      lapply(J.bp, function(expr_row) {
        do.call(cbind, lapply(
          expr_row, function(expr) eval(expr, envir = env)
        ))
      })
    );
  };
  
  solve.system <- function() {
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    # Возмущения для нулевых значений и двойственных переменных
    x.adj <- x(); # adj = adjusted
    from <- 1; to <- n.pers();
    
    for (v in core()$metadata$vars) {
      x.adj[from:to] <- pmax(x.adj[from:to], 1e2 * eps);
      from <- from + n.pers(); to <- to + n.pers();
      
      if (v$nonneg) {
        x.adj[from:to] <- rep(1e-1, n.pers());
        from <- from + n.pers(); to <- to + n.pers();
      }
    }
    
    from <- 1; to <- n.pers();
    x(x.adj);
    
    # Для проверки условия остановки 
    metric <- 0;
    n.steps <- 0;
    
    # Для хранения разбиения столбца-решения на переменные
    vars <- lapply(vars.index$vars, 
                   function(v) numeric());
    names(vars) <- vars.index$vars;
    
    withProgress(
      message = "Пересчёт решения",
      value = 0, {
        repeat {
          n.steps <- n.steps + 1;
          
          # Разбиение столбца-решения на переменные
          from <- 1; to <- n.pers();
          for (i in seq_along(vars.index$vars)) {
            vars[[i]] <- x()[from:to];
            from <- from + n.pers();
            to <- to + n.pers();
          }
          
          # Пересчёт вектора ошибок
          calc.err(vars = vars);
          
          # Пересчёт матрицы Якоби
          calc.J(vars = vars);
          
          # Пересчёт составляющих системы уравнений
          rhs <- t(w * J) %*% err;
          JJ <- t(J) %*% (w * J);
          
          # Предыдущее решение для сопоставления
          x.prev <- x();
          
          # Пересчёт решения
          x.step <- NULL; lambda.adj <- lambda;
          while (is.null(x.step)) {
            x.step <- tryCatch(
              solve(JJ + lambda.adj * .sparseDiagonal(ncol(JJ)), rhs), 
              error = function(e) NULL
            );
            
            lambda.adj <- 1e1 * lambda.adj;
          }
          
          x(x() - x.step);
          
          # Критерий останова: решение перестало меняться
          metric <- max(abs(x() - x.prev) * w.metric());
          
          incProgress(amount = NULL, 
                      detail = paste0(
                        "[", n.steps, "] ",  
                        "Критерий: ", metric)
                      );
          
          if (metric <= eps || n.steps >= max.steps) break;
        }
    });
    
    if (is.na(metric) || metric > eps) {
      showNotification(
        ui = "Ошибка поиска решения: алгоритм не сошёлся",
        closeButton = T, 
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
    
    # Ограничение нетривиальных значений сопряжённых переменных
    x.adj <- x();
    from <- 1; to <- n.pers();
    
    for (v in core()$metadata$vars) {
      from <- from + n.pers(); to <- to + n.pers();
      
      if (v$nonneg) {
        x.adj[from:to] <- pmin(x.adj[from:to], 1e-1);
        from <- from + n.pers(); to <- to + n.pers();
      }
    }
    
    x(x.adj);
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
      id <- paste("tables", config$frame, sep = ".");
      
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
      id <- paste("plots", src, sep = ".");
      colour <- config$colour;
      
      output[[id]] <- renderPlot({
        par(mfrow = c(ceiling(ncol(x.df()$main) / 2), 2), 
            mar = c(3, 3, 2, 1));
        
        for (i in seq_along(x.df()[[src]])[-1]) {
          limits <- range(range(x.df.bl()[[src]][[i]]), 
                          range(x.df()[[src]][[i]]));
          limits <- c(limits[1], max(limits[2], 1e2*eps));
          
          plot(x.df.bl()[[src]]$t, x.df.bl()[[src]][[i]], 
               main = names(x.df()[[src]])[i], ylim = limits, 
               type = "l", lty = "dashed", 
               lwd = 1.5, col = colour);
          
          lines(x.df()$main$t, x.df()$main[[i]], 
                ylim = limits, 
                type = "l", lty = "solid", 
                lwd = 2.0, col = colour);
        }
        
        par(mfrow = c(1, 1));
      }, 
      
      res = 96,
      
      height = function() {
        if (ncol(x.df()[[src]]) == 1) return(0);
        width <- session$clientData[[paste("output", id, "width", sep = "_")]] / 2;
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
  #     id <- paste("plots", src, sep = ".");
  #     colour <- config$colour;
  #     
  #     output[[id]] <- renderPlotly({
  #       shiny::req(x.df.bl(), nrow(x.df.bl()[[src]] > 0));
  #       
  #       plots <- lapply(
  #         seq_along(x.df()[[src]])[-1], 
  #         function(i) {
  #           plot_ly() |> 
  #             add_lines(
  #               x = x.df.bl()[[src]]$t, y = x.df.bl()[[src]][[i]], 
  #               line = list(dash = "dash", color = "orange", width = 1.5), 
  #               name = "Опорное решение", showlegend = (i == 2)
  #             ) |> 
  #             add_lines(
  #               x = x.df()[[src]]$t, y = x.df()[[src]][[i]], 
  #               line = list(color = "orange", width = 2.0), 
  #               name = "Текущее решение", showlegend = (i == 2)
  #             ) |> 
  #             layout(
  #               annotations = list(
  #                 list(text = names(x.df()[[src]])[i], 
  #                      xref = "paper", yref = "paper", 
  #                      x = 0.5, y = 1.05, showarrow = F)
  #               )
  #             );
  #         });
  #       
  #       subplot(plots, nrows = ceiling(length(plots) / 2), 
  #               shareX = T, titleY = T) |> 
  #         layout(height = 350 * ceiling(length(plots) / 2));
  #     });
  #   })
  # );
  
  output$dials <- renderUI({
    message(paste0(strrep("  ", depth), "Called ", 
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);
    
    panels <- lapply(seq_len(n.vars()), function(i) {
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
                                    ticks = F)), 
              class = "compact-row"
            )
          })
        )
      });
    
    return(do.call(accordion, c(list(open = F), panels)));
  });
  
  # Отражение изменений решения на странице
  observe({
    shiny::req(x.df());
    for (i in seq_len(n.vars())) {
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
    
    shiny::req(get.w(), get.trg());
  }) |> debounce(250);
  
  observeEvent(solver.launcher(), {
    message(paste0(strrep("  ", depth), "Called ",
                   deparse(sys.call(0)[[1]])));
    depth <<- depth + 1;
    on.exit(depth <<- depth - 1);

    solve.system();
  }, ignoreInit = T);
};
