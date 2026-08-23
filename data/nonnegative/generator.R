library(yyjsonr)

get_directory <- function() {
  args <- commandArgs(trailingOnly = F);
  
  if (args[1] == "RStudio") {
    frames <- sys.frames();
    frames.files <- unlist(lapply(
      frames, 
      function(x) {
        if (!is.null(x$ofile)) return(x$ofile);
        if ("fileName" %in% ls(x)) return(x$fileName);
        return(NULL)
      }
    ));
    
    if (length(frames.files)) {
      return(dirname(frames.files[length(frames.files)]));
    }
  }
  
  return(NULL);
};

read.metadata <- function() {
  setwd(get_directory());
  print(getwd());
  
  tryCatch(
    yyjsonr::read_json_file("./metadata.json"), 
	error = function(e) NULL
  )
};

#--- Метаданные модели ---
write.metadata <- function(metadata = NULL) {
  if (is.null(metadata)) metadata <- list();
  
  if (is.null(metadata$Id)) metadata$Id <- ids::uuid();
  metadata$name <- "Гармонический осциллятор (простой пример)";
  
  # Уравнения
  metadata$eqs <- list();
  
  metadata$eqs$x_init <- list();
  metadata$eqs$x_init$formula <- "x[1] - x0[1]";
  
  metadata$eqs$xp_init <- list();
  metadata$eqs$xp_init$formula <- "xp[1] - xp0[1]";
  
  metadata$eqs$xpp <- list();
  metadata$eqs$xpp$formula <- "dd(xp) + adj(q*x)";
  
  metadata$eqs$x_xp_aux <- list();
  metadata$eqs$x_xp_aux$formula <- "dd(x) - adj(xp)";
  
  # Настройка неравеств
  metadata$neqs <- list();
  
  # Переменные
  metadata$endogenous <- list();
  
  metadata$endogenous$x <- list();
  metadata$endogenous$x$name <- "x";
  metadata$endogenous$x$unit <- "безразм.";
  
  metadata$endogenous$xp <- list();
  metadata$endogenous$xp$name <- "x'";
  metadata$endogenous$xp$unit <- "безразм.";
  
  # Внешние переменные
  metadata$exogenous <- list();
  
  metadata$exogenous$x0 <- list();
  metadata$exogenous$x0$name <- "x(0)";
  metadata$exogenous$x0$unit <- "безразм.";
  
  metadata$exogenous$xp0 <- list();
  metadata$exogenous$xp0$name <- "x'(0)";
  metadata$exogenous$xp0$unit <- "безразм.";
  
  metadata$exogenous$q <- list();
  metadata$exogenous$q$name <- "q";
  metadata$exogenous$q$unit <- "безразм.";
  
  # Настройка параметров
  metadata$parameters <- list();
  
  yyjsonr::write_json_file(
    x = metadata, filename = "./metadata.json", 
    auto_unbox = T, pretty = T);
};

metadata <- read.metadata();
write.metadata(metadata = metadata);
