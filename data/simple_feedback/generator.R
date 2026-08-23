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
  metadata$name <- "Простая задача с обратной связью";
  
  # Настройка уравнений
  metadata$eqs <- list();
  
  metadata$eqs$x_def <- list();
  metadata$eqs$x_def$formula <- "x - sin(qx * t) - 5";
  
  metadata$eqs$y_def <- list();
  metadata$eqs$y_def$formula <- "y - qy * h";
  
  metadata$eqs$h_init <- list();
  metadata$eqs$h_init$formula <- "h[1] - h0[1]";
  
  metadata$eqs$h_def <- list();
  metadata$eqs$h_def$formula <- "dd(h) - x[-T] + y[-T]";
  
  # Настройка неравеств
  metadata$neqs <- list();
  
  # Настройка внутренних переменных
  metadata$endogenous <- list();
  
  metadata$endogenous$x <- list();
  metadata$endogenous$x$name <- "x0";
  metadata$endogenous$x$unit <- "безразм.";
  
  metadata$endogenous$y <- list();
  metadata$endogenous$y$name <- "xp0";
  metadata$endogenous$y$unit <- "безразм.";
  
  # Настройка внешних переменных
  metadata$exogenous <- list();
  
  metadata$exogenous$qx <- list();
  metadata$exogenous$qx$name <- "qx";
  metadata$exogenous$qx$unit <- "безразм.";
  
  metadata$exogenous$qy <- list();
  metadata$exogenous$qy$name <- "qy";
  metadata$exogenous$qy$unit <- "безразм.";
  
  metadata$exogenous$h0 <- list();
  metadata$exogenous$h0$name <- "h0";
  metadata$exogenous$h0$unit <- "безразм.";
  
  # Настройка параметров
  metadata$parameters <- list();
  
  yyjsonr::write_json_file(
    x = metadata, filename = "./metadata.json", 
    auto_unbox = T, pretty = T);
};

metadata <- read.metadata();
write.metadata(metadata = metadata);
