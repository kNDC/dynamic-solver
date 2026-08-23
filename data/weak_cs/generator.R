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
  metadata$name <- "Пример нестрогого условия дополняющей нежёсткости";
  
  # Уравнения
  metadata$eqs <- list();
  
  # Настройка неравеств
  metadata$neqs <- list();
  
  # Переменные
  metadata$endogenous <- list();
  
  metadata$endogenous$x <- list();
  metadata$endogenous$x$name <- "x";
  metadata$endogenous$x$nonneg <- TRUE;
  metadata$endogenous$x$initialiser <- TRUE;
  metadata$endogenous$x$unit <- "безразм.";
  
  # Внешние переменные
  metadata$exogenous <- list();
  
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
