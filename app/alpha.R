IDA_parse <- function(dat) {
  require(tidyverse)
  sample_names <- dat %>% filter(is.na(.[[1]]) & !is.na(.[[2]])) %>% .[[2]]
  sample_breaks <- c(0, which(is.na(dat[,1] & is.na(dat[,2]))))
  fn_obj <- list()
  if (length(sample_breaks) > 1){
    for (i in 2:length(sample_breaks)){
      fn_obj[i-1] <- list(dat %>% slice((sample_breaks[i-1]+2):(sample_breaks[i]-1)) %>% as.list())
      fn_obj[[i-1]][[2]] <- as.numeric(fn_obj[[i-1]][[2]])
      attr(fn_obj[[i-1]], "spec") <- NULL
    }
  } else {
    fn_obj[[1]] <- dat[-1,]
    fn_obj[[1]][[2]] <- as.numeric(fn_obj[[1]][[2]])
  }
  names(fn_obj) <- sample_names
  fn_obj <- fn_obj[order(names(fn_obj))]
  return(fn_obj)
}

IDA_background <- function(parsed_dat) {
  midpoint <- (max(parsed_dat)-min(parsed_dat))/2
  interest_region <- parsed_dat[parsed_dat < midpoint]
  dens <- density(interest_region)
  out <- dens$x[which(dens$y == max(dens$y))]
  out <- round(out,0)
  return(out)
}

IDA_stable <- function(dat, buffer = 10, tolerance = 0.1, expansion = 20) {
  require(tidyverse)
  dat <- as.data.frame(dat)
  interest_region <- dat %>% filter(.[[2]] > ((max(.[[2]])-min(.[[2]]))/2))
  d <- dim(interest_region)[1]
  interest_region <- interest_region[-c(d:(d-buffer)),]
  interest_region <- interest_region[-c(1:buffer),]
  start_point <- round(dim(interest_region)[1]/2)
  backward <- expansion
  forward <- expansion
  temp <- interest_region[c((start_point-backward):(start_point+forward)),]
  sd <- sd(temp$Ratio)
  
  # Go forward in time
  exit_check <- TRUE
  difference <- 0
  while (difference <= tolerance & exit_check) {
    temp <- interest_region[c((start_point-backward):(start_point+forward)),]
    difference <- (sd(temp$Ratio)-sd)/sd
    # Check ratio and 
    if (difference <= tolerance) {
      forward <- forward + expansion
    } else {
      exit_check <- !exit_check
    }
    if (forward + expansion > start_point) {
      if (forward + expansion < start_point + expansion) {
        forward <- dim(interest_region)[1] - start_point
      } else {
        forward <- forward - expansion
        exit_check <- !exit_check
      }
    }
  }
  temp <- interest_region[c((start_point-backward):(start_point+forward)),]
  
  # Go backward in time
  exit_check <- TRUE
  difference <- 0
  while (difference <= tolerance & exit_check) {
    temp <- interest_region[c((start_point-backward):(start_point+forward)),]
    difference <- (sd(temp$Ratio)-sd)/sd
    if (difference <= tolerance) {backward <- backward + expansion}
    if (backward + expansion > start_point) {
      if (backward + expansion < start_point + expansion) {
        backward <- start_point - 1
      } else {
        backward <- backward - expansion
        exit_check <- !exit_check
      }
    }
  }
  temp <- interest_region[c((start_point-backward):(start_point+forward)),]
  return(c(min(temp$Time), max(temp$Time)))
}

IDA_summary <- function(dat, stable.start, stable.end) {
  require(tidyverse)
  if (!is.data.frame(dat)) dat <- as.data.frame(dat)
  temp <- dat %>% filter(Time >= stable.start & Time <= stable.end)
  out <- data.frame(Mean = mean(temp$Ratio),
              StDev = sd(temp$Ratio),
              RSD = sd(temp$Ratio)/mean(temp$Ratio)*100)
  return(out)
}

IDA_calc <- function(dat, buffer = 10, tolerance = 0.1, expansion = 10, draw_stable_bounds) {
  require(tidyverse)
  m2 <- IDA_background(dat[[2]])
  m3 <- IDA_background(dat[[3]])
  dat2 <- data.frame(
    dat[[1]]/1000,
    dat[[2]] - m2,
    dat[[3]] - m3
  )
  names(dat2) <- names(dat)
  dat2 <- dat2 %>% mutate(Ratio=.[[2]]/.[[3]])
  t <- IDA_stable(dat2, buffer, tolerance, expansion)
  mostattributes(dat2) <- c(attributes(dat2),
                            "iso1.background" = m2,
                            "iso2.background" = m3,
                            "stable.start" = t[1],
                            "stable.end" = t[2],
                            "proc.buffer" = buffer,
                            "proc.tolerance" = tolerance,
                            "proc.expansion" = expansion)
  return(dat2)
}

IDA_graph <- function(dat, ind) {
  rat <- dat[[3]][[ind]]
  raw <- dat[[4]][[ind]]
  sam <- names(dat[[3]])[ind]
  iso <- attributes(dat)$isotopes
  return(IDA_graphs(raw, rat, sam, iso, draw_stable_bounds))
}

IDA_graphs <- function(dat_raw, dat_rat, sample_run, isotopes, draw_stable_bounds) {
  require(tidyverse)
  require(scales)
  stable_start <- attr(dat_rat, "stable.start")
  stable_end <- attr(dat_rat, "stable.end")
  isotope_label <- paste(strsplit(isotopes, " ")[[1]],
                    c("(black) ", "", "(blue)"),
                    collapse="")
  
  # Create signal plot
  signal <- dat_raw %>% as.data.frame %>%
    ggplot(aes(x = Time/1000))+
    geom_hline(yintercept=attr(dat_rat, "iso1.background"), colour='black')+
    geom_hline(yintercept=attr(dat_rat, "iso2.background"), colour='blue')+
    geom_line(aes(y=.[[2]]), colour='black')+
    geom_line(aes(y=.[[3]]), colour='blue')+
    labs(x = "Time (s)",
         y = "Signal Intensity",
         title = sample_run,
         subtitle = isotope_label)+
    theme_classic()+
    scale_y_continuous(labels=scientific_format(digits=2))
  if (draw_stable_bounds){
    signal <- signal+
    geom_vline(xintercept=stable_start, colour='red')+
    geom_vline(xintercept=stable_end, colour='red')
  }
  
  # Get functional RSD range limits
  y_range <- dat_rat %>% 
    as.data.frame() %>%
    filter(Time >= stable_start & Time <= stable_end) %>%
    pull(Ratio)
  y_range <- c(mean(y_range)-10*sd(y_range),
               mean(y_range)+10*sd(y_range))
  
  # Create ratio plot
  ratios <- dat_rat %>% 
    ggplot(aes(x = Time, y = Ratio))+
    geom_line()+
    labs(x = "Time (s)",
         y = paste("Signal Ratio", isotopes))+
    coord_cartesian(ylim=y_range)+
    theme_classic()+
    scale_y_continuous(labels=scientific_format(digits=2))
  if (draw_stable_bounds){
    ratios <- ratios+
      geom_vline(xintercept=stable_start, colour='red')+
      geom_vline(xintercept=stable_end, colour='red')+
      labs(caption = paste0(
        "\nSelected stability region is marked by red lines (",
        stable_start, "-", stable_end,
        " s)."))
  }
  return(list(Signal = signal, Ratio = ratios))
}

IDA_quality <- function(dat) {
  blanks <- c(grep("Blank", dat$Sample), grep("blank", dat$Sample))
  xbar <- mean(dat$RSD[-blanks])
  sd1 <- sd(dat$RSD[-blanks])
  n <- length(dat$Sample)
  samples_over <- dat[-blanks,] %>% filter(RSD > 5) %>% pull(Sample) %>% length()
  dat$Sample <- factor(dat$Sample, levels=rev(levels(dat$Sample)))
  out <-  ggplot(dat, aes(y=Sample, x=RSD))+
    geom_point(colour='white', alpha=0)
  if (n > 1) {
    out <- out+
      annotate("rect", xmin=xbar-2*sd1, xmax=xbar-sd1, ymin=-Inf, ymax=Inf, fill='blue', alpha=0.1)+
      annotate("rect", xmax=xbar+2*sd1, xmin=xbar+sd1, ymax=Inf, ymin=-Inf, fill='blue', alpha=0.1)+
      annotate("rect", xmax=xbar+sd1, xmin=xbar-sd1, ymax=Inf, ymin=-Inf, fill='green', alpha=0.2)+
      geom_vline(xintercept=xbar, colour='green')+
      geom_vline(xintercept=5, colour='red')
  }
  out <- out +
      geom_point(data=dat[-blanks,],
                 aes(y=Sample, x=RSD),
                 colour='black',
                 pch=16)+
      geom_point(data=dat[blanks,],
                 aes(y=Sample, x=RSD),
                 colour='blue',
                 pch=1)+
      geom_text(data=dat[blanks,],
                aes(label=Sample, y=Sample, x=RSD),
                colour='blue',
                hjust=1.1,
                vjust=0.5)+
      geom_point(data=dat[-blanks,] %>% filter(RSD > 5),
                 aes(y=Sample, x=RSD),
                 colour='red',
                 pch=16)+
      geom_text(data=dat[-blanks,] %>% filter(RSD > 5),
                aes(label=Sample, y=Sample, x=RSD),
                colour='red',
                hjust=-0.1,
                vjust=0.5)+
      labs(x="Relative Standard Deviation (%)",
           #title = "Quality Overview for this batch",
           subtitle = paste(samples_over, "samples over 5% RSD in this set.", sep=" "),
           caption = paste("",
                           "Any sample labeled 'blank' is displayed as open blue circles.",
                           "The green region is -/+1sd of measurements.",
                           "The blue region is 1sd-2sd of measurements.",
                           "The red line represents the 5% RSD threshold.",
                           sep="\n"))+
      theme_classic()
  return(out)
}

IDA <- function(file_name, buffer = 10, tolerance = 0.1, expansion = 10, draw_stable_bounds = TRUE) {
  require(tidyverse)
  dat <- read_csv(file_name)
  dat <- IDA_parse(dat)
  
  out <- as.data.frame(matrix(nrow=1, ncol=3))
  names(out) <- c("Mean", "StDev", "RSD")
  info <- as.data.frame(matrix(nrow=1, ncol=7))
  names(info) <- c("Isotope1", "Isotope2", "StableStart", "StableEnd", "Buffer", "Tolerance", "Expansion")
  full <- list()
  graphs <- list()
  
  isotopes <- paste(names(dat[[1]][2]), " / ", names(dat[[1]][3]), sep="")
  for (i in 1:length(dat)) {
    temp <- IDA_calc(dat[[i]], buffer, tolerance, expansion)
    out <- rbind(out,
                 IDA_summary(temp,
                             attr(temp, "stable.start"),
                             attr(temp, "stable.end")
                             )
                 )
    info <- rbind(info,
                  data.frame(Isotope1 = attr(temp, "iso1.background"),
                             Isotope2 = attr(temp, "iso2.background"),
                             StableStart = attr(temp, "stable.start"),
                             StableEnd = attr(temp, "stable.end"),
                             Buffer = attr(temp, "proc.buffer"),
                             Tolerance = attr(temp, "proc.tolerance"),
                             Expansion = attr(temp, "proc.expansion")
                             )
                  )
    full[[i]] <- as.list(temp)
    graphs[[i]] <- as.list(IDA_graphs(dat[[i]], temp, names(dat)[i], isotopes, draw_stable_bounds))
  }
  names(full) <- names(dat)
  names(graphs) <- names(dat)
  
  samples <- as.data.frame(strsplit(names(dat), "    "), row.names=NULL, stringsAsFactors=FALSE)
  names(samples) <- c()
  samples <- as.data.frame(t(samples))
  names(samples) <- c("Sample", "DateTime", "Run")
  out <- out[-1,]
  info <- info[-1,]
  names(info) <- c("Isotope1 Baseline", "Isotope2 Baseline", "Stable Time Start", "Stable Time End", "Buffer Size", "Tolerance", "Expansion Size")
  out <- cbind(samples, out)
  info <- cbind(samples, info)
  row.names(out) <- c()
  row.names(info) <- c()
  
  out <- list(Eval = out,
              Info = info,
              Processed = full,
              Raw = dat,
              Graphs = graphs,
              Quality = IDA_quality(out))
  
  mostattributes(out) <- c(attributes(out),
                           "isotopes" = isotopes)
  return(out)
}