#' ggplot2 extension for a MF object
#'
#' \code{ggMF}: the \code{\link[ggplot2]{ggplot}} extension for \code{MF} object to plot the correlation between species diversity and multi-functionality
#'
#' @param output the output from \code{MF_single} or \code{MF_multiple}. Note that, the column \code{Species.diversity} or the three columns \code{Species_gamma}, \code{Species_alpha} and \code{Species_beta} must be included in it, respectively.
#' @param by_group name of the column to be plot in difference colors by group. Default is \code{NULL}.\cr
#' For \code{fit} is selected to be linear mixed model, you must offer the \code{by_group} argument.
#' @param facets_scale Are scales shared across all facets (the default, \code{"fixed"}), or do they vary across rows (\code{"free_x"}), columns (\code{"free_y"}), or both rows and columns (\code{"free"})?
#' @param fit method of the fitted line. Select \code{fit = "lm"} for the linear model, or \code{fit = "LMM.intercept"}, \code{fit = "LMM.slope"} and \code{fit = "LMM.both"} for the linear mixed model with random effect 'intercept', 'slope' and 'both intercept and slope', respectively. Default is \code{fit = "LMM.intercept"}
#' @param text what text information would show in plots? Select \code{text = "Slope"} to see estimated of slopes, or \code{text = "R.squared"} to see model performance.Default is \code{text = "Slope"}.
#'
#' @import tidyverse
#' @import dplyr
#' @import broom
#' @import patchwork
#' @import ggpubr
#' @import grid
#'
#' @return For \code{MF_single} output, return a figure for multi-functionality with uncorrelated and correlated functions. For \code{MF_multiple} output, return a list of uncorrelated and correlated figures.
#'
#' @examples
#'
#' ## single ecosystem
#' data("Europe_Forest")
#' data("Europe_Forest_species")
#' Europe_Forest_function <- Europe_Forest[,4:29]
#' rownames(Europe_Forest_function) <- Europe_Forest$plotid
#' output1 <- MF_single(func_data = Europe_Forest_function, species_data = Europe_Forest_species)
#'
#' ## Display fitted line of linear mixed model with random effect 'intercept'.
#' ggMF(output1, facets_scale = 'fixed', fit = "LMM.intercept")
#'
#'
#' ## multiple ecosystems
#' Europe_Forest_function$Country <- Europe_Forest$Country
#' output2 <- MF_multiple(func_data = Europe_Forest_function,
#' species_data = Europe_Forest_species, by_group = "Country")
#'
#' ## Display fitted line of linear model
#' ggMF(output2, by_group = "Country", facets_scale = 'fixed', fit = "LMM.intercept")
#'
#'
#' @export

ggMF <- function(output, by_group = NULL, facets_scale = 'fixed', fit = "LMM.intercept", text = "Slope"){

  if(!(fit %in% c("lm","LMM.intercept","LMM.slope","LMM.both")))
    stop("Error: the argument fit should be `lm`, `LMM.intercept`, `LMM.slope` or `LMM.both`.")
  else if(fit != "lm" & is.null(by_group))
    stop("Error: linear mixed model must contain the argument by_group.")
  else if(fit != "lm")
    rand_eff <- strsplit(fit,split = "[.]") %>% unlist() %>% .[2]

  if(!is.null(by_group)){
    if(length(by_group)!=1) stop("Error: The number of the group variable should not be more than 1.")
    else if(!(by_group %in% names(output))) stop("Error: The group variable is not included in the given output data.")

  }

  if(!(text %in% c("Slope","R.squared")))
    stop("Error: the argument text should be `Slope` or `R.squared`.")

  stdPalette <- c("#FF88C2","gray55","purple2", "darkorange", "#00AAAA", "blue")
  #"steelblue1"
  #cbPalette <- c("#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

  if(!all(c("Type","Order.q","qMF","Species.diversity") %in% names(output)))
    stop("Error: For ggMF function, you should offer the MF_single output data including at least the columns `Type`, `Order.q`, `qMF` and `Species.diversity`, \n
            or the MF_multiple output data including at least the columns `Type`, `Scale`, `Order.q`, `qMF` and `Species.diversity`.")

  # For MF_single output
  else if(!("Scale" %in% names(output))){
    output$Type <- factor(output$Type, levels = c("Uncorrelated", "Correlated"))

    if(is.null(by_group)){
      lm_data <- output %>% group_by(Type, Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .)))
      lm_data <- mutate(lm_data,
                        Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Insignificant slope"),
                                            levels = c("Significant slope (P < 0.05)", "Insignificant slope")),
                        group="Linear model"
      ) %>%
        filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))

      output <- suppressMessages(output %>% dplyr::left_join(lm_data))

      plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF))+
        facet_grid(Type ~ Order.q, scales = facets_scale) +
        geom_point(size=0.7)+
        geom_smooth(aes(lty = Significance, col=group), method = "lm", se = F, size=1.2, formula = y ~ x)+
        geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=paste0("Slope = ", round(estimate, 4)), col=group),
                  hjust= -0.1, vjust= 2, size=3, key_glyph = draw_key_path)+
        scale_color_manual(values = "red")+
        theme_bw() +
        theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,0,0,0), legend.title = element_blank())+
        guides(
          linetype = guide_legend(override.aes = list(col = "#000000",size=0.7))
        )+
        labs(x = "Species diversity", y = "Multi-functionality")


    }
    else{
      output$group <- output %>% dplyr::select(by_group) %>% unlist()

      if(fit == "lm") lm_all <- output %>% group_by(Type, Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>%
          mutate(group="Linear model")
      else lm_all <- output %>% group_by(Type, Order.q) %>% do(Lmm_fit(.,r_effect = rand_eff)) %>% suppressMessages %>%
          mutate(group="Linear mixed model")

      lm_overall <- mutate(lm_all,
                           Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Insignificant slope"),
                                               levels = c("Significant slope (P < 0.05)", "Insignificant slope"))
      ) %>% dplyr::select(-c(std.error, statistic, p.value)) %>%
        pivot_wider(names_from = term,values_from = estimate) %>%
        group_by(Type,Order.q,group) %>%
        summarise(Significance=Significance[!is.na(Species.diversity)],
                  Intercept=sum(`(Intercept)`, na.rm = T),
                  estimate=sum(Species.diversity, na.rm=T)) %>% suppressMessages


      if(fit %in% c("lm","LMM.intercept")){
        lm_data <- output %>% group_by(Type, Order.q, group) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>%
          mutate(Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Insignificant slope"),
                                     levels = c("Significant slope (P < 0.05)", "Insignificant slope"))
          ) %>%
          filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))

        output <- suppressMessages(output %>% dplyr::full_join(lm_data))
        lm_data <- lm_data %>% ungroup() %>%
          dplyr::add_row(lm_overall %>% dplyr::select(-Intercept)) %>%
          mutate(group=factor(group,levels = unique(group)))

        plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF,col = group))+
          facet_grid(Type ~ Order.q, scales = facets_scale) +
          geom_point(size=0.9,alpha=0.2)+
          geom_smooth(aes(lty = Significance), method = "lm", se = F, size=0.8, formula = y ~ x)+
          geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.3, col="red",key_glyph = draw_key_path)
      }
      else{
        lm_data <- output %>% group_by(Type, Order.q) %>%
          do(Lmm_fit(.,r_effect = rand_eff,each_group = T)) %>% suppressMessages
        rangeinfo <- output %>% group_by(Type,Order.q,group) %>%
          summarise(x0 = min(Species.diversity),x1=max(Species.diversity))

        output <- suppressMessages(output %>% dplyr::full_join(left_join(lm_data,rangeinfo)))
        lm_data <- lm_data %>% ungroup() %>% select(-Intercept) %>%
          rename(estimate = Slope) %>%
          dplyr::add_row(lm_overall %>% dplyr::select(-c(Significance,Intercept))) %>%
          mutate(group=factor(group,levels = unique(group)))

        plot_output <- ggplot(data = output, aes(x = Species.diversity, y = qMF,col = group))+
          facet_grid(Type ~ Order.q, scales = facets_scale) +
          geom_point(size=0.9,alpha=0.2)+
          geom_segment(aes(x=x0,xend=x1,y=Intercept+Slope*x0,yend=Intercept+Slope*x1),size=0.5)+
          geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.3, col="red",key_glyph = draw_key_path)
      }


      lm_data <- lm_data %>% group_by(Type,Order.q,group) %>%
        summarise(Label=ifelse(round(estimate,2)==round(estimate,1),
                               paste0(round(estimate,2),"0"),
                               as.character(round(estimate,2))) %>%
                    paste("Slope = ",.,sep = ifelse(round(estimate,2)<0,""," "))) %>%
        mutate(h=rep(c(-0.1,-1.5),(length(unique(group))+1)/2)[1:length(unique(group))],
               v=rep(2+(0:((length(unique(group))+1)/2-1))*1.5,each=2)[1:length(unique(group))]) %>% suppressMessages

      output$group <- factor(output$group,levels = levels(lm_data$group))
      col_manual <- c(stdPalette[1:(length(levels(output$group))-1)],"red") %>%
        `names<-`(levels(output$group))


      if(text == "Slope"){
        plot_output <- plot_output +
          geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3,key_glyph = draw_key_path)+
          scale_colour_manual(values = col_manual) +
          theme_bw() +
          guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "#000000",size=0.6)),
                 col = guide_legend(title = by_group))+
          theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,0,0,0))+
          labs(x = "Species diversity", y = "Multi-functionality")
      }
      else{
        lm_text <- lm_all %>% select(Type, Order.q,R2_C,R2_M) %>%
          pivot_longer(cols = c(R2_C,R2_M),names_to = "R2_type",values_to = "Label") %>%
          mutate(Label=paste(ifelse(R2_type=="R2_C","R[C]^2==","R[M]^2=="),
                             round(Label,3),sep=" "),
                 h=ifelse(R2_type=="R2_C",-0.3,-1.8),
                 v=2)

        plot_output <- plot_output +
          geom_text(data = lm_text, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3,key_glyph = draw_key_path, parse = T,col="red")+
          scale_colour_manual(values = col_manual) +
          theme_bw() +
          guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "red",size=0.6)),
                 col = guide_legend(title = by_group))+
          theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,0,0,0))+
          labs(x = "Species diversity", y = "Multi-functionality")
      }
    }
  }
  else{
    multiple_plots <- function(data){

      scale_plot <- function(out, xlab, ylab, digit=2, h_j=-1.5){
        if(is.null(by_group)){
          lm_data <- out %>% group_by(Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .)))
          lm_data <- mutate(lm_data,
                            Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Insignificant slope"),
                                                levels = c("Significant slope (P < 0.05)", "Insignificant slope")),
                            group="Linear model"
          ) %>%
            filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))
          out <- suppressMessages(out %>% dplyr::left_join(lm_data))
          plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF))+
            facet_grid(~Order.q, scales = facets_scale) +
            geom_point(size=0.7,alpha=0.2)+
            geom_smooth(aes(lty = Significance, col=group), method = "lm", se = F, size=1.2, formula = y ~ x)+
            geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=paste0("Slope = ", round(estimate, 4)), col=group),
                      hjust= -0.1, vjust= 2, size=3, key_glyph = draw_key_path)+
            scale_color_manual(values = "red")+
            theme_bw() +
            theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,0,0,0), legend.title = element_blank())+
            guides(
              linetype = guide_legend(override.aes = list(col = "#000000",size=0.7))
            )+
            labs(x = xlab, y = ylab)


        }
        else{
          out$group <- out %>% dplyr::select(by_group) %>% unlist()

          if(fit == "lm") lm_all <- out %>% group_by(Order.q) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>%
              mutate(group="Linear model")
          else lm_all <- out %>% group_by(Order.q) %>% do(Lmm_fit(.,r_effect = rand_eff)) %>% suppressMessages %>%
              mutate(group="Linear mixed model")

          lm_overall <- mutate(lm_all,
                               Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Insignificant slope"),
                                                   levels = c("Significant slope (P < 0.05)", "Insignificant slope"))
          ) %>% dplyr::select(-c(std.error, statistic, p.value)) %>%
            pivot_wider(names_from = term,values_from = estimate) %>%
            group_by(Order.q,group) %>%
            summarise(Significance=Significance[!is.na(Species.diversity)],
                      Intercept=sum(`(Intercept)`, na.rm = T),
                      estimate=sum(Species.diversity, na.rm=T)) %>% suppressMessages


          if(fit %in% c("lm","LMM.intercept")){
            lm_data <- out %>% group_by(Order.q, group) %>% do(broom::tidy(lm(qMF ~ Species.diversity, .))) %>%
              mutate(Significance=factor(ifelse(p.value<0.05, "Significant slope (P < 0.05)", "Insignificant slope"),
                                         levels = c("Significant slope (P < 0.05)", "Insignificant slope"))
              ) %>%
              filter(term=="Species.diversity") %>% dplyr::select(-c(term, std.error, statistic, p.value))

            out <- suppressMessages(out %>% dplyr::full_join(lm_data))
            lm_data <- lm_data %>% ungroup() %>%
              dplyr::add_row(lm_overall %>% dplyr::select(-Intercept)) %>%
              mutate(group=factor(group,levels = unique(group)))

            plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF,col = group))+
              facet_grid( ~ Order.q, scales = facets_scale) +
              geom_point(size=0.7,alpha=0.05)+
              geom_smooth(aes(lty = Significance), method = "lm", se = F, size=0.8, formula = y ~ x)+
              geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.3, col="red",key_glyph = draw_key_path)
          }
          else{
            lm_data <- out %>% group_by(Order.q) %>%
              do(Lmm_fit(.,r_effect = rand_eff,each_group = T)) %>% suppressMessages
            rangeinfo <- out %>% group_by(Order.q,group) %>%
              summarise(x0 = min(Species.diversity),x1=max(Species.diversity))

            out <- suppressMessages(out %>% dplyr::full_join(left_join(lm_data,rangeinfo)))
            lm_data <- lm_data %>% ungroup() %>% select(-Intercept) %>%
              rename(estimate = Slope) %>%
              dplyr::add_row(lm_overall %>% dplyr::select(-c(Significance,Intercept))) %>%
              mutate(group=factor(group,levels = unique(group)))

            plot_output <- ggplot(data = out, aes(x = Species.diversity, y = qMF,col = group))+
              facet_grid(~Order.q, scales = facets_scale) +
              geom_point(size=0.7,alpha=0.05)+
              geom_segment(aes(x=x0,xend=x1,y=Intercept+Slope*x0,yend=Intercept+Slope*x1),size=0.5)+
              geom_abline(data = lm_overall, aes(slope=estimate, intercept=Intercept, lty = Significance), size=1.3, col="red",key_glyph = draw_key_path)
          }


          lm_data <- lm_data %>% group_by(Order.q,group) %>%
            summarise(Label=ifelse(round(estimate,digit)==0,paste0("0.",paste0(rep(0,digit),collapse = "")),
                                   ifelse(round(estimate,digit)==round(estimate,digit-1),
                                          paste0(round(estimate,digit),"0"),
                                          as.character(round(estimate,digit)))) %>%
                        paste("Slope = ",.,sep = ifelse(round(estimate,digit)<0,""," "))) %>%
            mutate(h=rep(c(-0.1,h_j),(length(unique(group))+1)/2)[1:length(unique(group))],
                   v=rep(2+(0:((length(unique(group))+1)/2-1))*1.5,each=2)[1:length(unique(group))]) %>% suppressMessages

          out$group <- factor(out$group,levels = levels(lm_data$group))
          col_manual <- c(stdPalette[1:(length(levels(out$group))-1)],"red") %>%
            `names<-`(levels(out$group))


          if(text == "Slope"){
            plot_output <- plot_output +
              geom_text(data = lm_data, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3,key_glyph = draw_key_path)+
              scale_colour_manual(values = col_manual) +
              theme_bw() +
              guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "#000000",size=0.6)),
                     col = guide_legend(title = by_group))+
              theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,0,0,0))+
              labs(x = xlab, y = ylab)
          }
          else{
            lm_text <- lm_all %>% select(Order.q,R2_C,R2_M) %>%
              pivot_longer(cols = c(R2_C,R2_M),names_to = "R2_type",values_to = "Label") %>%
              mutate(Label=paste(ifelse(R2_type=="R2_C","R[C]^2==","R[M]^2=="),
                                 round(Label,3),sep=" "),
                     h=ifelse(R2_type=="R2_C",-0.3,-1.8),
                     v=2)

            plot_output <- plot_output +
              geom_text(data = lm_text, aes(x = -Inf, y = Inf, label=Label, hjust= h, vjust= v), size=3,key_glyph = draw_key_path, parse = T,col="red")+
              scale_colour_manual(values = col_manual) +
              theme_bw() +
              guides(linetype = guide_legend(title = "",order = 1,override.aes = list(col = "red",size=0.6)),
                     col = guide_legend(title = by_group))+
              theme(legend.position = "bottom", legend.box = "vertical", legend.margin=margin(-6,0,0,0))+
              labs(x = xlab, y = ylab)
          }

          plot_output
        }
      }


      gamma_data <- data %>% filter(Scale=="Gamma")
      gamma <- scale_plot(gamma_data, "Gamma species diversity", "Gamma MF")

      alpha_data <- data %>% filter(Scale=="Alpha")
      alpha <- scale_plot(alpha_data, "Alpha species diversity", "Alpha MF")

      beta_data0 <- data %>% filter(Scale=="Beta", Order.q == "q = 0")
      beta0 <- scale_plot(beta_data0, NULL, "Beta MF", 3,h_j=-1.25)

      beta_data12 <- data %>% filter(Scale=="Beta", Order.q != "q = 0")
      beta12 <- scale_plot(beta_data12, NULL, NULL,h_j=-1.25)


      beta <- ggpubr::ggarrange(beta0,beta12,nrow = 1,widths=c(1.5,2.4),legend = "none") %>%
        ggpubr::annotate_figure(bottom = textGrob("Beta species diversity \n", gp = gpar(cex = 0.92),hjust = 0.35))

      combine_plots <- ggpubr::ggarrange(gamma,alpha,beta,ncol = 1,heights = c(1,1,1.12),legend = "bottom",common.legend = T)

      out <- list("ALL"=combine_plots,
                  "Gamma"=gamma,
                  "Alpha"=alpha,
                  "Beta"=ggpubr::ggarrange(alpha,beta,ncol = 1,heights = c(0,1),legend = "bottom",common.legend = T))
      return(out)
    }

    plot_output <- lapply(unique(output$Type),function(x) multiple_plots(output %>% filter(Type==x))) %>%
      `names<-`(unique(output$Type))

  }
  return(plot_output)
}


# Get the summary table of the linear mixed model fitted.
#
# \code{LMM_fit} Get the summary table of the linear mixed model fitted.
#
# @param data a data.frame contains 'qMF', 'Species.diversity' and 'group'.
# @param r_effect decide the random effect is "intercept", "slope" or "both".
# @each_group
# @return a summary table for the fitted model.

Lmm_fit <- function(data, r_effect = "intercept", each_group = FALSE){
  if(r_effect == "intercept") x <- lmerTest::lmer(qMF~Species.diversity+(1|group),data) %>% suppressWarnings()
  else if(r_effect == "slope") x <- lmerTest::lmer(qMF~Species.diversity+(0+Species.diversity|group),data) %>% suppressWarnings()
  else if(r_effect == "both") x <- lmerTest::lmer(qMF~Species.diversity+(1+Species.diversity|group),data) %>% suppressWarnings()

  if(each_group){
    ret <- as_tibble(data.frame(group=rownames(coef(x)$group),coef(x)$group))
    colnames(ret) <- c("group","Intercept","Slope")
  }

  else{
    johnson_r2 <- function(model){
      X <- model.matrix(model)
      n <- nrow(X)
      Beta <- lme4::fixef(model)
      Sf <- var(X %*% Beta)
      Sigma.list <- lme4::VarCorr(model)
      Sl <-
        sum(
          sapply(Sigma.list,
                 function(Sigma)
                 {
                   Z <-X[,rownames(Sigma)]
                   sum(diag(Z %*% Sigma %*% t(Z)))/n
                 }))
      Se <- attr(Sigma.list, "sc")^2
      Sd <- 0
      total.var <- Sf + Sl + Se + Sd
      Rsq.m <- Sf / total.var
      Rsq.c <- (Sf + Sl) / total.var

      return(c(Rsq.c,Rsq.m))
    }
    r2_j <- johnson_r2(x)
    ret <- as_tibble(cbind(summary(x)$coefficients[,-3],"R2_C"=r2_j[1],"R2_M"=r2_j[2]), rownames = "term")
    colnames(ret) <- c("term", "estimate", "std.error","statistic","p.value","R2_C","R2_M")
  }

  ret
}
