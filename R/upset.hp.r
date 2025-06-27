#' @title Visualization of VP and HP Using UpSet Diagram
#'
#' @description Visualization of variation partitioning (VP) and hierarchical partitioning (HP) with unlimited number of predictor variables (or matrices of predictors) using UpSet matrix layout.
#'
#' @param vp A matrix, which contains the output of variation partitioning (i.e. commonality analysis) from \code{rdacca.hp},\code{glmm.hp},\code{gam.hp},and \code{phylolm.hp}.
#' @param hp A matrix, which contains the output of hierarchical partitioning from \code{rdacca.hp},\code{glmm.hp},\code{gam.hp},and \code{phylolm.hp}.
#' @param plot.hp The default is \code{TRUE}, which plots the individual effect for each predictor on left column diagram. If \code{FALSE}, compute and plot the sum of unique effect and common effect for each predictor.
#' @param order.part How the VP components in matrix layout should be ordered. Options include \code{"effect"} (order the intersections by their effects) or \code{"degree"} (sort by the number of predictors involved in the intersection), default is \code{"effect"}.
#' @param decreasing.part How the intersections in \code{order.part} should be ordered. Default is \code{TRUE}, \code{"effect"} is decreasing (from greatest to least) or \code{"degree"} is increasing (from least to greatest).
#' @param order.var The predictors in the matrix layout should be ordered by. Default is \code{TRUE}, which orders the predictors by their effects. IF \code{FALSE}, sort by the order of predictors in input data.
#' @param decreasing.var If \code{order.var=TRUE}, how the predictors should be ordered. Default is \code{TRUE}, from greatest to least.
#' @param cutoff Effects below \code{cutoff} will not be displayed, default is \code{-1}. Note: Negative effects due to adjustment of R-squared mean negligible contributions, but they are included in the computation of the total contribution of each predictor category.
#' @param nVar Number of components in VP to plot, default is \code{30}.
#' @param col.width Width of bars in column diagram, default is \code{0.6}.
#' @param pch.size Size of points in matrix diagram, default is \code{3}.
#' @param line.lwd Width of lines in matrix diagram, default is \code{0.5}.
#' @param show.effect Show the relative importance of predictors (unique, common, or individual effects) above bars, default is \code{TRUE}.
#' @param effect.cex Font size of the effects, default is \code{2.7}.
#' @param title.cex Font size of axis titles, default is \code{10}.
#' @param axis.cex Font size of axis labels, default is \code{8}.
#' @param height.ratio Ratio between matrix and top column diagram, default is \code{c(2, 1)}.
#' @param width.ratio Ratio between matrix and left column diagram, default is \code{c(1, 3)}.
#' @param col Character. Color palette name: "nature" (default), "science", "cell", "bw" (black-white), or "cvd" (color-blind friendly).

#'
#' @details upset.hp diagram is an extension of UpSet technique to  and is used to visualize the object of \code{rdacca.hp},\code{glmm.hp},\code{gam.hp},and \code{phylolm.hp} (Lai et al. 2022a,2022b,2023,2024; Liu et al. 2023). The matrix layout enables the effective representation of relative importance of predictors, such as the unique effects and common effects in VP, as well as additional summary statistics or individual effects in HP. upset.hp diagram could, in principle, allow visualization of any number of predictor variables or groups of predictor variables. But considering the interpretability of data, we would like to recommend that the number of predictors (or groups of predictors) no more than 7.
#'
#' @return \itemize{Returns a ggplot2.}
#'
#' @references Lai J., Zou Y., Zhang J., Peres-Neto P. (2022) Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package. Methods in Ecology and Evolution, 13:782-788.
#' @references Lai J.,Zou Y., Zhang S.,Zhang X.,Mao L.(2022)glmm.hp: an R package for computing individual effect of predictors in generalized linear mixed models.Journal of Plant Ecology,15(6):1302-1307<DOI:10.1093/jpe/rtac096>
#' @references Lai J.,Zhu W., Cui D.,Mao L.(2023)Extension of the glmm.hp package to Zero-Inflated generalized linear mixed models and multiple regression.Journal of Plant Ecology,16(6):rtad038<DOI:10.1093/jpe/rtad038>
#' @references Liu Y., Yu X., Yu Y., et al. (2023) Application of "rdacca. hp" R package in ecological data analysis: case and progress. Chinese Journal of Plant Ecology, 27:134-144.
#' @references Lai J.,Tang J., Li T., Zhang A.,Mao L.(2024)Evaluating the relative importance of predictors in Generalized Additive Models using the gam.hp R package.Plant Diversity,46(4):542-546<DOI:10.1016/j.pld.2024.06.002>
#' @references Lai J.,He Y., Hou M., Zhang A.,Wang G., Mao L.(2025)Evaluating the relative importance of phylogeny and predictors in Phylogenetic Generalized Linear Models using the phylolm.hp R package. Plant Diversity. https://doi.org/10.1016/j.pld.2025.06.003

#' @export
#' @examples
#' library(glmm.hp)
#' #upset for glmm.hp() in lm()
#' m2<-lm(mpg~wt+carb+cyl,mtcars)
#' vp <- glmm.hp(m2,commonality=TRUE)$commonality.analysis
#' hp <- glmm.hp(m2)$hierarchical.partitioning
#' upset.hp(vp, hp, col = "cvd")



upset.hp <- function(vp, hp, plot.hp = TRUE, order.part = "effect", decreasing.part = TRUE, 
                     order.var = TRUE, decreasing.var = TRUE, cutoff = -1, nVar = 30, 
                     col.width = 0.6, pch.size = 3, line.lwd = 0.5, show.effect = TRUE, 
                     effect.cex = 2.7, title.cex = 10, axis.cex = 8, height.ratio = c(2, 1), 
                     width.ratio = c(1, 3), col = "nature") {
  
  # Check required packages
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required but not installed.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required but not installed.")
  }
  
  # Define color palettes
  color_palettes <- list(
    nature = list(
      pos_bars = c("#4E79A7", "#A0CBE8", "#76B7B2", "#59A14F", "#8CD17D"),
      neg_bars = c("#FFBE7D", "#F28E2B", "#E15759", "#FF9D9A", "#B07AA1"),
      point_active = "#2CA02C",
      point_inactive = "#D3D3D3",
      tile_even = "#F7F7F7",
      tile_odd = "#FFFFFF",
      line_color = "#636363",
      effect_text = "#000000",
      axis_text = "#000000"
    ),
    science = list(
      pos_bars = c("#1F77B4", "#8C564B", "#E377C2", "#7F7F7F", "#BCBD22"),
      neg_bars = c("#FF7F0E", "#D62728", "#9467BD", "#17BECF", "#2CA02C"),
      point_active = "#D62728",
      point_inactive = "#C7C7C7",
      tile_even = "#F0F0F0",
      tile_odd = "#FFFFFF",
      line_color = "#525252",
      effect_text = "#000000",
      axis_text = "#000000"
    ),
    cell = list(
      pos_bars = c("#3E7DCC", "#63B5F7", "#A5D8FF", "#4DBBD5", "#00A087"),
      neg_bars = c("#E64B35", "#F39B7F", "#FEB24C", "#FFD700", "#91D1C2"),
      point_active = "#00A087",
      point_inactive = "#DFDFDF",
      tile_even = "#F5F5F5",
      tile_odd = "#FFFFFF",
      line_color = "#737373",
      effect_text = "#000000",
      axis_text = "#000000"
    ),
    bw = list(
      pos_bars = rep("#666666", 5),
      neg_bars = rep("#CCCCCC", 5),
      point_active = "#000000",
      point_inactive = "#FFFFFF",
      tile_even = "#F0F0F0",
      tile_odd = "#FFFFFF",
      line_color = "#333333",
      effect_text = "#000000",
      axis_text = "#000000"
    ),
    cvd = list(
      pos_bars = c("#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677"),
      neg_bars = c("#AA3377", "#BBBBBB", "#EE8866", "#FFAABB", "#AAAA00"),
      point_active = "#004488",
      point_inactive = "#DDDDDD",
      tile_even = "#F7F7F7",
      tile_odd = "#FFFFFF",
      line_color = "#333333",
      effect_text = "#000000",
      axis_text = "#000000"
    )
  )
  
  # Validate color palette
  if (!col %in% names(color_palettes)) {
    warning(paste("Color palette '", col, "' not found. Using 'nature' instead. Available options:",
                  paste(names(color_palettes), collapse = ", ")))
    col <- "nature"
  }
  colors <- color_palettes[[col]]
  
  # Data validation
  if (missing(vp) || missing(hp) || nrow(vp) == 0 || nrow(hp) == 0) {
    stop("Input data 'vp' or 'hp' is missing or empty.")
  }
  
  # Process VP data
  Constrained <- 100 * tail(vp[, 1], 1)
  Var.part <- as.data.frame(vp[-nrow(vp), , drop = FALSE])
  Var.part$Var <- gsub("^\\s+|\\s+$", "", 
                       gsub("and ", "", 
                            gsub("Common to ", "", 
                                 gsub("Unique to ", "", rownames(Var.part)))))
  Var.part$Fractions <- 100 * Var.part[, 1]
  
  # Process HP data
  Hier.part <- as.data.frame(hp)
  Hier.part$Var <- rownames(Hier.part)
  Hier.part$Individual <- 100 * Hier.part[, 1]
  
  # Filter and order VP components
  if (nrow(Var.part) > 0) {
    Var.part$inter <- lengths(strsplit(as.character(Var.part$Var), ", "))
    Var.part$valid <- ifelse(Var.part$Fractions <= 0, "neg", "pos")
    Var.part <- Var.part[Var.part$Fractions >= 100 * cutoff, , drop = FALSE]
    
    if (nrow(Var.part) == 0) {
      warning("No variables passed the cutoff threshold.")
    } else {
      if (order.part == "effect") {
        Var.part <- Var.part[order(Var.part$Fractions, decreasing = decreasing.part), ]
      } else if (order.part == "degree") {
        Var.part <- Var.part[order(Var.part$inter, Var.part$Fractions, 
                                   decreasing = c(!decreasing.part, TRUE)), ]
      }
      if (nrow(Var.part) > nVar) Var.part <- Var.part[1:nVar, ]
      Var.part$Var <- factor(Var.part$Var, levels = Var.part$Var)
    }
  }
  
  # Create VP plot
  if (nrow(Var.part) > 0) {
    bar_colors <- ifelse(Var.part$valid == "pos",
                         rep(colors$pos_bars, length.out = nrow(Var.part)),
                         rep(colors$neg_bars, length.out = nrow(Var.part)))
    names(bar_colors) <- as.character(Var.part$Var)
    
    p.vp <- ggplot2::ggplot(Var.part, ggplot2::aes(x = Var, y = Fractions, fill = Var)) +
      ggplot2::geom_col(width = col.width) +
      ggplot2::scale_fill_manual(values = bar_colors) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),
        axis.line.y = ggplot2::element_line(color = colors$axis_text),
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_text(color = colors$axis_text, size = axis.cex),
        axis.ticks.x = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_line(color = colors$axis_text),
        axis.title = ggplot2::element_text(color = colors$axis_text, size = title.cex),
        plot.title = ggplot2::element_text(hjust = 0.5, size = title.cex),
        legend.position = "none"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(mult = c(ifelse(min(Var.part$Fractions) < 0, 0.1, 0), 0.1))
      ) +
      ggplot2::labs(
        y = "Fractions (%)", 
        x = "", 
        title = bquote(R^2 == .(Constrained)*"%" ~ "  Residual" ~ .(100 - Constrained)*"%")
      )
    
    if (show.effect) {
      p.vp <- p.vp + 
        ggplot2::geom_text(
          ggplot2::aes(label = round(Fractions, 1), 
                       vjust = ifelse(Fractions >= 0, -0.2, 1.2)),
          color = colors$effect_text, 
          size = effect.cex
        )
    }
    
    p.vp <- p.vp + 
      ggplot2::geom_hline(yintercept = 0, color = colors$axis_text)
  } else {
    p.vp <- ggplot2::ggplot() + 
      ggplot2::annotate(
        "text", x = 1, y = 1, 
        label = "No variables meet the cutoff criteria", 
        size = 6
      ) +
      ggplot2::theme_void()
  }
  
  # Create explanatory variables plot
  if (nrow(Var.part) > 0 && nrow(Hier.part) > 0) {
    Fractions <- sapply(rownames(Hier.part), function(i) {
      sum(Var.part[grep(i, Var.part$Var), "Fractions"])
    })
    Var.exp <- data.frame(Var = rownames(Hier.part), Fractions = Fractions)
    Var.exp$valid <- ifelse(Var.exp$Fractions <= 0, "neg", "pos")
    Var.exp <- Var.exp[Var.exp$Fractions >= 100 * cutoff, , drop = FALSE]
    
    if (order.var) {
      Var.exp <- Var.exp[order(Var.exp$Fractions, decreasing = !decreasing.var), ]
    }
    Var.exp$Var <- factor(Var.exp$Var, levels = Var.exp$Var)
    
    if (nrow(Var.exp) > 0) {
      exp_colors <- ifelse(Var.exp$valid == "pos",
                           rep(colors$pos_bars, length.out = nrow(Var.exp)),
                           rep(colors$neg_bars, length.out = nrow(Var.exp)))
      names(exp_colors) <- Var.exp$Var
      
      p.exp <- ggplot2::ggplot(Var.exp, ggplot2::aes(x = Var, y = Fractions, fill = Var)) +
        ggplot2::geom_col(width = col.width) +
        ggplot2::scale_fill_manual(values = exp_colors) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.line.x = ggplot2::element_line(color = colors$axis_text),
          axis.text.x = ggplot2::element_text(
            color = colors$axis_text, 
            size = axis.cex, 
            angle = 45, 
            hjust = 1
          ),
          axis.text.y = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_line(color = colors$axis_text),
          axis.ticks.y = ggplot2::element_blank(),
          axis.title = ggplot2::element_text(color = colors$axis_text, size = title.cex),
          legend.position = "none"
        ) +
        ggplot2::coord_flip() +
        ggplot2::scale_y_reverse(
          expand = ggplot2::expansion(mult = c(0.3, ifelse(min(Var.exp$Fractions) < 0, 0.3, 0)))
        ) +
        ggplot2::labs(y = "Fractions (%)", x = NULL)
      
      if (show.effect) {
        p.exp <- p.exp + 
          ggplot2::geom_text(
            ggplot2::aes(label = round(Fractions, 1), 
                         hjust = ifelse(Fractions >= 0, 1.2, -0.2)),
            color = colors$effect_text, 
            size = effect.cex
          )
      }
      
      p.exp <- p.exp + 
        ggplot2::geom_hline(yintercept = 0, color = colors$axis_text)
    } else {
      p.exp <- ggplot2::ggplot() + 
        ggplot2::annotate(
          "text", x = 1, y = 1, 
          label = "No explanatory variables meet the cutoff", 
          size = 6
        ) +
        ggplot2::theme_void()
    }
  } else {
    p.exp <- ggplot2::ggplot() + 
      ggplot2::annotate(
        "text", x = 1, y = 1, 
        label = "No data to display", 
        size = 6
      ) +
      ggplot2::theme_void()
  }
  
  # Create HP plot if requested
  if (plot.hp) {
    if (nrow(Hier.part) > 0) {
      Hier.part <- Hier.part[Hier.part$Individual >= 100 * cutoff, , drop = FALSE]
      
      if (order.var) {
        Hier.part <- Hier.part[order(Hier.part$Individual, decreasing = !decreasing.var), ]
      }
      Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
      Hier.part$valid <- ifelse(Hier.part$Individual <= 0, "neg", "pos")
      
      if (nrow(Hier.part) > 0) {
        hp_colors <- ifelse(Hier.part$valid == "pos",
                            rep(colors$pos_bars, length.out = nrow(Hier.part)),
                            rep(colors$neg_bars, length.out = nrow(Hier.part)))
        names(hp_colors) <- Hier.part$Var
        
        p.hp <- ggplot2::ggplot(Hier.part, ggplot2::aes(x = Var, y = Individual, fill = Var)) +
          ggplot2::geom_col(width = col.width) +
          ggplot2::scale_fill_manual(values = hp_colors) +
          ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            panel.background = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(color = colors$axis_text),
            axis.text.x = ggplot2::element_text(
              color = colors$axis_text, 
              size = axis.cex, 
              angle = 45, 
              hjust = 1
            ),
            axis.text.y = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_line(color = colors$axis_text),
            axis.ticks.y = ggplot2::element_blank(),
            axis.title = ggplot2::element_text(color = colors$axis_text, size = title.cex),
            legend.position = "none"
          ) +
          ggplot2::coord_flip() +
          ggplot2::scale_y_reverse(
            expand = ggplot2::expansion(mult = c(0.3, ifelse(min(Hier.part$Individual) < 0, 0.3, 0)))
          ) +
          ggplot2::labs(y = "Individual (%)", x = NULL)
        
        if (show.effect) {
          p.hp <- p.hp + 
            ggplot2::geom_text(
              ggplot2::aes(label = round(Individual, 1), 
                           hjust = ifelse(Individual >= 0, 1.2, -0.2)),
              color = colors$effect_text, 
              size = effect.cex
            )
        }
        
        p.hp <- p.hp + 
          ggplot2::geom_hline(yintercept = 0, color = colors$axis_text)
      } else {
        p.hp <- ggplot2::ggplot() + 
          ggplot2::annotate(
            "text", x = 1, y = 1, 
            label = "No hierarchical partitioning results meet the cutoff", 
            size = 6
          ) +
          ggplot2::theme_void()
      }
    } else {
      p.hp <- ggplot2::ggplot() + 
        ggplot2::annotate(
          "text", x = 1, y = 1, 
          label = "No hierarchical partitioning data available", 
          size = 6
        ) +
        ggplot2::theme_void()
    }
  }
  
  # Create matrix plot
  if (nrow(Var.part) > 0) {
    panel <- NULL
    if (plot.hp && exists("p.hp") && nrow(Hier.part) > 0) {
      Var <- Hier.part
    } else if (nrow(Var.exp) > 0) {
      Var <- Var.exp
    } else {
      Var <- data.frame(Var = character(0))
    }
    
    if (nrow(Var) > 0) {
      # Create all possible combinations
      panel <- expand.grid(X1 = Var$Var, X2 = Var.part$Var, stringsAsFactors = FALSE)
      panel$X3 <- "0"
      
      # Mark active connections
      for (i in 1:nrow(Var.part)) {
        i_var <- as.character(Var.part[i, "Var"])
        vars_in_component <- unlist(strsplit(i_var, ", "))
        panel$X3[panel$X1 %in% vars_in_component & panel$X2 == i_var] <- "1"
      }
      
      # Add row coloring
      panel$X4 <- ifelse(as.numeric(factor(panel$X1, levels = levels(panel$X1))) %% 2 == 0, 
                         "even", "odd")
      
      # Create matrix plot
      p.panel <- ggplot2::ggplot(panel, ggplot2::aes(x = X2, y = X1, color = X3, fill = X4)) +
        ggplot2::geom_tile(color = NA) +
        ggplot2::geom_point(size = pch.size) +
        ggplot2::scale_fill_manual(values = c(
          "even" = colors$tile_even, 
          "odd" = colors$tile_odd
        )) +
        ggplot2::scale_color_manual(values = c(
          "0" = colors$point_inactive, 
          "1" = colors$point_active
        )) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          panel.background = ggplot2::element_blank(),
          axis.text = ggplot2::element_text(color = colors$axis_text, size = title.cex),
          axis.ticks = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          legend.position = "none"
        ) +
        ggplot2::scale_x_discrete(expand = c(0, 0)) +
        ggplot2::scale_y_discrete(expand = c(0, 0)) +
        ggplot2::labs(y = NULL, x = NULL)
      
      # Add connecting lines
      for (i in levels(panel$X2)) {
        active_rows <- panel[panel$X2 == i & panel$X3 == "1", "X1"]
        if (length(active_rows) > 1) {
          active_rows <- sort(active_rows)
          p.panel <- p.panel + 
            ggplot2::annotate(
              "segment",
              x = i, xend = i,
              y = active_rows[1], yend = active_rows[length(active_rows)],
              color = colors$line_color,
              size = line.lwd
            )
        }
      }
    } else {
      p.panel <- ggplot2::ggplot() + 
        ggplot2::annotate(
          "text", x = 1, y = 1, 
          label = "No variables to display in the matrix", 
          size = 6
        ) +
        ggplot2::theme_void()
    }
  } else {
    p.panel <- ggplot2::ggplot() + 
      ggplot2::annotate(
        "text", x = 1, y = 1, 
        label = "Matrix plot disabled", 
        size = 6
      ) +
      ggplot2::theme_void()
  }
  
  # Image Layout Strategy
  design <- c(
    # Top plot (VP) - 对齐到右下方
    patchwork::area(1, 2, 2, 3),
    # Matrix plot (center)
    patchwork::area(3, 2, 4, 3),
    # Left plot (HP or Var.exp)
    patchwork::area(3, 1, 4, 2)
  )
  
  # Calculate the height and width of the plot based on proportions
  plot_heights <- c(height.ratio[1], 0.2, height.ratio[2]) 
  plot_widths <- c(width.ratio[1], width.ratio[2])
  
  # Ensure patchwork is loaded and Add spacing using patchwork::plot_spacer()
  if (plot.hp) {
    final_plot <- (patchwork::plot_spacer() | p.vp) / 
      (p.hp | p.panel) + 
      patchwork::plot_layout(
        heights = plot_heights,
        widths = plot_widths,
        design = design
      )
  } else {
    final_plot <- (patchwork::plot_spacer() | p.vp) / 
      (p.exp | p.panel) + 
      patchwork::plot_layout(
        heights = plot_heights,
        widths = plot_widths,
        design = design
      )
  }
  
  return(final_plot)
}
