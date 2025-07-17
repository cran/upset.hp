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
  
  colors <- color_palettes[[col]]
  
  process_vp_data <- function(vp) {
    if (nrow(vp) < 2) {
      warning("Variance partitioning data has insufficient rows.")
      return(NULL)
    }
    
    Constrained <- 100 * tail(vp[, 1], 1)
    Var.part <- as.data.frame(vp[-nrow(vp), , drop = FALSE])
    Var.part$Var <- gsub("^\\s+|\\s+$", "", 
                         gsub("and ", "", 
                              gsub("Common to ", "", 
                                   gsub("Unique to ", "", rownames(Var.part)))))
    Var.part$Fractions <- 100 * Var.part[, 1]
    Var.part$inter <- lengths(strsplit(as.character(Var.part$Var), ", "))
    Var.part$valid <- ifelse(Var.part$Fractions <= 0, "neg", "pos")
    Var.part <- Var.part[Var.part$Fractions >= 100 * cutoff, , drop = FALSE]
    
    if (nrow(Var.part) == 0) {
      warning("No variance components meet the cutoff criteria.")
      return(NULL)
    }
    
    # 排序逻辑
    if (order.part == "effect") {
      Var.part <- Var.part[order(Var.part$Fractions, decreasing = decreasing.part), ]
    } else if (order.part == "degree") {
      Var.part <- Var.part[order(Var.part$inter, Var.part$Fractions, 
                                 decreasing = c(!decreasing.part, TRUE)), ]
    }
    
    if (nrow(Var.part) > nVar) {
      Var.part <- Var.part[1:nVar, ]
    }
    
    Var.part$Var <- factor(Var.part$Var, levels = Var.part$Var)
    return(list(data = Var.part, constrained = Constrained))
  }
  
  # 3.2 处理HP数据
  process_hp_data <- function(hp) {
    if (nrow(hp) == 0) {
      warning("Hierarchical partitioning data is empty.")
      return(NULL)
    }
    
    Hier.part <- as.data.frame(hp)
    Hier.part$Var <- rownames(Hier.part)
    Hier.part$Individual <- 100 * Hier.part[, 1]
    Hier.part <- Hier.part[Hier.part$Individual >= 100 * cutoff, , drop = FALSE]
    
    if (nrow(Hier.part) == 0) {
      warning("No hierarchical components meet the cutoff criteria.")
      return(NULL)
    }
    
    if (order.var) {
      Hier.part <- Hier.part[order(Hier.part$Individual, decreasing = !decreasing.var), ]
    }
    
    Hier.part$Var <- factor(Hier.part$Var, levels = Hier.part$Var)
    Hier.part$valid <- ifelse(Hier.part$Individual <= 0, "neg", "pos")
    return(Hier.part)
  }
  
  # 4. 执行数据处理 ==========================================================
  vp_result <- process_vp_data(vp)
  if (is.null(vp_result)) {
    return(ggplot2::ggplot() + 
             ggplot2::annotate("text", x = 0.5, y = 0.5, 
                               label = "No variance components meet the cutoff", 
                               size = 6) + 
             ggplot2::theme_void())
  }
  
  Var.part <- vp_result$data
  Constrained <- vp_result$constrained
  
  # 5. 创建图形 ==============================================================
  # 5.1 创建VP图
  create_vp_plot <- function(data, colors) {
    bar_colors <- ifelse(data$valid == "pos", 
                         rep(colors$pos_bars, length.out = nrow(data)),
                         rep(colors$neg_bars, length.out = nrow(data)))
    names(bar_colors) <- as.character(data$Var)
    
    p <- ggplot2::ggplot(data, ggplot2::aes(x = Var, y = Fractions, fill = Var)) +
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
        expand = ggplot2::expansion(mult = c(ifelse(min(data$Fractions) < 0, 0.1, 0), 0.1))
      ) +
      ggplot2::labs(
        y = "Fractions (%)", 
        x = "", 
        title = bquote(R^2 == .(Constrained)*"%" ~ "  Residual" ~ .(100 - Constrained)*"%")
      )
    
    if (show.effect) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = round(Fractions, 1), 
                     vjust = ifelse(Fractions >= 0, -0.2, 1.2)),
        color = colors$effect_text, 
        size = effect.cex
      )
    }
    
    p + ggplot2::geom_hline(yintercept = 0, color = colors$axis_text)
  }
  
  p.vp <- create_vp_plot(Var.part, colors)
  
  # 5.2 创建EXP图
  create_exp_plot <- function(vp_data, hp_data, colors) {
    Fractions <- sapply(rownames(hp_data), function(i) {
      sum(vp_data[grep(i, vp_data$Var), "Fractions"])
    })
    
    Var.exp <- data.frame(Var = rownames(hp_data), Fractions = Fractions)
    Var.exp$valid <- ifelse(Var.exp$Fractions <= 0, "neg", "pos")
    Var.exp <- Var.exp[Var.exp$Fractions >= 100 * cutoff, , drop = FALSE]
    
    if (nrow(Var.exp) == 0) {
      warning("No explanatory variables meet the cutoff.")
      return(NULL)
    }
    
    if (order.var) {
      Var.exp <- Var.exp[order(Var.exp$Fractions, decreasing = !decreasing.var), ]
    }
    
    Var.exp$Var <- factor(Var.exp$Var, levels = Var.exp$Var)
    
    exp_colors <- ifelse(Var.exp$valid == "pos", 
                         rep(colors$pos_bars, length.out = nrow(Var.exp)),
                         rep(colors$neg_bars, length.out = nrow(Var.exp)))
    names(exp_colors) <- Var.exp$Var
    
    p <- ggplot2::ggplot(Var.exp, ggplot2::aes(x = Var, y = Fractions, fill = Var)) +
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
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = round(Fractions, 1), 
                     hjust = ifelse(Fractions >= 0, 1.2, -0.2)),
        color = colors$effect_text, 
        size = effect.cex
      )
    }
    
    p + ggplot2::geom_hline(yintercept = 0, color = colors$axis_text)
  }
  
  # 5.3 创建HP图
  create_hp_plot <- function(hp_data, colors) {
    hp_colors <- ifelse(hp_data$valid == "pos", 
                        rep(colors$pos_bars, length.out = nrow(hp_data)), 
                        rep(colors$neg_bars, length.out = nrow(hp_data)))
    names(hp_colors) <- hp_data$Var
    
    p <- ggplot2::ggplot(hp_data, ggplot2::aes(x = Var, y = Individual, fill = Var)) +
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
        expand = ggplot2::expansion(mult = c(0.3, ifelse(min(hp_data$Individual) < 0, 0.3, 0)))
      ) +
      ggplot2::labs(y = "Individual (%)", x = NULL)
    
    if (show.effect) {
      p <- p + ggplot2::geom_text(
        ggplot2::aes(label = round(Individual, 1), 
                     hjust = ifelse(Individual >= 0, 1.2, -0.2)),
        color = colors$effect_text, 
        size = effect.cex
      )
    }
    
    p + ggplot2::geom_hline(yintercept = 0, color = colors$axis_text)
  }
  
  # 5.4 创建矩阵图
  create_matrix_plot <- function(vp_data, hp_data, colors) {
    panel <- expand.grid(
      X1 = levels(hp_data$Var),
      X2 = levels(vp_data$Var),
      stringsAsFactors = FALSE
    )
    panel$X3 <- "0"
    
    for (i in 1:nrow(vp_data)) {
      current_var <- as.character(vp_data[i, "Var"])
      vars_in_combo <- unlist(strsplit(current_var, ", "))
      panel[panel$X1 %in% vars_in_combo & panel$X2 == current_var, "X3"] <- "1"
    }
    
    panel$X4 <- ifelse(as.numeric(factor(panel$X1, levels = levels(hp_data$Var))) %% 2 == 0, 
                       "even", "odd")
    
    p <- ggplot2::ggplot(panel, ggplot2::aes(x = X2, y = X1, color = X3, fill = X4)) +
      ggplot2::geom_tile(color = NA) +
      ggplot2::geom_point(size = pch.size) +
      ggplot2::scale_fill_manual(values = c(even = colors$tile_even, odd = colors$tile_odd)) +
      ggplot2::scale_color_manual(values = c("0" = colors$point_inactive, "1" = colors$point_active)) +
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
    
    # 添加连接线
    for (i in levels(panel$X2)) {
      panel_i <- panel[panel$X2 == i & panel$X3 == "1", ]
      panel_i <- panel_i[order(panel_i$X1), ]
      if (nrow(panel_i) > 1) {
        p <- p + ggplot2::annotate(
          "segment",
          x = i, xend = i,
          y = panel_i[1, "X1"], yend = panel_i[nrow(panel_i), "X1"],
          color = colors$line_color,
          size = line.lwd
        )
      }
    }
    
    return(p)
  }
  
  # 6. 图形组合 ==============================================================
  # 6.1 获取HP或EXP数据
  if (plot.hp) {
    hp_data <- process_hp_data(hp)
    if (is.null(hp_data)) {
      warning("Cannot create HP plot - using EXP plot instead")
      plot.hp <- FALSE
    }
  }
  
  if (!plot.hp) {
    exp_data <- hp  # 使用原始HP数据作为EXP数据
    p.exp <- create_exp_plot(Var.part, exp_data, colors)
    if (is.null(p.exp)) {
      warning("Cannot create EXP plot - showing VP plot only")
      return(p.vp)
    }
  } else {
    p.hp <- create_hp_plot(hp_data, colors)
  }
  
  # 6.2 创建矩阵图
  if (plot.hp) {
    p.panel <- create_matrix_plot(Var.part, hp_data, colors)
  } else {
    p.panel <- create_matrix_plot(Var.part, exp_data, colors)
  }
  
  # 6.3 定义布局
  left_margin <- 0.8  # 控制矩阵图向左移动的量
  
  design <- c(
    # VP图 (顶部右侧)
    patchwork::area(
      t = 1, 
      l = 1 + left_margin + width.ratio[1],  # 从左侧margin后开始
      b = 1 + height.ratio[1], 
      r = 1 + left_margin + width.ratio[1] + width.ratio[2]
    ),
    
    # EXP/HP图 (底部右侧)
    patchwork::area(
      t = 1 + height.ratio[1] + 1, 
      l = 1 + left_margin + width.ratio[1],  # 对齐VP图的左侧
      b = 1 + height.ratio[1] + height.ratio[2], 
      r = 1 + left_margin + width.ratio[1] + width.ratio[2]
    ),
    
    # 矩阵图 (底部左侧) - 向左移动
    patchwork::area(
      t = 1 + height.ratio[1] + 1, 
      l = 1,  # 从最左侧开始
      b = 1 + height.ratio[1] + height.ratio[2], 
      r = 1 + left_margin  # 控制矩阵图宽度
    )
  )
  
  # 调整整体宽度比例
  adjusted_widths <- c(left_margin, width.ratio[1], width.ratio[2])
  
  # 7. 最终组合
  if (plot.hp) {
    final_plot <- p.vp + p.panel + p.hp + 
      patchwork::plot_layout(
        design = design,
        heights = c(height.ratio[1], 0.2, height.ratio[2]),
        widths = adjusted_widths
      )
  } else {
    final_plot <- p.vp + p.panel + p.exp + 
      patchwork::plot_layout(
        design = design,
        heights = c(height.ratio[1], 0.2, height.ratio[2]),
        widths = adjusted_widths
      )
  }
  
  return(final_plot)
}

