# =============================================================================
# Network Clinical VIF-Pruned Regression Analysis & Visualization
# =============================================================================
# 前置条件：data_merge 已在环境中准备好，直接调用
# 输出：
#   - 每个 outcome 对应一份 DOCX 报告（包含简单斜率检验与鲁棒性检验图）
#   - 汇总 CSV 文件（模型拟合、系数、显著结果、VIF 历史等）
#   - 显著效应的可视化图片（PNG）和鲁棒性检验图片（PNG）
# 核心逻辑：总体模型 F 检验的 p 值需在对应家族内通过 FDR 校正 (<0.05)，才进行效应绘图
# =============================================================================

# -----------------------------------------------------------------------------
# 0. 加载依赖包
# -----------------------------------------------------------------------------
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)
library(broom)
library(flextable)
library(officer)
library(tibble)
library(readr)
library(car)
library(visreg)
library(ggplot2)
# 科学分析验证包
library(emmeans)     # 用于计算事后简单斜率
library(performance) # 用于检验模型假设 (check_model)
library(see)         # performance 绘图依赖
library(patchwork)   # performance 绘图依赖

# -----------------------------------------------------------------------------
# 1. 全局配置
# -----------------------------------------------------------------------------

# 输出根目录
output_dir <- "D:/youyi_fucha/code/network_clinical_vif_models"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# PA 组专用 outcomes（仅患者组回归）
pa_group_only <- c(
  "overeatting_time_perweek",
  "PA_time","DEBQ_33_constrain",
  "DEBQ_33_emotion",
  "DEBQ_33_external",
  "EDI_BN",
  "EAT_26",
  "BDI_21",
  "SAS_20",
  "EAT_Dieting",
  "EAT_Bulimia_and_Food_preoccupation",
  "EAT_Oral_control"
)

# 调节分析 outcomes（全样本主效应 及 全样本 × 组别交互）
interaction_var <- c(
  "DEBQ_33_constrain",
  "DEBQ_33_emotion",
  "DEBQ_33_external",
  "EDI_BN",
  "EAT_26",
  "BDI_21",
  "SAS_20",
  "EAT_Dieting",
  "EAT_Bulimia_and_Food_preoccupation",
  "EAT_Oral_control"
)

# 分组标签统一映射
patient_labels <- c("PA", "SUB", "bulimia", "patient")
control_labels <- c("CON", "HC", "control", "healthy")

# VIF 阈值与最少保留指标数
vif_threshold    <- 5
min_metric_keep  <- 1

# 图形输出子目录
plot_dir      <- file.path(output_dir, "significant_effect_plots")
plot_dir_pa   <- file.path(plot_dir, "PA_only_main_effects")
plot_dir_full <- file.path(plot_dir, "Full_sample_main_effects") 
plot_dir_int  <- file.path(plot_dir, "interaction_effects")
check_dir     <- file.path(output_dir, "model_assumption_checks") 

for (d in c(plot_dir, plot_dir_pa, plot_dir_full, plot_dir_int, check_dir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# =============================================================================
# 2. 辅助函数
# =============================================================================

# --- 2.1 统计输出格式化 ------------------------------------------------------

cleanPValue <- function(p_value) {
  case_when(
    is.na(p_value) ~ NA_character_,
    p_value < 0.001 ~ "< .001",
    TRUE ~ sprintf("%.3f", p_value)
  )
}

getSigMark <- function(p_value) {
  case_when(
    is.na(p_value) ~ "",
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    p_value < 0.10  ~ "†",
    TRUE            ~ ""
  )
}

formatCoefTable <- function(df_input) {
  df_input %>%
    mutate(
      estimate  = round(estimate,   3),
      std.error = round(std.error,  3),
      statistic = round(statistic,  3),
      conf.low  = round(conf.low,   3),
      conf.high = round(conf.high,  3),
      p         = cleanPValue(p.value),
      sig       = getSigMark(p.value),
      ci_95     = paste0("[", conf.low, ", ", conf.high, "]")
    )
}

extractModelFit <- function(model_fit, model_type, outcome_name, n_obs) {
  if (is.null(model_fit)) {
    return(tibble(
      model_type    = model_type,
      outcome       = outcome_name,
      n             = n_obs,
      r_squared     = NA_real_,
      adj_r_squared = NA_real_,
      statistic     = NA_real_,
      model_p       = NA_real_,
      aic           = NA_real_,
      bic           = NA_real_
    ))
  }
  fit_glance <- broom::glance(model_fit)
  tibble(
    model_type    = model_type,
    outcome       = outcome_name,
    n             = n_obs,
    r_squared     = fit_glance$r.squared,
    adj_r_squared = fit_glance$adj.r.squared,
    statistic     = fit_glance$statistic,
    model_p       = fit_glance$p.value,
    aic           = fit_glance$AIC,
    bic           = fit_glance$BIC
  )
}

# --- 2.2 数据与变量处理 -------------------------------------------------------

standardizeGroup <- function(group_vec) {
  case_when(
    group_vec %in% patient_labels ~ "Patient",
    group_vec %in% control_labels ~ "Control",
    TRUE ~ as.character(group_vec)
  )
}

makeSafeMetricName <- function(x) make.names(x)

# --- 2.3 flextable 主题 -------------------------------------------------------

makeFtTheme <- function(ft_obj, caption_text = NULL) {
  ft_obj <- ft_obj %>%
    bold(part = "header") %>%
    bg(part = "header", bg = "#2C5F8A") %>%
    color(part = "header", color = "white") %>%
    align(align = "center", part = "all") %>%
    fontsize(size = 10, part = "all") %>%
    font(fontname = "Times New Roman", part = "all") %>%
    border_outer(part = "all", border = fp_border(color = "#2C5F8A", width = 1.2)) %>%
    border_inner_h(part = "body", border = fp_border(color = "#BBBBBB", width = 0.5)) %>%
    set_table_properties(layout = "autofit", width = 1)
  
  if (!is.null(caption_text)) ft_obj <- ft_obj %>% set_caption(caption_text)
  ft_obj
}

# --- 2.4 绘图辅助 -------------------------------------------------------------

prettyMetricName  <- function(x) x %>% str_replace_all("\\.", "_") %>% str_replace_all("_", " ") %>% str_to_title()
prettyOutcomeName <- function(x) x %>% str_replace_all("_", " ") %>% str_to_title()
cleanFilename     <- function(x) x %>% str_replace_all("[^A-Za-z0-9_]+", "_") %>% str_replace_all("_+", "_") %>% str_replace_all("^_|_$", "")

# =============================================================================
# 3. 数据准备：长表转宽表
# =============================================================================

# 确认必需变量存在于 data_merge
required_vars <- unique(c(
  "participant", "group", "Metrics", "Degree", "BMI", "years",
  pa_group_only, interaction_var
))
missing_vars <- setdiff(required_vars, names(data_merge))
if (length(missing_vars) > 0) {
  stop(paste("data_merge 中缺少以下变量：", paste(missing_vars, collapse = ", ")))
}

# 获取显著指标列表（来自 FDR 校正后的显著性结果）
metric_sig <- unique(significant_metrics_fdr$Metrics)
data_merge$overeatting_time_perweek <- as.numeric(data_merge$overeatting_time_perweek )
# 长表：筛选显著指标，生成安全列名与标准化分组
data_metric_long <- data_merge %>%
  filter(Metrics %in% metric_sig) %>%
  mutate(
    metric_col = makeSafeMetricName(Metrics),
    group_std  = standardizeGroup(group)
  )

# 提取受试者基本信息（每人一行）
subject_base <- data_metric_long %>%
  select(participant, group, group_std, BMI, years,
         all_of(pa_group_only), all_of(interaction_var)) %>%
  distinct(participant, group, .keep_all = TRUE) %>%
  mutate(group_std = factor(group_std, levels = c("Control", "Patient")))

# 将 Degree 宽转：每个指标一列（多次测量取均值）
metric_wide <- data_metric_long %>%
  select(participant, metric_col, Degree) %>%
  group_by(participant, metric_col) %>%
  summarise(Degree = mean(Degree, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = metric_col, values_from = Degree)

# 合并为分析宽表
data_wide <- subject_base %>%
  left_join(metric_wide, by = "participant")

# 对所有指标列进行 z 标准化
metric_cols_safe      <- makeSafeMetricName(metric_sig)
existing_metric_cols  <- metric_cols_safe[metric_cols_safe %in% names(data_wide)]

data_wide <- data_wide %>%
  mutate(across(all_of(existing_metric_cols), ~ as.numeric(scale(.x))))

cat("宽表维度：", dim(data_wide), "\n")
cat("指标预测变量数量：", length(existing_metric_cols), "\n")

# =============================================================================
# 4. 公式构造器
# =============================================================================

makePaFormula <- function(outcome_name, metric_terms) {
  rhs <- paste(c(metric_terms, "BMI", "years"), collapse = " + ")
  as.formula(paste0(outcome_name, " ~ ", rhs))
}

makeFullSampleFormula <- function(outcome_name, metric_terms) {
  rhs <- paste(c(metric_terms, "group_std", "BMI", "years"), collapse = " + ")
  as.formula(paste0(outcome_name, " ~ ", rhs))
}

makeModerationFormula <- function(outcome_name, metric_terms) {
  interaction_terms <- paste0(metric_terms, ":group_std")
  rhs <- paste(c(metric_terms, "group_std", interaction_terms, "BMI", "years"), collapse = " + ")
  as.formula(paste0(outcome_name, " ~ ", rhs))
}

# =============================================================================
# 5. VIF 迭代剪枝
# =============================================================================

# 计算模型矩阵列级 VIF
computeColumnVIF <- function(model_fit) {
  if (is.null(model_fit)) return(tibble(term = character(), vif = numeric()))
  
  mm <- model.matrix(model_fit)
  if ("(Intercept)" %in% colnames(mm)) mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (ncol(mm) < 2) return(tibble(term = colnames(mm), vif = 1))
  
  vif_values <- sapply(seq_len(ncol(mm)), function(i) {
    y_col  <- mm[, i]
    x_cols <- mm[, -i, drop = FALSE]
    fit_i  <- lm(y_col ~ x_cols) 
    r2 <- summary(fit_i)$r.squared
    if (is.na(r2) || r2 >= 0.999999) return(Inf)
    1 / (1 - r2)
  })
  
  tibble(term = colnames(mm), vif = as.numeric(vif_values))
}

# 将列级 VIF 聚合为指标块级 VIF（主效应 + 交互项取最大值）
computeMetricBlockVIF <- function(vif_table, metric_terms, model_type) {
  if (nrow(vif_table) == 0) return(tibble(metric = character(), block_vif = numeric()))
  
  lapply(metric_terms, function(metric_name) {
    relevant_terms <- if (model_type %in% c("PA_only", "Full_Sample")) {
      metric_name
    } else {
      c(metric_name,
        paste0(metric_name, ":group_stdPatient"),
        paste0("group_stdPatient:", metric_name))
    }
    tmp <- vif_table %>% filter(term %in% relevant_terms)
    tibble(metric = metric_name,
           block_vif = ifelse(nrow(tmp) == 0, NA_real_, max(tmp$vif, na.rm = TRUE)))
  }) %>% bind_rows()
}

# 主迭代剪枝函数
runModelWithVIFPruning <- function(model_data, outcome_name, metric_terms,
                                   model_type = c("PA_only", "Moderation", "Full_Sample"),
                                   vif_threshold = 10, min_metric_keep = 1) {
  model_type        <- match.arg(model_type)
  current_metrics   <- metric_terms
  iteration_history <- list()
  removed_metrics   <- character()
  final_fit         <- NULL
  final_formula     <- NULL
  iter              <- 1
  
  repeat {
    if (length(current_metrics) < min_metric_keep) break
    
    model_formula <- if (model_type == "PA_only") {
      makePaFormula(outcome_name, current_metrics)
    } else if (model_type == "Full_Sample") {
      makeFullSampleFormula(outcome_name, current_metrics)
    } else {
      makeModerationFormula(outcome_name, current_metrics)
    }
    
    fit_obj <- lm(model_formula, data = model_data) 
    
    vif_col_tbl   <- computeColumnVIF(fit_obj)
    vif_block_tbl <- computeMetricBlockVIF(vif_col_tbl, current_metrics, model_type) %>%
      arrange(desc(block_vif))
    
    max_vif        <- if (nrow(vif_block_tbl) == 0) NA_real_ else vif_block_tbl$block_vif[1]
    remove_metric  <- NA_character_
    
    if (!is.na(max_vif) && max_vif > vif_threshold && length(current_metrics) > min_metric_keep) {
      remove_metric   <- vif_block_tbl$metric[1]
      current_metrics <- setdiff(current_metrics, remove_metric)
      removed_metrics <- c(removed_metrics, remove_metric)
    } else {
      final_fit     <- fit_obj
      final_formula <- model_formula
    }
    
    iteration_history[[iter]] <- vif_block_tbl %>%
      mutate(outcome = outcome_name, model_type = model_type,
             iteration = iter, threshold = vif_threshold,
             removed_metric = remove_metric)
    
    should_continue <- !is.na(remove_metric) &&
      max_vif > vif_threshold &&
      length(current_metrics) >= min_metric_keep
    
    if (should_continue) { iter <- iter + 1; next } else break
  }
  
  if (is.null(final_fit) && length(current_metrics) >= 1) {
    final_formula <- if (model_type == "PA_only") {
      makePaFormula(outcome_name, current_metrics)
    } else if (model_type == "Full_Sample") {
      makeFullSampleFormula(outcome_name, current_metrics)
    } else {
      makeModerationFormula(outcome_name, current_metrics)
    }
    final_fit <- lm(final_formula, data = model_data)
  }
  
  list(
    fit             = final_fit,
    formula         = final_formula,
    kept_metrics    = current_metrics,
    removed_metrics = removed_metrics,
    vif_history     = bind_rows(iteration_history)
  )
}

# =============================================================================
# 6. 单个 outcome 的建模主函数（加入事后检验及鲁棒性检验制图）
# =============================================================================

runSingleOutcomeModel <- function(outcome_name, model_type,
                                  data_wide, existing_metric_cols,
                                  vif_threshold = 5, min_metric_keep = 1) {
  # --- 按模型类型准备分析数据 ---
  if (model_type == "PA_only") {
    model_data <- data_wide %>%
      filter(group_std == "Patient") %>%
      select(all_of(c(outcome_name, "BMI", "years", existing_metric_cols))) %>%
      drop_na()
  } else { # Moderation 和 Full_Sample 都是全数据
    model_data <- data_wide %>%
      select(all_of(c(outcome_name, "group_std", "BMI", "years", existing_metric_cols))) %>%
      drop_na()
  }
  
  n_obs <- nrow(model_data)
  
  # --- 样本量不足时报错并停止 ---
  insufficient <- (model_type == "PA_only" && n_obs < 15) ||
    (model_type %in% c("Moderation", "Full_Sample") && (n_obs < 20 || length(unique(model_data$group_std)) < 2))
  
  if (insufficient) {
    stop(paste("数据样本量不足以运行模型：", outcome_name, "-", model_type))
  }
  
  # --- VIF 剪枝建模 ---
  vif_result <- runModelWithVIFPruning(
    model_data      = model_data,
    outcome_name    = outcome_name,
    metric_terms    = existing_metric_cols,
    model_type      = model_type,
    vif_threshold   = vif_threshold,
    min_metric_keep = min_metric_keep
  )
  
  fit_obj <- vif_result$fit
  
  # --- 鲁棒性检验 (Assumption Checking) ---
  cat("  -> 正在进行模型鲁棒性检验并保存图像: ", outcome_name, " (", model_type, ")\n", sep="")
  model_check_res <- performance::check_model(fit_obj)
  check_plot <- plot(model_check_res)
  check_filename <- file.path(check_dir, paste0("Check_", model_type, "_", cleanFilename(outcome_name), ".png"))
  ggsave(check_filename, plot = check_plot, width = 12, height = 10, dpi = 300, bg = "white")
  
  # --- 提取完整系数表 ---
  coef_full <- broom::tidy(fit_obj, conf.int = TRUE) %>%
    mutate(model_type = model_type, outcome = outcome_name)
  
  # --- 提取关键效应 ---
  key_effects <- if (model_type %in% c("PA_only", "Full_Sample")) {
    coef_full %>% filter(term %in% vif_result$kept_metrics)
  } else {
    int_terms_a <- paste0(vif_result$kept_metrics, ":group_stdPatient")
    int_terms_b <- paste0("group_stdPatient:", vif_result$kept_metrics)
    coef_full %>% filter(term %in% c(int_terms_a, int_terms_b))
  }
  
  # --- 事后检验 (Simple Slopes) ---
  simple_slopes_df <- tibble()
  if (model_type == "Moderation" && nrow(key_effects) > 0) {
    sig_ints <- key_effects %>% filter(p.value < 0.05)
    if (nrow(sig_ints) > 0) {
      sig_metrics <- unique(str_remove_all(sig_ints$term, ":group_stdPatient|group_stdPatient:"))
      slopes_list <- map(sig_metrics, function(m) {
        emt <- emtrends(fit_obj, specs = "group_std", var = m)
        res <- as.data.frame(test(emt))
        colnames(res)[1:2] <- c("Group", "Trend")
        res$Metric <- m
        res
      })
      simple_slopes_df <- bind_rows(slopes_list) %>%
        select(Metric, Group, Trend, SE = SE, df, t.ratio, p.value)
    }
  }
  
  list(
    outcome         = outcome_name,
    model_type      = model_type,
    n               = n_obs,
    fit             = fit_obj,
    formula         = vif_result$formula,
    kept_metrics    = vif_result$kept_metrics,
    removed_metrics = vif_result$removed_metrics,
    vif_history     = vif_result$vif_history,
    coef_full       = coef_full,
    key_effects     = key_effects,
    simple_slopes   = simple_slopes_df,
    check_plot_file = check_filename,
    fit_summary     = extractModelFit(fit_obj, model_type, outcome_name, n_obs)
  )
}

# =============================================================================
# 7. 批量运行全部模型 
# =============================================================================

all_outcomes <- c(
  pa_group_only, 
  interaction_var, 
  interaction_var
)
all_model_types <- c(
  rep("PA_only",     length(pa_group_only)),
  rep("Moderation",  length(interaction_var)),
  rep("Full_Sample", length(interaction_var))
)

cat("开始运行", length(all_outcomes), "个回归模型（含验证与制图）...\n")

analysis_results <- map2(
  .x = all_outcomes,
  .y = all_model_types,
  .f = ~ runSingleOutcomeModel(
    outcome_name         = .x,
    model_type           = .y,
    data_wide            = data_wide,
    existing_metric_cols = existing_metric_cols,
    vif_threshold        = vif_threshold,
    min_metric_keep      = min_metric_keep
  )
)

names(analysis_results) <- paste(all_outcomes, all_model_types, sep = "___")
cat("模型运行及鲁棒性检验完毕。\n\n")

# =============================================================================
# 8. 汇总结果与 FDR 校正（门控逻辑）
# =============================================================================

# 【核心更新】针对总体模型 F 检验的 p 值（model_p），在其各自模型类型内部进行 FDR 校正
all_fit_summary  <- bind_rows(map(analysis_results, "fit_summary")) %>%
  group_by(model_type) %>%
  mutate(model_p_fdr = p.adjust(model_p, method = "fdr")) %>%
  ungroup()

# 提取 F 检验过校正门槛（model_p_fdr < 0.05）的“合规模型”
valid_models <- all_fit_summary %>%
  filter(model_p_fdr < 0.05) %>%
  select(model_type, outcome)

all_coef_full    <- bind_rows(map(analysis_results, "coef_full"))
all_vif_history  <- bind_rows(map(analysis_results, "vif_history"))

all_key_effects <- bind_rows(map(analysis_results, "key_effects")) %>%
  formatCoefTable() %>%
  mutate(
    Metric = case_when(
      model_type %in% c("PA_only", "Full_Sample") ~ term,
      TRUE ~ term %>%
        str_remove(":group_stdPatient") %>%
        str_remove("group_stdPatient:")
    )
  ) %>%
  transmute(
    Model    = model_type,
    Outcome  = outcome,
    Metric   = Metric,
    Term     = term,
    Beta     = estimate,
    SE       = std.error,
    t        = statistic,
    p        = p,
    Sig      = sig,
    `95% CI` = ci_95,
    p_raw    = p.value
  ) %>%
  arrange(Model, Outcome, Metric)

model_metric_summary <- bind_rows(lapply(analysis_results, function(x) {
  tibble(
    Outcome         = x$outcome,
    Model           = x$model_type,
    N               = x$n,
    Kept_Metrics    = paste(x$kept_metrics,    collapse = "; "),
    Removed_Metrics = paste(x$removed_metrics, collapse = "; "),
    Final_Formula   = ifelse(
      length(x$formula) == 0 || all(is.na(x$formula)), NA_character_, deparse(x$formula)
    )
  )
}))

# --- 【核心更新】绘图门控：仅选取通过了总体 F 检验 FDR 校正（门控拦截）且自身效应也显著（p_raw < 0.05）的指标 ---
significant_results <- all_key_effects %>%
  filter(!is.na(p_raw), p_raw < 0.05) %>%
  inner_join(valid_models, by = c("Model" = "model_type", "Outcome" = "outcome")) %>%
  arrange(Model, Outcome, p_raw)

significant_pa_results <- significant_results %>% filter(Model == "PA_only")
significant_fs_results <- significant_results %>% filter(Model == "Full_Sample")
significant_interaction_results <- significant_results %>% filter(Model == "Moderation")

cat("F 检验通过校正且显著的 PA 主效应数量：",       nrow(significant_pa_results),          "\n")
cat("F 检验通过校正且显著的 全样本主效应数量：",    nrow(significant_fs_results),          "\n")
cat("F 检验通过校正且显著的 组别交互效应数量：",    nrow(significant_interaction_results), "\n\n")

# =============================================================================
# 9. 保存 CSV 汇总文件
# =============================================================================

write.csv(all_fit_summary,       file.path(output_dir, "All_model_fit_summary.csv"),          row.names = FALSE, fileEncoding = "UTF-8")
write.csv(all_coef_full,         file.path(output_dir, "All_full_coefficients.csv"),            row.names = FALSE, fileEncoding = "UTF-8")
write.csv(all_key_effects,       file.path(output_dir, "All_key_effects.csv"),                  row.names = FALSE, fileEncoding = "UTF-8")
write.csv(significant_results,   file.path(output_dir, "Significant_results_only.csv"),         row.names = FALSE, fileEncoding = "UTF-8")
write.csv(all_vif_history,       file.path(output_dir, "All_VIF_iteration_history.csv"),        row.names = FALSE, fileEncoding = "UTF-8")
write.csv(model_metric_summary,  file.path(output_dir, "Model_metric_retention_summary.csv"),   row.names = FALSE, fileEncoding = "UTF-8")

cat("CSV 文件保存完毕。\n\n")

# =============================================================================
# 10. 每个模型生成独立 DOCX 报告
# =============================================================================

buildDocxReport <- function(res_obj, vif_threshold, fit_summary_df) {
  
  outcome_name <- res_obj$outcome
  model_type   <- res_obj$model_type
  fit_obj      <- res_obj$fit
  
  # 获取该模型 FDR 校正后的 overall p 值
  current_fdr_p <- fit_summary_df %>%
    filter(model_type == res_obj$model_type, outcome == res_obj$outcome) %>%
    pull(model_p_fdr)
  
  # -- 各子表 --
  fit_tbl <- res_obj$fit_summary %>%
    mutate(`R²` = round(r_squared, 3), `Adj. R²` = round(adj_r_squared, 3),
           `F Statistic` = round(statistic, 3), 
           `Model p` = cleanPValue(model_p),
           `Model p (FDR)` = cleanPValue(current_fdr_p), # 【核心更新】插入 FDR P 值
           AIC = round(aic, 2), BIC = round(bic, 2)) %>%
    transmute(Model = model_type, Outcome = outcome, N = n,
              `R²`, `Adj. R²`, `F Statistic`, `Model p`, `Model p (FDR)`, AIC, BIC)
  
  coef_tbl <- formatCoefTable(res_obj$coef_full) %>%
    transmute(Term = term, Beta = estimate, SE = std.error,
              t = statistic, p = p, Sig = sig, `95% CI` = ci_95)
  
  key_tbl <- formatCoefTable(res_obj$key_effects) %>%
    mutate(Metric = case_when(
      model_type %in% c("PA_only", "Full_Sample") ~ term,
      TRUE ~ term %>% str_remove(":group_stdPatient") %>% str_remove("group_stdPatient:")
    )) %>%
    transmute(Metric, Term = term, Beta = estimate, SE = std.error,
              t = statistic, p = p, Sig = sig, `95% CI` = ci_95)
  
  vif_tbl <- if (nrow(res_obj$vif_history) == 0) {
    tibble(Iteration = NA_integer_, Metric = NA_character_, `Block VIF` = NA_real_,
           Threshold = NA_real_, `Removed This Round` = NA_character_)
  } else {
    res_obj$vif_history %>%
      mutate(block_vif = round(block_vif, 3)) %>%
      transmute(Iteration = iteration, Metric = metric, `Block VIF` = block_vif,
                Threshold = threshold, `Removed This Round` = removed_metric)
  }
  
  metric_tbl <- tibble(
    Item  = c("Kept Metrics", "Removed Metrics", "Final Formula"),
    Value = c(
      paste(res_obj$kept_metrics,    collapse = "; "),
      paste(res_obj$removed_metrics, collapse = "; "),
      ifelse(length(res_obj$formula) == 0 || all(is.na(res_obj$formula)),
             NA_character_, deparse(res_obj$formula))
    )
  )
  
  narrative_lines <- c(
    paste0("Outcome: ", outcome_name),
    paste0("Model type: ", model_type),
    paste0("Sample size used: ", res_obj$n),
    paste0("Initial number of metric predictors: ", length(existing_metric_cols)),
    paste0("Final number of retained metric predictors: ", length(res_obj$kept_metrics)),
    paste0("Number of removed metrics due to VIF pruning: ", length(res_obj$removed_metrics)),
    paste0("VIF threshold: ", vif_threshold),
    "",
    case_when(
      model_type == "PA_only" ~ "Interpretation focus: metric main effects within the patient group.",
      model_type == "Full_Sample" ~ "Interpretation focus: metric main effects in the pooled full sample, controlling for group differences.",
      model_type == "Moderation" ~ "Interpretation focus: metric × group interaction effects in the full sample."
    )
  )
  
  # -- 构建文档 --
  doc_obj <- read_docx() %>%
    body_add_par(paste0("Regression Report: ", outcome_name, " (", model_type, ")"), style = "heading 1") %>%
    body_add_par(format(Sys.time(), "Generated on %Y-%m-%d %H:%M:%S"), style = "Normal") %>%
    body_add_par("1. Analysis Summary", style = "heading 2")
  
  for (txt in narrative_lines) doc_obj <- body_add_par(doc_obj, txt, style = "Normal")
  
  doc_obj <- doc_obj %>%
    body_add_par("2. Model Fit", style = "heading 2") %>%
    body_add_flextable(flextable(fit_tbl) %>% makeFtTheme("Table 1. Model fit summary")) %>%
    body_add_par("3. Predictor Retention", style = "heading 2") %>%
    body_add_flextable(flextable(metric_tbl) %>% makeFtTheme("Table 2. Predictor retention summary")) %>%
    body_add_par("4. VIF Pruning History", style = "heading 2") %>%
    body_add_flextable(flextable(vif_tbl) %>% makeFtTheme("Table 3. VIF pruning history")) %>%
    body_add_par("5. Key Effects", style = "heading 2") %>%
    body_add_flextable(flextable(key_tbl) %>% makeFtTheme(
      case_when(
        model_type == "PA_only" ~ "Table 4. Key metric main effects (Patient)",
        model_type == "Full_Sample" ~ "Table 4. Key metric main effects (Full Sample)",
        model_type == "Moderation" ~ "Table 4. Key metric × group interaction effects"
      ))) %>%
    body_add_par("6. Full Coefficient Table", style = "heading 2") %>%
    body_add_flextable(flextable(coef_tbl) %>% makeFtTheme("Table 5. Full coefficient table"))
  
  # --- 插入简单斜率结果 ---
  if (model_type == "Moderation" && nrow(res_obj$simple_slopes) > 0) {
    slope_tbl <- res_obj$simple_slopes %>%
      mutate(Trend = round(Trend, 3), SE = round(SE, 3), t.ratio = round(t.ratio, 3),
             p.value = cleanPValue(p.value), Sig = getSigMark(p.value)) %>%
      rename(`p-value` = p.value)
    
    doc_obj <- doc_obj %>%
      body_add_par("7. Post-Hoc Simple Slopes Analysis", style = "heading 2") %>%
      body_add_par("Simple slopes calculated for metrics with significant interaction effects.", style = "Normal") %>%
      body_add_flextable(flextable(slope_tbl) %>% makeFtTheme("Table 6. Simple slopes by Group"))
  }
  
  # --- 插入鲁棒性检验图 ---
  if (file.exists(res_obj$check_plot_file)) {
    doc_obj <- doc_obj %>%
      body_add_par("Model Assumption Checks", style = "heading 2") %>%
      body_add_img(src = res_obj$check_plot_file, width = 6.5, height = 5.5)
  }
  
  out_path <- file.path(output_dir, paste0("Regression_Report_", model_type, "_", outcome_name, ".docx"))
  print(doc_obj, target = out_path)
  invisible(out_path)
}

cat("开始生成 DOCX 报告...\n")
walk(analysis_results, ~ buildDocxReport(.x, vif_threshold, all_fit_summary))
cat("DOCX 报告生成完毕。\n\n")

# =============================================================================
# 11. 绘图函数
# =============================================================================

getModelData <- function(outcome_name, model_type, kept_metrics, data_wide) {
  if (model_type == "PA_only") {
    data_wide %>%
      filter(group_std == "Patient") %>%
      select(all_of(c(outcome_name, "BMI", "years", kept_metrics))) %>%
      drop_na()
  } else {
    data_wide %>%
      select(all_of(c(outcome_name, "group_std", "BMI", "years", kept_metrics))) %>%
      drop_na()
  }
}

# 主效应绘图通用函数
plotMainEffect <- function(outcome_name, metric_name, fit_obj, data_plot, out_file, title_sub) {
  vis_obj <- visreg(fit_obj, data = data_plot, xvar = metric_name,
                    gg = TRUE, partial = FALSE, rug = FALSE, overlay = FALSE)
  
  p <- vis_obj +
    geom_point(data = data_plot,
               aes_string(x = metric_name, y = outcome_name),
               inherit.aes = FALSE, color = "#619bc1", alpha = 0.65, size = 4) +
    labs(
      title    = paste0(prettyOutcomeName(outcome_name), " ~ ", prettyMetricName(metric_name)),
      subtitle = title_sub,
      x        = prettyMetricName(metric_name),
      y        = prettyOutcomeName(outcome_name),
      caption  = "Line shows fitted association after covariate adjustment (BMI, years)."
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption  = element_text(size = 10, colour = "gray40"))
  
  ggsave(out_file, plot = p, width = 6.8, height = 5.2, dpi = 300, bg = "white")
  invisible(p)
}

plotInteractionEffect <- function(outcome_name, metric_name, fit_obj, data_plot, out_file) {
  data_plot <- data_plot %>%
    mutate(group_std = factor(group_std, levels = c("Control", "Patient")))
  
  vis_obj <- visreg(fit_obj, data = data_plot, xvar = metric_name,
                    by = "group_std", overlay = TRUE,
                    gg = TRUE, partial = FALSE, rug = FALSE)
  
  p <- vis_obj +
    geom_point(data = data_plot,
               aes_string(x = metric_name, y = outcome_name, color = "group_std"),
               inherit.aes = FALSE, alpha = 0.60, size = 3.0) +
    scale_color_manual(
      values = c("Control" = "#e88936", "Patient" = "#619bc1"),
      labels = c("Control" = "CON",     "Patient" = "Patient"),
      name   = "Group"
    ) +
    scale_fill_manual(
      values = c("Control" = alpha("#e88936", 0.2), "Patient" = alpha("#619bc1", 0.2)),
      guide  = "none"
    ) +
    labs(
      title    = paste0(prettyOutcomeName(outcome_name), " ~ ", prettyMetricName(metric_name), " × Group"),
      subtitle = "Interaction plot from final VIF-pruned model",
      x        = prettyMetricName(metric_name),
      y        = prettyOutcomeName(outcome_name),
      caption  = "Lines show fitted slopes after covariate adjustment (BMI, years)."
    ) +
    theme_classic(base_size = 13) +
    theme(plot.title    = element_text(face = "bold", hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          plot.caption  = element_text(size = 10, colour = "gray40"),
          legend.position = "top")
  
  ggsave(out_file, plot = p, width = 7.2, height = 5.4, dpi = 300, bg = "white")
  invisible(p)
}

# =============================================================================
# 12. 批量绘图（基于通过了模型FDR校正的显著结果）
# =============================================================================

# -- 12.1 PA 主效应图 ----------------------------------------------------------
cat("开始绘制显著 PA 主效应图...\n")
pa_plot_log <- if (nrow(significant_pa_results) > 0) {
  map_dfr(seq_len(nrow(significant_pa_results)), function(i) {
    outcome_name <- significant_pa_results$Outcome[i]
    metric_name  <- significant_pa_results$Metric[i]
    list_key     <- paste0(outcome_name, "___PA_only")
    res_i        <- analysis_results[[list_key]]
    
    data_plot <- getModelData(outcome_name, "PA_only", res_i$kept_metrics, data_wide)
    out_file  <- file.path(plot_dir_pa, paste0("PAonly_", cleanFilename(outcome_name), "_", cleanFilename(metric_name), ".png"))
    
    plotMainEffect(outcome_name, metric_name, res_i$fit, data_plot, out_file, "Patient group only")
    tibble(Plot_Type = "PA_only_main_effect", Outcome = outcome_name, Metric = metric_name, File = out_file)
  })
} else { tibble() }

# -- 12.2 全样本主效应图 -------------------------------------------------------
cat("开始绘制显著 全样本主效应图...\n")
fs_plot_log <- if (nrow(significant_fs_results) > 0) {
  map_dfr(seq_len(nrow(significant_fs_results)), function(i) {
    outcome_name <- significant_fs_results$Outcome[i]
    metric_name  <- significant_fs_results$Metric[i]
    list_key     <- paste0(outcome_name, "___Full_Sample")
    res_i        <- analysis_results[[list_key]]
    
    data_plot <- getModelData(outcome_name, "Full_Sample", res_i$kept_metrics, data_wide)
    out_file  <- file.path(plot_dir_full, paste0("FullSample_", cleanFilename(outcome_name), "_", cleanFilename(metric_name), ".png"))
    
    plotMainEffect(outcome_name, metric_name, res_i$fit, data_plot, out_file, "Pooled Full Sample (Controlling for Group)")
    tibble(Plot_Type = "Full_Sample_main_effect", Outcome = outcome_name, Metric = metric_name, File = out_file)
  })
} else { tibble() }

# -- 12.3 交互效应图 -----------------------------------------------------------
cat("开始绘制显著交互效应图...\n")
interaction_plot_log <- if (nrow(significant_interaction_results) > 0) {
  map_dfr(seq_len(nrow(significant_interaction_results)), function(i) {
    outcome_name <- significant_interaction_results$Outcome[i]
    metric_name  <- significant_interaction_results$Metric[i]
    list_key     <- paste0(outcome_name, "___Moderation")
    res_i        <- analysis_results[[list_key]]
    
    data_plot <- getModelData(outcome_name, "Moderation", res_i$kept_metrics, data_wide)
    out_file  <- file.path(plot_dir_int, paste0("Interaction_", cleanFilename(outcome_name), "_", cleanFilename(metric_name), ".png"))
    
    plotInteractionEffect(outcome_name, metric_name, res_i$fit, data_plot, out_file)
    tibble(Plot_Type = "Metric_Group_interaction", Outcome = outcome_name, Metric = metric_name, File = out_file)
  })
} else { tibble() }

# 保存绘图日志
all_plot_log <- bind_rows(pa_plot_log, fs_plot_log, interaction_plot_log)
if (nrow(all_plot_log) > 0) {
  write.csv(all_plot_log, file.path(plot_dir, "Significant_effect_plot_log.csv"), row.names = FALSE, fileEncoding = "UTF-8")
}

# =============================================================================
# 13. 最终控制台摘要
# =============================================================================
cat("\n============================================================\n")
cat("VIF 剪枝多指标回归分析（附总体 FDR 门控、事后检验与鲁棒性验证）全流程完成。\n")
cat("输出目录：", output_dir, "\n\n")
cat("生成文件：\n")
cat("  DOCX 报告：",      length(analysis_results), "份（内附假设检验图及简单斜率结果）\n")
cat("  所有模型鲁棒检验：", check_dir, "\n")
cat("============================================================\n")