# ==============================
# TRIX Shiny App (Bodrum & Milas, 2013-2022)
# ==============================
options(shiny.maxRequestSize = 30*1024^2)  # 30 MB

# ----- Libraries -----
suppressPackageStartupMessages({
  library(shiny)
  library(shinythemes)
  library(shinycssloaders)
  library(DT)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(ggplot2)
  library(scales)
  library(broom)
  library(janitor)
  library(zoo)
  library(zyp)
  library(Kendall)
  library(trend)
  library(fixest)
  library(zip)
})

# ----- Helpers (pure R; English only) -----
clean_numeric <- function(x){
  x <- as.character(x)
  x <- str_replace_all(x, "[\\u00A0\\s]", "")
  x <- str_replace(x, ",", ".")
  suppressWarnings(as.numeric(x))
}

to_long_trix <- function(df, dataset_label){
  df <- janitor::clean_names(df)
  # best-effort: first column = station id if "no" does not exist
  if (!"no" %in% names(df)) names(df)[1] <- "no"
  df <- df %>% mutate(across(everything(), as.character))
  # year columns like x2013..x2022 or 2013..2022
  year_cols <- grep("^(x)?20[0-2][0-9]$", names(df), value = TRUE)
  if (length(year_cols) == 0) stop("No year columns detected in the Excel sheet.")
  out <- df %>%
    rename(station = no) %>%
    pivot_longer(cols = all_of(year_cols),
                 names_to = "year_chr", values_to = "TRIX_raw") %>%
    mutate(
      year_chr = gsub("^x", "", year_chr),
      year = suppressWarnings(as.integer(year_chr)),
      TRIX = clean_numeric(TRIX_raw),
      dataset = dataset_label
    ) %>% filter(!is.na(year))
  out
}

mk_sen_yuepilon <- function(y, x){
  ok <- which(!is.na(x) & !is.na(y))
  x <- x[ok]; y <- y[ok]
  if(length(unique(x)) < 3 || length(y) < 3){
    return(tibble(tau = NA_real_, p_value = NA_real_, sen_slope = NA_real_, intercept = NA_real_))
  }
  fit <- tryCatch(zyp::zyp.trend.vector(y = y, x = x, method = "yuepilon"), error = function(e) NULL)
  if (is.null(fit)){
    return(tibble(tau = NA_real_, p_value = NA_real_, sen_slope = NA_real_, intercept = NA_real_))
  }
  tau_val <- tryCatch(unname(fit$tau), error = function(e) NA_real_)
  p_val   <- tryCatch(unname(fit$sig), error = function(e) NA_real_)
  slope   <- tryCatch(unname(fit$trend), error = function(e) NA_real_)
  intercept <- tryCatch(unname(fit$intercept), error = function(e) NA_real_)
  tibble(tau = tau_val, p_value = p_val, sen_slope = slope, intercept = intercept)
}

mk_sen_simple <- function(y, x){
  ok <- which(!is.na(x) & !is.na(y))
  x <- as.numeric(x[ok]); y <- as.numeric(y[ok])
  if(length(unique(x)) < 3){
    return(tibble(tau = NA_real_, p_value = NA_real_, sen_slope = NA_real_))
  }
  mk  <- Kendall::MannKendall(y)
  sen <- trend::sens.slope(ts(y, start = min(x), frequency = 1))
  tibble(tau = unname(mk$tau[1]),
         p_value = unname(mk$sl[1]),
         sen_slope = unname(sen$estimates[1]))
}

fit_its <- function(df, break_year = 2020){
  min_y <- min(df$year, na.rm = TRUE)
  d <- df %>%
    mutate(
      time = year - min_y,
      post = as.integer(year >= break_year),
      time_after = pmax(0, time - (break_year - min_y))
    )
  m <- fixest::feols(TRIX ~ time + post + time_after | station, data = d, cluster = ~station)
  broom::tidy(m, conf.int = TRUE) %>%
    mutate(model = "ITS", break_year = break_year)
}

fit_did <- function(df, treat_cut = "tertile", pre_end = 2018, post_start = 2020){
  base <- df %>%
    filter(year <= pre_end) %>%
    group_by(station, dataset) %>%
    summarize(pre_mean = mean(TRIX, na.rm = TRUE), .groups = "drop")
  if (treat_cut == "tertile"){
    thr <- quantile(base$pre_mean, probs = c(2/3), na.rm = TRUE)
    base <- base %>% mutate(treat = as.integer(pre_mean >= thr[1]))
  } else {
    thr <- median(base$pre_mean, na.rm = TRUE)
    base <- base %>% mutate(treat = as.integer(pre_mean >= thr))
  }
  dd <- df %>%
    left_join(base %>% select(station, dataset, treat), by = c("station","dataset")) %>%
    mutate(post = as.integer(year >= post_start))
  m <- fixest::feols(TRIX ~ treat:post | station + year, data = dd, cluster = ~station)
  broom::tidy(m, conf.int = TRUE) %>%
    mutate(model = "DiD", post_start = post_start)
}

plot_box_all <- function(dat){
  df_box <- dat %>% filter(!is.na(TRIX)) %>% mutate(year = as.integer(year))
  ggplot(df_box, aes(x = factor(year), y = TRIX)) +
    geom_boxplot(outlier.alpha = 0.5) +
    facet_wrap(~ dataset, ncol = 1, scales = "free_y") +
    labs(x = "Year", y = "TRIX", title = "Annual TRIX Distributions by Station") +
    theme_minimal(base_size = 12)
}

plot_forest <- function(trend_station_simple){
  trend_station_simple %>%
    group_by(dataset) %>%
    arrange(dataset, desc(sen_slope), .by_group = TRUE) %>%
    mutate(station_f = factor(station, levels = unique(station))) %>%
    ggplot(aes(x = sen_slope, y = station_f)) +
    geom_vline(xintercept = 0, linetype = 2) +
    geom_point() +
    facet_wrap(~dataset, ncol = 1, scales = "free_y") +
    labs(x = "Sen's slope (TRIX per year)", y = "Station ID",
         title = "Station-level TRIX Trends (Mann-Kendall + Sen's slope)",
         subtitle = "Point = slope estimate; dashed line = 0 (no trend)") +
    theme_minimal(base_size = 12)
}

plot_trend_dataset <- function(dat){
  df_mean <- dat %>% group_by(dataset, year) %>% summarize(TRIX = mean(TRIX, na.rm = TRUE), .groups = "drop")
  ggplot(df_mean, aes(year, TRIX, color = dataset)) +
    geom_point() +
    geom_smooth(method = "lm", se = TRUE) +
    scale_x_continuous(breaks = 2013:2022) +
    labs(x = "Year", y = "Mean TRIX", title = "Bodrum vs Milas: Annual Mean TRIX Trends") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
}

plot_ga <- function(dat){
  df_mean <- dat %>% group_by(dataset, year) %>% summarize(TRIX = mean(TRIX, na.rm = TRUE), .groups = "drop")
  ggplot(df_mean, aes(x = year, y = TRIX, color = dataset)) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 4, alpha = 0.08) +
    geom_hline(yintercept = 4, linetype = 2) +
    geom_point(size = 2) +
    geom_smooth(method = "lm", se = TRUE) +
    geom_vline(xintercept = 2020, linetype = 3) +
    scale_x_continuous(breaks = 2013:2022) +
    labs(x = "Year", y = "Mean TRIX", color = NULL) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom",
          plot.margin = margin(5.5, 8, 5.5, 8),
          panel.grid.minor = element_blank())
}

# ----- UI -----
ui <- navbarPage(
  title = "TRIX Analysis - Bodrum & Milas (2013-2022)",
  theme = shinytheme("flatly"),
  
  tabPanel("Data & Settings",
           sidebarLayout(
             sidebarPanel(
               h4("Inputs"),
               fileInput("bodrum_file", "Upload Bodrum Excel (.xlsx)", accept = c(".xlsx")),
               fileInput("milas_file",  "Upload Milas Excel (.xlsx)",  accept = c(".xlsx")),
               checkboxInput("use_default", "Use default files in working directory (DATA_BODRUM.xlsx, DATA_MILAS.xlsx)", TRUE),
               hr(),
               h4("Analysis options"),
               radioButtons("trend_method", "Trend method",
                            c("Simple (MK + Sen)" = "simple",
                              "YuePilon (zyp)"    = "yuepilon"),
                            selected = "simple"),
               numericInput("break_year", "ITS break year", value = 2020, min = 2015, max = 2022, step = 1),
               numericInput("pre_end", "DiD: pre period end", value = 2018, min = 2013, max = 2019, step = 1),
               numericInput("post_start", "DiD: post period start", value = 2020, min = 2020, max = 2022, step = 1),
               selectInput("treat_cut", "DiD: treatment cut",
                           choices = c("tertile","median"), selected = "tertile"),
               hr(),
               actionButton("run_btn", "Run analysis", class = "btn-primary btn-block")
             ),
             mainPanel(
               h4("Data preview"),
               verbatimTextOutput("data_info"),
               DTOutput("preview") %>% withSpinner(),
               tags$hr(),
               p("Notes: If no files are uploaded, the app will try to use DATA_BODRUM.xlsx and DATA_MILAS.xlsx in the working directory.")
             )
           )
  ),
  
  tabPanel("Summaries",
           fluidRow(
             column(12, h4("Annual summary by dataset and year")),
             column(12, DTOutput("tbl_yearly") %>% withSpinner())
           ),
           hr(),
           fluidRow(
             column(12, h4("Change 2013 -> 2022")),
             column(12, DTOutput("tbl_delta") %>% withSpinner())
           )
  ),
  
  tabPanel("Trends",
           fluidRow(
             column(6,
                    h4("Station-level trends"),
                    DTOutput("tbl_trend_station") %>% withSpinner()
             ),
             column(6,
                    h4("Dataset-level trend (summary)"),
                    DTOutput("tbl_trend_dataset") %>% withSpinner()
             )
           ),
           hr(),
           fluidRow(
             column(6, plotOutput("plt_forest", height = 450) %>% withSpinner()),
             column(6, plotOutput("plt_trend_ds", height = 450) %>% withSpinner())
           )
  ),
  
  tabPanel("ITS & DiD",
           fluidRow(
             column(6, h4("ITS (break @ year)"), DTOutput("tbl_its") %>% withSpinner()),
             column(6, h4("DiD (treat x post)"), DTOutput("tbl_did") %>% withSpinner())
           )
  ),
  
  tabPanel("Plots",
           fluidRow(
             column(12, plotOutput("plt_box_all", height = 550) %>% withSpinner())
           )
  ),
  
  tabPanel("Graphical Abstract",
           fluidRow(
             column(12, plotOutput("plt_ga", height = 450) %>% withSpinner()),
             column(4, downloadButton("dl_ga_png", "Download GA (PNG)")),
             column(4, downloadButton("dl_ga_svg", "Download GA (SVG)"))
           )
  ),
  
  tabPanel("Export",
           h4("Download outputs"),
           fluidRow(
             column(4, downloadButton("dl_csv_zip", "Download all CSV outputs (ZIP)")),
             column(4, downloadButton("dl_fig_zip", "Download all figures (ZIP)"))
           )
  )
)

# ----- Server -----
server <- function(input, output, session){
  
  # Read inputs (Excel) -------------------------------------------------------
  data_raw <- reactive({
    req(input$run_btn) # run on button
    isolate({
      # decide sources
      bodrum_path <- if (!is.null(input$bodrum_file$datapath)) input$bodrum_file$datapath else "DATA_BODRUM.xlsx"
      milas_path  <- if (!is.null(input$milas_file$datapath))  input$milas_file$datapath  else "DATA_MILAS.xlsx"
      if (!input$use_default && (is.null(input$bodrum_file) || is.null(input$milas_file))){
        validate(need(FALSE, "Please upload both Excel files or enable 'Use default files'."))
      }
      validate(need(file.exists(bodrum_path), "Bodrum Excel file not found."),
               need(file.exists(milas_path),  "Milas Excel file not found."))
      list(
        bodrum = readxl::read_excel(bodrum_path),
        milas  = readxl::read_excel(milas_path)
      )
    })
  })
  
  dat_all <- reactive({
    dfs <- data_raw()
    bodrum_long <- to_long_trix(dfs$bodrum, "Bodrum")
    milas_long  <- to_long_trix(dfs$milas,  "Milas")
    dat <- dplyr::bind_rows(bodrum_long, milas_long) %>%
      mutate(station = as.integer(station)) %>%
      arrange(dataset, station, year)
    validate(need(nrow(dat) > 0, "No data after reshaping. Check Excel formats."))
    dat
  })
  
  output$data_info <- renderText({
    req(dat_all())
    tbl <- dat_all()
    paste0("Rows: ", nrow(tbl), " | Cols: ", ncol(tbl),
           " | Datasets: ", paste(unique(tbl$dataset), collapse = ", "))
  })
  
  output$preview <- DT::renderDT({
    req(dat_all())
    DT::datatable(head(dat_all(), 20), options = list(scrollX = TRUE, pageLength = 20))
  })
  
  
  # Summaries -----------------------------------------------------------------
  yearly <- reactive({
    dat_all() %>%
      group_by(dataset, year) %>%
      summarize(
        n = sum(!is.na(TRIX)),
        mean = mean(TRIX, na.rm = TRUE),
        sd = sd(TRIX, na.rm = TRUE),
        median = median(TRIX, na.rm = TRUE),
        min = min(TRIX, na.rm = TRUE),
        max = max(TRIX, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  delta_13_22 <- reactive({
    out <- dat_all() %>%
      filter(year %in% c(2013, 2022)) %>%
      group_by(dataset, year) %>%
      summarize(mean = mean(TRIX, na.rm = TRUE), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = year, values_from = mean, names_prefix = "mean_")
    validate(need(all(c("mean_2013","mean_2022") %in% names(out)), "Delta table requires years 2013 and 2022."))
    out %>% mutate(delta_2022_2013 = mean_2022 - mean_2013)
  })
  
  
  output$tbl_yearly <- renderDT({
    DT::datatable(yearly() %>% mutate(across(where(is.numeric), ~round(., 2))),
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  output$tbl_delta <- renderDT({
    DT::datatable(delta_13_22() %>% mutate(across(where(is.numeric), ~round(., 2))),
                  options = list(dom = 't', scrollX = TRUE, pageLength = 10))
  })
  
  # Trends --------------------------------------------------------------------
  trend_station_simple <- reactive({
    dat_all() %>%
      group_by(dataset, station) %>%
      summarize(res = list(mk_sen_simple(TRIX, year)), .groups = "drop") %>%
      unnest(res) %>%
      mutate(p_adj_fdr = p.adjust(p_value, method = "fdr"))
  })
  
  trend_station_yue <- reactive({
    dat_all() %>%
      group_by(dataset, station) %>%
      summarize(res = list(mk_sen_yuepilon(TRIX, year)), .groups = "drop") %>%
      unnest(res) %>%
      mutate(p_adj_fdr = p.adjust(p_value, method = "fdr"))
  })
  
  trend_station_out <- reactive({
    if (input$trend_method == "yuepilon") trend_station_yue() else trend_station_simple()
  })
  
  trend_dataset <- reactive({
    agg <- dat_all() %>%
      group_by(dataset, year) %>%
      summarize(TRIX = mean(TRIX, na.rm = TRUE), .groups = "drop")
    if (input$trend_method == "yuepilon") {
      agg %>% group_by(dataset) %>% summarize(mk = list(mk_sen_yuepilon(TRIX, year)), .groups="drop") %>% unnest(mk)
    } else {
      agg %>% group_by(dataset) %>% summarize(mk = list(mk_sen_simple(TRIX, year)), .groups="drop") %>% unnest(mk)
    }
  })
  
  
  output$tbl_trend_station <- renderDT({
    DT::datatable(trend_station_out() %>%
                    mutate(across(where(is.numeric), ~round(., 4))),
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  output$tbl_trend_dataset <- renderDT({
    DT::datatable(trend_dataset() %>%
                    mutate(across(where(is.numeric), ~round(., 4))),
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  
  output$plt_forest <- renderPlot({
    plot_forest(trend_station_out())
  })
  output$plt_trend_ds <- renderPlot({
    plot_trend_dataset(dat_all())
  })
  
  # ITS & DiD -----------------------------------------------------------------
  its <- reactive({
    dat_all() %>% group_by(dataset) %>%
      group_modify(~ fit_its(.x, break_year = input$break_year)) %>%
      ungroup()
  })
  
  did <- reactive({
    dat_all() %>% group_split(dataset) %>% map_dfr(~ {
      fit_did(.x, treat_cut = input$treat_cut, pre_end = input$pre_end, post_start = input$post_start) %>%
        mutate(dataset = unique(.x$dataset))
    })
  })
  
  output$tbl_its <- renderDT({
    DT::datatable(its() %>% mutate(across(where(is.numeric), ~round(., 4))),
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  output$tbl_did <- renderDT({
    DT::datatable(did() %>% mutate(across(where(is.numeric), ~round(., 4))),
                  options = list(scrollX = TRUE, pageLength = 10))
  })
  
  # Plots ---------------------------------------------------------------------
  output$plt_box_all <- renderPlot({
    plot_box_all(dat_all())
  })
  
  # Graphical Abstract --------------------------------------------------------
  output$plt_ga <- renderPlot({
    plot_ga(dat_all())
  })
  
  output$dl_ga_png <- downloadHandler(
    filename = function(){ "graphical_abstract_trix.png" },
    content = function(file){
      g <- plot_ga(dat_all())
      ggsave(file, g, width = 8, height = 5, dpi = 300)
    }
  )
  output$dl_ga_svg <- downloadHandler(
    filename = function(){ "graphical_abstract_trix.svg" },
    content = function(file){
      g <- plot_ga(dat_all())
      ggsave(file, g, width = 8, height = 5, dpi = 300, device = "svg", bg = "white")
    }
  )
  
  # Export ZIPs ---------------------------------------------------------------
  output$dl_csv_zip <- downloadHandler(
    filename = function(){ "trix_outputs_csv.zip" },
    content = function(file){
      tmpdir <- tempdir()
      f1 <- file.path(tmpdir,"trix_yearly_summary.csv")
      f2 <- file.path(tmpdir,"trix_change_2013_2022.csv")
      f3 <- file.path(tmpdir,"trix_trend_station.csv")
      f4 <- file.path(tmpdir,"trix_trend_dataset.csv")
      f5 <- file.path(tmpdir,"trix_ITS.csv")
      f6 <- file.path(tmpdir,"trix_DiD.csv")
      readr::write_csv(yearly(), f1)
      readr::write_csv(delta_13_22(), f2)
      readr::write_csv(trend_station_out(), f3)
      readr::write_csv(trend_dataset(), f4)
      readr::write_csv(its(), f5)
      readr::write_csv(did(), f6)
      zip::zipr(zipfile = file, files = c(f1,f2,f3,f4,f5,f6), recurse = FALSE)
    }
  )
  
  output$dl_fig_zip <- downloadHandler(
    filename = function(){ "trix_figures.zip" },
    content = function(file){
      tmpdir <- tempdir()
      p1 <- file.path(tmpdir,"figure_box_yearly_all.png")
      p2 <- file.path(tmpdir,"figure_forest_trend_all.png")
      p3 <- file.path(tmpdir,"figure_trend_dataset_fixed.png")
      p4 <- file.path(tmpdir,"graphical_abstract_trix.png")
      ggsave(p1, plot_box_all(dat_all()), width = 9, height = 7, dpi = 300)
      ggsave(p2, plot_forest(trend_station_out()), width = 8, height = 10, dpi = 300)
      ggsave(p3, plot_trend_dataset(dat_all()), width = 8, height = 5, dpi = 300)
      ggsave(p4, plot_ga(dat_all()), width = 8, height = 5, dpi = 300)
      zip::zipr(zipfile = file, files = c(p1,p2,p3,p4), recurse = FALSE)
    }
  )
  
}

shinyApp(ui, server)
