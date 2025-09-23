library(shiny)
library(songbirdR)
library(shinyjs)
library(DT)
library(plotly)
library(tensorflow) 


ui <- fluidPage(
  useShinyjs(),
  div(style = "padding-top: 15px;",
      sidebarLayout(
        sidebarPanel(width = 3,  #Model setting and plot setting panel
                     tabsetPanel(
                       tabPanel("Model Settings",
                                div(style = "max-height: 875px; overflow-y: auto; overflow-x: hidden; box-sizing: border-box; width: 100%;",
                                    textInput("conda_env", "Conda environment name"),
                                    actionButton("load_env", "Load Conda Environment"),
                                    htmlOutput("env_message_ui"),
                                    tags$hr(),
                                    fileInput("otu_table", "Upload OTU Table (CSV)", accept = ".csv"),
                                    fileInput("metadata", "Upload Metadata (CSV)", accept = ".csv"),
                                    uiOutput("formula_ui"),
                                    checkboxInput("use_cluster", "Clustered resampling", value = FALSE),
                                    conditionalPanel(
                                      condition = "input.use_cluster == true",
                                      uiOutput("cluster_var_ui")
                                    ),
                                    uiOutput("random_effects_ui"),
                                    tags$small("Note: Categorical variables should be encoded as character or factor in your metadata. Numeric values like 1/2/3 will be treated as continuous unless converted. Delete variables with backspace"),
                                    uiOutput("reference_ui"),
                                    fluidRow(
                                      column(6,
                                             numericInput("differential_prior", "Differential Prior", value = 1, min = 1e-3, step = 0.1),
                                             numericInput("clipnorm", "Clipnorm", value = 10, min = 1, step = 1),
                                             numericInput("re_lambda", "Random-effect λ (L2)", value = 1e-3, min = 1e-6, step = 1e-4),
                                             numericInput("learning_rate", "Learning Rate", value = 0.001, min = 1e-8, step = 1e-4),
                                             numericInput("batch_size", "Batch Size", value = 5, min = 1, step = 1),
                                             numericInput("num_test_samples", "Number of Test Samples", value = 5, min = 1, step = 1),
                                             numericInput("seed", "Random Seed", value = 0),
                                             selectInput("parallel_mode", "Parallel mode",
                                                         choices = c("no","multicore","snow"),
                                                         selected = "no"),
                                             numericInput("n_cores", "Number of Cores", value = 1, min = 1, step = 1)
                                      ),
                                      column(6,
                                             numericInput("epochs", "Max Epochs per Resample", value = 1000, min = 2, step = 1),
                                             numericInput("n_boot", "Number of Bootstraps", value = 10, min = 0, step = 1),
                                             numericInput("n_perms", "Number of Permutations", value = 10, min = 0, step = 1),
                                             checkboxInput("early_stopping", "Early stopping", value = TRUE),
                                             numericInput("patience", "Patience (evals)", value = 20, min = 1, step = 1),
                                             numericInput("min_epochs", "Min epochs before stop", value = 200, min = 1, step = 1),
                                             numericInput("eval_every", "Eval every N steps", value = 20, min = 1, step = 1),
                                             numericInput("tol", "Improvement tolerance", value = 1e-4, step = 1e-5),
                                             checkboxInput("lr_plateau", "LR decay on plateau", value = TRUE),
                                             numericInput("lr_decay_factor", "LR decay factor", value = 0.5, min = 0.01, max = 1, step = 0.05),
                                             numericInput("lr_min", "LR minimum", value = 1e-6, min = 1e-8, step = 1e-6)
                                      )
                                    ),
                                    actionButton("choose_tb_dir", "Choose TensorBoard log directory"),
                                    div(style="margin-top:6px"),
                                    textOutput("tb_dir_text"),
                                    
                                    
                                    actionButton("run_model", "Run Regression"),
                                    htmlOutput("regression_message_ui"),
                                    
                                    div(style="margin-top:6px"),
                                    actionButton("open_tb", "Open TensorBoard"),
                                    tags$hr(),
                                )
                       ),#rank plot settings
                       tabPanel("Plot Settings",
                                uiOutput("coef_selector"),
                                textInput("highlight_otus", "Highlight OTUs (comma-separated)", value = ""),
                                numericInput("ylim_min", "Y-axis Min", value = -3),
                                numericInput("ylim_max", "Y-axis Max", value = 3),
                                numericInput("pval_threshold", "P-value Threshold", value = 0.05),
                                actionButton("run_plot", "Generate Plot")
                       )
                     )
        ),
        mainPanel(width = 9, # rank plot
                  fluidRow(
                    column(8,
                           uiOutput("plot_or_message")
                    ),
                    column(4, # csv browse and download, add to numerator/denominator
                           div(style = "height:280px; overflow-y:auto; border:1px solid #ccc; padding:5px; box-sizing:border-box;",
                               uiOutput("csv_browse_box")
                           ),
                           selectInput("csv_type", "Select CSV Type to Browse", choices = c("beta_mean", "beta_sd", "beta_pval", "beta_fdr")),
                           downloadButton("download_csv_dynamic", "Download CSV"),
                           downloadButton("download_rds", "Download whole RDS"),
                           div(style = "margin-top: 10px;",
                               actionButton("add_to_num", "Add to Numerator"),
                               actionButton("add_to_den", "Add to Denominator", style = "margin-left: 8px;")
                           )
                           
                    )
                  ),
                  div(style="height:6px"), # Log-ratio calc 
                  div(style = "height:540px; border:1px solid #ccc; padding:10px; box-sizing:border-box;",
                      fluidRow(
                        column(3,
                               textInput("numerator_input", "Numerator OTUs (comma-separated)", ""),
                               textInput("denominator_input", "Denominator OTUs (comma-separated)", "")
                        ),
                        column(3,
                               uiOutput("logratio_variable_selector"),
                               selectInput("logratio_mode", "Analysis Type", choices = c("categorical", "continuous")),
                               actionButton("run_logratio_plot", "Visualize Log-Ratio"),
                               downloadButton("download_logratios", "Download Log-Ratios"),
                               htmlOutput("logratio_message_ui")
                        ),
                        column(5,
                               uiOutput("logratio_plot_or_message")
                        )
                      ),
                      tags$hr(), # selbal calc
                      fluidRow(
                        column(4, uiOutput("selbal_variable_selector")),
                        column(2, checkboxInput("selbal_categorical", "Categorical Variable", value = TRUE)),
                        column(3, numericInput("selbal_nfold", "Number of folds", value = 5, min = 2)),
                        column(3, numericInput("selbal_niter", "Number of iterations", value = 10, min = 1))
                      ),
                      fluidRow(
                        column(9, actionButton("run_selbal", "Get Balances")),
                        column(3, htmlOutput("selbal_message_ui"))
                      ),
                      uiOutput("selbal_result_or_message")
                  )
        )
      )
  )
)



server <- function(input, output, session) {
  
  `%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

  #status and result containers
  betas_result        <- reactiveVal(NULL)
  metadata_data       <- reactiveVal(NULL)
  regression_status   <- reactiveVal("")
  env_status          <- reactiveVal("")
  logratio_status     <- reactiveVal("")
  selbal_status       <- reactiveVal("")
  tb_dir              <- reactiveVal("")
  logratios_result    <- reactiveVal(NULL)    
  selbal_result_text  <- reactiveVal(NULL)     
  
  #metadata upload
  observeEvent(input$metadata, {
    metadata <- read.csv(input$metadata$datapath, row.names = 1)
    metadata_data(metadata)
    updateSelectInput(session, "formula",
                      choices = names(metadata),
                      selected = names(metadata)[1])
  })
  
  #dynamic UI (formula and reference levels)
  output$formula_ui <- renderUI({
    req(metadata_data())
    selectizeInput("formula", "Model Formula Columns",
                   choices  = names(metadata_data()),
                   selected = names(metadata_data())[1],
                   multiple = TRUE)
  })
  output$cluster_var_ui <- renderUI({
    req(metadata_data())
    selectInput("cluster_var", "Cluster (bootstrap block) variable",
                choices = names(metadata_data()),
                selected = names(metadata_data())[1],
                multiple = FALSE)
  })
  # Random effects (batch/subject) variable picker
  output$random_effects_ui <- renderUI({
    req(metadata_data())
    selectizeInput(
      "random_effects_vars",
      "Random effects (batch/subject) variables",
      choices  = names(metadata_data()),
      selected = NULL,            # start empty (checkbox can auto-fill)
      multiple = TRUE,
      options = list(placeholder = 'Select one or more (e.g., Batch)')
    )
  })
  
  output$reference_ui <- renderUI({
    req(input$formula)
    lapply(input$formula, function(var) {
      var_data <- metadata_data()[[var]]
      if (is.factor(var_data) || is.character(var_data)) {
        selectInput(paste0("ref_", var),
                    paste("Reference level for", var),
                    choices = unique(var_data))
      }
    })
  })
  
  #environment message
  output$env_message_ui <- renderUI({ HTML(env_status()) })
  
  observeEvent(input$load_env, {
    env_status("<span style='color: orange; font-style: italic;'>Loading Conda environment...</span>")
    session$onFlushed(function() {
      later::later(function() {
        result <- tryCatch({
          library(reticulate)
          use_condaenv(isolate(input$conda_env), required = TRUE)
          Sys.setenv(RETICULATE_PYTHON = py_config()$python)
          paste0("<span style='color: #1E90FF; font-style: italic;'>Conda environment ",
                 isolate(input$conda_env), " loaded successfully.</span>")
        }, error = function(e) {
          paste0("<span style='color: red; font-style: italic;'>Failed to load environment: ",
                 isolate(input$conda_env), "<br>Error: ", e$message, "</span>")
        })
        env_status(result)  
      }, delay = 0.1)
    }, once = TRUE)
  })
  
  #regression message
  output$regression_message_ui <- renderUI({ HTML(regression_status()) })
  
  observeEvent(input$run_model, {
    # show status
    regression_status("<span style='color: orange; font-style: italic;'>Running regression...</span>")
    
    req(input$otu_table, input$metadata)
    otu_path     <- input$otu_table$datapath
    meta_path    <- input$metadata$datapath
    form_terms   <- input$formula %||% character(0)
    diff_prior   <- input$differential_prior
    lr           <- input$learning_rate
    clipnorm     <- input$clipnorm
    bsz          <- input$batch_size
    seed         <- input$seed
    ntest        <- input$num_test_samples
    epochs       <- input$epochs
    nboot        <- input$n_boot
    ncores       <- input$n_cores
    cluster_var  <- input$cluster_var
    parallel_sel <- input$parallel_mode %||% "auto"
    nperms       <- input$n_perms
    re_lambda    <- input$re_lambda
    early_stop   <- isTRUE(input$early_stopping)
    pat          <- input$patience
    min_ep       <- input$min_epochs
    eval_every   <- input$eval_every
    tol          <- input$tol
    lr_plateau   <- isTRUE(input$lr_plateau)
    lr_factor    <- input$lr_decay_factor
    lr_minimum   <- input$lr_min
    
    # resampling mode
    resample_mode <- if (isTRUE(input$use_cluster)) "cluster" else "row"
    if (resample_mode == "cluster") req(cluster_var)
    
    # capture reference levels if present 
    ref_levels <- NULL
    if (length(form_terms)) {
      ref_levels <- lapply(form_terms, function(var) input[[paste0("ref_", var)]])
      names(ref_levels) <- form_terms
    }
    
    # RE input
    re_vars_raw <- input$random_effects_vars %||% character(0)
    re_vars_raw <- re_vars_raw[nzchar(re_vars_raw)]
    random_effects_param <- if (length(re_vars_raw)) re_vars_raw else NULL
    
    # capture TensorBoard dir
    tb_path <- isolate(tb_dir())  
    
    later::later(function() {
      tryCatch({
        otu_table <- read.csv(otu_path,  row.names = 1, check.names = FALSE)
        metadata  <- read.csv(meta_path, row.names = 1)
        
        # align rows by sample ids (robust)
        if (!all(rownames(otu_table) %in% rownames(metadata)) ||
            !all(rownames(metadata) %in% rownames(otu_table))) {
          common <- intersect(rownames(otu_table), rownames(metadata))
          otu_table <- otu_table[common, , drop = FALSE]
          metadata  <- metadata [common, , drop = FALSE]
        }
        stopifnot(all(rownames(otu_table) == rownames(metadata)))
        
        # apply reference levels to fixed-effect variables
        if (length(form_terms)) {
          for (var in form_terms) {
            lvl <- ref_levels[[var]]
            if (!is.null(lvl)) {
              metadata[[var]] <- factor(
                metadata[[var]],
                levels = c(lvl, setdiff(unique(metadata[[var]]), lvl))
              )
            }
          }
        }
        
        # build fixed effect formula object
        formula_obj <- if (length(form_terms)) {
          as.formula(paste("~", paste(form_terms, collapse = "+")))
        } else {
          as.formula("~ 1")
        }
        
        betas <- run_regression(
          otu_table        = otu_table,
          metadata         = metadata,
          formula          = formula_obj,
          cluster_var      = if (resample_mode == "cluster") cluster_var else NULL,
          random_effects   = random_effects_param,
          reference_levels = ref_levels,
          differential_prior = diff_prior,
          learning_rate      = lr,
          clipnorm           = clipnorm,
          batch_size         = bsz,
          seed               = seed,
          num_test_samples   = ntest,
          epochs             = epochs,
          n_boot             = nboot,
          n_perms            = nperms,        
          n_cores            = ncores,
          parallel           = parallel_sel,
          re_lambda          = re_lambda,     
          resample           = resample_mode,
          early_stopping     = early_stop,  
          patience           = pat,
          min_epochs         = min_ep,
          eval_every         = eval_every,
          tol                = tol,
          lr_decay_on_plateau= lr_plateau,
          lr_decay_factor    = lr_factor,
          lr_min             = lr_minimum,
          tb_logdir          = if (nzchar(tb_path)) tb_path else NULL
        )
        
        
        betas_result(betas)
        regression_status("<span style='color: #1E90FF; font-style: italic;'>Regression done.</span>")
      }, error = function(e) {
        regression_status(paste0(
          "<span style='color: red; font-style: italic;'>Regression failed: ",
          e$message, "</span>"
        ))
      })
    }, delay = 0.05)
  })
  
  
  
  observeEvent(input$choose_tb_dir, {
    # Native folder dialog 
    path <- tryCatch({
      if (.Platform$OS.type == "windows") {
        utils::choose.dir(caption = "Select TensorBoard log directory")
      } else {
        if (!requireNamespace("tcltk", quietly = TRUE)) stop("Package 'tcltk' is required on non-Windows for folder selection.")
        tcltk::tk_choose.dir(caption = "Select TensorBoard log directory")
      }
    }, error = function(e) NA_character_)
    
    if (!is.na(path) && nzchar(path)) {
      tb_dir(normalizePath(path, winslash = "/"))
    }
  })
  
  output$tb_dir_text <- renderText({
    p <- tb_dir()
    if (nzchar(p)) paste("Logdir:", p) else "Logdir: (none – TensorBoard logging disabled)"
  })
  
  observeEvent(input$open_tb, {
    req(nzchar(tb_dir()))
    tryCatch({
      tensorflow::tensorboard(tb_dir())  
    }, error = function(e) {
      showNotification(paste("Failed to open TensorBoard:", e$message), type = "error")
    })
  })
  
  
  
  #rank plot
  output$coef_selector <- renderUI({
    req(betas_result())
    selectInput("coef_name", "Coefficient to Plot",
                choices = colnames(betas_result()$beta_mean))
  })
  
  output$songbird_plotly <- renderPlotly({
    req(betas_result(), input$coef_name, input$ylim_min, input$ylim_max, input$pval_threshold)
    
    res <- betas_result()
    df  <- coef_df() 
    
    arr <- res$beta_array
    j   <- match(input$coef_name, dimnames(arr)[[2]])
    stopifnot(!is.na(j))
    
    A <- arr[, j, , drop = FALSE]
    K <- dim(A)[1]; B <- dim(A)[3]
    A <- array(A, dim = c(K, B))  
    
    qlo <- 0.025; qhi <- 0.975
    lo  <- apply(A, 1, function(v) { v <- v[is.finite(v)]; if (length(v) >= 2) quantile(v, qlo, names = FALSE) else NA_real_ })
    hi  <- apply(A, 1, function(v) { v <- v[is.finite(v)]; if (length(v) >= 2) quantile(v, qhi, names = FALSE) else NA_real_ })
    
    # align to df order via Feature
    feats <- rownames(res$beta_mean)
    lo_map <- setNames(lo, feats)
    hi_map <- setNames(hi, feats)
    df$ci_lo <- unname(lo_map[df$Feature])
    df$ci_hi <- unname(hi_map[df$Feature])
    
    # error bar lengths
    df$err_up <- pmax(0, df$ci_hi - df$beta_mean)
    df$err_dn <- pmax(0, df$beta_mean - df$ci_lo)
    
    #split groups
    sig_flag <- is.finite(df$p_value) & (df$p_value <= input$pval_threshold)
    hl_raw <- trimws(unlist(strsplit(input$highlight_otus, ",")))
    hl <- unique(hl_raw[hl_raw != ""])
    is_hl <- if (length(hl)) df$Feature %in% hl else rep(FALSE, nrow(df))
    
    df_hl   <- df[ is_hl, , drop = FALSE]
    df_sig  <- df[!is_hl &  sig_flag, , drop = FALSE]
    df_nsig <- df[!is_hl & !sig_flag, , drop = FALSE]
    
    hover <- function(d) {
      if (!nrow(d)) return(character(0))
      paste0(
        "OTU: ", d$Feature,
        "<br>Rank: ", d$rank,
        "<br>Beta mean: ", sprintf("%.4f", d$beta_mean),
        "<br>CI: [", ifelse(is.finite(d$ci_lo), sprintf("%.4f", d$ci_lo), "NA"),
        ", ", ifelse(is.finite(d$ci_hi), sprintf("%.4f", d$ci_hi), "NA"), "]",
        "<br>P-value: ", ifelse(is.finite(d$p_value), sprintf("%.4g", d$p_value), "NA")
      )
    }
    
    p <- plot_ly(source = "songbird")
    
    if (nrow(df_nsig)) {
      p <- add_trace(
        p, data = df_nsig,
        x = ~rank, y = ~beta_mean,
        type = "scatter", mode = "markers",
        text = hover(df_nsig), hoverinfo = "text",
        name = paste0("p > ", format(input$pval_threshold)),
        marker = list(size = 6, opacity = 0.35),
        error_y = list(type = "data",
                       array = df_nsig$err_up,
                       arrayminus = df_nsig$err_dn,
                       thickness = 1.2,
                       width = 0)
      )
    }
    
    if (nrow(df_sig)) {
      p <- add_trace(
        p, data = df_sig,
        x = ~rank, y = ~beta_mean,
        type = "scatter", mode = "markers",
        text = hover(df_sig), hoverinfo = "text",
        name = paste0("p \u2264 ", format(input$pval_threshold)),
        marker = list(size = 7),
        error_y = list(type = "data",
                       array = df_sig$err_up,
                       arrayminus = df_sig$err_dn,
                       thickness = 1.2,
                       width = 0)
      )
    }
    
    if (nrow(df_hl)) {
      p <- add_trace(
        p, data = df_hl,
        x = ~rank, y = ~beta_mean,
        type = "scatter", mode = "markers",
        text = hover(df_hl), hoverinfo = "text",
        name = "Highlighted OTUs",
        marker = list(size = 11, symbol = "diamond-open", line = list(width = 1.5)),
        error_y = list(type = "data",
                       array = df_hl$err_up,
                       arrayminus = df_hl$err_dn,
                       thickness = 1.5,
                       width = 0)
      )
    }
    
    p %>%
      layout(
        xaxis = list(title = "Rank"),
        yaxis = list(title = "CLR Beta Coefficient",
                     range = c(input$ylim_min, input$ylim_max)),
        legend = list(orientation = "h", y = -0.2)
      )
  })
  
  
  
  
  # csv browse interaction with rank plot
  csv_proxy <- DT::dataTableProxy("csv_preview")
  observeEvent(plotly::event_data("plotly_click", source = "songbird"), {
    click <- plotly::event_data("plotly_click", source = "songbird")
    req(click)
    
    clicked_rank <- as.integer(round(click$x))
    df_coef <- coef_df()
    req(clicked_rank >= 1, clicked_rank <= nrow(df_coef))
    
    clicked_feature <- df_coef$Feature[clicked_rank]
    
    df_csv <- csv_current_df()  
    idx <- which(df_csv$Feature == clicked_feature)
    if (length(idx)) {
      DT::selectRows(csv_proxy, idx)
      DT::dataTableAjax(session, "csv_preview")
    }
  })
  

  
  coef_df <- reactive({
    req(betas_result(), input$coef_name)
    res <- betas_result()
    bm  <- res$beta_mean[, input$coef_name, drop = TRUE]
    pv  <- res$beta_pval[, input$coef_name, drop = TRUE]
    
    ord <- order(bm, decreasing = FALSE)
    data.frame(
      Feature   = rownames(res$beta_mean)[ord],
      beta_mean = as.numeric(bm[ord]),
      p_value   = as.numeric(pv[ord]),
      rank      = seq_along(ord),
      stringsAsFactors = FALSE
    )
  })
  
  
  output$plot_or_message <- renderUI({
    if (is.null(betas_result()) || input$run_plot == 0) {
      div("Rank plot appears here",
          style = "text-align:center; font-style:italic; font-size: 24px; padding-top: 100px;")
    } else {
      plotlyOutput("songbird_plotly")
    }
  })
  
  #CSV browse box
  output$csv_browse_box <- renderUI({
    if (is.null(betas_result())) {
      div("Run regression to browse",
          style = "text-align:center; font-style:italic; font-size: 20px; padding-top: 100px;")
    } else {
      DT::DTOutput("csv_preview")
    }
  })
  #selection input for logratios
  csv_current_df <- reactive({
    req(betas_result())
    result <- betas_result()
    data_to_show <- switch(input$csv_type,
                           beta_mean = result$beta_mean,
                           beta_sd   = result$beta_sd,
                           beta_pval = result$beta_pval,
                           beta_fdr  = result$beta_fdr)
    df <- data.frame(Feature = rownames(data_to_show),
                     data_to_show, check.names = FALSE)
    rownames(df) <- NULL
    df
  })
  
  
  output$csv_preview <- DT::renderDT({
    df <- csv_current_df()
    DT::datatable(
      df,
      options = list(
        pageLength = nrow(df),
        paging     = FALSE,
        dom        = 'ft',
        autoWidth  = TRUE
      ),
      rownames  = FALSE,
      selection = "multiple"
    )
  })
  
  #appenders to numerator and denominator inputs from CSV selection
  get_selected_features <- reactive({
    idx <- input$csv_preview_rows_selected
    if (length(idx) == 0) return(character(0))
    csv_current_df()$Feature[idx]
  })
  
  append_to_input <- function(input_id, new_items) {
    if (length(new_items) == 0) return(invisible(NULL))
    cur <- isolate(input[[input_id]])
    cur_items <- if (is.null(cur) || cur == "") character(0) else trimws(unlist(strsplit(cur, ",")))
    cur_items <- cur_items[cur_items != ""]
    combined <- unique(c(cur_items, new_items))
    updateTextInput(session, input_id, value = paste(combined, collapse = ", "))
  }
  
  observeEvent(input$add_to_num, {
    append_to_input("numerator_input", get_selected_features())
  })
  
  observeEvent(input$add_to_den, {
    append_to_input("denominator_input", get_selected_features())
  })
  
  
  #log-ratio UI and plotting
  output$logratio_variable_selector <- renderUI({
    req(metadata_data())
    selectInput("logratio_variable", "Metadata Variable",
                choices = colnames(metadata_data()))
  })
  
  output$logratio_message_ui <- renderUI({ HTML(logratio_status()) })
  
  output$logratio_plot_or_message <- renderUI({
    if (is.null(logratios_result()) || input$run_logratio_plot == 0) {
      div("Log-ratio plot appears here",
          style = "text-align:center; font-style:italic; font-size: 18px; padding-top: 100px;")
    } else {
      plotOutput("logratio_plot", height = "260px")
    }
  })
  
  observeEvent(input$run_logratio_plot, {
    logratio_status("<span style='color: orange; font-style: italic;'>Calculating log-ratios...</span>")
    session$onFlushed(function() {
      later::later(function() {
        tryCatch({
          otu_table  <- read.csv(isolate(input$otu_table)$datapath,  row.names = 1, check.names = FALSE)
          metadata   <- read.csv(isolate(input$metadata)$datapath,  row.names = 1)
          numerator  <- trimws(unlist(strsplit(isolate(input$numerator_input), ",")))
          denominator<- trimws(unlist(strsplit(isolate(input$denominator_input), ",")))
          mode       <- isolate(input$logratio_mode)
          variable   <- isolate(input$logratio_variable)
          
          logratios <- compute_logratios(otu_table, numerator, denominator)
          logratios_result(logratios)
          
          output$logratio_plot <- renderPlot({
            if (mode == "categorical") {
              plot_logratios_by_group(logratios, metadata, variable)
            } else {
              plot_logratios_by_continuous(logratios, metadata, variable)
            }
          })
          
          logratio_status("<span style='color: #1E90FF; font-style: italic;'>Log-ratio plot generated.</span>")
        }, error = function(e) {
          logratio_status(paste0("<span style='color: red; font-style: italic;'>Failed to generate log-ratios: ",
                                 e$message, "</span>"))
        })
      }, delay = 0.1)
    }, once = TRUE)
  })
  
  output$download_logratios <- downloadHandler(
    filename = function() { "logratios.csv" },
    content  = function(file) {
      write.csv(logratios_result(), file)
    }
  )
  
  #Selbal
  output$selbal_variable_selector <- renderUI({
    req(metadata_data())
    selectInput("selbal_variable", "Metadata Variable",
                choices = names(metadata_data()))
  })
  
  output$selbal_message_ui <- renderUI({ HTML(selbal_status()) })
  
  observeEvent(input$run_selbal, {
    selbal_status("<span style='color: orange; font-style: italic;'>Running selbal...</span>")
    session$onFlushed(function() {
      later::later(function() {
        tryCatch({
          otu_table <- read.csv(isolate(input$otu_table)$datapath, row.names = 1, check.names = FALSE)
          metadata  <- read.csv(isolate(input$metadata)$datapath, row.names = 1)
          
          selrez <- run_selbal_balance(
            otu_table    = otu_table,
            metadata     = metadata,
            variable_name= isolate(input$selbal_variable),
            n_fold       = isolate(input$selbal_nfold),
            n_iter       = isolate(input$selbal_niter),
            categorical  = isolate(input$selbal_categorical)
          )
          
          balance <- selrez$global.balance
          num <- paste(balance$Taxa[balance$Group == "NUM"], collapse = ", ")
          den <- paste(balance$Taxa[balance$Group == "DEN"], collapse = ", ")
          
          selbal_result_text(paste0("Numerator OTUs: ", num, "\nDenominator OTUs: ", den))
          selbal_status("<span style='color: #1E90FF; font-style: italic;'>Balances calculated.</span>")
        }, error = function(e) {
          selbal_result_text(paste("Error:", e$message))
          selbal_status("<span style='color: red; font-style: italic;'>Selbal failed.</span>")
        })
      }, delay = 0.1)
    }, once = TRUE)
  })
  
  output$selbal_result_or_message <- renderUI({
    if (is.null(selbal_result_text())) {
      div("Balances appear here",
          style = "text-align:center; font-style:italic; font-size: 16px; padding-top: 20px;")
    } else {
      verbatimTextOutput("selbal_result")
    }
  })
  
  output$selbal_result <- renderText({ selbal_result_text() })
  
  #downloads
  output$download_rds <- downloadHandler(
    filename = function() { "regression_result.rds" },
    content  = function(file) { saveRDS(betas_result(), file) }
  )
  
  output$download_csv_dynamic <- downloadHandler(
    filename = function() { paste0(input$csv_type, ".csv") },
    content  = function(file) {
      result <- betas_result()
      data_to_write <- switch(input$csv_type,
                              beta_mean = result$beta_mean,
                              beta_sd   = result$beta_sd,
                              beta_pval = result$beta_pval,
                              beta_fdr  = result$beta_fdr
      )
      write.csv(data_to_write, file)
    }
  )
}


shinyApp(ui = ui, server = server)
