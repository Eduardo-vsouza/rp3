library(shiny)
library(DT)
library(ggplot2)
library(plotly)

args <- commandArgs(trailingOnly = TRUE)
paths_to_dfs <- strsplit(args[1], ",")[[1]]
pdf_dir <- args[2]

df_list <- lapply(paths_to_dfs, function(path) {
  read.table(path, sep='\t', header=TRUE, stringsAsFactors = FALSE)
})

df_names <- sapply(paths_to_dfs, basename)

ui <- fluidPage(
  titlePanel("Dynamic Data Frame Viewer with Filters"),
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_df", "Select Data Frame", choices = paths_to_dfs),
      uiOutput("filter_ui"),
      actionButton("apply_filters", "Apply Filters", class = "btn-primary"),
      actionButton("clear_filters", "Clear Filters", class = "btn-secondary"),
      actionButton("plot_volcano", "Plot Volcano Plot", class = "btn-success")
    ),
    mainPanel(
      DTOutput("data_table"),
      plotlyOutput("volcano_plot")
    )
  )
)

server <- function(input, output, session) {
  current_df <- reactive({
    if (is.null(input$selected_df)) {
      return(NULL)
    }
    df_index <- match(input$selected_df, paths_to_dfs)
    if (!is.na(df_index)) {
      df_list[[df_index]]
    } else {
      NULL
    }
  })
  
  output$filter_ui <- renderUI({
    df <- current_df()
    if (is.null(df) || ncol(df) == 0) {
      return(NULL)
    }
    
    filters <- lapply(names(df), function(col) {
      if (is.numeric(df[[col]])) {
        sliderInput(paste0("slider_", col), sprintf("Filter %s:", col),
                    min = min(df[[col]], na.rm = TRUE), max = max(df[[col]], na.rm = TRUE),
                    value = c(min(df[[col]], na.rm = TRUE), max(df[[col]], na.rm = TRUE)))
      } else {
        textInput(paste0("text_", col), sprintf("Filter %s:", col), value = "")
      }
    })
    do.call(tagList, filters)
  })
  
  observeEvent(input$apply_filters, {
    df <- current_df()
    if (is.null(df)) return()
    
    for (col in names(df)) {
      input_id_slider <- paste0("slider_", col)
      input_id_text <- paste0("text_", col)
      
      if (!is.null(input[[input_id_slider]]) && length(input[[input_id_slider]]) == 2) {
        range <- input[[input_id_slider]]
        df <- df[df[[col]] >= range[1] & df[[col]] <= range[2], ]
      }
      if (!is.null(input[[input_id_text]]) && nchar(input[[input_id_text]]) > 0) {
        text <- input[[input_id_text]]
        df <- df[grepl(text, df[[col]], ignore.case = TRUE), ]
      }
    }
    output$data_table <- renderDataTable({
      datatable(df)
    })
  })
  
  observeEvent(input$clear_filters, {
    output$data_table <- renderDataTable({
      datatable(current_df())
    })
  })
  
  observeEvent(input$plot_volcano, {
    df <- current_df()  # Assuming you have a reactive called filtered_data that applies the filters
    req(df)
    
    # Add a new column for color categorization based on your conditions
    df$ColorCategory <- ifelse(abs(df$`Protein.Log2.Fold.Change`) > 1 & df$`False.Discovery.Rate` < 0.05, "red",
                               ifelse(abs(df$`Protein.Log2.Fold.Change`) > 1, "green",
                                      ifelse(df$`False.Discovery.Rate` < 0.05, "blue", "grey")))
    
    output$volcano_plot <- renderPlotly({
      p <- ggplot(df, aes(x = `Protein.Log2.Fold.Change`, y = -log10(`False.Discovery.Rate`), text = Protein.Group)) +
        geom_point(aes(color = ColorCategory), size = 3) +
        geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +  # Log2 Fold-Change cutoff lines
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +  # FDR cutoff line
        scale_color_manual(values = c("red" = "red", "green" = "green", "blue" = "blue", "grey" = "grey")) +
        labs(x = "Log2 Fold Change", y = "-log10 FDR", title = "Volcano Plot") +
        theme_minimal()
      ggplotly(p, tooltip = "text")
    })
  })
  
}

shinyApp(ui = ui, server = server)
