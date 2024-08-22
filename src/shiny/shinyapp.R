library(shiny)
library(DT)
library(ggplot2)
library(plotly)

# Arguments from command line
args <- commandArgs(trailingOnly = TRUE)
paths_to_dfs <- strsplit(args[1], ",")[[1]]  # Comma-separated paths to data frames
image_columns <- strsplit(args[2], ",")[[1]]  # Comma-separated image column names
image_dirs <- strsplit(args[3], ",")[[1]]  # Comma-separated directories for each image column
df_aliases <- strsplit(args[4], ",")[[1]]  # Comma-separated aliases for the data frames

# Load data frames
df_list <- lapply(paths_to_dfs, function(path) {
  read.table(path, sep='\t', header=TRUE, stringsAsFactors = FALSE)
})

# Check that df_aliases is the same length as paths_to_dfs
if (length(df_aliases) != length(paths_to_dfs)) {
  stop("The number of aliases provided must match the number of data frames.")
}

ui <- fluidPage(
  titlePanel("Data Frame Viewer with Filtering, Volcano Plot, and Image Comparison"),
  sidebarLayout(
    sidebarPanel(
      selectInput("selected_df", "Select Data Frame", choices = df_aliases),
      uiOutput("column_selector"),
      actionButton("apply_filters", "Apply Filters", class = "btn-primary"),
      actionButton("clear_filters", "Clear Filters", class = "btn-secondary"),
      actionButton("plot_volcano", "Plot Volcano Plot", class = "btn-success")
    ),
    mainPanel(
      DTOutput("data_table"),
      plotlyOutput("volcano_plot"),
      fluidRow(
        uiOutput("image_output_panels")  # This will render multiple image panels side by side
      )
    )
  )
)

server <- function(input, output, session) {
  # Register each image directory as a resource path with a unique alias
  for (i in seq_along(image_dirs)) {
    addResourcePath(paste0("image_dir_", i), image_dirs[i])
  }
  
  # Reactive expression to get the currently selected data frame
  current_df <- reactive({
    req(input$selected_df)
    # Find the index of the selected alias in df_aliases
    df_index <- match(input$selected_df, df_aliases)
    if (!is.na(df_index) && df_index > 0 && df_index <= length(df_list)) {
      df <- df_list[[df_index]]
      return(df)
    } else {
      return(NULL)
    }
  })
  
  # UI for column selector
  output$column_selector <- renderUI({
    df <- current_df()
    if (is.null(df)) return(NULL)
    checkboxGroupInput("show_columns", "Show/Hide Columns", choices = names(df), selected = names(df))
  })
  
  # Reactive data for filtering and showing selected columns
  filtered_data <- reactive({
    df <- current_df()
    if (is.null(df)) return(NULL)
    
    # Only show selected columns
    if (!is.null(input$show_columns)) {
      df <- df[, input$show_columns, drop = FALSE]
    }
    
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
    df
  })
  
  # Render the data table with clickable "Show Image" links
  output$data_table <- renderDT({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    
    # Ensure image columns are associated with their respective directories
    for (i in seq_along(image_columns)) {
      col <- image_columns[i]
      
      if (col %in% names(df)) {
        # Construct relative image path using the resource path alias and show "Show Image" link
        df[[col]] <- ifelse(!is.na(df[[col]]) & grepl("\\.png$", df[[col]]),
                            {
                              relative_img_path <- file.path(paste0("image_dir_", i), basename(df[[col]]))
                              paste0('<a href="#" onclick="Shiny.setInputValue(\'image_click_', i, '\', \'', relative_img_path, '\')">Show Image</a>')
                            },
                            df[[col]])
      }
    }
    
    datatable(df, escape = FALSE, filter = "top", options = list(autoWidth = TRUE))
  }, server = TRUE)
  
  # Dynamically create UI outputs for each image column
  output$image_output_panels <- renderUI({
    img_panels <- lapply(seq_along(image_columns), function(i) {
      column(
        width = 6,
        h4(paste("Image from", image_columns[i])),
        uiOutput(paste0("image_output_", i))
      )
    })
    do.call(fluidRow, img_panels)
  })
  
  # Handle image display for each image column independently
  for (i in seq_along(image_columns)) {
    local({
      column_index <- i
      observeEvent(input[[paste0("image_click_", column_index)]], {
        relative_img_path <- input[[paste0("image_click_", column_index)]]
        print(paste("Image clicked for column", image_columns[column_index], ":", relative_img_path))  # Debugging
        
        output[[paste0("image_output_", column_index)]] <- renderUI({
          tags$img(src = relative_img_path, width = "100%", height = "auto")
        })
      })
    })
  }
  
  # Event to handle volcano plot rendering
  observeEvent(input$plot_volcano, {
    df <- filtered_data()
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
  
  # Event to clear filters
  observeEvent(input$clear_filters, {
    df <- current_df()
    if (is.null(df)) return()
    lapply(names(df), function(col) {
      if (is.numeric(df[[col]])) {
        updateSliderInput(session, paste0("slider_", col), value = c(min(df[[col]], na.rm = TRUE), max(df[[col]], na.rm = TRUE)))
      }
      if (exists(paste0("text_", col), where = input)) {
        updateTextInput(session, paste0("text_", col), value = "")
      }
    })
  })
}

shinyApp(ui = ui, server = server)
