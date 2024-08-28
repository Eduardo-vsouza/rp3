library(shiny)
library(DT)
library(ggplot2)
library(plotly)
library(shinyBS)  # For accordion panels

# Arguments from command line
args <- commandArgs(trailingOnly = TRUE)
paths_to_dfs <- strsplit(args[1], ",")[[1]]  # Comma-separated paths to data frames
image_columns <- strsplit(args[2], ",")[[1]]  # Comma-separated image column names
image_dirs <- strsplit(args[3], ",")[[1]]  # Comma-separated directories for each image column
df_aliases <- strsplit(args[4], ",")[[1]]  # Comma-separated aliases for the data frames
pdf_columns <- strsplit(args[5], ",")[[1]]  # Comma-separated PDF column names
pdf_dirs <- strsplit(args[6], ",")[[1]]
# Load data frames

# Load data frames
df_list <- lapply(paths_to_dfs, function(path) {
  read.table(path, sep='\t', header=TRUE, stringsAsFactors = FALSE)
})

# Check that df_aliases is the same length as paths_to_dfs
if (length(df_aliases) != length(paths_to_dfs)) {
  stop("The number of aliases provided must match the number of data frames.")
}

ui <- fluidPage(
  titlePanel("RpS"),
  
  # Custom CSS for resizable images and PDFs
  tags$head(
    tags$style(HTML("
      .resizable-panel {
        overflow: auto;
        resize: both;
        border: 1px solid #ddd;
        padding: 10px;
      }
      .data-table-container {
        height: 400px;
        overflow: auto;
        resize: vertical;
      }
      .zoom-container {
        text-align: center;
        margin-bottom: 10px;
      }
      .zoomable-image {
        max-width: 100%;
        height: auto;
        transition: transform 0.25s ease;
      }
      .resizable-media {
        resize: both;
        overflow: auto;
        border: 1px solid #ddd;
        padding: 10px;
        width: 100%;
        height: auto;
      }
    "))
  ),
  
  fluidRow(
    column(
      width = 2,
      selectInput("selected_df", "Select Data Frame", choices = df_aliases)
    ),
    column(
      width = 10,
      div(
        style = "display: flex; flex-wrap: wrap;",
        uiOutput("column_selector")
      )
    )
  ),
  
  fluidRow(
    column(
      width = 12,
      actionButton("apply_filters", "Apply Filters", class = "btn-primary"),
      actionButton("clear_filters", "Clear Filters", class = "btn-secondary"),
      actionButton("plot_volcano", "Erupt Volcano", class = "btn-success"),
      downloadButton("download_filtered_df", "Download Filtered Data", class = "btn-info") # New button to download
    )
  ),
  
  fluidRow(
    column(
      width = 12,
      bsCollapse(id = "accordion_left", multiple = TRUE,
                 bsCollapsePanel(
                   "Data Table", 
                   div(DTOutput("data_table"), class = "resizable-panel data-table-container"), 
                   style = "info"
                 ),
                 bsCollapsePanel(
                   "Volcano Plot", 
                   div(plotlyOutput("volcano_plot"), class = "resizable-panel"), 
                   style = "info"
                 )
      )
    ),
    column(
      width = 12,
      bsCollapse(id = "accordion_right", multiple = TRUE,
                 bsCollapsePanel(
                   "PGContext", 
                   fluidRow(uiOutput("image_output_panels")), 
                   style = "info"
                 ),
                 bsCollapsePanel(
                   "MSA for homologs across the genome", 
                   fluidRow(uiOutput("pdf_output_panels")), 
                   style = "info"
                 )
      )
    )
  )
)

server <- function(input, output, session) {
  # Register each image directory as a resource path with a unique alias
  for (i in seq_along(image_dirs)) {
    addResourcePath(paste0("image_dir_", i), image_dirs[i])
  }
  
  # Register each PDF directory as a resource path with a unique alias
  for (i in seq_along(pdf_dirs)) {
    addResourcePath(paste0("pdf_dir_", i), pdf_dirs[i])
  }
  
  # Reactive expression to get the currently selected data frame
  current_df <- reactive({
    req(input$selected_df)
    df_index <- match(input$selected_df, df_aliases)
    if (!is.na(df_index) && df_index > 0 && df_index <= length(df_list)) {
      df <- df_list[[df_index]]
      return(df)
    } else {
      return(NULL)
    }
  })
  
  # Update checkbox group for column selector
  output$column_selector <- renderUI({
    df <- current_df()
    if (is.null(df)) return(NULL)
    
    # Divide the checkbox inputs into 4 columns
    column_names <- names(df)
    column_groups <- split(column_names, ceiling(seq_along(column_names) / 4))
    
    tagList(
      lapply(seq_along(column_groups), function(i) {
        column(3, checkboxGroupInput(paste0("show_columns_group_", i), NULL, choices = column_groups[[i]], selected = column_groups[[i]]))
      })
    )
  })
  
  # Reactive data for filtering and showing selected columns
  filtered_data <- reactive({
    df <- current_df()
    if (is.null(df)) return(NULL)
    
    # Combine selected columns from all groups
    column_groups <- split(names(df), ceiling(seq_along(names(df)) / 4))
    selected_columns <- unlist(lapply(seq_along(column_groups), function(i) {
      input[[paste0("show_columns_group_", i)]]
    }), use.names = FALSE)
    
    # Only show selected columns
    if (!is.null(selected_columns)) {
      df <- df[, selected_columns, drop = FALSE]
    }
    
    df
  })
  
  # Render the data table with clickable "Show Image" and "Show PDF" links
  output$data_table <- renderDT({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    
    # Ensure image columns are associated with their respective directories
    for (i in seq_along(image_columns)) {
      col <- image_columns[i]
      if (col %in% names(df)) {
        df[[col]] <- ifelse(!is.na(df[[col]]) & grepl("\\.png$", df[[col]]),
                            paste0('<a href="#" onclick="Shiny.setInputValue(\'image_click_', i, '\', \'', file.path(paste0("image_dir_", i), basename(df[[col]])), '\')">Show Image</a>'),
                            df[[col]])
      }
    }
    
    # Ensure PDF columns are associated with their respective directories
    for (i in seq_along(pdf_columns)) {
      col <- pdf_columns[i]
      if (col %in% names(df)) {
        df[[col]] <- ifelse(!is.na(df[[col]]) & grepl("\\.pdf$", df[[col]]),
                            paste0('<a href="#" onclick="Shiny.setInputValue(\'pdf_click_', i, '\', \'', file.path(paste0("pdf_dir_", i), basename(df[[col]])), '\')">Show PDF</a>'),
                            df[[col]])
      }
    }
    
    datatable(df, escape = FALSE, filter = "top", options = list(autoWidth = TRUE))
  }, server = TRUE)
  
  # Download handler for the filtered data frame
  output$download_filtered_df <- downloadHandler(
    filename = function() {
      paste0(input$selected_df, "_filtered.csv")  # Use .tsv extension for tab-separated files
    },
    content = function(file) {
      df <- filtered_data()
      if (!is.null(df)) {
        write.table(df, file, sep = "\t", row.names = FALSE, quote = FALSE)  # Specify tab separator
      }
    }
  )
  
  
  # Dynamically create UI outputs for each image column
  output$image_output_panels <- renderUI({
    img_panels <- lapply(seq_along(image_columns), function(i) {
      column(
        width = 6,
        h4(paste("Image from", image_columns[i])),
        div(class = "zoom-container",
            actionButton(paste0("zoom_in_", i), "Zoom In", class = "btn-info"),
            actionButton(paste0("zoom_out_", i), "Zoom Out", class = "btn-info")
        ),
        uiOutput(paste0("image_output_", i), class = "resizable-media")
      )
    })
    do.call(fluidRow, img_panels)
  })
  
  # Dynamically create UI outputs for each PDF column
  output$pdf_output_panels <- renderUI({
    pdf_panels <- lapply(seq_along(pdf_columns), function(i) {
      column(
        width = 6,
        h4(paste("PDF from", pdf_columns[i])),
        uiOutput(paste0("pdf_output_", i), class = "resizable-media")
      )
    })
    do.call(fluidRow, pdf_panels)
  })
  
  # Handle image display for each image column independently
  for (i in seq_along(image_columns)) {
    local({
      column_index <- i
      zoom_level <- reactiveVal(1)
      
      observeEvent(input[[paste0("image_click_", column_index)]], {
        relative_img_path <- input[[paste0("image_click_", column_index)]]
        output[[paste0("image_output_", column_index)]] <- renderUI({
          tags$img(src = relative_img_path, class = "zoomable-image", style = paste0("transform: scale(", zoom_level(), ");"))
        })
      })
      
      observeEvent(input[[paste0("zoom_in_", column_index)]], {
        zoom_level(zoom_level() * 1.2)
        output[[paste0("image_output_", column_index)]] <- renderUI({
          tags$img(src = input[[paste0("image_click_", column_index)]], class = "zoomable-image", style = paste0("transform: scale(", zoom_level(), ");"))
        })
      })
      
      observeEvent(input[[paste0("zoom_out_", column_index)]], {
        zoom_level(zoom_level() / 1.2)
        output[[paste0("image_output_", column_index)]] <- renderUI({
          tags$img(src = input[[paste0("image_click_", column_index)]], class = "zoomable-image", style = paste0("transform: scale(", zoom_level(), ");"))
        })
      })
    })
  }
  
  # Handle PDF display for each PDF column independently
  for (i in seq_along(pdf_columns)) {
    local({
      column_index <- i
      observeEvent(input[[paste0("pdf_click_", column_index)]], {
        relative_pdf_path <- input[[paste0("pdf_click_", column_index)]]
        output[[paste0("pdf_output_", column_index)]] <- renderUI({
          tags$iframe(src = relative_pdf_path, class = "resizable-media", width = "100%", height = "600px")
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