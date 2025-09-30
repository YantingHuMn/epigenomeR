calc_ht_size_heatmap <- function(mat, cell_width = 0.5, cell_height = 0.5) {
  if (is.null(rownames(mat)) || is.null(colnames(mat))) {
    stop("Matrix must have rownames and colnames!")
  }

  # Measure label sizes
  row_labels <- rownames(mat)
  col_labels <- colnames(mat)

  # Create a temporary viewport to measure text
  grid.newpage()
  pushViewport(viewport())

  label_width <- max(stringWidth(row_labels))
  label_height <- max(stringHeight(col_labels))

  label_width_inch <- convertWidth(label_width, unitTo = "inches", valueOnly = TRUE)
  label_height_inch <- convertHeight(label_height, unitTo = "inches", valueOnly = TRUE)

  # Pop viewport
  popViewport()

  # Calculate dynamic width and height
  heatmap_width <- length(col_labels) * cell_width + label_width_inch + 1  # add padding
  heatmap_height <- length(row_labels) * cell_height + label_height_inch + 1
  return(c(heatmap_width, heatmap_height))
}
