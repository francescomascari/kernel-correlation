fancy_png <- function(
    plot,
    out_path,
    width = 12.5,
    height = 11.67,
    dpi = 300) {
  # ---------------------------------------------------------------------------
  # Render a plot to a high-quality PNG using the tikz device and LaTeX backend.
  # ---------------------------------------------------------------------------
  #
  # Arguments:
  #   plot     : ggplot or base R plot object – the plot to render.
  #
  #   out_path : character – path to save the final PNG file.
  #
  #   width    : numeric – width of the plot in inches (default 12.5).
  #
  #   height   : numeric – height of the plot in inches (default 11.67).
  #
  #   dpi      : numeric – resolution in dots per inch for the PNG (default 300).
  #
  # Returns:
  #   None – saves a PNG file at the specified `out_path`.

  # create a new temporary working directory
  # to save the intermediate rendering steps (.tex and .pdf)
  # of the final plot.
  tmpdir <- tempfile("tikz_work")
  dir.create(tmpdir)

  # inizialize the .tex and .pdf files
  texfile <- file.path(tmpdir, "tmp_plot.tex")
  pdffile <- file.path(tmpdir, "tmp_plot.pdf")

  # open tikz device in the temporary working directory
  tikz(texfile, standAlone = TRUE, width = width, height = height)

  # print the plot and close it
  print(plot)
  dev.off()

  # compile TeX -> PDF
  tinytex::latexmk(texfile, engine = "pdftex")

  # convert PDF -> PNG
  image_write(image_read_pdf(pdffile, density = dpi), out_path, format = "png")

  # remove the temporary folder and all its content
  unlink(tmpdir, recursive = TRUE)
}
