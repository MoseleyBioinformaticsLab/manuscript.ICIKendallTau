render_docx = function(input_rmd, output_docx, other_input = NULL){
  rmarkdown::render(input_rmd, output_file = output_docx)
}
