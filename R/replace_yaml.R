replace_yaml = function(alt_dep, rmd_doc, yaml_doc,
                        out_file = "docs/ici_kt_manuscript_html.Rmd")
{
  # alt_dep = tar_read(manuscript)
  # rmd_doc = "docs/ici_kt_manuscript.Rmd"
  # yaml_doc = "docs/html_document_opts.yaml"
  # out_file = "docs/ici_kt_manuscript_html.Rmd"
  rmd_contents = readLines(rmd_doc)
  yaml_locs = grep("---", rmd_contents)
  yaml_range = seq(yaml_locs[1], yaml_locs[2], 1)
  rmd_noyaml = rmd_contents[-yaml_range]
  
  new_yaml = readLines(yaml_doc)
  new_yaml_delim = c("---", new_yaml, "---")
  rmd_newyaml = c(new_yaml_delim, rmd_noyaml)
  cat(rmd_newyaml, file = out_file, sep = "\n", append = FALSE)
  return(out_file)
}

copy_figure = function(manuscript_dep, from_loc, to_loc)
{
  fs::file_copy(from_loc, to_loc)
  to_loc
}
