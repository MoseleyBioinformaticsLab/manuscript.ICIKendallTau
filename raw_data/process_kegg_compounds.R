process_compound_file = function(in_file)
{
  # in_file = "data/kegg_compounds/C00043.txt"
  in_compound = readLines(in_file)
  start_field = grepl("^[A-Z]", in_compound)
  split_entries = strsplit(in_compound, "[ ]{2, }")
  
  section_starters = purrr::map_lgl(split_entries, \(x){nchar(x[1]) > 0})
  use_starters = which(section_starters)
  compound_sections = vector("list", length(use_starters))
  
  isection = 1
  while (isection <= length(use_starters)) {
    jsection = isection + 1
    if (jsection > length(use_starters)) {
      use_range = seq(use_starters[isection], length(split_entries))
    } else {
      use_range = seq(use_starters[isection], use_starters[jsection] - 1)
    }
    
    section_entry = split_entries[use_range]
    section_id = section_entry[[1]][1]
    names(compound_sections)[isection] = section_id 
    section_entry[[1]][1] = ""
    
    if (section_id %in% "REACTION") {
      new_entry = strsplit(unlist(section_entry), " ", fixed = TRUE) |> unlist() |> unique()
      new_entry = new_entry[nchar(new_entry) > 0]
      section_entry = new_entry
    } else if (section_id %in% c("PATHWAY", "MODULE", "NETWORK")) {
      new_entry = purrr::map(section_entry, \(in_entry){tibble::tibble(id = in_entry[2],
                                                                       description = in_entry[3])}) |>
        dplyr::bind_rows()
      section_entry = new_entry
    } else {
      section_entry = unique(unlist(section_entry))
      section_entry = section_entry[nchar(section_entry) > 0]
    }
    compound_sections[[isection]] = section_entry
    isection = isection + 1
  }
  compound_sections
}

compound_files = dir("data/kegg_compounds", pattern = "txt$", full.names = TRUE)

compound_data = purrr::map(compound_files, process_compound_file)

compounds_from_filename = gsub(".txt", "", basename(compound_files))
compounds_from_data = purrr::map_chr(compound_data, \(x){x$ENTRY[1]})

all.equal(compounds_from_data, compounds_from_filename)

names(compound_data) = compounds_from_data
saveRDS(compound_data, file = "data/kegg_compound_info.rds")

all_meta = purrr::map(compound_data, \(in_compound){
  dplyr::bind_rows(in_compound[c("PATHWAY", "MODULE", "NETWORK")])
}) |>
  dplyr::bind_rows() |>
  dplyr::distinct() |>
  dplyr::mutate(type = dplyr::case_when(
    grepl("^map", id) ~ "pathway",
    grepl("^nt", id) ~ "network",
    grepl("^M", id) ~ "module"
  ))

compound_2_pathway = purrr::map(compound_data, \(in_compound){
  if (!is.null(in_compound$PATHWAY)) {
    out_res = in_compound$PATHWAY |>
      dplyr::select(-description) |>
      dplyr::mutate(compound = in_compound$ENTRY[1])
    return(out_res)
  } else {
    return(NULL)
  }
  
}) |>
  dplyr::bind_rows() |>
  dplyr::distinct()

compound_2_network = purrr::map(compound_data, \(in_compound){
  if (!is.null(in_compound$NETWORK)) {
    out_res = in_compound$NETWORK |>
      dplyr::select(-description) |>
      dplyr::mutate(compound = in_compound$ENTRY[1])
    return(out_res)
  } else {
    return(NULL)
  }
  
}) |>
  dplyr::bind_rows() |>
  dplyr::distinct()

compound_2_module = purrr::map(compound_data, \(in_compound){
  if (!is.null(in_compound$MODULE)) {
    out_res = in_compound$MODULE |>
      dplyr::select(-description) |>
      dplyr::mutate(compound = in_compound$ENTRY[1])
    return(out_res)
  } else {
    return(NULL)
  }
  
}) |>
  dplyr::bind_rows() |>
  dplyr::distinct()

compound_mappings = list(meta = all_meta,
                         network = compound_2_network,
                         module = compound_2_module,
                         pathway = compound_2_pathway)
saveRDS(compound_mappings, file = "data/kegg_compound_mapping.rds")
