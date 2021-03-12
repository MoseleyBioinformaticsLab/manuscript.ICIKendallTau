##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title
##' @param transcript_data
##' @return
##' @author rmflight
##' @export
zero_to_na <- function(transcript_data) {

  transcript_na = transcript_data
  transcript_na[transcript_data == 0] = NA
  transcript_na

}
