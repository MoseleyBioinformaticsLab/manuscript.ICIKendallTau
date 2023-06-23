resize_images = function(figure_path, options){
  new_path = dn_modify_path(figure_path, options)
  tmp_image = image_read(new_path)
  tmp_info = image_info(tmp_image)
  tmp_height = tmp_info$height
  tmp_width = tmp_info$width
  new_width = round(tmp_width * 0.73)
  new_height = round(tmp_height * 0.73)
  out_size = paste0(new_width, "x", new_height)
  out_image = image_scale(tmp_image, out_size)
  image_write(out_image, path = new_path, format = "png")
  new_path
}
