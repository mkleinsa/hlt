
logo = function() {
  df = data.frame(x = c("h", "l", "t"), y = c(1, 0.8, 0.7))
  p <- ggdotchart(df, x = "x", y = "y",
                  color = "x",
                  shape = 15,
                  palette = "jco",
                  add = "segments",                   
                  rotate = TRUE,               
                  dot.size = 4,                  
                  label = df$x,  
                  sort = "desc",
                  font.label = list(color = "white", size = 4,
                                    vjust = 0.5),             
                  ggtheme = theme_void()) + 
    theme_transparent() + theme(legend.position = "none") 
  p
  
  sticker(p, package="", p_size=4.5, s_x=0.9, s_y=1, s_width=1.7, s_height=1.3,
          p_x = 1.1, p_y = 0.9,
          url = "https://cran.r-project.org/package=hlt", u_color = "white", u_size = 1,
          h_fill="black", h_color="grey",
          filename="man/figures/logo.png")
}

