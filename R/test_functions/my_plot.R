my_plot <- function(mat){
  x_names <- rownames(mat)
  ggplot(data = mat) + scale_x_continuous(trans='log10') + geom_pointrange(aes(x = N,y = Mean,ymin = Q1,ymax = Q3),size=0.2, color="blue", fill="white") + geom_path(aes(x = N,y = True,group = 1),color = "red")
}