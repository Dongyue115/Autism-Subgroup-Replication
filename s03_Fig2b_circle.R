library(R.matlab)
library(readxl)
library(igraph)
library(ggraph)
library(scales)  # For controlling color scales


orderT <- read_excel("rrb_ccT_0.001.xlsx")
name <- orderT$Group_Buch
net_order <- orderT$network_order

fc <- readMat("rrb_circle_0.001.mat")$fc.dis

rownames(fc) <- name
colnames(fc) <- name

gg <- graph_from_adjacency_matrix(fc, mode = "undirected", weighted = TRUE)

ldfr <- length(name)
angles <- 90 - 360 * 0:(ldfr-1) / ldfr
hjust <- ifelse(angles < -90, 1, 0)
angle <- ifelse(angles < -90, angles + 180, angles)

V(gg)$net_order <- net_order

# Define specific colors for specific net_order values
specific_colors <- c(
  "1" = rgb(0.2902, 0.2039, 0.5882),
  "2" = rgb(0.3137, 0.3569, 0.6667),
  "3" = rgb(0.0235, 0.5804, 0.8235),
  "4" = rgb(0.2863, 0.7412, 0.5961),
  "5" = rgb(0.9373, 0.9216, 0.4),
  "6" = rgb(0.9686, 0.6863, 0.2157),
  "7" = rgb(0.9529, 0.3961, 0.5686),
  "8" = rgb(0.7176, 0.2039, 0.4275),
  "9" = rgb(0.4078, 0.1294, 0.3961),
  "10" = rgb(0.2549, 0.1294, 0.3922)
)

custom_edge_colors <- c(rgb(0.4392,0.8039,0.8667),rgb(0.2667,0.4706,0.7373),rgb(0.2275,0.3255,0.6431),
                        rgb(0.1843,0.2549,0.6039),rgb(0.1569,0.1647,0.4549),rgb(0,0,0),
                        rgb(0.498,0.0784,0.0863),rgb(0.749,0.1255,0.1451),rgb(0.9294,0.1333,0.1412),
                        rgb(0.9608,0.502,0.1255),rgb(0.9647,0.9216,0.0863))


color_positions <- seq(-0.4, 0.4, length.out = length(custom_edge_colors))

gg_filtered <- delete_edges(gg, E(gg)[weight > -0.3 & weight < 0.3]) #not display edge within -0.3~0.3

a <- ggraph(gg_filtered, layout = 'linear', circular = TRUE) +
   geom_edge_arc(aes(color = weight), show.legend = TRUE, edge_width = 1.5) +
  geom_node_point(aes(x = x*1.05, y = y*1.05, color = factor(net_order)), size = 12) +
  geom_node_text(aes(x = x*1.14, y = y*1.14, label = name, angle = angle, hjust = hjust), size = 8) +  # Node labels
  scale_color_manual(values = specific_colors) +
  scale_edge_color_stepsn(colors = custom_edge_colors,        
                          values = rescale(color_positions),  
                          name = "Edge Weight",              
                          breaks = seq(-0.4, 0.4, by = 0.6/11),
                          limits = c(-0.4, 0.4),   
                          oob = squish) +
  theme_void() +
  expand_limits(x = c(-1.8, 1.8), y = c(-1.8, 1.8))# Adjust circle size
                                        

ggsave("FC_loading_fig/rrb39_0.001_filter0.3.png", plot = a, width = 10, height = 8, dpi = 300)
