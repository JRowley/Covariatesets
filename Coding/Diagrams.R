# Script to compute the x and y coordinates of sets and to draw them onto a graph.
Boundingcoordinates <- function(first.u, first.l, second.u, second.l){
  
  Bounds <- second.u - first.l
  Bounds <- c(Bounds, second.u - first.u)
  Bounds <- c(Bounds, second.l - first.u)
  Bounds <- c(Bounds, second.l - first.l)
  
  Summary <- data.frame(x = c(first.l, first.u, first.u, first.l), y = Bounds)
  
  return(Summary)
}
data <- read.csv("U:/My work/RStudio working directory/PUMS80M/Intersection bounds.csv")
Bound <- NULL
for(j in 1:max(data$Class)){
  
  i = 1 + 4*(j-1)
  Bound <- rbind(Bound, Boundingcoordinates(data[i,1],data[i+1,1],data[i+2,1],data[i+3,1]))
  
}
rm(i)
rm(j)
Bound <- cbind(Bound, Class = data$Class)
Hull <- NULL
for(j in 1:max(data$Class)){
  
  index = chull(matrix(c(Bound[Bound$Class == j,1], Bound[Bound$Class == j,2]), nrow = 4, ncol = 2))
  Hull <- rbind(Hull, Bound[4 * (j - 1) + index,])
  rm(index)
  
}
rm(j)
for(j in 1:max(data$Class)){
  
  name <- paste("Bound.", j, sep = "")
  assign(name, Hull[Hull$Class == j, 1:3])
  rm(name)
  
}
rm(j)

# Composite.point <- Boundingcoordinates(0.718232486626734, 0.533521253083058, 0.476392823418319, 0.476392823418319)
# Composite.correct <- Boundingcoordinates(0.734414500000000, 0.526428100000000, 0.497670000000000, 0.460281800000000)
# Multi2nd.point <- Boundingcoordinates(0.732695536477967, 0.528653340513828, 0.476392823418319, 0.476392823418319)
# Multi2nd.correct <- Boundingcoordinates(0.734414500000000, 0.526428100000000, 0.497670000000000, 0.460281800000000)
# Genderneg.point <- Boundigncoordinates(0.722931792053168, 0.532568467801628, 0.523891587028043, 0.186017546608178)
# Genderneg.correct <- Boundingcoordinates(0.725358600000000, 0.529430400000000, 0.527021400000000, 0.183937900000000)
# Genderpos.point <- Boundingcoordinates(0.523891587028043, 0.370923035733972, 0.777048409838634, 0.532568467801628)
# Genderpos.correct <- Boundingcoordinates(0.527021400000000, 0.368293800000000, 0.779269800000000, 0.529430400000000)
# Boysneg.point <- Boundingcoordinates(0.732095603926590, 0.529945582586428, 0.523402862498884, 0.179635195048651)
# Boysneg.correct <- Boundingcoordinates(0.734379300000000, 0.527362600000000, 0.527716300000000, 0.176360500000000)
# Boyspos.point <- Boundingcoordinates(0.523402862498884, 0.358151941954759, 0.780611182194186, 0.529945582586428)
# Boyspos.correct <- Boundingcoordinates(0.527716300000000, 0.355986600000000, 0.783674600000000, 0.527362600000000)
# Girlsneg.point <- Boundingcoordinates(0.729805157357061, 0.529396317350470, 0.524435032921947, 0.193114515435265)
# Girlsneg.correct <- Boundingcoordinates(0.731767900000000, 0.526858700000000, 0.528983600000000, 0.190091800000000)
# Girlspos.point <- Boundingcoordinates(0.524435032921947, 0.361524756961608, 0.773086722032889, 0.529396317350470)
# Girlspos.correct <- Boundingcoordinates(0.528983600000000, 0.359398100000000, 0.776316400000000, 0.526858700000000)

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.2, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_line(data = Bound.1, aes(x = x, y=y), colour = "#000000") +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Composite instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Composite instrument.pdf")

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.4, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_line(data = Bound.3, aes(x = x, y=y), colour = "#000000") +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Multiple instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Multiple instrument.pdf")

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.6, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.5, aes(x = x, y=y), fill = "#000000", alpha = 1) +
  geom_polygon(data = Bound.8, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.7, aes(x = x, y=y), fill = "#000000", alpha = 1) +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Gender instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Gender instrument.pdf")

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.10, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.9, aes(x = x, y=y), fill = "#000000", alpha = 1) +
  geom_polygon(data = Bound.12, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.11, aes(x = x, y=y), fill = "#000000", alpha = 1) +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Boys instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Boys instrument.pdf")

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.14, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.13, aes(x = x, y=y), fill = "#000000", alpha = 1) +
  geom_polygon(data = Bound.16, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.15, aes(x = x, y=y), fill = "#000000", alpha = 1) +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Girls instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Girls instrument.pdf")

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.10, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.14, aes(x = x, y=y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.12, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.16, aes(x = x, y=y), fill = "#000000", alpha = 0.5) +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Boys and girls instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Boys and girls instrument.pdf")

p <- ggplot(Hull, aes(x = x, y = y))
q <- p + theme_bw() + 
  theme(text = element_text(family = "CM Roman"),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.ticks = element_blank()) +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = 'black'),
        panel.background = element_blank(),
        axis.line = element_line(colour = 'black')) +
  xlab(expression(rho[0])) +
  ylab(expression(rho[1])) + 
  geom_polygon(data = Bound.4, aes(x = x, y=y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.6, aes(x = x, y = y), fill = "#000000", alpha = 0.5) +
  geom_polygon(data = Bound.8, aes(x = x, y=y), fill = "#000000", alpha = 0.5) +
  geom_hline(yintercept = 0) + 
  xlim(0.3, 0.8) +
  ylim(-0.6, 0.6)
ggsave("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Gender and multiple instrument.pdf", q)
embed_fonts("U:/My work/RStudio working directory/Confidence regions for AngEv98/Pictures/Gender and multiple instrument.pdf")
