M <- readRDS('Data/Markers.rds')

M[1:5,1:5]

M <- M[rownames(M) %in% Gkeep_YT_22_23, ]

pcs = prcomp(M, scale = T, center = T)
str(pcs)

pc.vars <- data.frame(PC.num = 1:length(pcs$sdev),
                      PC= colnames(pcs$x),
                      var=(pcs$sdev)^2/sum((pcs$sdev)^2))

head(pc.vars)

ggplot(pc.vars, 
       aes(x=PC.num,
           y=var))+
  geom_point()+
  geom_line()+
  labs(title="Scree plot of PC variances of marker data")


check_prefixes <- c('Kaskaskia','07','12','16','17','18MSFRS', 'US16', 'US17')
data.for.PC.plot <- data.frame(pcs$x) |>
  mutate(Gen=row.names(pcs$x)) |>
  left_join(
    data.frame(Gen=row.names(M)) |>
      mutate(group = sub("-.*", "", row.names(M)))
    ) |>
  mutate(
    group = if_else(group %in% check_prefixes, 'Z_checks', group)
  )


data.for.PC.plot |>
  select(Gen, group, PC1, PC2) |>
  group_by(group) |>
  summarise(nGen = n_distinct(Gen), .groups = "drop") |>
  glimpse()


ggplot(data.for.PC.plot |> filter(group%in%c('IL2021','IL16LCSDH')),
       aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group), size = 2, alpha = 0.8) + 
  scale_color_viridis_d(option = "H") +
  theme_minimal()

group_summary <- data.for.PC.plot |>
  select(Gen, group, PC1, PC2) |>
  group_by(group) |>
  summarise(nGen = n_distinct(Gen), .groups = "drop")

data_plot_with_summary <- data.for.PC.plot |>
  left_join(group_summary, by = "group")

group_centroids <- data_plot_with_summary |>
  group_by(group) |>
  summarise(PC1 = median(PC1), PC2 = median(PC2), nGen = first(nGen), .groups = "drop")

ggplot(data_plot_with_summary, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group), size = 2, alpha = 0.8) +
  geom_text(
    data = group_centroids,
    aes(label = paste0(group, "\n(n=", nGen, ")")),
    size = 3, fontface = "bold", vjust = -0.5
  ) +
  scale_color_viridis_d(option = "H") +
  theme_minimal()

library(ggrepel)

ggplot(data_plot_with_summary, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = group), size = 1.8, alpha = 0.7) +
  geom_label_repel(
    data = group_centroids,
    aes(label = paste0(group, "\n(n=", nGen, ")")),
    size = 3,
    fontface = "bold",
    box.padding = 0.4,
    max.overlaps = Inf,
    segment.color = "grey30"
  ) +
  scale_color_viridis_d(option = "H") +
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "right",
    panel.grid = element_blank()
  )

save.image('Data/PopStructure.RData')
