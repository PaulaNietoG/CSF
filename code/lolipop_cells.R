df <- data.frame(
  "Patient" = c("P01", "P02", "P03", "P04", "P05", "P06", "P07", "P08", "P09", "P10", "P11", "P12", "P13", "P14", "P10", "P15", "P16", "P12", "P16", "P07", "P12", "P12"),
  "Sample" = c("S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "S21", "S22"),
  "Cells" = c(3238, 3305, 6626, 8616, 7438, 4692, 9547, 3798, 1795, 993, 9118, 210, 6106, 2931, 1118, 4585, 7270, 52, 46, 6448, 165, 47)
)

df$patient_sample <- paste(df$Patient, df$Sample)

# Create the lollipop plot
ggplot(df, aes(x = patient_sample, y = Cells, color = Patient)) +
  geom_segment(aes(xend = patient_sample, yend = 0)) +
  geom_point(size = 5) +
  labs(title = "Number of Cells per Sample",
       x = "Sample",
       y = "Cells") +
  theme_minimal() +
  theme(text = element_text(size = 16)) +
  scale_color_manual(values = as.vector(pals::polychrome())) +
  coord_flip()
