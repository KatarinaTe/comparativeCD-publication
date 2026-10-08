library(tidyverse)
library(googlesheets4)

gs4_auth(email = "elinor@broadinstitute.org")

url <- "<Google Sheet link removed>"

nc2 <- read_tsv("name_conversion_v2.tsv", show_col_types = FALSE) %>%
  mutate(include = as.logical(include))

# Create the new sheet if it doesn't exist, then write
sheet_name <- "name conversion"
existing <- sheet_names(url)
if (!sheet_name %in% existing) {
  sheet_add(url, sheet = sheet_name)
}

range_clear(ss = url, sheet = sheet_name, reformat = FALSE)
range_write(nc2, ss = url, sheet = sheet_name, range = "A1", col_names = TRUE, reformat = FALSE)

message("Done. Wrote ", nrow(nc2), " rows to '", sheet_name, "'")
message("Big Five breakdown (included rows):")
nc2 %>%
  filter(include) %>%
  filter(short %in% c("Neuroticism","Extraversion","Openness","Conscientiousness","Agreeableness","Other temperament")) %>%
  count(short, sort = TRUE) %>%
  print()
