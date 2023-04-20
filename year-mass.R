# Header ------------------------------------------------------------------

# This is a script to take Michigan DNR trawl data and calculate yearly species biomass estimates
# Weighed fish are used to estimate the mass of unweighed fish of the same species
# Relative abundances are estimated using fish weights and counts
# Smoothing plots are used to visualize trends over time at various scales
# An NMDS plot is made to visualize trends in fish community composition through time

# Author: Justin Meyer
# Date: 2023-03-30
# Email: meyer443@purdue.edu
# github data repository: https://github.com/meyer443/DSB_finalproject

# Set directory -----------------------------------------------------------
rm(list=ls())
setwd("C:/Users/justi/OneDrive/Desktop/commdata")

# Load required packages --------------------------------------------------
library(tidyverse)
library(readxl)
library(ggplot2)
library(GGally)
library(vegan)
library(RVAideMemoire)

# create function to use later
`%notin%` <- Negate(`%in%`)

# load data ---------------------------------------------------------------
card1a <- read_xlsx("./SagBay2009-2021.xlsx", sheet = 1)
card1b <- read_xlsx("./SagBay1970-2008.xlsx", sheet = 1)
card2a <- read_xlsx("./SagBay2009-2021.xlsx", sheet = 2)
card2b <- read_xlsx("./SagBay1970-2008.xlsx", sheet = 2)
card3a <- read_xlsx("./SagBay2009-2021.xlsx", sheet = 3)
card3b <- read_xlsx("./SagBay1970-2008.xlsx", sheet = 3) # warnings not about data we care about
card4a <- read_xlsx("./SagBay2009-2021.xlsx", sheet = 5)
card4b <- read_xlsx("./SagBay1970-2008.xlsx", sheet = 5)
card5a <- read_xlsx("./SagBay2009-2021.xlsx", sheet = 6)
card5b <- read_xlsx("./SagBay1970-2008.xlsx", sheet = 6)
sppcodes <- read_xlsx("./SpeciesCodes.xlsx")

# append cards
card1 <- card1a %>% bind_rows(card1b)
card2 <- card2a %>% bind_rows(card2b)
card3 <- card3a %>% bind_rows(card3b)
card4 <- card4a %>% bind_rows(card4b)
card5 <- card5a %>% bind_rows(card5b)

# remove bad tows and blanks
exclusion <- card1 %>% filter(CPE == 2 | CPE %in% NA) %>% select(SYear, Idnum)
card1 <- card1 %>% filter(CPE == 1)

card2 <- card2 %>% anti_join(exclusion)
card3 <- card3 %>% anti_join(exclusion)
card4 <- card4 %>% anti_join(exclusion)
card5 <- card5 %>% anti_join(exclusion)

# Card 5 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# card 5
# isolate age classes, differentiating between weighed and unweighed fish

# convert all unknown ages (0+) to code 5
card5 <- 
  card5 %>%
  mutate(Age = replace(Age, Age == 0, 5)) %>% 
  mutate(Age = replace(Age, Age == 7, 5)) %>% 
  mutate(Age = replace(Age, Age %in% NA, 5))

# Age 0 fish --------------------------------------------------------------

# for 'AGE == 1' age 0 fish with weights
card5.0 <- 
  card5 %>% 
  group_by(SYear, Species) %>% 
  filter(Age == 1) %>% 
  filter(Totwt != 'NA' &
           Totwt != 0) %>% 
  summarise(Age = mean(Age),
            biomass = sum(Totwt),
            numAge = sum(Totnum)) %>% 
  filter(numAge %notin% NA)

# age 0 fish without weights
card5.0.noweight <- 
  card5 %>% 
  group_by(SYear, Species) %>% 
  filter(Age == 1) %>% 
  filter(Totwt == 'NA' |
           Totwt == 0 | 
           Totwt %in% NA) %>% 
  summarise(Age = mean(Age),
            numAge = sum(Totnum)) %>% 
  filter(numAge %notin% NA)

# Adult fish --------------------------------------------------------------

# for 'AGE = 2-4', for age 1, 1+, and 2+ fish with weights
card5.1 <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 2 |
           Age == 3 |
           Age == 4) %>% 
  filter(Totwt != 'NA' &
           Totwt != 0) %>% 
  summarise(biomass = sum(Totwt),
            numAge = sum(Totnum))

# for 'AGE = 2-4' without weights
card5.1.noweight <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 2 |
           Age == 3 |
           Age == 4) %>% 
  filter(Totwt == 'NA' |
           Totwt == 0) %>% 
  summarise(numAge = sum(Totnum))

# Unaged fish -------------------------------------------------------------

# for AGE == 5, age 0+ fish with weights
card5.2 <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 5) %>%
  filter(Totwt != 'NA' &
           Totwt != 0) %>% 
  summarise(biomass = sum(Totwt),
            numAge = sum(Totnum))

# for AGE ==  5 without weights
card5.2.noweight <- 
  card5 %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Age == 5) %>% 
  filter(Totwt == 'NA' |
           Totwt == 0) %>% 
  summarise(numAge = sum(Totnum)) %>% 
  filter(numAge %notin% NA)

# Card 4 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# card 4
# isolate age classes, differentiating between weighed and unweighed fish

# convert all unknown age codes to 5
card4 <- 
  card4 %>%
  mutate(AGE = replace(AGE, AGE == 0, 5))

# there is 1 row where the weight and count = 0, so I am removing it since I can't be sure the fish was actually caught
card4 <- 
  card4 %>%
  filter(NUMBER != 0)

# Age 0 fish --------------------------------------------------------------

# for 'AGE == 1' age 0 fish with weights
card4.0 <- 
  card4 %>% 
  group_by(SYear, SPECIES) %>% 
  filter(AGE == 1) %>% 
  filter(Weight != 'NA' &
           Weight != 0) %>% 
  summarise(AGE = mean(AGE),
            biomass = sum(Weight),
            numAge = sum(NUMBER))

# age 0 fish without weights
card4.0.noweight <- 
  card4 %>% 
  group_by(SYear, SPECIES) %>% 
  filter(AGE == 1) %>% 
  filter(Weight == 'NA' | 
           Weight == 0) %>% 
  summarise(AGE = mean(AGE),
            numAge = sum(NUMBER))

# Adult fish --------------------------------------------------------------

# for 'AGE = 2-4', for age 1, 1+, and 2+ fish with weights
card4.1 <- 
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 2 |
           AGE == 3 |
           AGE == 4) %>% 
  filter(Weight != 'NA' &
           Weight != 0) %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(NUMBER))

# for 'AGE = 2-4' weighout weights
card4.1.noweight <- 
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 2 |
           AGE == 3 |
           AGE == 4) %>% 
  filter(Weight == 'NA' |
           Weight == 0) %>% 
  summarise(numAge = sum(NUMBER))

# Unaged fish -------------------------------------------------------------

# for 'AGE == 5', age 0+ fish with weights
card4.2 <- 
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 5) %>%
  filter(Weight != 'NA' &
           Weight != 0) %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(NUMBER))

# for 'AGE ==  5' without weights
card4.2.noweight <-
  card4 %>% 
  group_by(SYear, SPECIES, AGE) %>% 
  filter(AGE == 5) %>% 
  filter(Weight == 'NA' |
           Weight == 0) %>% 
  summarise(numAge = sum(NUMBER))

# Card 2 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# card 2
# isolate age classes, differentiating between weighed and unweighed fish
# all aged fish have weights

# edit ages so 0, 7, and NA are 5 (age 0+/unknown)
card2 <- 
  card2 %>% 
  mutate(Age = replace(Age, Age == 0, 5)) %>%
  mutate(Age = replace(Age, Age == 7, 5)) %>%
  mutate(Age = replace(Age, Age %in% NA, 5))


# age 0 with weights
card2.0 <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 1) %>% 
  group_by(SYear, Species) %>% 
  filter(Weight != 0) %>% 
  summarise(Age = mean(Age),
            biomass = sum(Weight),
            numAge = sum(Totnum))

# age 0 fish without weights
card2.0.noweight <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 1) %>% 
  group_by(SYear, Species) %>% 
  filter(Weight == 0 |
           Weight == 'NA') %>% 
  summarise(Age = mean(Age),
            numAge = sum(Totnum))

# Adult fish --------------------------------------------------------------

# age 1, 2, 3 with weights
card2.1 <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 2 |
           Age == 3 | 
           Age == 4) %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Weight != 0) %>% 
  summarise(biomass = sum(Weight),
            numAge = sum(Totnum))

# adult fish of unknown weight
card2.1.noweight <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 2 |
           Age == 3 | 
           Age == 4) %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Weight == 0 |
           Weight == 'NA') %>% 
  summarise(numAge = sum(Totnum))

# Unaged fish -------------------------------------------------------------

# unknown ages with weights
card2.2 <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 5) %>% 
  group_by(SYear, Species, Age) %>% 
  filter(Weight != 0) %>% 
  summarise(Age = mean(Age),
            biomass = sum(Weight),
            numAge = sum(Totnum))

# unknown ages no weights
card2.2.noweight <- 
  card2 %>% 
  select(-Maturity, -Sex, -'Inch Group', -'cm Group', -N) %>% 
  distinct() %>% 
  filter(Age == 5) %>% 
  filter(Totnum != 0) %>% 
  filter(Weight == 0 |
           Weight == 'NA') %>% 
  group_by(SYear, Species, Age) %>% 
  summarise(numAge = sum(Totnum))

# Card 3 ------------------------------------------------------------------

# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # #                 Summarize Biomass and Unweighed Fish                # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #
# # # = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = # # #

# Calculate yearly species biomass and counts -----------------------------

# card 3
# isolate age classes, differentiating between weighed and unweighed fish

# age = actual age
# convert ages to match other cards
card3 <- 
  card3 %>% 
  mutate(Age = replace(Age, Age > 4 & Age < 99, 4)) %>%
  mutate(Age = replace(Age, Age == -9, 5)) %>% 
  mutate(Age = replace(Age, Age == 99, 5)) %>% 
  mutate(Age = replace(Age, Age == 3, 4)) %>%
  mutate(Age = replace(Age, Age == 2, 3)) %>% 
  mutate(Age = replace(Age, Age == 1, 2)) %>%
  mutate(Age = replace(Age, Age == 0, 1)) %>% 
  mutate(Age = replace(Age, Age %in% NA, 5))

# All weighed fish --------------------------------------------------------

# all fish with weights
card3biomass <- 
  card3 %>% 
  filter(Weight != 'NA' & 
           Weight != 0) %>% 
  group_by(SYear, Species, Age) %>% 
  summarise(biomass = sum(Weight),
            numAge = n())

# pull out age classes
card3.0 <-
  card3biomass %>%
  filter(Age == 1)

card3.1 <- 
  card3biomass %>%
  filter(Age != 1 & Age != 5)

card3.2 <- 
  card3biomass %>%
  filter(Age == 5)

# Unweighed fish ----------------------------------------------------------

# all fish without weights
card3.noweight <- 
  card3 %>% 
  filter(Weight == "NA" |
           Weight == 0) %>%
  group_by(SYear, Species, Age) %>%
  summarise(numAge = n())

# pull out age classes
card3.0.noweight <- 
  card3.noweight %>% 
  filter(Age == 1)

card3.1.noweight <- 
  card3.noweight %>%
  filter(Age != 1 &
           Age != 5)

card3.2.noweight <- 
  card3.noweight %>% 
  filter(Age == 5)

# Biomass estimates -------------------------------------------------------

card4.0 <- 
  card4.0 %>%
  rename(Age = AGE, 
         Species = SPECIES)

card4.1 <- 
  card4.1 %>%
  rename(Age = AGE, 
         Species = SPECIES)
card4.2 <- 
  card4.2 %>%
  rename(Age = AGE, 
         Species = SPECIES)

card4.0.noweight <- 
  card4.0.noweight %>%
  rename(Age = AGE, 
         Species = SPECIES)

card4.1.noweight <- 
  card4.1.noweight %>%
  rename(Age = AGE, 
         Species = SPECIES)

card4.2.noweight <- 
  card4.2.noweight %>%
  rename(Age = AGE, 
         Species = SPECIES)

# tables ------------------------------------------------------------------

age0.weights <- 
  card5.0 %>% 
  bind_rows(card4.0) %>%
  bind_rows(card3.0) %>% 
  bind_rows(card2.0)

age1.weights <-
  card5.1 %>% 
  bind_rows(card4.1) %>%
  bind_rows(card3.1) %>%
  bind_rows(card2.1)

no.age.weights <- 
  card5.2 %>% 
  bind_rows(card4.2) %>%
  bind_rows(card3.2) %>%
  bind_rows(card2.2)

age0.noweights <- 
  card5.0.noweight %>% 
  bind_rows(card4.0.noweight) %>% 
  bind_rows(card3.0.noweight) %>%
  bind_rows(card2.0.noweight)

age1.noweights <- 
  card5.1.noweight %>% 
  bind_rows(card4.1.noweight) %>% 
  bind_rows(card3.1.noweight) %>%
  bind_rows(card2.1.noweight)

no.age.noweights <- 
  card5.2.noweight %>% 
  bind_rows(card4.2.noweight) %>% 
  bind_rows(card3.2.noweight) %>%
  bind_rows(card2.2.noweight)

# merge weighted and unweighed --------------------------------------------

all.weights <- 
  age0.weights %>%
  bind_rows(age1.weights)

# combine same species/age from different cards
all.weights1 <- 
  all.weights %>% 
  group_by(SYear, Species, Age) %>%
  summarise(biomass = sum(biomass),
            numAge = sum(numAge))

all.noweights <-
  age0.noweights %>% 
  bind_rows(age1.noweights)

# combine same species/age from different cards
all.noweights1 <-
  all.noweights %>% 
  group_by(SYear, Species, Age) %>% 
  summarise(numAge = sum(numAge))

# calculate weight at age -------------------------------------------------

# break by decade

all.weights.decades <- 
  all.weights1 %>% 
  mutate(decade = ifelse(SYear < 1980, '70s',
                         ifelse(SYear > 1979 & SYear < 1990, '80s',
                                ifelse(SYear > 1989 & SYear < 2000, '90s',
                                       ifelse(SYear > 1999 & SYear < 2010, '00s',
                                              ifelse(SYear > 2009 & SYear, '10s', SYear))))))

# year-specific
avg.weight <- 
  all.weights.decades %>% 
  group_by(SYear, Species, Age) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

# decade-specific
avg.weight.decade <-
  all.weights.decades %>% 
  group_by(decade, Species, Age) %>% 
  summarise(avgweight = sum(biomass/sum(numAge)))

# join avg weight to unaged/weighed fish
no.weight <-
  all.noweights1 %>% 
  left_join(avg.weight) %>% 
  mutate(biomass = numAge*avgweight)

no.weight.years <- no.weight %>% filter(biomass %notin% NA) %>% select(-avgweight)

no.weight.fail <- no.weight %>% filter(biomass %in% NA) %>% select(-avgweight, -biomass)

# add decade column to no.weight.fail
no.weight.fail <-
  no.weight.fail %>% 
  mutate(decade = ifelse(SYear < 1980, '70s',
                         ifelse(SYear > 1979 & SYear < 1990, '80s',
                                ifelse(SYear > 1989 & SYear < 2000, '90s',
                                       ifelse(SYear > 1999 & SYear < 2010, '00s',
                                              ifelse(SYear > 2009 & SYear, '10s', SYear))))))

# join decade averages to year-specific fails
no.weight.dec <-
  no.weight.fail %>% 
  left_join(avg.weight.decade) %>% 
  mutate(biomass = numAge*avgweight)

no.weight.decades <- no.weight.dec %>% filter(biomass %notin% NA) %>% select(-decade, -avgweight)
no.weight.decades.fail <- no.weight.dec %>% filter(biomass %in% NA) %>% select(-decade, -avgweight, -biomass)

# species-specific for year-specific and decade-specific failures
# species-specific
avg.weight.sp <- 
  all.weights.decades %>% 
  group_by(Species, Age) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

no.weight.sp <-
  no.weight.decades.fail %>% 
  left_join(avg.weight.sp) %>% 
  mutate(biomass = numAge*avgweight)

no.weight.sp %>% print(n = 39)

no.weight.species <- no.weight.sp %>% filter(biomass %notin% NA) %>% select(-avgweight)
no.weight.sp.fail <- no.weight.sp %>% filter(biomass %in% NA) %>% select(-avgweight, -biomass)

# bind newly weighted species together, then to weighted dataset
new.weights <-
  no.weight.years %>% 
  bind_rows(no.weight.decades) %>% 
  bind_rows(no.weight.species)

fish.weights <-
  new.weights %>% 
  bind_rows(all.weights1)

# now have fish.weights and no.weight.sp.fail

# unaged fish -------------------------------------------------------------

no.age.weights
no.age.noweights

# use age proportions of fish.weights to calculate unaged fish weights
# Unaged fish -------------------------------------------------------------

# age structure
age.proportion <- 
  fish.weights %>%
  group_by(Species, Age) %>%
  summarise(count = sum(numAge)) %>% 
  ungroup() %>%
  group_by(Species) %>% 
  mutate(prop = count/sum(count)) %>% 
  select(-count)

# species-specific
avg.weight.all <- 
  all.weights.decades %>% 
  group_by(Species, Age) %>% 
  summarise(avgweight = sum(biomass)/sum(numAge))

# join with allweights
propweight <- 
  avg.weight.all %>%
  left_join(age.proportion) %>%
  group_by(Species) %>%
  summarise(avgweight = sum(prop*avgweight))

# fill in weights
no.age.newweight <- 
  no.age.noweights %>%
  left_join(propweight) %>%
  mutate(biomass = numAge*avgweight) %>% 
  select(-avgweight)

no.age.weights.all <-
  no.age.newweight %>% filter(biomass %notin% NA) %>% bind_rows(no.age.weights)
no.age.weights.fail <-
  no.age.newweight %>% filter(biomass %in% NA) %>% select(-biomass)

# now have no.age.weights.all, fish.weights, no.age.weights.fail, and no.weight.sp.fail
fish.weights %>% bind_rows(no.age.weights.all)
no.weight.sp.fail %>% bind_rows(no.age.weights.fail)

# ignore unaged fish we can't calculate from this dataset for now (they could be any weight, right?)
# unweighed gobies remaining are from card 5, use avg weight of age 0 and age 2 fish to calculate age 1 weight

goby.weight <- 
  card5 %>%
  filter(Species == 906 & Totwt != 0.000) %>%
  group_by(Age) %>% 
  summarise(avg = sum(Totwt)/sum(Totnum)) %>%
  filter(Age == 1 | Age == 2) %>%
  summarise(avg2 = sum(avg)/2)
# 0.000525

# do the same thing for the bluntnose minnow
minnow.weight <-
  card5 %>% 
  filter(Species == 513 & Totwt != 0.000) %>% 
  group_by(Age) %>% 
  summarise(avg = sum(Totwt)/sum(Totnum)) %>% 
  filter(Age == 1 | Age == 3) %>% 
  summarise(avg2 = sum(avg)/2)
# 0.00526

# add weights to no.weight.sp.fail
no.weight.sp.fail <- no.weight.sp.fail %>% mutate(biomass = ifelse(Species == 513, 0.00526,
                                                                   ifelse(Species == 906, 0.000525,99))) %>% 
  filter(biomass != 99)

# merge dataset to one object
fish.dataset <-
  fish.weights %>% 
  bind_rows(no.age.weights.all) %>% 
  bind_rows(no.weight.sp.fail)

# Standardize effort ------------------------------------------------------

# remove unidentified fish, dreissenid mussels, and the invasive sea lamprey since they are not informative of commmunity productivity
# join species info from 'sppcodes' object
df.main <- 
  fish.dataset %>% 
  left_join(sppcodes) %>%
  filter(Species != 999 &
           Species != 998 &
           Species != 920 &
           Species != 5) %>% 
  mutate(Func_group = str_replace(Func_group, "_", " "))

# convert species id code to character from double
df.main$Species <- as.character(df.main$Species)

# === need to standardize per trawl and log-transform biomass === #

# create trawls/year tibble from card1 then join it to dataframe
trawl <- 
  card1 %>%
  group_by(SYear) %>%
  summarise(towtime = sum(Towtime))

# join
all.trawl <- 
  df.main %>%
  left_join(trawl)

# write_csv(all.trawl, "catch_trawl.csv")

# create new biomass/trawl column
fishpertrawl <-
  all.trawl %>%
  mutate(BioPerTrawl = biomass/towtime)
# fishpertrawl will be used for nmds later

# write_csv(fishpertrawl, "fishpertrawl.csv")

# fishpertrawl <- read_csv("fishpertrawl.csv")

# Bay-wide ----------------------------------------------------------------

# get yearly bay-wide values
yearly <- 
  fishpertrawl %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

# log-transform biomass since fish growth is not linear
forgraph <- 
  yearly %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)

# low sampling effort in 1970 and only yellow perch reported
forgraph <- forgraph %>% filter(SYear != 1970)

# bay-wide biomass/trawl graph
ggplot(forgraph) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5, color = "#CFB53B") +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute") +
  scale_y_continuous(expand =c(0,0)) +
  theme_bw(base_size = 22)
# ggsave(filename = "biomass.jpg", device='jpg', dpi=700)

# not log transformed
ggplot(yearly) +
  aes(x = SYear, y = bioper) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5, color = "#CFB53B") +
  xlab("Year") +
  ylab(expression(paste("Biomass (kg \u00D7 minute"^"-1"*")"))) + 
  ggtitle("Fish biomass per tow minute") +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6), expand=c(0.005,0.005)) + 
  scale_x_continuous(limits = c(1970,2022), expand=c(0.005,0.005)) +
  theme_bw(base_size = 22)
# ggsave(filename = "biomass_untr.jpg", device='jpg', dpi=700)

# Functional group graph ------------------------------------------

# biomass per trawl grouped by feeding functional group
feedingpertrawl <-   
  fishpertrawl %>%
  group_by(SYear, Func_group) %>%
  summarise(biomass = sum(biomass),
            bioper = sum(BioPerTrawl)) %>% 
  filter(Func_group != 'Parasite')
# write_csv(feedingpertrawl, "FuncEffort.csv")

# log-transform bio
# get rid of parasite data, not interested
feedgraph <- 
  feedingpertrawl %>% 
  mutate(Log_Biomass = log10(biomass)+3,
         Log_Biomass_per_Trawl = log10(bioper)+3)

# graph of feeding functional group biomass per trawl
# I added a color-blind friendly color palette after reading the report rubric
# cbbPalette and scale_colour_manual are the only changes from the figure in my class presentation

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#CFB53B", "#0072B2", "#D55E00", "#CC79A7")

ggplot(feedgraph) +
  aes(x = SYear, y = Log_Biomass_per_Trawl, group = Func_group) +
  # geom_point(aes(shape = Func_group, colour = Func_group), size = 2.5) +
  geom_smooth(aes(color = Func_group), se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg \u00D7 minute"^"-1"*")+3"))) +
  ggtitle("Feeding functional group biomass per tow minute") +
  labs(color="Functional Group") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=cbbPalette) +
  scale_y_continuous(expand = c(0,0))
# ggsave(filename = "func.jpg", device='jpg', dpi=700)

# not log transformed
ggplot(feedingpertrawl) +
  aes(x = SYear, y = bioper, group = Func_group) +
  # geom_point(aes(shape = Func_group, colour = Func_group), size = 2.5) +
  geom_smooth(aes(color = Func_group), se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Biomass (kg \u00D7 minute"^"-1"*")"))) +
  ggtitle("Feeding functional group biomass per tow minute") +
  scale_colour_manual(name="", values=cbbPalette) +
  theme_bw(base_size = 22) +
  theme(legend.justification=c(0,1), legend.position=c(0,1.1), legend.background = element_blank(),
        legend.box.background = element_blank(),
        legend.key = element_blank()) +
  scale_y_continuous(breaks=c(0,1,2,3,4,5,6), expand=c(0.005,0.005)) + 
  scale_x_continuous(limits = c(1970,2022), expand=c(0.005,0.005)) +
  guides(color=guide_legend(override.aes=list(fill=NA)))
# ggsave(filename = "func_untr.jpg", device='jpg', dpi=700)

# Pairs plot --------------------------------------------------------------

# create object to compare groups to one another
func.compare <- 
  feedgraph %>% 
  group_by(SYear) %>% 
  select(-bioper, -Log_Biomass, -biomass) %>%
  pivot_wider(names_from="Func_group", values_from="Log_Biomass_per_Trawl") %>%
  filter(SYear != 1970) %>% 
  ungroup() %>% 
  rename('Ben invert' = 'Benthic invertivore', 'Gen invert' = 'General invertivore') %>% 
  select(-SYear)

# plot pairwise correlations
ggpairs(func.compare) +
  ggtitle("Functional Group Biomass Correlations") +
  theme_bw(base_size = 14)
# ggsave(filename = "pairsplot.jpg", device='jpg', dpi=700)

# not log transformed
func.compare2 <- 
  feedingpertrawl %>% 
  group_by(SYear) %>% 
  select(-biomass) %>%
  pivot_wider(names_from="Func_group", values_from="bioper") %>% 
  ungroup() %>% 
  rename('Ben invert' = 'Benthic invertivore', 'Gen invert' = 'General invertivore') %>% 
  select(-SYear)

# plot pairwise correlations
ggpairs(func.compare2) +
  ggtitle("Functional Group Biomass Correlations") +
  theme_bw(base_size = 14)
# ggsave(filename = "pairsplot_untr.jpg", device='jpg', dpi=700)


# Calculate catch-per-unit-effort -----------------------------------------

# cpue all fish
cpue <- 
  fishpertrawl %>%
  ungroup() %>%
  group_by(SYear) %>%
  summarise(cpue = sum(numAge)/towtime) %>%
  unique()

# log-transform
cpue <- 
  cpue %>% 
  mutate(logcpue = log10(cpue)+1)

# plot
ggplot(cpue) +
  aes(x = SYear, y = logcpue) +
  geom_point(size = 2.5) +
  geom_smooth(linewidth = 1.5, se = TRUE) +
  xlab("Year") +
  ylab(expression(paste("log(CPUE (no. fish \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Count per tow minute") +
  theme_bw(base_size = 16)
# ggsave(filename = "cpueplot.jpg", device='jpg', dpi=700)

# not log transformed
ggplot(cpue) +
  aes(x = SYear, y = cpue) +
  geom_point(size = 2.5) +
  geom_smooth(linewidth = 1.5, se = TRUE) +
  xlab("Year") +
  ylab(expression(paste("log(CPUE (no. fish \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Count per tow minute") +
  scale_y_continuous(expand=c(0.005,0.005)) + 
  scale_x_continuous(limits = c(1970,2022), expand=c(0.005,0.005)) +
  theme_bw(base_size = 16)
# ggsave(filename = "cpueplot_untr.jpg", device='jpg', dpi=700)

# NMDS --------------------------------------------------------------------

# create decade column for nmds grouping
fishpertrawl <- 
  fishpertrawl %>% 
  mutate(decade = ifelse(SYear < 1980, '70s',
                         ifelse(SYear > 1979 & SYear < 1990, '80s',
                                ifelse(SYear > 1989 & SYear < 2000, '90s',
                                       ifelse(SYear > 1999 & SYear < 2010, '00s',
                                              ifelse(SYear > 2009 & SYear < 2020, '10s',
                                                     ifelse (SYear > 2019, '20s', SYear)))))),
         Func_group = str_replace(Func_group, "_", " "))

# convert data types
fishpertrawl$Species <- as.character(fishpertrawl$Species)
fishpertrawl$SYear <- as.character(fishpertrawl$SYear)

# Species nmds ------------------------------------------------------------

# in order to help nmds converge, only keep species present in at least 15 years
# sum cpue/tow for each species-year combination and pivot wider for matrix conversion
fishpertrawl %>% group_by(SYear, Species) %>% summarise(count = n()) %>% ungroup() %>% group_by(Species) %>% summarise(count = n()) %>% filter(count >= 15) %>% print(n = 27)

fish <- 
  fishpertrawl %>% 
  filter(Species == 106 |
           Species == 108 |
           Species == 109 |
           Species == 113 |
           Species == 119 |
           Species == 131 |
           Species == 132 | 
           Species == 133 |
           Species == 134 | 
           Species == 203 |
           Species == 401 |
           Species == 402 | 
           Species == 403 |
           Species == 405 |
           Species == 413 | 
           Species == 504 | 
           Species == 508 |
           Species == 510 |
           Species == 511 | 
           Species == 601 | 
           Species == 604 |
           Species == 706 | 
           Species == 707 |
           Species == 801 |
           Species == 803 |
           Species == 906) %>% 
  group_by(SYear, Common_name, decade) %>% 
  summarise(BioPerTrawl = sum(BioPerTrawl)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Common_name", values_from = BioPerTrawl)

fish2 <- 
  fishpertrawl %>% 
  filter(Common_name %notin% NA) %>% 
  group_by(SYear, Common_name, decade) %>% 
  summarise(BioPerTrawl = sum(BioPerTrawl)) %>% 
  ungroup() %>% 
  pivot_wider(names_from = "Common_name", values_from = BioPerTrawl)

# NAs are for years where certain species were not caught, change this to 0 for relative abundance estimates in nmds
# object created will be used for nmds plot creation
fish <- 
  fish %>%
  replace(is.na(.), 0)

fish2 <- 
  fish2 %>%
  replace(is.na(.), 0)

# add 1 before log transforming
fish3 <- 
  fish %>% 
  pivot_longer(cols = 3:28, names_to = 'name', values_to = 'bio') %>% 
  mutate(bio = bio + 1) %>% 
  mutate(logbio = log10(bio)) %>% 
  select(-bio) %>% 
  filter(SYear != 1970) %>% 
  pivot_wider(names_from = 'name', values_from = logbio)  

# select only species data
com <- 
  fish %>%
  select(-SYear, -decade) %>% 
  as.matrix()

com2 <- 
  fish2 %>%
  select(-SYear, -decade) %>% 
  as.matrix()

com3 <-
  fish3 %>% 
  select(-SYear, -decade) %>% 
  as.matrix()

# set variable for nmds
set.seed(123)

# run Bray-Curtis nmds on matrix
nmds <- metaMDS(com, distance = 'euclid')
nmds2 <- metaMDS(com2, distance = 'euclid')
nmds3 <- metaMDS(com3, distance = 'euclid')

# generate nmds plot with ellipses and labels
x <- -1.5:1.5
y <- -2:1

# common species
# tiff("ellipse_species.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
plot(x,y, type = "n", main='Species NMDS', xlab='Axis 1', ylab='Axis 2')
ordiellipse(nmds, fish$decade, draw = "polygon", label = TRUE)
orditorp(nmds,display="species",col="red",air=1,cex=1)
# dev.off()

# generate nmds plot with ellipses and labels
# all species
# tiff("ellipse_all_species.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
plot(nmds2, type = "n", main='Every Species NMDS', xlab='', ylab='')
orditorp(nmds2,display="species",col="black")
ordiellipse(nmds2, fish2$decade, draw = "polygon", col="#CFB53B", label = TRUE)
title(ylab = "Axis 2", line = 2)   
title(xlab = "Axis 1", line = 2) 
# dev.off()

# generate nmds plot with arrows
x <- -2:3
y <- -3:2
# all species
# tiff("arrows_all_species.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
plot(x,y, type = "n", main='Every Species NMDS', xlab='', ylab='', cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
# points(nmds2, pch = 20)
# orditorp(nmds2, display = "species", labels = "n", pch = 18)
ordiarrows(nmds2, fish2$decade, lwd=5, col=cbbPalette, label=TRUE, cex=1.1)
title(ylab = "Axis 2", line = 2, cex.lab=1.5)   
title(xlab = "Axis 1", line = 2, cex.lab=1.5)
points(x = 0.5, y = 0.2, pch = 13, cex=1.5)
points(x = 0, y = -0.4, pch = 13, cex=1.5)
points(x = 0.6, y = 0.4, pch = 13, cex = 1.5)
points(x = 0.4, y = 0.5, pch = 13, cex=1.5)
points(x=-0.9, y = -0.4, pch = 13, cex = 1.5)
points(x = 0.2, y = 0.2, pch = 13, cex = 1.5)
points(x = 1.2, y = 0, pch = 13, cex = 1.5)
# dev.off()

# generate nmds plot with transformed data
plot(nmds3)
x <- -1.5:1.5
y <- -1.5:1.5
# all species
# tiff("arrows_all_species.tiff", width = 9.85, height = 6.32, units = 'in', res = 300)
plot(x,y, type = "n", main='Every Species NMDS', xlab='', ylab='')
points(nmds3, pch = 19)
ordiarrows(nmds3, fish2$decade, lwd=3, col=cbbPalette, label=TRUE, cex=1.1)
title(ylab = "Axis 2", line = 2)   
title(xlab = "Axis 1", line = 2) 
# dev.off()

# PERMANOVA
com.dist <- vegdist(com, mothod='euclid')
set.seed(36)
com.div <- adonis2(com.dist~decade, data=fish, permutations = 999, method = 'euclid', stata='PLOT')
com.div

# plot distance from centroid (visualize permanova)
dispersion <- betadisper(com.dist, group=fish$decade)
permutest(dispersion)
plot(dispersion, hull=FALSE, ellipse=TRUE)

# pairwise permanova by decades
pairwise.perm.manova(com.dist, fish$decade, nperm=999, p.method='bonferroni', F=TRUE, R2=TRUE)

# simper of decades
simper <- with(fish, simper(com, fish$decade))
summary(simper)

# dput(summary(simper)$'70s_80s', 'table.txt')

# PERMANOVA
com.dist2 <- vegdist(com2, mothod='euclid')
set.seed(36)
com.div2 <- adonis2(com.dist2~decade, data=fish2, permutations = 999, method = 'euclid', stata='PLOT')
com.div2

# plot distance from centroid (visualize permanova)
dispersion2 <- betadisper(com.dist2, group=fish2$decade)
permutest(dispersion2)
plot(dispersion2, hull=FALSE, ellipse=TRUE)

# pairwise permanova by decades
# significant difference in dispersion found, pairwise test not appropriate
pairwise.perm.manova(com.dist2, fish2$decade, nperm=999, p.method='bonferroni', F=TRUE, R2=TRUE)

# simper of decades
simper2 <- with(fish2, simper(com2, fish2$decade))
summary(simper2)

# dput(summary(simper)$'70s_80s', 'table_all.txt')

# PERMANOVA
com.dist3 <- vegdist(com3, mothod='euclid')
set.seed(36)
com.div3 <- adonis2(com.dist3~decade, data=fish3, permutations = 999, method = 'euclid', stata='PLOT')
com.div3

# plot distance from centroid (visualize permanova)
dispersion3 <- betadisper(com.dist3, group=fish3$decade)
permutest(dispersion3)
plot(dispersion3, hull=FALSE, ellipse=TRUE)

str(dispersion3)

# pairwise permanova by decades
# actually shouldn't do this because we rejected the null from dispersion test
pairwise.perm.manova(com.dist3, fish3$decade, nperm=999, p.method='bonferroni', F=TRUE, R2=TRUE)

# simper of decades
simper3 <- with(fish3, simper(com3, fish3$decade))
summary(simper2)

# Functional graphs -------------------------------------------------------

# piscivores
pisc <- fishpertrawl %>% filter(Func_group == 'Piscivore') %>% filter(Common_name != 'Longnose Gar' & Common_name != 'Rock Bass' & Common_name != 'Smallmouth Bass' & Common_name != 'Largemouth Bass') %>% group_by(SYear, Common_name) %>% summarise(bio = log10(sum(BioPerTrawl))+5)

cbbPalette.pisc <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00")

ggplot(pisc) +
  aes(x = SYear, y = bio, group = Common_name) +
  geom_line(aes(color = Common_name), linewidth = 1.5) +
  geom_point(aes(colour = Common_name), size = 1.5) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+5"))) +
  ggtitle("Piscivore biomass per tow minute") +
  labs(color="Species") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=cbbPalette.pisc)
# ggsave(filename = "piscivore.jpg", device='jpg', dpi=700)

# walleye, white bass, northern pike

# general invertivores
gen <- fishpertrawl %>% filter(Func_group == 'General invertivore') %>% filter(Common_name != 'Bluegill' & Common_name != 'Ninespine Stickleback' & Common_name != 'Pumpkinseed') %>% group_by(SYear, Common_name) %>% summarise(bio = log10(sum(BioPerTrawl))+5)

ggplot(gen) +
  aes(x = SYear, y = bio, group = Common_name) +
  geom_line(aes(color = Common_name), linewidth = 1.5) +
  geom_point(aes(colour = Common_name), size = 1.5) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+5"))) +
  ggtitle("General invertivore biomass per tow minute") +
  labs(color="Species") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=cbbPalette)
# ggsave(filename = "gen_invert.jpg", device='jpg', dpi=700)

# yellow perch, trout-perch, spottail shiner, mimic shiner more recently

# benthic invertivores
ben <- fishpertrawl %>% filter(Func_group == 'Benthic invertivore') %>% filter(Common_name != 'Unidentified Chubs' & Common_name != 'Unidentified Redhorse') %>% group_by(SYear, Common_name) %>% summarise(bio = log10(sum(BioPerTrawl))+5)

ggplot(ben) +
  aes(x = SYear, y = bio, group = Common_name) +
  geom_line(aes(color = Common_name), linewidth = 1.5) +
  geom_point(aes(colour = Common_name), size = 2.5) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+5"))) +
  ggtitle("Benthic invertivore biomass per tow minute") +
  labs(color="Species") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom")
  # scale_colour_manual(values=cbbPalette)
# ggsave(filename = "ben_invert.jpg", device='jpg', dpi=700)

# white perch, white sucker, freshwater drum, round goby

# omnivores
omn <- fishpertrawl %>% filter(Func_group == 'Omnivore') %>% filter(Common_name != 'Unidentified Minnow' & Common_name != 'Goldfish x Carp Hybrid') %>% group_by(SYear, Common_name) %>% summarise(bio = log10(sum(BioPerTrawl))+6)

ggplot(omn) +
  aes(x = SYear, y = bio, group = Common_name) +
  geom_line(aes(color = Common_name), linewidth = 1.5) +
  geom_point(aes(colour = Common_name), size = 1.5) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+6"))) +
  ggtitle("Omnivore biomass per tow minute") +
  labs(color="Species") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") 
# scale_colour_manual(values=cbbPalette)
# ggsave(filename = "omnivore.jpg", device='jpg', dpi=700)

# carp, channel catfish

# planktivores
plank <- fishpertrawl %>% filter(Func_group == 'Planktivore') %>% group_by(SYear, Common_name) %>% summarise(bio = log10(sum(BioPerTrawl))+5)

cbbPalette.plank <- c("#000000", "#E69F99", "#56B4E9", "#009E73", "#F0E442", "#D55E00")

ggplot(plank) +
  aes(x = SYear, y = bio, group = Common_name) +
  geom_line(aes(color = Common_name), linewidth = 1.5) +
  geom_point(aes(colour = Common_name), size = 1.5) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+5"))) +
  ggtitle("Planktivore biomass per tow minute") +
  labs(color="Species") +
  theme_bw(base_size = 15) +
  theme(legend.position = "bottom") +
  scale_colour_manual(values=cbbPalette.plank)
# ggsave(filename = "planktivore.jpg", device='jpg', dpi=700)

# alewife, rainbow smelt, gizzard shad

# Bay-wide w/o walleye ----------------------------------------------------
nowall <- fishpertrawl %>% filter(Common_name != 'Walleye') %>% group_by(SYear) %>%
  summarise(biomass = log10(sum(BioPerTrawl))) %>% filter(SYear != 1970)

ggplot(nowall) +
  aes(x = SYear, y = biomass) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg " %*% " minute"^"-1"*")"))) +
  ggtitle("Fish biomass per tow minute (no walleye)") +
  theme_bw(base_size = 16)
# ggsave(filename = "biomass_nowalleye.jpg", device='jpg', dpi=700)
# same trend as with walleye, fish stocking not influential in bay-wide trend

early <- 
  fishpertrawl %>%
  filter(SYear < 1995 & SYear != 1970) %>% 
  group_by(SYear) %>%
  summarise(biomass = log10(sum(BioPerTrawl)))

ggplot(early) +
  aes(x = SYear, y = biomass) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg " %*% " minute"^"-1"*")"))) +
  ggtitle("Pre-1995") +
  theme_bw(base_size = 16)
# ggsave(filename = "early_trend.jpg", device='jpg', dpi=700)

late <- 
  fishpertrawl %>%
  filter(SYear > 2005) %>%
  group_by(SYear) %>%
  summarise(biomass = log10(sum(BioPerTrawl)))

ggplot(late) +
  aes(x = SYear, y = biomass) +
  geom_point(size = 2) +
  geom_smooth(method = 'lm', se = TRUE, linewidth = 1.5) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg " %*% " minute"^"-1"*")"))) +
  ggtitle("Post-2005") +
  theme_bw(base_size = 16)
# ggsave(filename = "late_trend.jpg", device='jpg', dpi=700)

earl <- lm(biomass ~ SYear, data = early)
noche <- lm(biomass ~ SYear, data = late)

compare.coeff <- function(b1, se1, b2, se2){
  return((b1-b2)/sqrt(se1^2+se2^2))
}
b1 <- summary(earl)$coefficients[2,1]
b2 <- summary(noche)$coefficients[2,1]
se1 <- summary(earl)$coefficients[2,2]
se2 <- summary(noche)$coefficients[2,2]

p_value <- 2*pnorm(-abs(compare.coeff(b1,se1,b2,se2)))
p_value

cor(x = early$SYear, y = early$biomass)^2
cor(x = late$SYear, y = late$biomass)

# Species-specific plots --------------------------------------------------

# read in fishpertrawl
fishpertrawl <- read_csv("fishpertrawl.csv")

# fill in for species you want to remove
spp <- fishpertrawl %>% filter(Common_name != 'Alewife' & Common_name != 'Walleye' & Common_name != 'Carp')
yoy <- fishpertrawl %>% filter(Age != 2) #?

yearly.spp <- 
  spp %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

# log-transform biomass since fish growth is not linear
forgraph.spp <- 
  yearly.spp %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)

# low sampling effort in 1970 and only yellow perch reported
forgraph.spp <- forgraph.spp %>% filter(SYear != 1970)

# bay-wide biomass/trawl graph
ggplot(forgraph.spp) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (alewife, carp, and walleye removed)")
  
ggsave(filename = "carp_ale_wall.jpg", device='jpg', dpi=700)

yearly.yoy <- 
  yoy %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

# log-transform biomass since fish growth is not linear
forgraph.yoy <- 
  yearly.yoy %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)

# low sampling effort in 1970 and only yellow perch reported
forgraph.yoy <- forgraph.yoy %>% filter(SYear != 1970)

# bay-wide biomass/trawl graph
ggplot(forgraph.yoy) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (Age-1 removed)")

# ggsave(filename = "age1.jpg", device='jpg', dpi=700)

alecarp <- fishpertrawl %>% filter(Common_name == 'Carp' | Common_name == 'Alewife' | Common_name == 'Walleye')

yearly.alecarp <- 
  alecarp %>%
  group_by(SYear, Common_name) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

# log-transform biomass since fish growth is not linear
forgraph.alecarp <- 
  yearly.alecarp %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)

# low sampling effort in 1970 and only yellow perch reported
forgraph.alecarp <- forgraph.alecarp %>% filter(SYear != 1970)

cbbPalette.alecarp <- c("#000000", "#D55E00", "#56B4E9")

# bay-wide biomass/trawl graph
ggplot(forgraph.alecarp) +
  aes(x = SYear, y = Log_Biomass_per_Trawl, group = Common_name, color = Common_name) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("Log biomass (kg \u00D7 minute"^"-1"*")"))) +
  ggtitle("Fish biomass per tow minute (Age-0 removed)") +
  scale_colour_manual(values=cbbPalette.alecarp) +
  guides(color = guide_legend(title = "Species")) +
  theme(legend.position = 'bottom')

# ggsave(filename = "wallalecarp.jpg", device='jpg', dpi=700)

# rank abundance
cpue.all <- 
  fishpertrawl %>% 
  group_by(SYear, Common_name) %>%
  summarise(cpue = sum(numAge)/towtime) %>%
  unique()

# log-transform
cpue.all <-
  cpue.all %>% 
  mutate(logcpue = log10(cpue)+3) %>% arrange(logcpue)

cpue.log <- 
  cpue.all %>% 
  group_by(SYear) %>%
  filter(Common_name != 'NA' & logcpue != -Inf) %>%
  summarise(catch = (sum(logcpue)))


# plot
ggplot(cpue.log) +
  aes(x = reorder(Common_name, -avg), y = avg, fill=factor(ifelse(Common_name == "Carp" | Common_name == 'Alewife' | Common_name == 'Walleye',"Carp","Others"))) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(name = "Species", values=c("red","grey50")) +
  ggtitle("Average species CPUE") +
  ylab(expression(paste("Log CPUE (no. fish \u00D7 minute"^"-1"*")"))) +
  xlab("Species") +
  theme_bw(base_size = 16) +
  theme(legend.position = 'bottom')

# ggsave(filename = "counts.jpg", device='jpg', dpi=700)

# without functional ------------------------------------------------------

# piscivores
no.pisc <- 
  fishpertrawl %>% 
  filter(Func_group != 'Piscivore') %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

no.pisc.log <- 
  no.pisc %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)


forgraph.pisc <- no.pisc.log %>% filter(SYear != 1970)

ggplot(forgraph.pisc) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (piscivores removed)")

# ggsave(filename = "nopisc.jpg", device='jpg', dpi=700)

# omnivores
no.omni <- 
  fishpertrawl %>% 
  filter(Func_group != 'Omnivore') %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

no.omni.log <- 
  no.omni %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)


forgraph.omni <- no.omni.log %>% filter(SYear != 1970)

ggplot(forgraph.omni) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (omnivores removed)")

# ggsave(filename = "noomni.jpg", device='jpg', dpi=700)

# general invertivores
no.gen <- 
  fishpertrawl %>% 
  filter(Func_group != 'General invertivore') %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

no.gen.log <- 
  no.gen %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)


forgraph.gen <- no.gen.log %>% filter(SYear != 1970)

ggplot(forgraph.gen) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (general invertivores removed)")

# ggsave(filename = "nogen.jpg", device='jpg', dpi=700)

# benthic invertivores
no.ben <- 
  fishpertrawl %>% 
  filter(Func_group != 'Benthic invertivore') %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

no.ben.log <- 
  no.ben %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)


forgraph.ben <- no.ben.log %>% filter(SYear != 1970)

ggplot(forgraph.omni) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (benthic invertivores removed)")

# ggsave(filename = "noben.jpg", device='jpg', dpi=700)

# planktivores
no.plank <- 
  fishpertrawl %>% 
  filter(Func_group != 'Planktivore') %>%
  group_by(SYear) %>%
  summarise(biomass = sum(biomass), 
            bioper = sum(BioPerTrawl))

no.plank.log <- 
  no.plank %>% 
  mutate(Log_Biomass = log10(biomass)+1,
         Log_Biomass_per_Trawl = log10(bioper)+1)


forgraph.plank <- no.plank.log %>% filter(SYear != 1970)

ggplot(forgraph.plank) +
  aes(x = SYear, y = Log_Biomass_per_Trawl) +
  geom_point(size = 2) +
  geom_smooth(se = TRUE, linewidth = 1.5) +
  theme_bw(base_size = 16) +
  xlab("Year") +
  ylab(expression(paste("log(biomass (kg \u00D7 minute"^"-1"*"))+1"))) +
  ggtitle("Fish biomass per tow minute (planktivores removed)")

# ggsave(filename = "noplank.jpg", device='jpg', dpi=700)
