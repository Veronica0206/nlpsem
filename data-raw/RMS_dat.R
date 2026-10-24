# Derivation of the RMS_dat example data set from the public-use ECLS-K:2011 files.
# The three SAS files (ecls_k_2011_assessment.sas7bdat, demo_ses.sas7bdat, ecls_k_e2.sas7bdat) are not
# redistributed with the package; obtain them from https://nces.ed.gov/ecls/dataproducts.asp and run this
# script from the folder that contains them. It writes RMS_dat.RData (copy to data/) and two figures.
# Kept verbatim from the original preparation script apart from this header.

library(dplyr)
library(tidyr)
library(ggplot2)

Access0 <- sas7bdat::read.sas7bdat(file = "ecls_k_2011_assessment.sas7bdat", debug = F)
names(Access0)

Info0 <- sas7bdat::read.sas7bdat(file = "demo_ses.sas7bdat", debug = F)
names(Info0)

Skill0 <- sas7bdat::read.sas7bdat(file = "ecls_k_e2.sas7bdat", debug = F)
names(Skill0)

dat0 <- Access0[, c(1, 40:74, 120, 123:132)] %>%
  left_join(Info0[, c(1, 62:105)], by = "CHILDID") %>%
  left_join(Skill0[, c(1, 40, 50, 60, 70, 80, 90, 98)], by = "CHILDID") %>% arrange(CHILDID)

#### 1: ID; 2-10: Reading IRT score; 11-19: Math IRT score; 20-27: Science IRT score;
#### 28:36: age (in month); 37: sex; 38: race;
#### 39: baseline locale type; 48/53: baseline parent 1/2 education;
#### 58: baseline income (measured in fall); 64: baseline school type
#### 92: Teacher report approaches to learning
#### 93: Teacher report self-control
#### 94: Teacher report interpersonal
#### 95: Teacher report external prob behaviors
#### 96: Teacher report internal prob behaviors
#### 97: Teacher report attention focus
#### 98: Teacher report inhibitory control
dat <- dat0[, c(1, 2:27, 28:36, 37:39, 48, 53, 58, 64, 92:98)]
dat[dat == "NaN" | dat == -9 | dat == -8 | dat == -7 | dat == -1] <- NA
dim(dat)
summary(dat)

three_TIC0 <- na.omit(dat)
dim(three_TIC0)

set.seed(20191029)
#set.seed(20221011)
selected <- sample(unique(three_TIC0$CHILDID), 500, replace = F)
three_TIC <- three_TIC0[three_TIC0$CHILDID %in% selected, ]
names(three_TIC) <- c("ID", paste0("R", 1:9), paste0("M", 1:9), paste0("S", 2:9), paste0("T", 1:9),
                      "SEX", "RACE", "LOCALE", "P1EDU", "P2EDU", "INCOME", "SCHOOL_TYPE",
                      "Approach_to_Learning", "Self_control", "Interpersonal", "External_prob_Behavior",
                      "Internal_prob_Behavior", "Attention_focus", "Inhibitory_Ctrl")
three_TIC$EDU <- apply(three_TIC[, 40:41], 1, max)
three_TIC <- three_TIC[, -c(40:41)]

RMS_dat <- three_TIC
save(RMS_dat, file = "RMS_dat.RData")

### Plot
set.seed(19850326)
selected <- sample(unique(three_TIC$ID), 30, replace = F)

Month_long <- gather(data = three_TIC[, c(1, 28:36)],
                     key = Wave, value = Month, T1:T9)
Month_long$Wave <- substr(Month_long$Wave, 2, 2)
RSCALK5_long <- gather(data = three_TIC[, c(1:10)],
                       key = Wave, value = RSCALK5, R1:R9)
RSCALK5_long$Wave <- substr(RSCALK5_long$Wave, 2, 2)
RSCALK5_dat <- Month_long %>%
  left_join(RSCALK5_long, by = c("ID", "Wave")) %>% arrange(ID)
head(RSCALK5_dat)
RSCALK5_dat0 <- RSCALK5_dat[RSCALK5_dat$ID %in% selected, ]
RSCALK5_dat0$Year <- RSCALK5_dat0$Month/12

reading_IRT <- ggplot(data = RSCALK5_dat0, aes(x = Year, y = RSCALK5, group = ID)) +
  geom_point(aes(group = ID), color = "grey", size = 1.2) +
  geom_line(aes(group = ID), color = "grey") +
  #geom_smooth(aes(group = 1), size = 3, col = "light blue", se = F) +
  labs(x = "Age (in Year)", y = "Reading Scores from Grade K to Grade 5") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 1),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.placement = "outside",
        legend.position = "bottom") +
  scale_x_continuous(breaks = c(4, 5, 6, 7, 8, 9, 10, 11, 12))
ggsave("reading_IRT.pdf", reading_IRT, device = "pdf", height = 4, width = 6)

Month_long <- gather(data = three_TIC[, c(1, 28:36)],
                     key = Wave, value = Month, T1:T9)
Month_long$Wave <- substr(Month_long$Wave, 2, 2)
MSCALK5_long <- gather(data = three_TIC[, c(1, 11:19)],
                       key = Wave, value = MSCALK5, M1:M9)
MSCALK5_long$Wave <- substr(MSCALK5_long$Wave, 2, 2)
MSCALK5_dat <- Month_long %>%
  left_join(MSCALK5_long, by = c("ID", "Wave")) %>% arrange(ID)
head(MSCALK5_dat)
MSCALK5_dat0 <- MSCALK5_dat[MSCALK5_dat$ID %in% selected, ]
MSCALK5_dat0$Year <- MSCALK5_dat0$Month/12

math_IRT <- ggplot(data = MSCALK5_dat0, aes(x = Year, y = MSCALK5, group = ID)) +
  geom_point(aes(group = ID), color = "grey", size = 1.2) +
  geom_line(aes(group = ID), color = "grey") +
  #geom_smooth(aes(group = 1), size = 3, col = "light blue", se = F) +
  labs(x = "Age (in Year)", y = "Mathematics Scores from Grade K to Grade 5") +
  theme_bw() +
  theme(strip.text.y = element_text(size = 1),
        strip.background = element_rect(color = "white", fill = "white"),
        strip.placement = "outside",
        legend.position = "bottom") +
  scale_x_continuous(breaks = c(4, 5, 6, 7, 8, 9, 10, 11, 12))
ggsave("math_IRT.pdf", math_IRT, device = "pdf", height = 4, width = 6)
