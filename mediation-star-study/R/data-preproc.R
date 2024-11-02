library(haven)
library(dplyr)
students <- haven::read_spss("./data/STAR_Students.sav")

dat <- students %>% 
  filter(FLAGSG1 == 1,
         FLAGSGK == 1,
         FLAGSG2 == 1,
         FLAGSG3 == 1,
         flagg4 == 1,
         flagg8 == 1) %>%
  select(
    # id/grouping variables
    stdntid, gkschid,
    # pretreatment variables
    gender, race, birthmonth, birthyear, gkfreelunch, gksurban,
    # randomized treatment
    gkclasstype, 
    # Grade K posttreatment variables
    gkabsent, gkpresent, 
    gktreadss, gktmathss, gktlistss, gkwordskillss,
    gkmotivraw, gkselfconcraw,
    # Grade 1 posttreatment variables
    g1absent, g1present,
    g1treadss, g1tmathss, g1tlistss, g1wordskillss,
    g1readbsraw, g1mathbsraw, g1readbsobjpct, g1mathbsobjpct,
    g1motivraw, g1selfconcraw,
    # Grade 2 posttreatment variables
    g2treadss, g2tmathss, g2tlistss, g2wordskillss, 
    g2readbsraw, g2mathbsraw, g2readbsobjpct, g2mathbsobjpct,
    g2motivraw, g2selfconcraw,
    # Grade 3 posttreatment variables
    g3present, g3absent,
    g3treadss, g3tmathss, g3tlangss, g3tlistss, g3socialsciss, g3spellss,
    g3vocabss, g3mathcomputss, g3mathnumconcss, g3mathapplss,
    g3wordskillss, g3readbsraw, g3mathbsraw, g3readbsobjpct, g3mathbsobjpct,
    g3motivraw, g3selfconcraw,
    # Grade 4 posttreatment variables
    # g4present, g4absent,
    g4treadss, g4tmathss, g4tlangss, g4socialsciss, g4sciencess, g4tbattss,
    g4readcomprehss, g4spellss, g4vocabss,
    g4langexpss, g4langmechss, g4studyskillss, g4mathbsobjraw,
    g4mathcomputss, g4mathconcapplss, # same test as g8mathco, g8math_a at 4th gr
    # Grade 8 outcomes
    g8mathco, g8math_a
    # Grade 8 possible predictors
    # g8surban
    ) %>% 
  mutate(
    yearsto85 = 1985 + (8/12) - (birthyear + (birthmonth/12)),
    yearsto85 = round(yearsto85, 1),
    gkabsent = gkabsent/(gkabsent + gkpresent),
         g1absent = g1absent/(g1absent + g1present),
         g3absent = g3absent/(g3absent + g3present),
         gkpresent = 1-gkabsent,
         g1present = 1-g1absent,
         g3present = 1-g3absent
  ) %>% 
  select(-birthmonth, -birthyear)

md <- dat %>% mice::md.pattern(plot=FALSE)
md.df <- data.frame(num = as.numeric(rownames(md)), md)
md.df <- md.df %>% arrange(desc(num))
# View(md.df)

dat <- dat %>% na.omit
dim(dat)

# Objects for Analysis

X_df <- select(dat, gender, race, yearsto85, gkfreelunch, gksurban) %>% 
  mutate(genderM = 1*(gender == 1),
         racewhite = 1*(race == 1),
         gkfreelunch = 1*(gkfreelunch == 1),
         gksurban = factor(gksurban, levels = attr(dat$gksurban, "labels"),
                           labels = names(attr(dat$gksurban, "labels")))
  ) %>% 
  select(genderM, racewhite, gkfreelunch, gksurban)
# levels(X_df$gksurban) <- attr(dat$gksurban, "labels")
Xmat <- model.matrix(~., data=X_df)[,-1]

d_vec <- pull(dat, gkclasstype)
d_bin <- 1*(d_vec == 1) # binary == 1 iff small class size
M_df <- dat %>% select( 
  # Grade K posttreatment variables
  gkabsent, # gkpresent, REMOVE since pct absent already captured
  gktreadss, gktmathss, gktlistss, gkwordskillss,
  gkmotivraw, gkselfconcraw,
  # Grade 1 posttreatment variables
  g1absent, # g1present, REMOVE since pct absent already captured
  g1treadss, g1tmathss, g1tlistss, g1wordskillss,
  g1readbsraw, g1mathbsraw, g1readbsobjpct, g1mathbsobjpct,
  g1motivraw, g1selfconcraw,
  # Grade 2 posttreatment variables
  g2treadss, g2tmathss, g2tlistss, g2wordskillss, 
  g2readbsraw, g2mathbsraw, g2readbsobjpct, g2mathbsobjpct,
  g2motivraw, g2selfconcraw,
  # Grade 3 posttreatment variables
  g3absent, # g3present, REMOVE since pct absent already captured
  g3treadss, g3tmathss, g3tlangss, g3tlistss, g3socialsciss, g3spellss,
  g3vocabss, g3mathcomputss, g3mathnumconcss, g3mathapplss,
  g3wordskillss, g3readbsraw, g3mathbsraw, g3readbsobjpct, g3mathbsobjpct,
  g3motivraw, g3selfconcraw,
  # Grade 4 posttreatment variables
  # g4present, g4absent,
  # g4treadss, SINGULAR 
  g4tmathss, g4tlangss, g4socialsciss, g4sciencess, g4tbattss,
  g4readcomprehss, g4spellss, g4vocabss,
  g4langexpss, g4langmechss, g4studyskillss, g4mathbsobjraw,
  g4mathcomputss, g4mathconcapplss
)
Y1 <- pull(dat, g8mathco)
Y2 <- pull(dat, g8math_a)
