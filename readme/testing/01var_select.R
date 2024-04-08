devtools::document()
devtools::install()
# devtools::load_all()
library(data.table)

extract_dt <- readRDS("../inst/extdata/2019extracted_dt.rds")
# write.csv(extract_dt, "../inst/extdata/2019extracted_dt.csv")
names(extract_dt)
hist(extract_dt[, h_canopy_ge0])
extract_dt = extract_dt[n_canopy_total > 5 & n_ground > 5]

y = "h_canopy_ge0"
xcols = names(extract_dt)[-(1:10)]
prepend_class(extract_dt, "icesat2.atl03_atl08_seg_dt")
set.seed(178724)
res <- var_select_local(extract_dt, y, xcols, nboots = 10, nTrees = 100, train_split = 0.7, delta = 0.0000000001, spacing = 30)
cols = strsplit(rev(res$property)[2], ", ")[[1]]
# cols = c("red_diss_mean", "blue_savg_min", "num_stdDev", "swir1_diss")

# spacing = 30
# train_split = 0.7
# sample_i <- sample(extract_dt, spacedSampling(0.9999999, spacing / degree_to_meter_factor))
# sample_ii <- sample(sample_i, stratifiedSampling(20, y, breaks = c(0, 5, 10, 15, 20, 999)))[, I := .I]
# train_i <- sample(sample_ii, randomSampling(train_split))

# x_i <- sample_ii[train_i$I, .SD, .SDcols = "blue"]
# y_i <- sample_ii[train_i$I][[y]]
# validation_sample <- sample_ii[-train_i$I]

library(randomForest)
set.seed(178724)
colsToUse <- c(y, cols)


library(VSURF)

x = extract_dt[, .SD, .SDcols = -(1:10)]
y = extract_dt[[y]]
var_sel <- VSURF::VSURF(x, y, parallel = TRUE)

degree_to_meter_factor <- 111139
summary(extract_dt[, h_canopy_ge0])

sample2 <- sample(extract_dt, stratifiedSampling(124, variable = y, breaks = c(0, 5, 10, 15, 20, 999)))
extract_dt[, I := .I]
train <- sample(extract_dt, randomSampling(0.7))
test <- extract_dt[-train$I]



rf <- randomForest(data = train[, .SD, .SDcols = colsToUse], h_canopy_ge0 ~ ., mtry = 1, ntree = 1000)
forest_string <- build_forest(rf)
gee_forest <- ee$Classifier$decisionTreeEnsemble(forest_string)

predicted <- predict(rf, test)

png("here.png")
lims <- c(0, max(extract_dt[[y]]))
stats_model(test[[y]], predicted, xlim = lims, ylim = lims, unit = "m")
dev.off()


