# Independent base-R implementation of invented paired mechanics, version 0.1.0.
# Only reads the newly generated synthetic resampling schedule; never MIMIC.
options(digits=17)
validate_rows <- function(rows) {
  stopifnot(is.data.frame(rows), nrow(rows)>0,
            all(rows$scope=="SYNTHETIC"), all(rows$invented),
            !anyDuplicated(rows$item_id), all(nzchar(rows$item_id)), all(nzchar(rows$cluster_id)),
            all(is.na(rows$reference) | rows$reference %in% c(0,1)),
            !anyNA(rows$candidate), !anyNA(rows$comparator),
            all(rows$candidate %in% c(0,1)), all(rows$comparator %in% c(0,1)))
}
stratum <- function(rows, reference=1) {
  validate_rows(rows); stopifnot(length(reference)==1, !is.na(reference), reference %in% c(0,1))
  selected <- rows[!is.na(rows$reference) & rows$reference==reference,,drop=FALSE]
  stopifnot(nrow(selected)>0); selected
}
paired_summary <- function(rows, reference=1) {
  selected <- stratum(rows,reference)
  gain <- sum(selected$candidate==1 & selected$comparator==0)
  loss <- sum(selected$candidate==0 & selected$comparator==1)
  data.frame(reference_stratum=reference,n_total=nrow(rows),
             n_unknown_reference=sum(is.na(rows$reference)),n_stratum=nrow(selected),
             n_clusters=length(unique(selected$cluster_id)),
             candidate_rate=mean(selected$candidate),comparator_rate=mean(selected$comparator),
             paired_difference=(gain-loss)/nrow(selected),gains=gain,losses=loss)
}
cluster_bootstrap <- function(rows, plan, reference=1) {
  selected <- stratum(rows,reference)
  groups <- split(selected,selected$cluster_id)
  stopifnot(length(groups)>=2,length(plan)>=2)
  values <- vapply(plan,function(draws) {
    stopifnot(length(draws)==length(groups),all(draws %in% names(groups)))
    # Preserve repeated sampled clusters and every paired item within them.
    sampled <- do.call(rbind,lapply(draws,function(key) groups[[key]]))
    mean(sampled$candidate-sampled$comparator)
  },numeric(1))
  list(values=values,lower=unname(quantile(values,.025,type=7)),
       upper=unname(quantile(values,.975,type=7)))
}
exact_discordance <- function(rows,reference=1) {
  selected <- stratum(rows,reference)
  stopifnot(nrow(selected)==length(unique(selected$cluster_id)))
  value <- paired_summary(rows,reference)
  total <- value$gains+value$losses
  if(total==0) return(1)
  binom.test(value$gains,total,p=.5,alternative="two.sided")$p.value
}
expect_error <- function(expr) {
  stopifnot(inherits(tryCatch({force(expr);NULL},error=function(e)e),"error"))
}

args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==1,dir.exists(args[[1]]))
out <- args[[1]]
outputs <- file.path(out,c("r_summary.csv","r_bootstrap_values.csv","r_mechanics.json","r_session.txt"))
stopifnot(!any(file.exists(outputs)))
# Constructed independently of Python code and unrelated to clinical forms.
rows <- data.frame(scope=rep("SYNTHETIC",8),invented=rep(TRUE,8),
                   item_id=paste0("i",1:8),cluster_id=c("c1","c1","c2","c3","c3","c4","c4","c5"),
                   reference=c(1,1,1,1,0,0,0,NA),candidate=c(1,1,0,1,1,0,0,1),
                   comparator=c(0,1,1,0,0,1,0,1),stringsAsFactors=FALSE)
pos <- paired_summary(rows,1); neg <- paired_summary(rows,0)
stopifnot(pos$n_stratum==4,pos$n_clusters==3,pos$paired_difference==.25,
          pos$candidate_rate==.75,pos$comparator_rate==.5,pos$gains==2,pos$losses==1,
          pos$n_unknown_reference==1,neg$n_stratum==3,neg$paired_difference==0)
hand <- cluster_bootstrap(rows,list(c("c1","c1","c2"),rep("c3",3)))
stopifnot(isTRUE(all.equal(hand$values,c(.2,1))),abs(hand$lower-.22)<1e-12,abs(hand$upper-.98)<1e-12)
expect_error(exact_discordance(rows))
independent <- rows[1:4,]; independent$cluster_id <- paste0("ind",1:4)
independent$candidate <- 1; independent$comparator <- 0
stopifnot(exact_discordance(independent)==.125)
bad <- rows; bad$candidate[1] <- NA; expect_error(paired_summary(bad))
bad <- rows; bad$reference <- NA; expect_error(paired_summary(bad))
expect_error(cluster_bootstrap(rows,list(rep("unknown",3),rep("unknown",3))))
schedule <- read.csv(file.path(out,"synthetic_resampling_plan.csv"),stringsAsFactors=FALSE)
stopifnot(all(schedule$scope=="SYNTHETIC"),nrow(schedule)==6000,
          identical(sort(unique(schedule$replicate)),1:2000))
plan <- lapply(1:2000,function(i) {
  group <- schedule[schedule$replicate==i,,drop=FALSE]
  stopifnot(identical(group$draw,1:3)); group$cluster_id
})
boot <- cluster_bootstrap(rows,plan)
write.csv(rbind(pos,neg),outputs[[1]],row.names=FALSE)
write.csv(data.frame(replicate=seq_along(boot$values),value=boot$values),outputs[[2]],row.names=FALSE)
writeLines(sprintf('{"scope":"SYNTHETIC_MECHANICS_ONLY","assertions":"PASS","bootstrap_lower":%.17g,"bootstrap_upper":%.17g,"exact_independent_p":%.17g}',boot$lower,boot$upper,exact_discordance(independent)),outputs[[3]])
capture.output(sessionInfo(),file=outputs[[4]])
cat("R_INVENTED_STATISTICAL_MECHANICS=PASS\nNO_CLINICAL_INFERENCE\n")
