library(ggplot2)
library(reshape2)
library(dplyr)
library(simulateGP)
library(dplyr)

###########

###########

make_families <- function(af, nfam) {
	nsnp <- length(af)
	dads <- matrix(0, nfam, nsnp)
	mums <- matrix(0, nfam, nsnp)
	sibs1 <- matrix(0, nfam, nsnp)
	sibs2 <- matrix(0, nfam, nsnp)
	ibd <- matrix(0, nfam, nsnp)
	ibs <- matrix(0, nfam, nsnp)
	for(i in 1:nsnp)
	{
		dad1 <- rbinom(nfam, 1, af[i]) + 1
		dad2 <- (rbinom(nfam, 1, af[i]) + 1) * -1
		mum1 <- rbinom(nfam, 1, af[i]) + 1
		mum2 <- (rbinom(nfam, 1, af[i]) + 1) * -1

		dadindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		dadh <- rep(NA, nfam)
		dadh[dadindex] <- dad1[dadindex]
		dadh[!dadindex] <- dad2[!dadindex]

		mumindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		mumh <- rep(NA, nfam)
		mumh[mumindex] <- mum1[mumindex]
		mumh[!mumindex] <- mum2[!mumindex]

		sib1 <- cbind(dadh, mumh)

		dadindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		dadh <- rep(NA, nfam)
		dadh[dadindex] <- dad1[dadindex]
		dadh[!dadindex] <- dad2[!dadindex]

		mumindex <- sample(c(TRUE, FALSE), nfam, replace=TRUE)
		mumh <- rep(NA, nfam)
		mumh[mumindex] <- mum1[mumindex]
		mumh[!mumindex] <- mum2[!mumindex]

		sib2 <- cbind(dadh, mumh)

		ibd[,i] <- (as.numeric(sib1[,1] == sib2[,1]) + as.numeric(sib1[,2] == sib2[,2])) / 2


		sibs1[,i] <- rowSums(abs(sib1) - 1)
		sibs2[,i] <- rowSums(abs(sib2) - 1)
		dads[,i] <- dad1 - 1 + abs(dad2) - 1
		mums[,i] <- mum1 - 1 + abs(mum2) - 1

		# l[[i]] <- (sum(sib1[,1] == sib2[,1]) / nsnp + sum(sib1[,2] == sib2[,2]) / nsnp) / 2

	}

	# This may not be correct - getting some really large values
	ibs <- scale(sibs1) * scale(sibs2)

	# Just count how many alleles are in common
	ibs_unw <- abs(abs(sibs1 - sibs2) - 2) / 2

	return(list(dads=dads, mums=mums, sibs1=sibs1, sibs2=sibs2, ibd=ibd, ibs=ibs, ibs_unw=ibs_unw))
}



makePhen <- function(effs, indep, vy=1, vx=rep(1, length(effs)), my=0) {
	if(is.null(dim(indep))) indep <- cbind(indep)
	stopifnot(ncol(indep) == length(effs))
	stopifnot(length(vx) == length(effs))
	cors <- effs * vx / sqrt(vx) / sqrt(vy)
	stopifnot(sum(cors^2) <= 1)
	cors <- c(cors, sqrt(1-sum(cors^2)))
	indep <- t(t(scale(cbind(indep, rnorm(nrow(indep))))) * cors * c(vx, 1))
	y <- drop(scale(rowSums(indep)) * sqrt(vy)) + my
	return(y)
}

chooseEffects <- function(nsnp, totvar, sqrt=TRUE) {
	eff <- rnorm(nsnp)
	aeff <- abs(eff)
	sc <- sum(aeff) / totvar
	out <- eff / sc
	if(sqrt)
	{
		out <- sqrt(abs(out)) * sign(out)
	}
	return(out)
}

make_phenotypes <- function(fam, eff_gx, eff_xy, vx, vy, mx, my) {
	a <- lapply(fam, function(g)
	{
		u <- rnorm(nrow(g))
		x <- makePhen(c(eff_gx), cbind(g), vy=vx, my=mx)
		y <- makePhen(c(eff_xy), cbind(x), vy=vy, my=my)
		cc <- rbinom(length(y), 1, plogis(y))		
		return(data.frame(x=x, y=y, cc=cc))
	})
	names(a) <- names(fam)

	# generate sibling sex
	a$sibs1$sex <- rbinom(nrow(a$sibs1), 1, 0.5)
	a$sibs2$sex <- rbinom(nrow(a$sibs2), 1, 0.5)

	# only females can be cases
	a$sibs1$cco <- a$sibs1$cc
	a$sibs1$cc[a$sibs1$sex==0] <- 0
	a$sibs2$cc[a$sibs2$sex==0] <- 0

	return(a)
}

join_populations <- function(l) {
	dads <- do.call(rbind, lapply(l, function(x) x$dads))
	mums <- do.call(rbind, lapply(l, function(x) x$mums))
	sibs1 <- do.call(rbind, lapply(l, function(x) x$sibs1))
	sibs2 <- do.call(rbind, lapply(l, function(x) x$sibs2))
	ibd <- do.call(rbind, lapply(l, function(x) x$ibd))
	ibs <- do.call(rbind, lapply(l, function(x) x$ibs))
	return(list(dads=dads, mums=mums, sibs1=sibs1, sibs2=sibs2, ibd=ibd, ibs=ibs))
}

sample_populations <- function(l, n) {
	x <- nrow(l$dads)
	index <- sort(sample(1:x, n, replace=FALSE))
	l$dads <- l$dads[index,]
	l$mums <- l$mums[index,]
	l$sibs1 <- l$sibs1[index,]
	l$sibs2 <- l$sibs2[index,]
	l$ibd <- l$ibd[index,]
	l$ibs <- l$ibs[index,]
	l$ibs_unw <- l$ibs_unw[index,]
	return(l)
}


fastAssoc <- function(y, x) {
	index <- is.finite(y) & is.finite(x)
	n <- sum(index)
	y <- y[index]
	x <- x[index]
	vx <- var(x)
	bhat <- cov(y, x) / vx
	ahat <- mean(y) - bhat * mean(x)
	# fitted <- ahat + x * bhat
	# residuals <- y - fitted
	# SSR <- sum((residuals - mean(residuals))^2)
	# SSF <- sum((fitted - mean(fitted))^2)

	rsq <- (bhat * vx)^2 / (vx * var(y))
	fval <- rsq * (n-2) / (1-rsq)
	tval <- sqrt(fval)
	se <- abs(bhat / tval)

	# Fval <- (SSF) / (SSR/(n-2))
	# pval <- pf(Fval, 1, n-2, lowe=F)
	p <- pf(fval, 1, n-2, lowe=F)
	return(list(
		ahat=ahat, bhat=bhat, se=se, fval=fval, pval=p
	))
}

fastGLM <- function(cc, x) {
	require(fastglm)
	index <- is.finite(cc) & is.finite(x)
	n <- sum(index)
	cc <- cc[index]
	x <- cbind(rep(1, n), x[index])
	fit <- summary(fastglm(x, cc, family=binomial()))
	return(list(
		ahat=fit$coefficients[1,1],
		bhat=fit$coefficients[2,1],
		se=fit$coefficients[2,2],
		zval=fit$coefficients[2,3],
		pval=fit$coefficients[2,4]
	))
}

gwas <- function(y, g) {
	out <- matrix(0, ncol(g), 5)
	for(i in 1:ncol(g))
	{
		o <- fastAssoc(y, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

gwasGLM <- function(cc, g) {
	out <- matrix(0, ncol(g), 5)
	for(i in 1:ncol(g))
	{
		o <- fastGLM(cc, g[,i])
		out[i, ] <- unlist(o)
	}
	out <- as.data.frame(out)
	names(out) <- names(o)
	return(out)
}

do_mr_standard <- function(x, y, g) {
	gwasx <- gwas(x, g)
	gwasy <- gwas(y, g)
    ind <- gwasx$pval < 5e-8
	out <- mr_ivw(gwasx$bhat[ind], gwasy$bhat[ind], gwasx$se[ind], gwasy$se[ind])
	return(out)
}

do_mr_glm <- function(x, cc, g) {
	gwasx <- gwas(x, g)
	gwascc <- gwasGLM(cc, g)
    ind <- gwasx$pval < 5e-8
	out <- mr_ivw(gwasx$bhat[ind], gwascc$bhat[ind], gwasx$se[ind], gwascc$se[ind])
	return(out)
}

mr_ivw <- function(b_exp, b_out, se_exp, se_out, parameters = default_parameters()) {
  if (sum(!is.na(b_exp) & !is.na(b_out) & !is.na(se_exp) & !is.na(se_out)) < 2) {
    return(list(b = NA, se = NA, pval = NA, nsnp = NA))
  }

  ivw.res <- summary(stats::lm(b_out ~ -1 + b_exp, weights = 1 / se_out^2))
  b <- ivw.res$coef["b_exp", "Estimate"]
  se <- ivw.res$coef["b_exp", "Std. Error"] / min(1, ivw.res$sigma) #sigma is the residual standard error
  pval <- 2 * stats::pnorm(abs(b / se), lower.tail = FALSE)
  Q_df <- length(b_exp) - 1
  Q <- ivw.res$sigma^2 * Q_df
  Q_pval <- stats::pchisq(Q, Q_df, lower.tail = FALSE)
  # from formula phi =  Q/DF rearranged to to Q = phi*DF, where phi is sigma^2
  # Q.ivw<-sum((1/(se_out/b_exp)^2)*(b_out/b_exp-ivw.reg.beta)^2)
  return(list(
    b = b,
    se = se,
    pval = pval,
    nsnp = length(b_exp),
    Q = Q,
    Q_df = Q_df,
    Q_pval = Q_pval
  ))
}


library(parallel)
nsim <- 100
res <- mclapply(1:nsim, function(i) {
    fam <- make_families(af=seq(0.1, 0.5, length.out=50), nfam=10000)
    eff <- chooseEffects(50, totvar=0.4, sqrt=TRUE)
    phen <- make_phenotypes(fam, eff, -0.42, 1, 1, 0, 0)
	ind <- which(phen$sibs1$sex == 1)
    bind_rows(
        do_mr_standard(phen$sibs1$x, phen$sibs1$y, fam$sibs1) %>% as_tibble() %>% mutate(what="self"),
        do_mr_glm(phen$sibs1$x, phen$sibs2$cc, fam$sibs1) %>% as_tibble() %>% mutate(what="sibling"),
        do_mr_glm(phen$sibs1$x, phen$mums$cc, fam$sibs1) %>% as_tibble() %>% mutate(what="mother")
    )
}, mc.cores=10) %>% bind_rows()

res <- res %>% mutate(source="simulated")

saveRDS(res, "mr_proxy_simulation_results.rds")

obs_res <- tibble(
    what = c("self", "sibling", "mother"),
    b = log(c(0.66, 0.89, 0.84)),
    lower = log(c(0.52, 0.69, 0.67)),
    upper = log(c(0.84, 1.14, 1.06)),
    se = (c((log(0.84)-log(0.66))/1.96, (log(1.14)-log(0.89))/1.96, (log(1.06)-log(0.84))/1.96))
) %>% mutate(source="empirical")
obs_res
ggplot(bind_rows(res, obs_res), aes(x=what, y=exp(b), color=source)) +
    geom_point(position=position_dodge(width=0.5), size=3) +
    geom_errorbar(aes(ymin=exp(lower), ymax=exp(upper)), width=0, position=position_dodge(width=0.5)) +
    labs(x="Individual used for outcome GWAS", y="MR Estimate (odds ratio)", colour="") +
    geom_hline(yintercept=1, linetype="dashed")
ggsave("mr_proxy_analysis.pdf", width=6, height=6)



