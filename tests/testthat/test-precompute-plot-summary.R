make_fixture2 <- function(n = 60, p1 = 25, p2 = 20, K = 2L, seed = 11L) {
    set.seed(seed)
    groups <- rep(1:3, each = n / 3)
    X1 <- matrix(rnorm(n * p1), n, p1)
    X2 <- matrix(rnorm(n * p2), n, p2)
    miss1 <- seq(2, n, by = 6)
    miss2 <- seq(4, n, by = 6)
    X1_m <- X1
    X1_m[miss1, ] <- NA
    X2_m <- X2
    X2_m[miss2, ] <- NA
    anchor_idx <- setdiff(seq_len(n), union(miss1, miss2))
    list(
        X1 = X1_m, X2 = X2_m, groups = groups,
        anchor_idx = anchor_idx, K = K, n = n
    )
}

test_that("woven_precompute() Laplacians are reused identically by woven()", {
    d <- make_fixture2()
    X_list <- list(d$X1, d$X2)

    precomp <- woven_precompute(X_list, k_nn = 10L)
    expect_length(precomp, 2L)
    expect_true(all(vapply(precomp, function(L) is.matrix(L) || inherits(L, "Matrix"), logical(1))))
    expect_true(all(vapply(precomp, nrow, integer(1)) == d$n))

    fit_precomp <- woven(X_list,
        Y = d$groups, anchor_idx = d$anchor_idx,
        K = d$K, precomp = precomp, verbose = FALSE
    )
    fit_direct <- woven(X_list,
        Y = d$groups, anchor_idx = d$anchor_idx,
        K = d$K, verbose = FALSE
    )
    expect_equal(Z(fit_precomp), Z(fit_direct), tolerance = 1e-8)
    expect_equal(WList(fit_precomp)[[1]], WList(fit_direct)[[1]], tolerance = 1e-8)
})

test_that("a precomp with the wrong number of blocks errors clearly", {
    d <- make_fixture2()
    precomp_bad <- woven_precompute(list(d$X1, d$X2), k_nn = 10L)[1] # drop 1 block
    expect_error(
        woven(list(d$X1, d$X2),
            Y = d$groups, anchor_idx = d$anchor_idx,
            K = d$K, precomp = precomp_bad, verbose = FALSE
        ),
        regexp = "woven_precompute"
    )
})

test_that("woven_precompute() reused across two different anchor sets/fits", {
    d <- make_fixture2()
    X_list <- list(d$X1, d$X2)
    precomp <- woven_precompute(X_list, k_nn = 10L)

    half <- d$anchor_idx[seq_len(floor(length(d$anchor_idx) / 2))]
    fit_full <- woven(X_list, Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, precomp = precomp, verbose = FALSE)
    fit_half <- woven(X_list, Y = d$groups, anchor_idx = half, K = d$K, precomp = precomp, verbose = FALSE)

    expect_equal(nrow(Z(fit_full)), d$n)
    expect_equal(nrow(Z(fit_half)), d$n)
    expect_false(isTRUE(all.equal(Z(fit_full), Z(fit_half))))
})

test_that("summary() with no labels prints without error and returns NULL invisibly", {
    d <- make_fixture2()
    fit <- woven(list(d$X1, d$X2), Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, verbose = FALSE)
    expect_output(res <- summary(fit), "Pass labels = Y")
    expect_null(res)
})

test_that("summary() with labels returns an invisible named metrics vector", {
    d <- make_fixture2()
    fit <- woven(list(d$X1, d$X2), Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, verbose = FALSE)
    expect_output(res <- summary(fit, labels = d$groups), "Silhouette")
    expect_named(res, c("Silhouette", "Davies-Bouldin", "NMI", "ESS"))
    expect_true(all(is.finite(res)))
})

test_that("plot() returns a visible ggplot object, with and without labels", {
    skip_if_not_installed("ggplot2")
    d <- make_fixture2()
    fit <- woven(list(d$X1, d$X2), Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, verbose = FALSE)

    p1 <- plot(fit, labels = d$groups)
    expect_s3_class(p1, "ggplot")

    p2 <- plot(fit)
    expect_s3_class(p2, "ggplot")

    p3 <- plot(fit, labels = d$groups, highlight_anchors = FALSE)
    expect_s3_class(p3, "ggplot")
})

test_that("plot() errors clearly if requested dims exceed K", {
    skip_if_not_installed("ggplot2")
    d <- make_fixture2()
    fit <- woven(list(d$X1, d$X2), Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, verbose = FALSE)
    expect_error(plot(fit, labels = d$groups, dims = c(1L, 5L)))
})

test_that("woven_nystrom_error() returns a finite, non-degenerate value (regression test for the fit$Za_list bug)", {
    d <- make_fixture2()
    fit <- woven(list(d$X1, d$X2), Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, verbose = FALSE)
    err <- woven_nystrom_error(fit, X_list = list(d$X1, d$X2))
    expect_true(is.numeric(err))
    expect_length(err, 1L)
    expect_true(is.finite(err))
    expect_false(is.nan(err))
    expect_gte(err, 0)
})

test_that("woven_nystrom_error() respects n_loo and errors on non-woven input", {
    d <- make_fixture2()
    fit <- woven(list(d$X1, d$X2), Y = d$groups, anchor_idx = d$anchor_idx, K = d$K, verbose = FALSE)
    err <- woven_nystrom_error(fit, X_list = list(d$X1, d$X2), n_loo = 5L)
    expect_true(is.finite(err))
    expect_error(woven_nystrom_error(list(a = 1), X_list = list(d$X1, d$X2)), regexp = "woven object")
})
