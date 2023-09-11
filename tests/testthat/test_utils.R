context("Consenus")
test_that("Clade node counting functions agree with ground truth",
{
    swap_n <- function(x, i, j)
    {
        sapply(x, function(y) if(y==i) j else if(y==j) i else y)
    }

    tr1 <- read.tree(text=c("(((A:1,B:1,C:1,D:1):1, X:1), ((E:1, F:1):1, G:1, H:1):1, I:1):1;"))
    tr2 <- read.tree(text=c("((((A:1,B:1,C:1):1,D:1):1, X:1), ((E:1, F:1):1, (G:1, H:1):1):1, I:1):1;"))
    tr3 <- read.tree(text=c("(((((A:1,B:1):1,C:1):1,D:1):1, X:1), (((E:1, F:1):1, (G:1, H:1):1):1, I:1):1):1;"))

    gt1 <- c(0L, 1L, 3L)
    gt2 <- c(1L, 2L, 3L)
    gt3 <- c(2L, 2L, 4L)

    to_test <- list(c("C","B","A", "D"), c("F","E","G", "H"), c("A","B","D","C","G","E","F","H","I","X"))
    clade_mat <- as.matrix(do.call(rbind, lapply(to_test, function(x) clade_as_vec(x, tr1$tip.label))))

    res1 <- nodes_in_clades(tr1, clade_mat)
    res2 <- nodes_in_clades(tr2, clade_mat)
    res3 <- nodes_in_clades(tr3, clade_mat)

    expect_equal(res1, gt1)
    expect_equal(res2, gt2)
    expect_equal(res3, gt3)
})