coal2newick <- function(samp_t, n_samp, coal_t, n_merge, node_name_prefix="N")
{
    n_leaves <- sum(n_samp)
    samp_t_expanded <- expand_samp_times(samp_t, n_samp)

    tree_nodes <- seq(1, n_leaves)
    tree_nodes <- sapply(tree_nodes, function (x) return (paste0("S", x)))

    dates <- samp_t_expanded
    names(dates) <- tree_nodes

    extant_entries <- c(1)
    extant_times <- c(samp_t_expanded[1])
    
    coal_idx <- 1
    s_idx <- 1
    t <- samp_t_expanded[1]
    if (length(coal_t) > 0) 
    {
        while (coal_idx <= length(coal_t)) 
        {
            if (s_idx < length(samp_t_expanded) && (samp_t_expanded[s_idx + 1] < coal_t[coal_idx])) {
                s_idx <- s_idx + 1
                t <- samp_t_expanded[s_idx]
                extant_entries <- c(extant_entries, s_idx)
                extant_times <- c(extant_times, t)
            } else {
                t <- coal_t[coal_idx]
                k <- n_merge[coal_idx]
                
                to_merge_indices <- sample(1:length(extant_entries), k)
                to_merge_nodes <- extant_entries[to_merge_indices]

                br_lens <- sapply(to_merge_indices, function(i) t - extant_times[i])
                node_entries <- tree_nodes[to_merge_nodes]

                to_del <- to_merge_indices[1:(k-1)]
                ret <- to_merge_indices[k]
                
                node_content <- paste0(node_entries, ":", br_lens, collapse=",")

                tree_nodes[to_merge_nodes[k]] <- paste0("(", node_content, ")",node_name_prefix,coal_idx)
                extant_times[ret] <- t
                coal_idx <- coal_idx + 1

                extant_times <- extant_times[-to_del]
                extant_entries <- extant_entries[-to_del]
            }
        }
    }
    stopifnot(length(extant_entries)==1)
    tree_str <- paste0(tree_nodes[extant_entries[1]], ":0;")
    return(list(str=tree_str, dates=dates))
}

expand_samp_times <- function(samp_t, n_samp)
{
    return(unlist(lapply(1:length(samp_t), function(i) rep(samp_t[i], n_samp[i]))))
}