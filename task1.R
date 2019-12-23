rm(list=ls())
library(igraph)
require(ggplot2)
require(reshape2)

till_end <- FALSE
n <- 1000
max_t <- 100
beta <- 0.2
gamma <- 0.1
p0 <- 0.01

graphs <- list(
  full_graph = make_full_graph(n),
  er_graph = erdos.renyi.game(n, 0.01),
  ws_graph = watts.strogatz.game(1, n, 4, 0.05),
  ba_graph = barabasi.game(n, directed=FALSE)
)

get_stat <- function(graph) {
  i <- sum(V(graph)$labels=='i')
  r <- sum(V(graph)$labels=='r')
  s <- sum(V(graph)$labels=='s')
  print(i)
}

experiment <- function(graph, beta, gamma) {

  V(graph)$labels <- rep("s", n) # All are susceptible
  
  # Select initially infected ones
  inf.init.number <- round(p0*n)
  inf.init.nodes <- sample.int(n, size = inf.init.number)
  V(graph)$labels[inf.init.nodes] <- rep("i", inf.init.number)
  infected.trace <- c(inf.init.number)
  t <- 0
  while ((sum(V(graph)$labels=="r") | !till_end) != n & t < max_t) {
    get_stat(graph)
    t <- t + 1
    # Recover
    infected <- V(graph)[V(graph)$labels=='i']
    to_recover.sample <- runif(length(infected))
    to_recover <- infected[to_recover.sample < gamma]
    
    # Infect
    to_infect.all <- c()
    for (v in infected) {
      adjacent <- V(graph)[neighbors(graph, v)]
      if(length(adjacent) == 0 | sum(adjacent$labels=="s")==0) {
        next
      }
      
      adj.infectable <- adjacent[adjacent$labels=="s"]
      if(length(adj.infectable) == 0) {
        next
      }
      
      to_infect.sample <- runif(length(adj.infectable))
      if(sum(to_infect.sample < beta) == 0) {
        next
      }
      
      to_infect <- adj.infectable[to_infect.sample < beta]
      to_infect <- adj.infectable[to_infect.sample < beta]
      to_infect.all <- c(to_infect.all, to_infect)
    }
    
    # Perform timestep
    V(graph)[to_recover]$labels <- rep("r", length(to_recover))
    V(graph)[to_infect.all]$labels <- rep("i", length(to_infect.all))
    
    infected.trace <- c(infected.trace, sum(V(graph)$labels=="i"))
  }
  return(infected.trace)
}

all.results <- list()
i <- 1
for (graph in graphs) {
  print(graph$name)
  result <- experiment(graph, beta, gamma)
  all.results[[i]] <- result
  i <- i + 1
}

names(all.results) <- sapply(graphs, function(graph) graph$name[1])
all.results <- as.data.frame(all.results)
all.results$time <- seq.int(nrow(all.results))

df <- melt(all.results,  id.vars = 'time', variable.name = "graph")

plot <- ggplot(df, aes(time, value)) + 
  geom_line(aes(colour = graph)) + 
  ggtitle("Number of infected through time") +
  ylab("Infected") + 
  xlab("Time")

show(plot)

eigen_vals <- sapply(graphs, function(graph) spectrum(graph, which=list(pos="LA", howmany=1))$values)

#
# Task 2
#

format_data <- function(data, title, color, frac) {
    
    to_df <- as.data.frame(data)
    to_df$time <- seq.int(nrow(to_df))
    final_df <- melt(to_df,  id.vars = 'time', variable.name = "simulation")
    
    final.plot <- ggplot(final_df, aes(time, value)) + 
        geom_line(aes(colour = simulation)) + 
        ggtitle(title) +
        scale_color_manual(labels = c(title), values = c(color)) +
        ylab("Infected") + 
        xlab("Time") +
        theme(legend.position = "none") +
        annotate("text", x=max(final_df$time)-15, 
                 y=max(final_df$value)-15, label= paste("β/γ = ", frac), size=5)
    
    show(final.plot)
}

run_all <- function(graph_n, b, g, title, color) {
    
    beta <- b
    gamma <- g
    cat(graph_n$name, title,"β = ", beta, ", γ = ", gamma,"\n")
    result.epidemic <- experiment(graph_n, beta, gamma)
    format_data(result.epidemic, title, color, b/g)
}


for (graph in graphs) {
    if (graph$name == "Full graph"){
        run_all(graph, b=0.002, g=0.1,
                "Full Graph - Epidemic occurs", "red")
        run_all(graph, b=0.00005, g=0.1,
                "Full Graph - No Epidemic occurs", "blue")
    }
    if (graph$name == "Erdos renyi (gnp) graph"){
        run_all(graph, b=0.03, g=0.1,
                "Erdos-Renyi - Epidemic occurs", "red")
        run_all(graph, b=0.007, g=0.1,
                "Erdos-Renyi - No Epidemic occurs", "blue")
    }
    if (graph$name == "Watts-Strogatz random graph"){
        run_all(graph, b=0.07, g=0.1,
                "Watts-Strogatz - Epidemic occurs", "red")
        run_all(graph, b=0.005, g=0.1,
                "Watts-Strogatz - No Epidemic occurs", "blue")
    }
    if (graph$name == "Barabasi graph"){
        run_all(graph, b=0.5, g=0.1,
                "Barabasi - Epidemic occurs", "red")
        run_all(graph, b=0.005, g=0.1,
                "Barabasi - No Epidemic occurs", "blue")
    }
}
