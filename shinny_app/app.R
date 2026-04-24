#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
#

library(shiny)
library(ggplot2)



# Generate synthetic data if file doesn't exist for demonstration
if(!file.exists("aluminium.csv")) {
  set.seed(123)
  true_tau <- 0.005
  x_sim <- seq(0, 0.015, length.out = 100)
  y_sim <- ifelse(x_sim <= true_tau, 70000 * x_sim, 70000 * true_tau + 5000 * (x_sim - true_tau)) + rnorm(100, 0, 15)
  #write.csv(data.frame(strain = x_sim, stress = y_sim), "aluminium.csv", row.row = FALSE)
}

ui <- fluidPage(
  titlePanel("Bayesian Piecewise Linear Yield Point Analysis"),
  
  sidebarLayout(
    sidebarPanel(
      h4("Prior Settings"),
      sliderInput("prior_beta", "Prior Mean Slopes (beta1, beta2):", 
                  min = 10000, max = 100000, value = 30000, step = 5000),
      sliderInput("prior_b0", "Inverse-Gamma Scale (b0):", 
                  min = 1, max = 50, value = 13.12, step = 1),
      hr(),
      helpText("This app calculates the posterior distribution of the yield point (tau) and averages the regression lines over all possible breakpoints.")
    ),
    
    mainPanel(
      plotOutput("stressStrainPlot"),
      hr(),
      fluidRow(
        column(6, h4("Posterior of Yield Point"), plotOutput("tauDist")),
        column(6, h4("Numerical Summaries"), verbatimTextOutput("stats"))
      )
    )
  )
)

server <- function(input, output) {
  
  # Reactive calculation
  model_results <- reactive({
    df <- read.csv("aluminium.csv")
    x <- df$strain
    y <- df$stress
    n <- length(y)
    
    # Hyperparameters from inputs
    beta0 <- c(input$prior_beta, input$prior_beta)
    V0 <- diag(20000^2, 2)
    invV0 <- solve(V0)
    a0 <- 4
    b0 <- input$prior_b0
    
    taus <- sort(unique(x))
    taus <- taus[taus > quantile(x, 0.05) & taus < quantile(x, 0.95)]
    
    log_marg_lik <- numeric(length(taus))
    beta_means <- matrix(0, nrow = length(taus), ncol = 2)
    
    for (i in seq_along(taus)) {
      tau <- taus[i]
      Z <- cbind(pmin(x, tau), pmax(0, x - tau))
      Vn <- solve(invV0 + t(Z) %*% Z)
      betan <- Vn %*% (invV0 %*% beta0 + t(Z) %*% y)
      an <- a0 + n/2
      bn <- b0 + 0.5 * as.numeric(t(y)%*%y + t(beta0)%*%invV0%*%beta0 - t(betan)%*%solve(Vn)%*%betan)
      log_marg_lik[i] <- 0.5 * determinant(Vn)$modulus - an * log(bn)
      beta_means[i, ] <- betan
    }
    
    post_probs <- exp(log_marg_lik - max(log_marg_lik))
    post_probs <- post_probs / sum(post_probs)
    
    # Posterior Mean Curve
    x_grid <- seq(min(x), max(x), length.out = 200)
    y_post_mean <- numeric(length(x_grid))
    for (i in seq_along(taus)) {
      y_post_mean <- y_post_mean + post_probs[i] * (beta_means[i,1]*pmin(x_grid, taus[i]) + beta_means[i,2]*pmax(0, x_grid - taus[i]))
    }
    
    list(df=df, x_grid=x_grid, y_post_mean=y_post_mean, taus=taus, post_probs=post_probs)
  })
  
  output$stressStrainPlot <- renderPlot({
    res <- model_results()
    ggplot(res$df, aes(strain, stress)) + 
      geom_point(alpha = 0.4) +
      geom_line(data = data.frame(x=res$x_grid, y=res$y_post_mean), aes(x, y), color = "red", size = 1.2) +
      labs(title = "Fitted Posterior Mean Curve", subtitle = "Red line represents the Bayesian Model Average over all taus") +
      theme_minimal()
  })
  
  output$tauDist <- renderPlot({
    res <- model_results()
    ggplot(data.frame(tau = res$taus, prob = res$post_probs), aes(tau, prob)) +
      geom_col(fill = "steelblue") +
      labs(x = "Strain (tau)", y = "Posterior Probability") +
      theme_minimal()
  })
  
  output$stats <- renderPrint({
    res <- model_results()
    map_tau <- res$taus[which.max(res$post_probs)]
    mean_tau <- sum(res$taus * res$post_probs)
    cat("MAP Estimate of Yield Point:  ", map_tau, "\n")
    cat("Posterior Mean Yield Point: ", mean_tau, "\n")
  })
}

shinyApp(ui = ui, server = server)