
library('shiny')
library(shinyjs)

ui <- fluidPage(
  useShinyjs(),
  fluidRow(h1("Bias classification"), align="center"),
  fluidRow(wellPanel(
    fluidRow(
      column(12, sliderInput(inputId="n", 
                             label="N", 
                             value=1e3, 
                             min=100, 
                             max=1e5)
      )
    ),
    fluidRow(
      column(12, plotOutput(outputId="mainPlot", click="setX1", dblclick="setX2"))
    ),
    fluidRow(
      column(5, sliderInput("testN", "Test set size", value=100, min=100, max=1000)),
      column(4, selectInput("method", "Distribution generation: ", 
                            c("Standard" = "standard",
                              "Central limit theorem" = "clt",
                              "Polar coords" = "polar"))
      ),
      column(3, sliderInput("clt_k", "k", value=100, min=10, max=1000))
    ),
    fluidRow(
      column(6, 
          sliderInput("bandwidth1", "Bandwidth for kernel density estimation for X1", 
                      value = 1, min=0.01, max = 2)),
      column(6, 
          sliderInput("bandwidth2", "Bandwidth for kernel density estimation for X2", 
                      value = 1, min=0.01, max = 2))
    )
  )
  ),
  fluidRow(
    column(6, 
           wellPanel(
             fluidRow(h3("First class", align="center")),
             fluidRow(
               column(12, 
                      sliderInput(inputId="p1", label="P(1)", value=0.5, min=0.01, max=0.99))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="m11", label="m11", value=0, min=-10, max=10))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="m12", label="m12", value=0, min=-10, max=10))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="sd11", label="sd11", value=1, min=0.1, max=10))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="sd12", label="sd12", value=1, min=0.1, max=10))
             )
           )
    ),
    column(6, 
           wellPanel(
             fluidRow(h3("Second class", align="center")),
             fluidRow(
               column(12,
                      sliderInput(inputId="p2", label="P(2)", value=0.5, min=0.01, max=0.99))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="m21", label="m21", value=5, min=-10, max=10))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="m22", label="m22", value=5, min=-10, max=10))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="sd21", label="sd21", value=1, min=0.1, max=10))
             ),
             fluidRow(
               column(12, 
                      sliderInput(inputId="sd22", label="sd22", value=1, min=0.1, max=10))
             )
           )
    )
  ),
  fluidRow(wellPanel(
    fluidRow(
      column(10, h2("Estimation for densities")),
      column(2, checkboxInput("estimation_show", "Show? (not easy to render for sure)", value = FALSE))
      
    ),
    fluidRow(
      column(6, plotOutput(outputId="estimation_1")),
      column(6, plotOutput(outputId="estimation_2")) 
    ) 
  )
  )
)

normal_clt <- function(n, m=0, sd=1, k=100) {
  numbers = c()
  
  for (j in 1:n) {
    rndsum = 0
    for (i in 1:k) {
      rndsum = rndsum + runif(1) - 0.5
    }
    numbers = c(numbers, sd*sqrt(12/k)*rndsum + m)
  }
  
  return(numbers)
}

normal_polar = function(m1=0, sd1=1, m2=0, sd2=1) {
  s = 1
  while(s >= 1) {
    v1 = 2*runif(1) - 1
    v2 = 2*runif(1) - 1
    s = v1^2 + v2^2
  }
  
  x1 = v1*sqrt((-2)*log(s)/s)
  x2 = v2*sqrt((-2)*log(s)/s)
  
  return(c(sd1*x1 + m1, sd2*x2 + m2))
}

K = function(z) {
  if (abs(z) <= 1) return(1 - abs(z))
  return(0)
}
getMean = function(x) {
  mean_x = rep(0, ncol(x))
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      mean_x[j] = mean_x[j] + x[i, j]
    }
  }
  return(mean_x/nrow(x))
}
getSd = function(x) {
  m_x = getMean(x)
  sd_x = rep(0, ncol(x))
  for (i in 1:nrow(x)) {
    for (j in 1:ncol(x)) {
      sd_x[j] = sd_x[j] + (x[i, j] - m_x[j])^2
    }
  }
  return(sd_x/nrow(x))
}
getP = function(x, n) {return(nrow(x)/n)}
estimate = function(x, sample, kernel=K, h=1.06*getSd(sample)*nrow(sample)^(-1/5)) {
  if (length(x) != ncol(sample) || length(x) != length(h)) stop("wrong args")
  result = 0
  for (i in 1:nrow(sample)) {
    kmult = 1
    for (j in 1:ncol(sample)) {
      kmult = kmult*kernel((x[j] - sample[i, j])/h[j])/h[j]
    }
    result = result + kmult
  }
  result = result/nrow(sample)
  return(result)
}
getRandPoint = function(dist_func, m1=0, sd1=1, m2=0, sd2=1, k=100) {
  if (identical(dist_func, normal_polar)) {
    return(dist_func(m1, sd1, m2, sd2))
  } else if (identical(dist_func, normal_clt)) {
    return(c(dist_func(1, m1, sd1, k=k), dist_func(1, m2, sd2, k=k)))
  } else {
    return(c(dist_func(1, m1, sd1), dist_func(1, m2, sd2)))
  }
}

server <- function(input, output, session) {
  values <- reactiveValues()
  observeEvent(input$p1, {
    updateSliderInput(session, "p2", value=1-floor(input$p1*100)/100)
  })
  observeEvent(input$p2, {
    updateSliderInput(session, "p1", value=1-floor(input$p2*100)/100)
  })
  observeEvent(input$method, {
    if(input$method == "clt") enable("clt_k")
    else disable("clt_k")
  })
  
  output$mainPlot <- renderPlot({
    method = c("standard" = rnorm,
               "clt" = normal_clt,
               "polar" = normal_polar)
    
    h=c(input$bandwidth1, input$bandwidth2)
    func = method[[input$method]]
    n = input$n
    p1 = input$p1
    p2 = 1 - p1
    
    m11 = input$m11
    m12 = input$m12
    m21 = input$m21
    m22 = input$m22
    sd11 = input$sd11
    sd12 = input$sd12
    sd21 = input$sd21
    sd22 = input$sd22
    
    x = c()
    y = c()
    for (i in 1:n) {
      if (runif(1) <= p1) x = c(x, getRandPoint(func, m11, sd11, m12, sd12))
      else y = c(y, getRandPoint(func, m21, sd21, m22, sd22))
    }
    x = matrix(x, ncol = 2, byrow = TRUE)
    y = matrix(y, ncol = 2, byrow = TRUE)
    n1 = nrow(x)
    n2 = nrow(y)
    plot(x, ylim = c(min(x[, 2], y[, 2]), max(x[, 2], y[, 2])), 
         xlim = c(min(x[, 1], y[, 1]), max(x[, 1], y[, 1])), pch = 1)
    points(y, pch = 2)
    
    ## render estimations
    renderEstimations = input$estimation_show
    if (renderEstimations) {
      est_plot_size = 100
      output$estimation_1 <- renderPlot({
        test = seq(min(x[, 1], x[, 2]), max(x[, 1], x[, 2]), length=est_plot_size)
        est = c()
        for (i in 1:est_plot_size) {
          for (j in 1:est_plot_size) 
            est = c(est, estimate(c(test[i], test[j]), x, h=h))
          if (i %% 10 == 0)
            print(paste0(i/est_plot_size*100, "% rendered of 1 estimation plot"))
        }
        est = matrix(est, nrow = est_plot_size)
        image(test, test, est, xlab = "X1", ylab = "X2", main="FIRST")
      })
      output$estimation_2 <- renderPlot({
        test = seq(min(y[, 1], y[, 2]), max(y[, 1], y[, 2]), length=est_plot_size)
        est = c()
        for (i in 1:est_plot_size) {
          for (j in 1:est_plot_size) 
            est = c(est, estimate(c(test[i], test[j]), y, h=h))
          if (i %% 10 == 0)
            print(paste0(i/est_plot_size*100, "% rendered of 2 estimation plot"))
        }
        est = matrix(est, nrow = est_plot_size)
        image(test, test, est, xlab = "X1", ylab = "X2", main="SECOND")
      })
    }
    
    test_len = input$testN
    test_points = c()
    spread_coef = 1.5
    class = c()
    for (i in 1:test_len) {
      k = if (input$clt_k) input$clt_k else 100
      if (runif(1) <= p1) {
        test_points = c(test_points, getRandPoint(func, m11, sd11*spread_coef, m12, sd12*spread_coef, k=k))
        class = c(class, 1)
      } else {
        test_points = c(test_points, getRandPoint(func, m21, sd21*spread_coef, m22, sd22*spread_coef, k=k))
        class = c(class, 2)
      }
    }
    # count error
    estimation = c()
    error = 0
    unclear = 0
    for (i in 1:test_len) {
      est_x = estimate(c(test_points[i*2 - 1], test_points[i*2]), x, h=h)*p1
      est_y = estimate(c(test_points[i*2 - 1], test_points[i*2]), y, h=h)*p2
      
      col = "green"
      
      if (est_x - est_y < 0) {
        estimation = c(estimation, 2)
      } else if (est_x - est_y > 0) {
        estimation = c(estimation, 1)
      } else {
        estimation = c(estimation, 0)
        col = "gray"
      }
      
      if (estimation[i] == 0) {
        unclear = unclear + 1
      } else if (estimation[i] != class[i]) {
        error = error + 1
        col = "red"
      }
      print(paste0(i, " - ",estimation[i],":", col, ", ", est_x-est_y))
      points(test_points[i*2-1], test_points[i*2], pch=estimation[i]+15, col=col, cex=1.5)
    }
    print(paste(error/test_len*100, "%; Unclear ", unclear/test_len*100, "%"))
    legend("bottomright", 
           legend=paste0("Error ", error/test_len*100, "%; Unclear ", unclear/test_len*100, "%"))
  })
  
}

shinyApp(ui = ui, server = server)

