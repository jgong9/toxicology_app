#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
# Created by Joonho Gong
# Email : jgong9@ncsu.edu


## For Shiny and plots
library(shiny)
library(shinythemes)
library(DT)
library(latex2exp)
library(shinydashboard)
library(markdown)
library(shinyWidgets)

## For model fitting
library(readxl)
library(minpack.lm)
library(nlstools)


##### Read in the original data file
isotherm.meta<-read_excel("data/Glyphosate_Paraquat_master.xlsx",sheet="isotherm.meta")
isotherm.data <- read_excel("data/Glyphosate_Paraquat_master.xlsx",sheet="isotherm.data",skip=10, .name_repair = "minimal")
data_set <- isotherm.data
data_meta <- data.frame(isotherm.meta)

n_datasets <- ncol(isotherm.data) / 2

x_matrix <- matrix(0, n_datasets, 101)
n_vec <- vector("numeric", length=n_datasets)



##### Functions required
## 1. Langmuir model
negLogLikelihood <- function(par, data_x, data_y) {

  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(data_y, mean = (par[1] * data_x * exp(par[2])) / (1 +  exp(par[2]) *data_x) , sd = par[3])
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle <- function(Q, K_new, C, sigma, cov_mat){
  der_vec <- rep(NA, 3)
  der_vec[1] <- (exp(K_new) * C) / (1 + exp(K_new) * C)
  der_vec[2] <- (Q  * exp(K_new) * C) / (1 + exp(K_new) * C)^(2) 
  der_vec[3] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec  ) 
  
}


## 2. Freundlich model
negLogLikelihood_fre <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(data_y, mean = (par[1] * data_x^(1/par[2]) ), sd = par[3])
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

negLogLikelihood_fre_trans <- function(par, data_x, data_y) {
  
  
  if(par[3] < 0){
    output <- 10000000
  } else {
    likelihoods <- dnorm(data_y, mean = (exp(par[1]) * data_x^(1/par[2]) ), sd = par[3])
    log.likelihoods <- log(likelihoods)
    output <- -sum(log.likelihoods)
    return(output)
  }
} # end of negLogLikelihood

se_for_CI_mle_fre <- function(K_F, n, C, sigma, cov_mat){
  der_vec <- rep(NA, 3)
  der_vec[1] <- C^(1/n)
  der_vec[2] <- K_F * log(C) * C^(1/n) * (-1/n^(2))
  der_vec[3] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec ) 
  
}

se_for_CI_mle_fre_trans <- function(K_F_prime, n, C, sigma, cov_mat){
  der_vec <- rep(NA, 3)
  der_vec[1] <- exp(K_F_prime) * C^(1/n)
  der_vec[2] <- exp(K_F_prime) * log(C) * C^(1/n) * (-1/n^(2))
  der_vec[3] <- 0
  
  se <- sqrt( t(der_vec) %*% cov_mat %*% der_vec ) 
  
}
##### End of functions required







sidebar <- dashboardSidebar(
  sidebarMenu(id="tabs",
              menuItem("Models", tabName = "models", icon=icon("equals"), selected = T),
              menuItem("Initial values", tabName="initial", icon=icon("file-text-o"), selected=F),
              menuItem("Confidence envelopes", tabName = "conf", icon=icon("line-chart"), selected=F),
              menuItem("AIC", tabName = "aic", icon=icon("pencil"), selected=F)
            #  menuItem("Question", tabName = "ques", icon = icon("question"))
  )
)

body <- dashboardBody(
  tabItems(
    tabItem(tabName = "models",
            br(),
            fluidRow(
              column(
                p(
                  strong("1. Langmuir model"),
                  br(),
                  br(),
                  "Suppose we observed \\(n\\) pairs of \\((C_{W,i}, q_{i})_{i=1,\\dots,n}\\) of isotherm curve. Then, our model is ",
                  "$$q_{i} = \\frac{Q_{max} \\exp(K^{\\prime}_{d}) C_{W,i}}{ 1 + \\exp(K^{\\prime}_{d}) C_{i}} + \\varepsilon_{i},$$",
                  "where \\(\\varepsilon_{i}\\) is a measurement error distributed with \\(N(0, \\sigma^{2})\\) and \\(K_{d}\\) is reparameterized with \\(K^{\\prime}_{d} = \\log(K_{d})\\). In order to estimate those three parameters, \\(Q_{max}\\), \\(K^{\\prime}_{d}\\), and \\(\\sigma^{2}\\), we use the maximum likelihood estimation method based the likelihood function given as ",
                  "$$ L\\left( Q_{max}, K^{\\prime}_{d}, \\sigma^{2} \\right) = \\prod^{n}_{i=1} \\frac{1}{\\sqrt{2\\pi \\sigma^{2}}} \\exp \\left( -  \\frac{ \\left( q_{i} - \\frac{Q_{max} \\exp(K^{\\prime}_{d}) C_{W,i}}{ 1 + \\exp(K^{\\prime}_{d}) C_{W,i}} \\right)^{2} }{2\\sigma^{2}} \\right).$$",
                  "The 3 estimators which maximize \\(L\\left( Q_{max}, K^{\\prime}_{d}, \\sigma^{2} \\right)\\) are our maximum likelihood estimators (MLE).",
                  style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"
                ),
                br(),
                p(
                  strong("2. Freundlich model"),
                  br(),
                  br(),
                  "Similarly, given the same N pairs of observations \\((C_{W,i}, q_{i})_{i=1,\\dots,N}\\), we assume that ",
                  "$$ q_{i} = \\exp(K^{\\prime}_{F}) \\cdot C_{W,i}^{\\frac{1}{n}}+ \\varepsilon_{i},$$ ",
                  "where \\(K_{F} = \\exp(K^{\\prime}_{F})\\) and the distribution of the measurement error \\( \\varepsilon_{i}\\) is \\(N(0, \\sigma^{2})\\). As we have shown for Langmuir model, we set the likelihood function as ",
                  "$$ L\\left(K^{\\prime}_{F}, n, \\sigma^{2} \\right)=\\prod^{N}_{i=1} \\frac{1}{\\sqrt{2\\pi \\sigma^{2}}} \\exp \\left( - \\frac{ \\left( q_{i} - \\exp(K^{\\prime}_{F}) \\cdot C_{W,i}^{\\frac{1}{n}} \\right)^{2} }{2\\sigma^{2}} \\right),$$",
                  "and then attain the MLE that maximizes \\(L\\left(K^{\\prime}_{F}, n, \\sigma^{2} \\right)\\)",
                  style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"
                ),
                width=11),
            ),
    ), # end of tabItme 0
    tabItem(tabName = "initial",
              br(),

            
              fluidRow(
                column(
                  p(
                    strong("1. Langmuir model"),
                    br(),
                    br(),
                       "In order to get the initial values for \\(Q_{max}\\), \\(K_{d}\\), and \\(\\sigma^{2}\\), we transform them as 
                       \\(Q^{\\star}_{max} = Q_{max}^{-1}\\), \\(K^{\\star}_{d} = K_{d}^{-1}\\). Let \\(Y^{\\star}=q^{-1}\\) and \\(C^{\\star}_{W}= C^{-1}_{W}\\). Then, ",
                       br(),
                       br(),
                       "\\( q = \\frac{Q_{max} K_{d} C_{W}}{1 + K_{d} C_{W}} \\Longleftrightarrow q = \\frac{Q_{max} C_{W}}{K_{d}^{\\star} +  C_{W}}\\)",
                       br(),
                      "\\(\\Longrightarrow \\frac{1-q}{q} = \\frac{K_{d}^{\\star} + C_{W} - Q_{max}C_{W}}{Q_{max}C_{W}}\\)",
                      br(),
                      "\\(\\Longleftrightarrow \\frac{1-q}{q} = Q_{max}^{\\star} + K_{d}^{\\star} Q_{max}^{\\star} C_{W}^{\\star} - 1\\)",
                      br(),
                      "\\(\\Longrightarrow Y^{\\star} =  Q_{max}^{\\star} + K_{d}^{\\star} Q_{max}^{\\star} C_{W}^{\\star}\\).",
                      br(),
                    br(),
                      "Since we now have a linear model with an intercept \\(Q_{max}^{\\star}\\) and a slope \\(K_{d}^{\\star}Q_{max}^{\\star}\\), we can apply the ordinary least square (OLS) method to estimate both of them and get the initial values. 
                      Specifically, we set \\(\\hat{Q}_{max} = \\frac{1}{\\hat{Q}^{\\star}_{max}}\\) and \\(\\hat{K}_{d}=\\frac{\\hat{Q}^{\\star}_{max}}{\\widehat{K_{d}^{\\star} Q_{max}^{\\star}}}\\). For \\(\\sigma^{2}\\), we adopt the sum of squared errors \\(\\hat{\\sigma}^{2} = \\sum^{n}_{i=1} ( y^{\\star}_{i} - \\hat{y}^{\\star}_{i}  ). \\)", 
                      style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"
                    ),
                  br(),
                  p(
                    strong("2. Freundlich model"),
                    br(),
                    br(),
                    "Similarly, we reparameterize the Freundlich model with \\(q^{\\star}=\\log{q}\\), \\(C^{\\star}_{W}=\\log{C_{W}}\\), \\(K^{\\star}_{F} = \\log{K_{F}}\\), and \\(n^{\\star} = \\frac{1}{n}\\) . Then, ",
                    br(),
                    br(),
                    "\\( q^{\\star} = K^{\\star}_{F} + n^{\\star} C^{\\star}_{W}\\).",
                    br(),
                    br(),
                    "Since we now have a linear model with an intercept \\( K^{\\star}_{F}\\) and a slope \\(n^{\\star}\\), we can apply OLS again to estimate both of them and get the initial values. 
                      In particular, we set \\(\\hat{K}_{F} = exp(\\hat{K}^{\\star}_{F}) \\) and \\(\\hat{n}= \\frac{1}{\\hat{n}^{\\star} }\\). For \\(\\sigma^{2}\\), we adopt the sum of squared errors \\(\\hat{\\sigma}^{2} = \\sum^{n}_{i=1} ( q^{\\star}_{i} - \\hat{q}^{\\star}_{i}  ). \\)", 
                    style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"
                ),
                  width=11),
            ),
    ), # end of tabItme 1
    tabItem(tabName = "conf",
            br(),
            
            fluidRow(
              column(
                p(
                  strong("Delta method"),
                  br(),
                  br(),
                  "Let \\( \\mathbf{\\theta}_{0} = [Q_{max}, \\: K^{\\prime}_{d}, \\: \\sigma^{2} ]\\) and denote its maximum likelihood estimator (MLE) by \\( \\hat{\\mathbf{\\theta}}_{MLE}\\). Since the aymptotic distribution of MLE is multivariate normal, we know ",
                  "$$ \\sqrt{n}( \\hat{\\mathbf{\\theta}}_{MLE} - \\mathbf{\\theta}_{0} ) \\; \\overset{\\cdot}{\\sim}\\; N( \\mathbf{0}, \\mathbf{I}^{-1}(\\mathbf{\\theta}_{0}) ))$$",
                  "where \\(\\mathbf{I}\\) is the Fisher information matrix. Next, define a continuous function \\( g(\\mathbf{\\theta}_{0})=\\frac{Q_{max}exp(K^{\\prime}_{d}) C_{W} }{1+exp(K^{\\prime}_{d})C_{W} }\\). Since its partial derivaties exist, the delta method is applicable and 
                  we can get the \\((1-\\alpha)\\)% confidence envelops as below.",
                  "$$ g(\\hat{\\mathbf{\\theta}}_{MLE}) \\pm z_{\\frac{\\alpha}{2}} \\times \\sqrt{ \\mathbf{g}^{\\prime T}(\\hat{\\mathbf{\\theta}}_{MLE}) \\mathbf{H}^{-1}(\\hat{\\mathbf{\\theta}}_{MLE}) \\mathbf{g}^{\\prime}(\\hat{\\mathbf{\\theta}}_{MLE})  }, $$",
                  "where \\(\\mathbf{g}^{\\prime}\\) is derivative of the function \\(g\\), \\(\\mathbf{H}\\) is numerical hessian matrix, and \\(\\alpha\\) is a given significance level with corresponding critical value \\(z\\). ",
                  br(),
                  br(),
                  "Similarly, the confidence envelops for the Freundlich model can be obtained through the delta method with a continuous function \\(h(\\mathbf{\\theta})=exp(K_{F}) C^{\\frac{1}{n}}_{W}\\).",
                  style="text-align:justify;color:black;background-color:#D5F5E3;padding:15px;border-radius:10px"
                ),
                width=11),
            )

    ), # end of tabItem 2
    tabItem(tabName = "aic",
            fluidRow(
              column(
                br(),

                
                p(
                  strong("Akaike information criterion (AIC)"),
                  br(),
                  br(),
                  "Recall the models with Gaussian random error. Let k be the number of parameters in each model and \\(L\\)\ the likelihood function for the model. Then the AIC value of the model is defined as ",
                  "$$ \\text{AIC} = 2k - \\ln{L(\\hat{\\mathbf{\\theta}}_{MLE} )}.$$",
                  style="text-align:justify;color:black;background-color:#D8E8F1;padding:15px;border-radius:10px"
                ),
                width=11),
            )
    )
    
  ), # end of tabItems
) # end of dashboardBody

######################################################################################






#####################################################################################


ui <- tagList(
  tags$style("@import url(https://use.fontawesome.com/releases/v5.6.0/css/all.css);"),
  tags$head(HTML("<title>Toxicology</title> <link rel='icon' type=image/gif/png href='warning.png'>")),
  
  navbarPage(
        title = span("Toxin Analysis", style="font-weight: bold; font-weight: 900;"), 
        theme = shinytheme("flatly"),
       tabPanel("Model",
                style='padding-left:10px; padding-right:25px; padding-top:10px; padding-bottom:5px',
                
                fluidRow(
                  column(
                    br(),
                    div( img(src="acs-logo.svg",width="169px",height="54px"),  align= "center"),
                    p("For more information, please see ",em("Wang et al. (2019)"),
                      br(),
                      a(href="https://pubs.acs.org/doi/10.1021/acsomega.9b02051", "Here",target="_blank",style = "color: blue;"),style="text-align:center;color:black"),
                    br(),
                    br(),
                    br(),
                    div(img(height = 50, width = 140, src = "tam_logo2.png"),
                        strong("X"),
                        img(height = 50, width = 140, src = "ncstate_logo.png"),
                        align= "center"
                    ),
                    br(),
                    div("This is a product of   ", 
                        
                        a(actionButton(inputId = "email1", label = "Joonho Gong", 
                                       icon = icon("envelope", lib = "font-awesome"),
                                       style="color: #0303FC; background-color: #EBEDEF; border-color: #DCDCDC"
                        ),
                        href="mailto:jgong9@ncsu.edu"),
                        align= "center"
                    ),
                    width=3),
                  
                  
                  column(
                    h4(strong("Isotherm Models")),
                    p( strong("1. Langmuir model"),
                       withMathJax("$$ q = \\frac{Q_{max} K_{d} C_{W}}{1+K_{d} C_{W} },$$"),
                       "where \\(q\\) is toxin absorbed (mol/kg), \\(Q_{max}\\) is the maximum capacity (mol/kg), \\(K_{d}\\) is a distribution constant, and \\(C_{W}\\) is the toxin equilibirum conncentration.",style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                    br(),
                    
                    p( strong("2. Freundlich model"),
                       withMathJax("$$ q = K_{F} \\cdot C_{W}^{\\frac{1}{n}} $$"),
                       "where \\(q\\) is toxin absorbed (mol/kg), \\(K_{F}\\) is the Freundlich constant, \\(\\frac{1}{n}\\) is the measure of intensity, and \\(C_{W}\\) is the toxin equilibirum conncentration.",
                       style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                    br(),
                    
                    width=8)
                  
                )
       ), ## end of tab 1
        tabPanel("Plot",
                 sidebarLayout(
                     sidebarPanel(
                         selectizeInput('main', 'Toxicant', c("Glyphosate", "Paraquat")
                                        ),
                      
                         conditionalPanel(
                             condition = "input.main=='Glyphosate'",
                             sliderInput("k",
                                         "Data ID of Glyphosate:",
                                         min = 1,
                                         max = 8,
                                         value = 1,
                                         step=1,
                                         round =T),
                         ),
                         conditionalPanel(
                             condition = "input.main=='Paraquat'",
                             sliderInput("k2",
                                         "Data ID of Paraquat:",
                                         min = 9,
                                         max = 12,
                                         value = 9,
                                         step=1,
                                         round=T),
                         ),
                         
                         selectizeInput('Rfunction', 'Choose your R optimization function', c("nlm", "optim")
                                        ),
                         
                         conditionalPanel(
                           condition = "input.Rfunction=='optim'",
                           selectizeInput('method', 'Method of minimization', c( "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"  ))
                         ),

                        selectizeInput("conf_level", "Confidence level for MLE band", c("99%", "95%", "90%", "80%")),
                        
                         br(),
                         br(),
                         br(),
                         br(),
                         br()
                         
                     ),
                     mainPanel(
                                h3(strong("1. Langmuir Model")),
                                h4("1-1. Non-transformed OLS vs Transformed MLE"),
                                plotOutput("trajectoryPlot"),
                                withMathJax("For MLE, we assume \n
                                     $$ q = \\frac{Q_{max} \\cdot \\exp{K_{d}^{\\prime}} \\cdot C_{W}}{(1+ \\exp{K_{d}^{\\prime}} \\cdot C_{W})} + \\varepsilon, $$
                                     \n where \\(K_{d} = \\exp(K_{d}^{\\prime})\\) and \\( \\varepsilon \\sim N(0, \\sigma^{2}) \\).
                                     \n"),
                                br(),
                                br(),
                                h4("1-2. MLE summary of Langmuir"),
                                tableOutput("mle_table"),
                        
                                h4("1-3. AIC of Langmuir"),
                                verbatimTextOutput("aic_lan"),
                                
                                br(),
                                br(),
                                h3(strong("2. Freundlich vs Langmuir")),
                                h4("2-1. Trajectory plot"),
                                plotOutput("plot_fre"),
                                withMathJax("For MLE, we assume \n
                                     $$ q = \\exp{K^{\\prime}_{F}} \\cdot C_{W}^{\\frac{1}{n}}+ \\varepsilon, $$
                                     \n where \\(K_{F} = \\exp(K^{\\prime}_{F})\\) and \\( \\varepsilon \\sim N(0, \\sigma^{2}) \\).
                                     \n"),
                                br(),
                                br(),
                                h4("2-2. MLE summary of Freundlich"),
                                tableOutput("mle_table_fre"),
                                h4("2-3. AIC of Freundlich"),
                                verbatimTextOutput("aic_fre"),
                                
                                br(),
                                br()
                                
                                )  # end of mainpanel
                     ) # end of sidebar layout
        ), ## end of tab 2

   

        tabPanel("Data",
                 style='padding-left:50px; padding-right:50px; padding-top:40px; padding-bottom:5px',
                 
                 dataTableOutput("table")
        ), ## end of tab 3  
  
       
       
       tabPanel("Method",
                
              dashboardPage(
                        dashboardHeader(title = NULL, titleWidth = 0, disable =T
                                        ),
                        sidebar,
                        body
                      ),

              tags$style(type = "text/css", ".container-fluid {padding-left:0px;
                    padding-right:0px;}"),
              tags$style(type = "text/css", ".navbar {margin-bottom: .0px;padding-left:20px}")

       ), ## end of tab 4  

navbarMenu("Do-It-Yourself Tool",
      tabPanel("Your Data",
           style='padding-left:00px; padding-right:50px; padding-top:00px; padding-bottom:5px',
         sidebarLayout(
           sidebarPanel(

             fileInput("file1", "Upload your CSV file",
                       multiple = TRUE,
                       accept = c(
                                  ".csv"
                       )),

             
             # Input: Checkbox if file has header ----
             prettyRadioButtons(
               inputId = "header",
               label = "Header",
               thick = TRUE,
               choices = c(Yes = "Yes",
                           No = "No"
                           ),
               selected = "Yes",
               outline = T,
               shape = "round",
               animation = "pulse",
               fill = F,
               status = "primary"
             ),
             
             # Input: Select separator ----
             prettyRadioButtons(
               inputId = "sep",
               label = "Separator",
               thick = TRUE,
               choices = c(Comma = ",",
                           Semicolon = ";",
                           Tab = "\t"),
               selected = ",",
               outline = T,
               shape = "round",
               animation = "pulse",
               fill = F,
               status = "primary"
             ),
             
             # Input: Checkbox if file has the pair of zeors at first ----
             prettyRadioButtons(
               inputId = "zeros",
               label = "A pair of two zeros in the first row",
               thick = TRUE,
               choices = c(Yes = "Yes",
                           No= "No"
               ),
               selected = "Yes",
               outline = T,
               shape = "round",
               animation = "pulse",
               fill = F,
               status = "primary"
             ),
             
             # Input: Select quotes ----
             prettyRadioButtons(
               inputId = "quote",
               label = "Quote",
               thick = TRUE,
               choices = c(None = "",
                           "Double Quote" = '"',
                           "Single Quote" = "'"),
               selected = "",
               outline = T,
               shape = "round",
               animation = "pulse",
               fill = F,
               status = "primary"
             )
             ,
             
             fluidRow(
               column(
                 downloadButton("downloadData_ex", "Download an example file"),
                 width=12,align= "center")
             )
             
             
           ) # end of sidebaPanel
           ,
           # Main panel for displaying outputs ----
           mainPanel(
             hr(),
             # Output: Data file ----
             dataTableOutput("user_data")
             
           )
           
         ) # end of sidbarLayout
         

      ) ## end of Your Data
      ,
      
      tabPanel(
        "Your Model",
        sidebarLayout(
          sidebarPanel(
            
            selectizeInput('user_model', 'Choose your isothem model', c("Langmuir", "Freundlich")
            ),
            
            selectizeInput('user_Rfunction', 'Choose your R optimization function', c("nlm", "optim")
            ),
            
            conditionalPanel(
              condition = "input.user_Rfunction=='optim'",
              selectizeInput('user_method', 'Method of minimization', c( "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN"  ))
            ),
            
            selectizeInput("user_conf_level", "Confidence level for MLE band", c("99%", "95%", "90%", "80%")),
            
          ) # end of sidebarPanel
          ,
          mainPanel(
            h3(strong("Your Model")),
            h4("1. Plot"),
            plotOutput("user_trajectoryPlot"),
            hr(),
            h4("2. Estimates"),
            tableOutput("mle_table_user"),
          )  # end of mainpanel
        ) # end of sidebar layout
        
        
      ) ## end of Your Model
      ,
      tabPanel(
        "Download results",
        sidebarLayout(
          
          # Sidebar panel for inputs ----
          sidebarPanel(
            fluidRow(
            # Input: Choose dataset ----
            column(
            selectInput("down_dataset", "Choose the output file you want:",
                        choices = c("Summary", "Trajactory")),
            align = "center", width = 12
            ),
            # Button
            column(
              downloadButton("downloadData", "Download"),
              align = "center", width = 12
            )
            )
          ) # end of sidePanel
          ,
          
          # Main panel for displaying outputs ----
          mainPanel(
            br(),
            tableOutput("user_down")
            
          ) # end of mainPanel
          
        ) ## end of sidbarLayout
        
      ) ## end of Download
      
), ### end of navarMenu



      navbarMenu("More",
               tabPanel("Resources",
                        br(),
                        br(),
                        fluidRow(
                          column("",width=1),
                          column(
                          box(title = "Github", HTML("R Shiny code is available on github. <br/>  <br/>"), width=15,
                              column(12,align="center",
                                              actionButton(inputId='ab1', label=HTML("<img src='github_logo.png' height='40'/> &nbsp; Go to GitHub"), 
                                                         #         icon = icon("github "), 
                                                           style='padding:10px; font-size:160%',
                                                           
                                                                  onclick ="window.open('https://github.com/jgong9/toxicology_app', '_blank')")
                              )
                        , status = "primary",solidHeader = T, collapsible = TRUE,
                        ), width=5)
                        ,
                        
                          
                        column(
                              box(title = "Personal Shiny Server", HTML("This app is accessible through Shiny Server on Raspberry Pi. <br/>  <br/>"), 
                                  width=15,
                                  column(12, align="center",
                                  actionButton(inputId='rasp_button', label=HTML("<img src='raspberry_pi_logo.png' height='40'/> &nbsp; Visit Shiny Server"), 
                                            #   icon = icon("linux"), 
                                               style='padding:10px; font-size:160%',
                                               onclick ="window.open('http://174.99.0.141:9000/toxicology_app', '_blank')")
                                  
                                 ) , status = "danger",solidHeader = T, collapsible = TRUE,
                              ), width=5)
                  ),  
                  fluidRow(
                    column("",width=1),
                    column(
                      box(title = "Contact", HTML("Contact me for any questions or inquiries. <br/>  <br/>"), 
                          width=15,
                          column(12, align="center",
                                 a(actionButton(inputId = "email2", label = HTML("&nbsp; Email Joonho Gong"), 
                                                icon = icon("envelope", lib = "font-awesome"),
                                                style='padding:10px; font-size:160%')
                                 ,
                                 href="mailto:jgong9@ncsu.edu"
                                 )
                          ) , status = "warning",solidHeader = T, collapsible = TRUE,
                      ), width=5)
                  ,  
                  column(
                    box(title = "R Shiny", HTML("This app is available with Shinyapps.io <br/>  <br/>"), width=15,
                        column(12,align="center",
                               actionButton(inputId='ab1', label=HTML("<img src='rstudio_logo3.png' height='40'/> &nbsp; Visit Shinyapps.io"), 
                                            #         icon = icon("github "), 
                                            style='padding:10px; font-size:160%',
                                            
                                            onclick ="window.open('https://jgong9.shinyapps.io/chemistry_app/', '_blank')")
                        )
                        , status = "info",solidHeader = T, collapsible = TRUE,
                    ), width=5)
                  ),
                  style='padding-left:20px; padding-right:10px; padding-top:0px; padding-bottom:5px'
                ),
                tabPanel("About",
                         fluidPage(
                           includeMarkdown("README.md")
                         ),
                         style='padding-left:30px; padding-right:30px; padding-top:10px; padding-bottom:5px',

                )
         ) ## end of tab 4

    ) # end of navbarPage
) # end of tagList





server <- function(input, output) {
    

    inputVar <- reactive({
        if(input$main =="Glyphosate") {
            k_val <- input$k
        } else if(input$main =="Paraquat"){
            k_val <- input$k2
        }


        value <- c(k_val, input$main, 
                   input$conf_level, input$Rfunction, input$method,
                   input$model)
    })
    
    output$otherString <- renderText({
        value <- inputVar()

        paste0(value[2], ", ", value[1])
    })
 
  
    output$trajectoryPlot <- renderPlot({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      ### 1. Original fit
      fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
                           data=data_kth,
                           start=list(Qmax.mol.kg = initial_values[[1]],
                                      Kd.l.mol = initial_values[[2]]))
      
      est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
      max_y <- max( est_curve_origin )
      summary_1 <-summary(fit_nonBoot)
      
      mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      
      ### 3. MLE
      # 3-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                  sqrt(sigma2_initial))

  
      # 3-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      # 3-3 get estimates
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]

      # 3-4 prepare for plotting
      
      Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * data.matrix(data_kth[,1]) ) / (1 + exp(coef_mle[2]) * data.matrix(data_kth[,1]))
      res_mle <- data_kth[,2] - Y_original_pred_mle
      mse_3 <- (sum(res_mle^(2)) ) / (n_data_kth  )
      # compute mse
      
      
      fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_matrix[kth, ]) / (1 + exp(coef_mle[2]) * x_matrix[kth, ])
      # fitted values based on mle for plotting 
      
      
      # 3-5 MLE Confidence band 
      if( input$conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est <- solve(fit_mle$hessian)
      c_pred_vec <- x_matrix[kth, ]
      
      mle_band_upper <- rep(NA, length(c_pred_vec))
      mle_band_lower <- rep(NA, length(c_pred_vec))
      
      for(i in 1:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                sigma = sigma_est_mle, cov_mat=Cov_est)
        
        term <- conf_val * se_est
        mle_band_upper[i] <- fitted_values_mle[i] + term
        mle_band_lower[i] <- fitted_values_mle[i] - term
      }
      
      
      ### 4. Plot
      
      max_final <- max(c(max_y, Y_original_pred, unlist(data_kth[-1,2]), mle_band_upper ))
      
      # main plot
      plot(x_matrix[kth, ], est_curve_origin, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
           ylim = c(0, max_final ), lwd = 2, col= 'blue',
           main=paste("NonTrans-NonLinear-OLS vs Trans-NonLinear-MLE; ", value[2], ", Data ID: ", value[1]))
      lines(c(x_com_vec), c(Y_com_pred), col="green", lty=2, lwd=2)
      lines(x_matrix[kth, ], fitted_values_mle, col="red", lty=3, lwd=2)
      
      # Confidence band
      lines(c_pred_vec, mle_band_upper, col="black", lty=3, lwd=2)
      lines(c_pred_vec, mle_band_lower, col="black", lty=3, lwd=2)
      
      # observed data points
      points(data_set[,(2* kth -1):(2*kth)])
      
      # legend
      legend("bottomright",legend=c(
        paste("NonTrans+NonLinear+OLS : Qm=",
              signif(coefficients(fit_nonBoot)["Qmax.mol.kg"],6),
              ", Kd=",round(coefficients(fit_nonBoot)["Kd.l.mol"]),
              ", MSE=",format(mse_1, scientific = T),sep=""),
        paste("Trans+NonLinear+MLE : Qm=",
              signif(Q_est_mle,6),
              ", Kd=",round(K_est_mle),
              ", MSE=",format(mse_3, scientific = T),sep=""),
        paste("Commercial Product : Qm=",
              signif(Q_com,6),
              ", Kd=",round(K_com),
              ", MSE=",format(mse_0,scientific = T),sep=""),
        paste( input$conf_level, " Confidence band of MLE"),
        paste("Observed data point")
      ),
      bty='n',
      col=c("blue","red","green", "black", "black"),lwd=c(2, 2, 2, 2, NA),lty=c(1,2,3, 3,NA),
      density = c(0, 0, 0, 0, NA), fill = c("blue", "red","green", "black", "white"),
      border=c(NA, NA, NA, NA, NA),
      pch = c(NA, NA, NA, NA, 1),
      cex=1
      )
    })
    
    
    
    output$mle_table <- renderTable({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      ### 1. Original fit
      fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
                           data=data_kth,
                           start=list(Qmax.mol.kg = initial_values[[1]],
                                      Kd.l.mol = initial_values[[2]]))
      
      est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
      max_y <- max( est_curve_origin )
      summary_1 <-summary(fit_nonBoot)
      
      mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute rse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      ### 3. MLE
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
  
      
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       #    ndigit = 20,
                       #    method = "CG",
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      
      
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]

      
      Cov_est <- solve(fit_mle$hessian)
      std_error1 <- sqrt(Cov_est[1,1]  )
      std_error2 <- sqrt(Cov_est[2,2]  )
      std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
      std_error3 <- sqrt(Cov_est[3,3]  )
      
      if( input$conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if (input$conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
      SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
      upper <- c(coef_mle[1] + conf_val * std_error1,
                 coef_mle[2] + conf_val * std_error2, 
                 exp(coef_mle[2] + conf_val * std_error2),
                 # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                 coef_mle[3] + conf_val * std_error3 )
      lower <- c(coef_mle[1] - conf_val * std_error1,
                 coef_mle[2] - conf_val * std_error2, 
                 exp(coef_mle[2] - conf_val * std_error2),
                 # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                 coef_mle[3] - conf_val * std_error3 )
      
      

      table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T), Lower=format(lower, scientific = T))
      

    },
    bordered=T)
    
    
    output$aic_lan <- renderText({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      ### 1. Original fit
      fit_nonBoot <- nlsLM(q.mol.kg ~ Qmax.mol.kg * Kd.l.mol * Cw.mol.l / (1 + Kd.l.mol*Cw.mol.l),
                           data=data_kth,
                           start=list(Qmax.mol.kg = initial_values[[1]],
                                      Kd.l.mol = initial_values[[2]]))
      
      est_curve_origin <- predict(fit_nonBoot, list( Cw.mol.l=x_matrix[kth,] ) )
      max_y <- max( est_curve_origin )
      summary_1 <-summary(fit_nonBoot)
      
      mse_1 <- sum( (summary_1$residuals)^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute rse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      ### 3. MLE
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       #    ndigit = 20,
                       #    method = "CG",
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      
      
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]
      

      
      AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2]))
      
      paste("AIC of Langmuir model is ", AIC_lan)
          
      })
    
    

###################################################################
#### Model 2. Freundlich
    
    
    output$plot_fre <- renderPlot({
      value <- inputVar()
      
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      Q_com <- initial_values[[1]]
      K_com <- initial_values[[2]]
      
      x_com_vec <- x_matrix[kth,]
      Y_com_pred <- (Q_com * K_com * x_com_vec) / (1 + K_com * x_com_vec)
      
      mse_0 <- sum( (data_kth[,2] - Y_com_pred )^2 ) / n_data_kth
      
      
      ### 2. Linear model
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      
      ### 3. MLE
      # 3-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      # 3-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle_lan <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle_lan <- fit_mle_lan$par
        
      } else if(value[4] == "nlm"){
        fit_mle_lan <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle_lan <- fit_mle_lan$estimate
        
      }
      
      # 3-3 get estimates
      Q_est_mle_lan <- coef_mle_lan[1]
      K_est_mle_lan <- exp( coef_mle_lan[2] )
      sigma_est_mle_lan <- coef_mle_lan[3]
      
      # 3-4 prepare for plotting
      
      Y_original_pred_mle_lan <- ( coef_mle_lan[1] * exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]) ) / (1 + exp(coef_mle_lan[2]) * data.matrix(data_kth[,1]))
      res_mle_lan <- data_kth[,2] - Y_original_pred_mle_lan
      mse_lan <- (sum(res_mle_lan^(2)) ) / (n_data_kth  )
      # compute mse
      
      
      fitted_values_mle_lan <-  ( coef_mle_lan[1] * exp(coef_mle_lan[2]) *  x_matrix[kth, ]) / (1 + exp(coef_mle_lan[2]) * x_matrix[kth, ])
      # fitted values based on mle for plotting 
      
      
      # 3-5 MLE Confidence band 
      if( input$conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est_lan <- solve(fit_mle_lan$hessian)
      c_pred_vec <- x_matrix[kth, ]
      
      mle_band_upper_lan <- rep(NA, length(c_pred_vec))
      mle_band_lower_lan <- rep(NA, length(c_pred_vec))
      
      for(i in 1:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle(Q = Q_est_mle_lan, K_new=coef_mle_lan[2], C=c_now,
                                sigma = sigma_est_mle_lan, cov_mat=Cov_est_lan
                                )
        
        term <- conf_val * se_est
        mle_band_upper_lan[i] <- fitted_values_mle_lan[i] + term
        mle_band_lower_lan[i] <- fitted_values_mle_lan[i] - term
      }
      
    
      ### 2. Freundlich model
      ## 2-1 initial value
      Y_trans <-  log(data_kth[-1,2] )
      X_trans <- log(data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      
      
      res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      

      
        ## 2-2. MLE
        initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))


        # 2-2 fit a new model of reparameterization for Kd
        if(value[4] == "optim"){
          fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                           par = initial_value_mle,
                           hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
          )
          coef_mle_fre <- fit_mle_fre$par

        } else if(value[4] == "nlm"){
          fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                         p = initial_value_mle,
                         hessian = T,
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
          )

          coef_mle_fre <- fit_mle_fre$estimate

        }

        # 2-3 get estimates
        K_F_prime_est_mle_fre <- coef_mle_fre[1]
        K_F_est_mle_fre <- exp(coef_mle_fre[1])

        n_est_mle_fre <- coef_mle_fre[2]
        sigma_est_mle_fre <- coef_mle_fre[3]

        # 2-4 prepare for plotting
        Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * data.matrix(data_kth[,1])^(1/coef_mle_fre[2])
        res_mle_fre <- data_kth[,2] - Y_original_pred_mle_fre
        mse_fre <- (sum(res_mle_fre^(2)) ) / (n_data_kth  )
        # compute mse


        fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_matrix[kth, ]^(1/coef_mle_fre[2])
        # fitted values based on mle for plotting


        # 2-5 MLE Confidence band
        if( input$conf_level == "95%"){
          conf_val <- qnorm(0.975)
        } else if( input$conf_level == "90%"){
          conf_val <- qnorm(0.95)
        }  else if (input$conf_level == "80%"){
          conf_val<- qnorm(0.90)
        } else if(input$conf_level == "99%"){
          conf_val<- qnorm(0.995)
        }

        Cov_est_fre <- solve(fit_mle_fre$hessian)
        c_pred_vec <- x_matrix[kth, ]

        mle_band_upper_fre <- rep(NA, length(c_pred_vec))
        mle_band_lower_fre <- rep(NA, length(c_pred_vec))

        for(i in 2:length(c_pred_vec)){
          c_now <- c_pred_vec[i]
          se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)

          term <- conf_val * se_est
          mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
          mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
        }

        mle_band_upper_fre[1] <- 0
        mle_band_lower_fre[1] <- 0
        # since derivative[2] does not exist when C= 0

      ### 4. Plot
      
      max_final <- max(c(unlist(data_kth[-1,2]), mle_band_upper_fre, mle_band_upper_lan ))
      
      # main plot
      plot(x_matrix[kth, ], fitted_values_mle_fre, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
           ylim = c(0, max_final ), lwd = 2, col= 'blue',
           main=paste( value[2], ", Data ID: ", value[1]))
      
      
      # Confidence band
      lines(c_pred_vec, mle_band_upper_fre, col="blue", lty=3, lwd=2)
      lines(c_pred_vec, mle_band_lower_fre, col="blue", lty=3, lwd=2)
      
      
      lines(c_pred_vec, fitted_values_mle_lan, col="red", lty=1, lwd=2)
      
      lines(c_pred_vec, mle_band_upper_lan, col="red", lty=3, lwd=2)
      lines(c_pred_vec, mle_band_lower_lan, col="red", lty=3, lwd=2)
      lines(c(x_com_vec), c(Y_com_pred), col="green", lty=2, lwd=2)
      
      # observed data points
      points(data_set[,(2* kth -1):(2*kth)])
      
  
      legend("bottomright",legend=c(
        paste("Langmuir : Qm=",
              signif(Q_est_mle_lan,6),
              ", Kd=",round(K_est_mle_lan),
              ", MSE=",format(mse_lan,scientific = T),sep=""),
        paste("Freundlich : K_F=",
              signif(K_F_est_mle_fre),
              ", n=",signif(n_est_mle_fre),
              ", MSE=", format(mse_fre,scientific = T ),sep=""),
        paste( input$conf_level, " Confidence band of Langmuir MLE"),
        paste( input$conf_level, " Confidence band of Freundlich MLE"),
        paste("Commercial Product : Qm=",
              signif(Q_com,6),
              ", Kd=",round(K_com),
              ", MSE=",format(mse_0, scientific = T),sep=""),
        paste("Observed data point")
      ),
      bty='n',
      col=c("red","blue","red", "blue", "green", "black"),lwd=c(2, 2, 2, 2,2, NA),lty=c(1,1,3,3,2,NA),
      density = c(0, 0, 0, 0, 0, NA), fill = c("red", "blue","red", "blue", "green", "white"),
      border=c(NA, NA, NA, NA, NA, NA),
      pch = c(NA, NA, NA, NA, NA, 1),
      cex=1
      )
    })
    
    
    
    output$mle_table_fre <- renderTable({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product

      
      ### 2. Freundlich model
      ## 2-1 initial value
      Y_trans <-  log(data_kth[-1,2] )
      X_trans <- log(data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      
      res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      

      
      # 2-2. MLE
       initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))


       # 2-2 fit a new model of reparameterization for Kd
       if(value[4] == "optim"){
         fit_mle <- optim(fn = negLogLikelihood_fre_trans,
                          par = initial_value_mle,
                          hessian = T,
                          method = value[5],
                          data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
         )
         coef_mle <- fit_mle$par

       } else if(value[4] == "nlm"){
         fit_mle <- nlm(f = negLogLikelihood_fre_trans,
                        p = initial_value_mle,
                        hessian = T,
                        data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
         )

         coef_mle <- fit_mle$estimate

       }



       K_F_prime_est_mle_fre <- coef_mle[1]
       K_F_est_mle_fre <- exp(coef_mle[1])

       n_est_mle_fre <- coef_mle[2]
       sigma_est_mle <- coef_mle[3]


       Cov_est <- solve(fit_mle$hessian)
       std_error1 <- sqrt(Cov_est[1,1])
       std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle[1]))^(2))
       std_error2 <- sqrt(Cov_est[2,2])
       std_error3 <- sqrt(Cov_est[3,3])

       if( input$conf_level == "95%"){
         conf_val <- qnorm(0.975)
       } else if( input$conf_level == "90%"){
         conf_val <- qnorm(0.95)
       }  else if (input$conf_level == "80%"){
         conf_val<- qnorm(0.90)
       } else if (input$conf_level == "99%"){
         conf_val<- qnorm(0.995)
       }

       MLE <- c(coef_mle[1], exp(coef_mle[1]), coef_mle[2], coef_mle[3])
       SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
       upper <- c(coef_mle[1] + conf_val * std_error1, 
                  # exp(coef_mle[1]) + conf_val * std_error1_delta_method,
                  exp(coef_mle[1] + conf_val * std_error1),
                  coef_mle[2] + conf_val * std_error2, 
                  coef_mle[3] + conf_val * std_error3 )
       lower <- c(coef_mle[1] - conf_val * std_error1, 
                  # exp(coef_mle[1]) - conf_val * std_error1_delta_method,
                  exp(coef_mle[1] - conf_val * std_error1),
                  coef_mle[2] - conf_val * std_error2, 
                  coef_mle[3] - conf_val * std_error3 )


      table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                 MLE = format(MLE, scientific = T),
                                 se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T),
                                 Lower=format(lower, scientific = T))

    },
    bordered=T)
    
    
    output$aic_fre <- renderText({
      value <- inputVar()
      
      kth <- as.numeric(value[1])
      
      data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
      data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]
      
      n_data_kth <- dim(data_kth)[1]
      x_matrix[kth, ] <- (0:100)/100*max( data_kth[,1] )
      initial_values <- data_meta[ kth, c(9,10)]
      location_error_vec <- NA
      colnames(data_kth) <- c("Cw.mol.l", "q.mol.kg")
      
      
      ### 0. Commerical product
      
      
      ### 2. Freundlich model
      ## 2-1 initial value
      Y_trans <-  log(data_kth[-1,2] )
      X_trans <- log(data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      
      res_lm <- data_kth[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      
      
      
      # 2-2. MLE
      initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
      
      
      # 2-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle <- optim(fn = negLogLikelihood_fre_trans,
                         par = initial_value_mle,
                         hessian = T,
                         method = value[5],
                         data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle <- fit_mle$par
        
      } else if(value[4] == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood_fre_trans,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      
      K_F_prime_est_mle_fre <- coef_mle[1]
      K_F_est_mle_fre <- exp(coef_mle[1])
      
      n_est_mle_fre <- coef_mle[2]
      sigma_est_mle <- coef_mle[3]
      
      
      ### Langmuir
      Y_trans <-  (1 / data_kth[-1,2] )
      X_trans <- (1 / data_kth[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vec <- 1 / x_matrix[kth,-1]
      Y_trans_pred <- predict(fit_lm, list(Cw.mol.l = x_vec  ))
      Y_original_pred <- 1/ Y_trans_pred
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * data_kth[-1,1]
      
      res_lm <- data_kth[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_data_kth  )
      
      # 3-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_data_kth -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      # 3-2 fit a new model of reparameterization for Kd
      if(value[4] == "optim"){
        fit_mle_lan <- optim(fn = negLogLikelihood,
                             par = initial_value_mle,
                             hessian = T,
                             method = value[5],
                             data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        coef_mle_lan <- fit_mle_lan$par
        
      } else if(value[4] == "nlm"){
        fit_mle_lan <- nlm(f = negLogLikelihood,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2])
        )
        
        coef_mle_lan <- fit_mle_lan$estimate
        
      }

      
      
      
      AIC_fre <- 2 * length(coef_mle) + negLogLikelihood_fre_trans(coef_mle, data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2]))
      AIC_lan <- 2 * length(coef_mle_lan) + negLogLikelihood(coef_mle_lan, data_x = data.matrix(data_kth[,1]) , data_y = data.matrix(data_kth[,2]))
      if(AIC_lan < AIC_fre){
        AIC_winner <- "Langmuir"
      } else if(AIC_lan > AIC_fre){
        AIC_winner <- "Freundlich"
      } else { AIC_winner <- "even" }
      
      
      paste("AIC of Freundlich model is ", AIC_fre, "\nAIC of Langmuir model is ", AIC_lan, "\nWinner is ", AIC_winner)
      
    })

    
    output$table <- renderDataTable({
        value <- inputVar()

        kth <- as.numeric(value[1])

        data_kth_wNA <- data_set[,(2* kth -1):(2*kth) ]
        data_kth <- data_kth_wNA[!is.na(data_kth_wNA[,1]),]


        datatable(data_kth)
    })
    
    output$user_data <- renderDataTable({
      
      # input$file1 will be NULL initially. After the user selects
      # and uploads a file, head of that data file by default,
      # or all rows if selected, will be shown.
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      req(input$file1)
      
      df <- read.csv(input$file1$datapath,
                     # header = input$header,
                     # sep = input$sep,
                     # quote = input$quote)
                     header = input$header == "Yes",
                     sep = input$sep,
                     quote = input$quote)

        datatable(df)

      
    })
    
    
    
    
    
    output$user_trajectoryPlot <- renderPlot({
      value <- inputVar()
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      req(input$file1)
      
      user_data <- read.csv(input$file1$datapath,
                     # header = input$header,
                     # sep = input$sep,
                     # quote = input$quote)
                     header = input$header == "Yes",
                     sep = input$sep,
                     quote = input$quote)
      
      colnames(user_data) <- c("Cw.mol.l", "q.mol.kg")
      n_user_data <- dim(user_data)[1]
      
      
      if(input$user_model == "Langmuir"){
      
      if( input$zeros == "Yes" ){
        #### Langmuir
        ### 1. initial values
        Y_trans <-  (1 / user_data[-1,2] )
        X_trans <- (1 / user_data[-1,1] )
        
        data_trans <- data.frame(q.mol.kg = Y_trans, 
                                 Cw.mol.l = X_trans)
        
        fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
        coef_lm <- fit_lm$coefficients
        Q_est_lm <- 1 / coef_lm[1]
        K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
        
        x_vals <- (0:100)/100*max( user_data[,1] )

        # Compute mse
        fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
        res_lm <- user_data[-1,2] - 1/fitted_values_lm
        mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
        
        
        ### 2. MLE
        # 2-1 initial values from the previous linear model
        sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
        initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                sqrt(sigma2_initial))
        
        
        # 2-2 fit a new model of reparameterization for Kd
        if(input$user_Rfunction == "optim"){
          fit_mle <- optim(fn = negLogLikelihood,
                           par = initial_value_mle,
                           hessian = T,
                           method = input$user_method,
                           data_x = user_data[,1] , data_y = user_data[,2]
          )
          coef_mle <- fit_mle$par
          
        } else if( input$user_Rfunction == "nlm"){
          fit_mle <- nlm(f = negLogLikelihood,
                         p = initial_value_mle,
                         hessian = T,
                         data_x = user_data[,1] , data_y = user_data[,2]
          )
          
          coef_mle <- fit_mle$estimate
          
        }
        
        # 2-3 get estimates
        Q_est_mle <- coef_mle[1]
        K_est_mle <- exp( coef_mle[2] )
        sigma_est_mle <- coef_mle[3]
        
        # 2-4 prepare for plotting
        
        Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
        res_mle <- user_data[,2] - Y_original_pred_mle
        mse_3 <- (sum(res_mle^(2)) ) / (n_user_data  )
        # compute mse
        
        
        fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_vals) / (1 + exp(coef_mle[2]) * x_vals)
        # fitted values based on mle for plotting 
        
        
        # 2-5 MLE Confidence band 
        if( input$user_conf_level == "95%"){
          conf_val <- qnorm(0.975)
        } else if( input$user_conf_level == "90%"){
          conf_val <- qnorm(0.95)
        }  else if (input$user_conf_level == "80%"){
          conf_val<- qnorm(0.90)
        } else if(input$user_conf_level == "99%"){
          conf_val<- qnorm(0.995)
        }
        
        Cov_est <- solve(fit_mle$hessian)
        c_pred_vec <- x_vals
        
        mle_band_upper <- rep(NA, length(c_pred_vec))
        mle_band_lower <- rep(NA, length(c_pred_vec))
        
        for(i in 1:length(c_pred_vec)){
          c_now <- c_pred_vec[i]
          se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                  sigma = sigma_est_mle, cov_mat=Cov_est)
          
          term <- conf_val * se_est
          mle_band_upper[i] <- fitted_values_mle[i] + term
          mle_band_lower[i] <- fitted_values_mle[i] - term
        }
        
        
        ### 3. Plot
        
        max_final <- max(c(user_data[,2], mle_band_upper ))
        
        # main plot
        plot(x_vals, fitted_values_mle, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
             ylim = c(0, max_final ), lwd = 2, col= 'blue',
             main=paste( "Your Langmuir model") )

        # Confidence band
        lines(c_pred_vec, mle_band_upper, col="blue", lty=3, lwd=2)
        lines(c_pred_vec, mle_band_lower, col="blue", lty=3, lwd=2)
        
        # observed data points
        points(user_data)
        
        # legend
        legend("bottomright",legend=c(
          paste("MLE of Langmuir: Qm=",
                signif(Q_est_mle,6),
                ", Kd=",round(K_est_mle),
                ", MSE=",format(mse_3, scientific = T),sep=""),
          paste( input$user_conf_level, " Confidence band of MLE"),
          paste("Observed data point")
        ),
        bty='n',
        col=c("blue","blue", "black"),lwd=c(2, 2, NA),lty=c(1,2,NA),
        density = c(0, 0, NA), fill = c("blue", "blue", "white"),
        border=c(NA, NA, NA),
        pch = c(NA, NA, 1),
        cex=1
        )
      } else {
        
        ### 1. initial values
        Y_trans <-  (1 / user_data[,2] )
        X_trans <- (1 / user_data[,1] )
        
        data_trans <- data.frame(q.mol.kg = Y_trans, 
                                 Cw.mol.l = X_trans)
        
        fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
        coef_lm <- fit_lm$coefficients
        Q_est_lm <- 1 / coef_lm[1]
        K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
        
        x_vals <- (0:100)/100*max( user_data[,1] )
        
        # Compute mse
        fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[ ,1]
        res_lm <- user_data[ ,2] - 1/fitted_values_lm
        mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
        
        
        ### 2. MLE
        # 2-1 initial values from the previous linear model
        sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
        initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                sqrt(sigma2_initial))
        
        
        # 2-2 fit a new model of reparameterization for Kd
        if( input$user_Rfunction == "optim"){
          fit_mle <- optim(fn = negLogLikelihood,
                           par = initial_value_mle,
                           hessian = T,
                           method = input$user_method,
                           data_x = user_data[,1] , data_y = user_data[,2]
          )
          coef_mle <- fit_mle$par
          
        } else if( input$user_Rfunction == "nlm"){
          fit_mle <- nlm(f = negLogLikelihood,
                         p = initial_value_mle,
                         hessian = T,
                         data_x = user_data[,1] , data_y = user_data[,2]
          )
          
          coef_mle <- fit_mle$estimate
          
        }
        
        # 2-3 get estimates
        Q_est_mle <- coef_mle[1]
        K_est_mle <- exp( coef_mle[2] )
        sigma_est_mle <- coef_mle[3]
        
        # 2-4 prepare for plotting
        
        Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
        res_mle <- user_data[,2] - Y_original_pred_mle
        mse_3 <- (sum(res_mle^(2)) ) / (n_user_data  )
        # compute mse
        
        
        fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_vals) / (1 + exp(coef_mle[2]) * x_vals)
        # fitted values based on mle for plotting 
        
        
        # 2-5 MLE Confidence band 
        if( input$user_conf_level == "95%"){
          conf_val <- qnorm(0.975)
        } else if( input$user_conf_level == "90%"){
          conf_val <- qnorm(0.95)
        }  else if (input$user_conf_level == "80%"){
          conf_val<- qnorm(0.90)
        } else if(input$user_conf_level == "99%"){
          conf_val<- qnorm(0.995)
        }
        
        Cov_est <- solve(fit_mle$hessian)
        c_pred_vec <- x_vals
        
        mle_band_upper <- rep(NA, length(c_pred_vec))
        mle_band_lower <- rep(NA, length(c_pred_vec))
        
        for(i in 1:length(c_pred_vec)){
          c_now <- c_pred_vec[i]
          se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                  sigma = sigma_est_mle, cov_mat=Cov_est)
          
          term <- conf_val * se_est
          mle_band_upper[i] <- fitted_values_mle[i] + term
          mle_band_lower[i] <- fitted_values_mle[i] - term
        }
        
        
        ### 3. Plot
        
        max_final <- max(c(user_data[,2], mle_band_upper ))
        
        # main plot
        plot(x_vals, fitted_values_mle, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
             ylim = c(0, max_final ), lwd = 2, col= 'blue',
             main=paste( "Your Langmuir model") )
        
        # Confidence band
        lines(c_pred_vec, mle_band_upper, col="blue", lty=3, lwd=2)
        lines(c_pred_vec, mle_band_lower, col="blue", lty=3, lwd=2)
        
        # observed data points
        points(user_data)
        
        # legend
        legend("bottomright",legend=c(
          paste("MLE of Langmuir: Qm=",
                signif(Q_est_mle,6),
                ", Kd=",round(K_est_mle),
                ", MSE=",format(mse_3, scientific = T),sep=""),
          paste( input$user_conf_level, " Confidence band of MLE"),
          paste("Observed data point")
        ),
        bty='n',
        col=c("blue","blue", "black"),lwd=c(2, 2, NA),lty=c(1,2,NA),
        density = c(0, 0, NA), fill = c("blue", "blue", "white"),
        border=c(NA, NA, NA),
        pch = c(NA, NA, 1),
        cex=1
        )
      } # end of zeors == "No"
     
      } ## end of Langmuir
      else if( input$user_model == "Freundlich" ) {
        if( input$zeros == "Yes"){
          
          ### Freundlich model
          ## 1. initial values
          Y_trans <-  log(user_data[-1,2] )
          X_trans <- log(user_data[-1,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
          
          K_F_est_lm <- exp( coef_lm[1] )
          n_est_lm <- ( 1 / coef_lm[2] ) 
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          
          res_lm <- user_data[-1,2] - exp(fitted_values_lm)
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          
       
          ## 2. MLE
          initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
          
          
          # 2-1 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                                 par = initial_value_mle,
                                 hessian = T,
                                 method = input$user_method,
                                 data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle_fre <- fit_mle_fre$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                               p = initial_value_mle,
                               hessian = T,
                               data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle_fre <- fit_mle_fre$estimate
            
          }
          
          # 2-3 get estimates
          K_F_prime_est_mle_fre <- coef_mle_fre[1]
          K_F_est_mle_fre <- exp(coef_mle_fre[1])
          
          n_est_mle_fre <- coef_mle_fre[2]
          sigma_est_mle_fre <- coef_mle_fre[3]
          
          # 2-4 prepare for plotting
          Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          # compute mse
          
          
          fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_vals^(1/coef_mle_fre[2])
          # fitted values based on mle for plotting
          
          
          # 2-5 MLE Confidence band
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if(input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          Cov_est_fre <- solve(fit_mle_fre$hessian)
          c_pred_vec <- x_vals
          
          mle_band_upper_fre <- rep(NA, length(c_pred_vec))
          mle_band_lower_fre <- rep(NA, length(c_pred_vec))
          
          for(i in 2:length(c_pred_vec)){
            c_now <- c_pred_vec[i]
            se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)
            
            term <- conf_val * se_est
            mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
            mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
          }
          
          mle_band_upper_fre[1] <- 0
          mle_band_lower_fre[1] <- 0
          # since derivative[2] does not exist when C= 0
          
          ### 4. Plot
          
          max_final <- max(c(user_data[,2], mle_band_upper_fre ))
          
          # main plot
          plot(x_vals, fitted_values_mle_fre, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
               ylim = c(0, max_final ), lwd = 2, col= 'blue',
               main=paste( "Your Freundlich model" ))
          
          
          # Confidence band
          lines(c_pred_vec, mle_band_upper_fre, col="blue", lty=3, lwd=2)
          lines(c_pred_vec, mle_band_lower_fre, col="blue", lty=3, lwd=2)

          # observed data points
          points( user_data )
          
          
          legend("bottomright",legend=c(
            paste("Freundlich : K_F=",
                  signif(K_F_est_mle_fre),
                  ", n=",signif(n_est_mle_fre),
                  ", MSE=", format(mse_fre,scientific = T ),sep=""),
            paste( input$user_conf_level, " Confidence band of Freundlich MLE"),
            paste("Observed data point")
          ),
          bty='n',
          col=c("blue","blue", "black"),lwd=c(2, 2, NA),lty=c(1,3,NA),
          density = c(0, 0, NA), fill = c("blue", "blue", "white"),
          border=c(NA, NA, NA),
          pch = c(NA, NA, 1),
          cex=1
          )
        } else {
          
          ### Freundlich model
          ## 1. initial values
          Y_trans <-  log(user_data[,2] )
          X_trans <- log(user_data[,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[,1]
          
          K_F_est_lm <- exp( coef_lm[1] )
          n_est_lm <- ( 1 / coef_lm[2] ) 
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          
          res_lm <- user_data[,2] - exp(fitted_values_lm)
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          
          
          ## 2. MLE
          initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
          
          
          # 2-1 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                                 par = initial_value_mle,
                                 hessian = T,
                                 method = input$user_method,
                                 data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle_fre <- fit_mle_fre$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                               p = initial_value_mle,
                               hessian = T,
                               data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle_fre <- fit_mle_fre$estimate
            
          }
          
          # 2-3 get estimates
          K_F_prime_est_mle_fre <- coef_mle_fre[1]
          K_F_est_mle_fre <- exp(coef_mle_fre[1])
          
          n_est_mle_fre <- coef_mle_fre[2]
          sigma_est_mle_fre <- coef_mle_fre[3]
          
          # 2-4 prepare for plotting
          Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          # compute mse
          
          
          fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_vals^(1/coef_mle_fre[2])
          # fitted values based on mle for plotting
          
          
          # 2-5 MLE Confidence band
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if(input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          Cov_est_fre <- solve(fit_mle_fre$hessian)
          c_pred_vec <- x_vals
          
          mle_band_upper_fre <- rep(NA, length(c_pred_vec))
          mle_band_lower_fre <- rep(NA, length(c_pred_vec))
          
          for(i in 2:length(c_pred_vec)){
            c_now <- c_pred_vec[i]
            se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)
            
            term <- conf_val * se_est
            mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
            mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
          }
          
          mle_band_upper_fre[1] <- 0
          mle_band_lower_fre[1] <- 0
          # since derivative[2] does not exist when C= 0
          
          ### 4. Plot
          
          max_final <- max(c(user_data[,2], mle_band_upper_fre ))
          
          # main plot
          plot(x_vals, fitted_values_mle_fre, type='l', xlab = 'Cw (mol/l)', ylab = 'q (mol/kg)',
               ylim = c(0, max_final ), lwd = 2, col= 'blue',
               main=paste( "Your Freundlich model" ))
          
          
          # Confidence band
          lines(c_pred_vec, mle_band_upper_fre, col="blue", lty=3, lwd=2)
          lines(c_pred_vec, mle_band_lower_fre, col="blue", lty=3, lwd=2)
          
          # observed data points
          points( user_data )
          
          
          legend("bottomright",legend=c(
            paste("Freundlich : K_F=",
                  signif(K_F_est_mle_fre),
                  ", n=",signif(n_est_mle_fre),
                  ", MSE=", format(mse_fre,scientific = T ),sep=""),
            paste( input$user_conf_level, " Confidence band of Freundlich MLE"),
            paste("Observed data point")
          ),
          bty='n',
          col=c("blue","blue", "black"),lwd=c(2, 2, NA),lty=c(1,3,NA),
          density = c(0, 0, NA), fill = c("blue", "blue", "white"),
          border=c(NA, NA, NA),
          pch = c(NA, NA, 1),
          cex=1
          )
        } #end of zeros == "No"
      } ## end of Freundlich
    })
    
    
    
    
    
    output$mle_table_user <- renderTable({
      value <- inputVar()
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      req(input$file1)
      
      user_data <- read.csv(input$file1$datapath,
                            # header = input$header,
                            # sep = input$sep,
                            # quote = input$quote)
                            header = input$header == "Yes",
                            sep = input$sep,
                            quote = input$quote)
      
      colnames(user_data) <- c("Cw.mol.l", "q.mol.kg")
      n_user_data <- dim(user_data)[1]
      
      
      if(input$user_model == "Langmuir"){
        
        if( input$zeros == "Yes" ){
          #### Langmuir
          ### 1. initial values
          Y_trans <-  (1 / user_data[-1,2] )
          X_trans <- (1 / user_data[-1,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          Q_est_lm <- 1 / coef_lm[1]
          K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
          
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          # Compute mse
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
          res_lm <- user_data[-1,2] - 1/fitted_values_lm
          mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
          
          
          ### 2. MLE
          # 2-1 initial values from the previous linear model
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                  sqrt(sigma2_initial))
          
          
          # 2-2 fit a new model of reparameterization for Kd
          if(input$user_Rfunction == "optim"){
            fit_mle <- optim(fn = negLogLikelihood,
                             par = initial_value_mle,
                             hessian = T,
                             method = input$user_method,
                             data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle <- fit_mle$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle <- nlm(f = negLogLikelihood,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle <- fit_mle$estimate
            
          }
          
          ### MSE and AIC ###
          Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
          res_mle <- user_data[,2] - Y_original_pred_mle
          mse_lan <- (sum(res_mle^(2)) ) / (n_user_data  )
          
          AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x =  user_data[,1]  , data_y =  user_data[,2]  ) 
          ###
          
          
          Q_est_mle <- coef_mle[1]
          K_est_mle <- exp( coef_mle[2] )
          sigma_est_mle <- coef_mle[3]
          
          
          Cov_est <- solve(fit_mle$hessian)
          std_error1 <- sqrt(Cov_est[1,1]  )
          std_error2 <- sqrt(Cov_est[2,2]  )
          std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
          std_error3 <- sqrt(Cov_est[3,3]  )
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
          SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
          upper <- c(coef_mle[1] + conf_val * std_error1,
                     coef_mle[2] + conf_val * std_error2, 
                     exp(coef_mle[2] + conf_val * std_error2),
                     # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                     coef_mle[3] + conf_val * std_error3 )
          lower <- c(coef_mle[1] - conf_val * std_error1,
                     coef_mle[2] - conf_val * std_error2, 
                     exp(coef_mle[2] - conf_val * std_error2),
                     # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                     coef_mle[3] - conf_val * std_error3 )
          
          
          
          table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                     Upper=format(upper, scientific = T), Lower=format(lower, scientific = T),
                                     AIC = c(format(AIC_lan, scientific = T), rep("",3) ),
                                     MSE = c(format(mse_lan, scientific = T), rep("",3 ) )
                                     )
          
          
        } else {
          
          ### 1. initial values
          Y_trans <-  (1 / user_data[,2] )
          X_trans <- (1 / user_data[,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          Q_est_lm <- 1 / coef_lm[1]
          K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
          
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          # Compute mse
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[ ,1]
          res_lm <- user_data[ ,2] - 1/fitted_values_lm
          mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
          
          
          ### 2. MLE
          # 2-1 initial values from the previous linear model
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                  sqrt(sigma2_initial))
          
          
          # 2-2 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle <- optim(fn = negLogLikelihood,
                             par = initial_value_mle,
                             hessian = T,
                             method = input$user_method,
                             data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle <- fit_mle$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle <- nlm(f = negLogLikelihood,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle <- fit_mle$estimate
            
          }
          
          
          
          ### MSE and AIC ###
          Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
          res_mle <- user_data[,2] - Y_original_pred_mle
          mse_lan <- (sum(res_mle^(2)) ) / (n_user_data  )
          
          AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x =  user_data[,1]  , data_y =  user_data[,2]  ) 
          ###
          
          Q_est_mle <- coef_mle[1]
          K_est_mle <- exp( coef_mle[2] )
          sigma_est_mle <- coef_mle[3]
          
          
          Cov_est <- solve(fit_mle$hessian)
          std_error1 <- sqrt(Cov_est[1,1]  )
          std_error2 <- sqrt(Cov_est[2,2]  )
          std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
          std_error3 <- sqrt(Cov_est[3,3]  )
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
          SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
          upper <- c(coef_mle[1] + conf_val * std_error1,
                     coef_mle[2] + conf_val * std_error2, 
                     exp(coef_mle[2] + conf_val * std_error2),
                     # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                     coef_mle[3] + conf_val * std_error3 )
          lower <- c(coef_mle[1] - conf_val * std_error1,
                     coef_mle[2] - conf_val * std_error2, 
                     exp(coef_mle[2] - conf_val * std_error2),
                     # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                     coef_mle[3] - conf_val * std_error3 )
          
          
          table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                     Upper=format(upper, scientific = T), Lower=format(lower, scientific = T),
                                     AIC = c(format(AIC_lan, scientific = T), rep("",3) ),
                                     MSE = c(format(mse_lan, scientific = T), rep("",3 ) )
          )
        } # end of zeors == "No"
        
      } ## end of Langmuir
      else if( input$user_model == "Freundlich" ) {
        if( input$zeros == "Yes"){
          
          ### Freundlich model
          ## 1. initial values
          Y_trans <-  log(user_data[-1,2] )
          X_trans <- log(user_data[-1,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
          
          K_F_est_lm <- exp( coef_lm[1] )
          n_est_lm <- ( 1 / coef_lm[2] ) 
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          
          res_lm <- user_data[-1,2] - exp(fitted_values_lm)
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          
          
          ## 2. MLE
          initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
          
          
          # 2-1 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                                 par = initial_value_mle,
                                 hessian = T,
                                 method = input$user_method,
                                 data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle <- fit_mle_fre$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                               p = initial_value_mle,
                               hessian = T,
                               data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle <- fit_mle_fre$estimate
            
          }
          
          
          ### MSE and AIC ###
          Y_original_pred_mle_fre <-  exp(coef_mle[1]) * user_data[,1]^(1/coef_mle[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          
          AIC_fre <- 2 * length(coef_mle) + negLogLikelihood_fre_trans(coef_mle, data_x = user_data[,1] , data_y = user_data[,2])
          ###

          K_F_prime_est_mle_fre <- coef_mle[1]
          K_F_est_mle_fre <- exp(coef_mle[1])
          
          n_est_mle_fre <- coef_mle[2]
          sigma_est_mle <- coef_mle[3]
          
          
          Cov_est <- solve(fit_mle_fre$hessian)
          std_error1 <- sqrt(Cov_est[1,1])
          std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle[1]))^(2))
          std_error2 <- sqrt(Cov_est[2,2])
          std_error3 <- sqrt(Cov_est[3,3])
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle[1], exp(coef_mle[1]), coef_mle[2], coef_mle[3])
          SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
          upper <- c(coef_mle[1] + conf_val * std_error1, 
                     # exp(coef_mle[1]) + conf_val * std_error1_delta_method,
                     exp(coef_mle[1] + conf_val * std_error1),
                     coef_mle[2] + conf_val * std_error2, 
                     coef_mle[3] + conf_val * std_error3 )
          lower <- c(coef_mle[1] - conf_val * std_error1, 
                     # exp(coef_mle[1]) - conf_val * std_error1_delta_method,
                     exp(coef_mle[1] - conf_val * std_error1),
                     coef_mle[2] - conf_val * std_error2, 
                     coef_mle[3] - conf_val * std_error3 )
          
          
          table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                     MLE = format(MLE, scientific = T),
                                     se = format( SE, scientific = T),
                                     Upper=format(upper, scientific = T),
                                     Lower=format(lower, scientific = T),
                                     AIC = c(format(AIC_fre, scientific = T), rep("",3) ),
                                     MSE = c(format(mse_fre, scientific = T), rep("",3 ) )
                                     
                                     )
        } else {
          
          ### Freundlich model
          ## 1. initial values
          Y_trans <-  log(user_data[,2] )
          X_trans <- log(user_data[,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[,1]
          
          K_F_est_lm <- exp( coef_lm[1] )
          n_est_lm <- ( 1 / coef_lm[2] ) 
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          
          res_lm <- user_data[,2] - exp(fitted_values_lm)
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          
          
          ## 2. MLE
          initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
          
          
          # 2-1 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                                 par = initial_value_mle,
                                 hessian = T,
                                 method = input$user_method,
                                 data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle <- fit_mle_fre$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                               p = initial_value_mle,
                               hessian = T,
                               data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle <- fit_mle_fre$estimate
            
          }
          
          
          
          ### MSE and AIC ###
          Y_original_pred_mle_fre <-  exp(coef_mle[1]) * user_data[,1]^(1/coef_mle[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          
          AIC_fre <- 2 * length(coef_mle) + negLogLikelihood_fre_trans(coef_mle, data_x = user_data[,1] , data_y = user_data[,2])
          ###
          
          K_F_prime_est_mle_fre <- coef_mle[1]
          K_F_est_mle_fre <- exp(coef_mle[1])
          
          n_est_mle_fre <- coef_mle[2]
          sigma_est_mle <- coef_mle[3]
          
          
          Cov_est <- solve(fit_mle_fre$hessian)
          std_error1 <- sqrt(Cov_est[1,1])
          std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle[1]))^(2))
          std_error2 <- sqrt(Cov_est[2,2])
          std_error3 <- sqrt(Cov_est[3,3])
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle[1], exp(coef_mle[1]), coef_mle[2], coef_mle[3])
          SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
          upper <- c(coef_mle[1] + conf_val * std_error1, 
                     # exp(coef_mle[1]) + conf_val * std_error1_delta_method,
                     exp(coef_mle[1] + conf_val * std_error1),
                     coef_mle[2] + conf_val * std_error2, 
                     coef_mle[3] + conf_val * std_error3 )
          lower <- c(coef_mle[1] - conf_val * std_error1, 
                     # exp(coef_mle[1]) - conf_val * std_error1_delta_method,
                     exp(coef_mle[1] - conf_val * std_error1),
                     coef_mle[2] - conf_val * std_error2, 
                     coef_mle[3] - conf_val * std_error3 )
          
          AIC_fre <- 2 * length(coef_mle) + negLogLikelihood_fre_trans(coef_mle, data_x = user_data[,1] , data_y = user_data[,2])

          table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                     MLE = format(MLE, scientific = T),
                                     se = format( SE, scientific = T),
                                     Upper=format(upper, scientific = T),
                                     Lower=format(lower, scientific = T),
                                     AIC = c(format(AIC_fre, scientific = T), rep("",3) ),
                                     MSE = c(format(mse_fre, scientific = T), rep("",3 ) )
                                     
          )
        } #end of zeros == "No"
      } ## end of Freundlich
    },
    bordered = T
    )
    
    

    
    # Table of selected dataset ----
    output$user_down <- renderTable({
      value <- inputVar()
      
      file <- input$file1
      ext <- tools::file_ext(file$datapath)
      validate(need(ext == "csv", "Please upload a csv file"))
      
      req(input$file1)
      
      user_data <- read.csv(input$file1$datapath,
                            # header = input$header,
                            # sep = input$sep,
                            # quote = input$quote)
                            header = input$header == "Yes",
                            sep = input$sep,
                            quote = input$quote)
      
      colnames(user_data) <- c("Cw.mol.l", "q.mol.kg")
      n_user_data <- dim(user_data)[1]
      
      
      if(input$user_model == "Langmuir"){
        
        if( input$zeros == "Yes" ){
          #### Langmuir
          ### 1. initial values
          Y_trans <-  (1 / user_data[-1,2] )
          X_trans <- (1 / user_data[-1,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          Q_est_lm <- 1 / coef_lm[1]
          K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
          
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          # Compute mse
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
          res_lm <- user_data[-1,2] - 1/fitted_values_lm
          mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
          
          
          ### 2. MLE
          # 2-1 initial values from the previous linear model
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                  sqrt(sigma2_initial))
          
          
          # 2-2 fit a new model of reparameterization for Kd
          if(input$user_Rfunction == "optim"){
            fit_mle <- optim(fn = negLogLikelihood,
                             par = initial_value_mle,
                             hessian = T,
                             method = input$user_method,
                             data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle <- fit_mle$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle <- nlm(f = negLogLikelihood,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle <- fit_mle$estimate
            
          }
          
          ### MSE and AIC ###
          Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
          res_mle <- user_data[,2] - Y_original_pred_mle
          mse_lan <- (sum(res_mle^(2)) ) / (n_user_data  )
          
          AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x =  user_data[,1]  , data_y =  user_data[,2]  ) 
          ###
          
          
          # 2-3 get estimates
          Q_est_mle <- coef_mle[1]
          K_est_mle <- exp( coef_mle[2] )
          sigma_est_mle <- coef_mle[3]
          
          # 2-4 prepare for plotting
          
          Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
          res_mle <- user_data[,2] - Y_original_pred_mle
          mse_3 <- (sum(res_mle^(2)) ) / (n_user_data  )
          # compute mse
          
          
          fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_vals) / (1 + exp(coef_mle[2]) * x_vals)
          # fitted values based on mle for plotting 
          
          
          # 2-5 MLE Confidence band 
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if(input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          Cov_est <- solve(fit_mle$hessian)
          c_pred_vec <- x_vals
          
          mle_band_upper <- rep(NA, length(c_pred_vec))
          mle_band_lower <- rep(NA, length(c_pred_vec))
          
          for(i in 1:length(c_pred_vec)){
            c_now <- c_pred_vec[i]
            se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                    sigma = sigma_est_mle, cov_mat=Cov_est)
            
            term <- conf_val * se_est
            mle_band_upper[i] <- fitted_values_mle[i] + term
            mle_band_lower[i] <- fitted_values_mle[i] - term
          }
          
          
          fitted_values <- fitted_values_mle
          upper_bound <- mle_band_upper
          lower_bound <- mle_band_lower
          
          
          trajectory_result <- data.frame("Cw_values"=x_vals,
                                          "fitted_values"=fitted_values,
                                          "upper_bound"=upper_bound,
                                          "lower_bound"=lower_bound)
          
          
          ### MLE Table
          Q_est_mle <- coef_mle[1]
          K_est_mle <- exp( coef_mle[2] )
          sigma_est_mle <- coef_mle[3]
          
          
          Cov_est <- solve(fit_mle$hessian)
          std_error1 <- sqrt(Cov_est[1,1]  )
          std_error2 <- sqrt(Cov_est[2,2]  )
          std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
          std_error3 <- sqrt(Cov_est[3,3]  )
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
          SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
          upper <- c(coef_mle[1] + conf_val * std_error1,
                     coef_mle[2] + conf_val * std_error2, 
                     exp(coef_mle[2] + conf_val * std_error2),
                     # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                     coef_mle[3] + conf_val * std_error3 )
          lower <- c(coef_mle[1] - conf_val * std_error1,
                     coef_mle[2] - conf_val * std_error2, 
                     exp(coef_mle[2] - conf_val * std_error2),
                     # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                     coef_mle[3] - conf_val * std_error3 )
          
          
          
          table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                       Upper=format(upper, scientific = T), Lower=format(lower, scientific = T),
                                       AIC = c(format(AIC_lan, scientific = T), rep("",3) ),
                                       MSE = c(format(mse_lan, scientific = T), rep("",3 ) )
          )
          
        } else {
          
          ### 1. initial values
          Y_trans <-  (1 / user_data[,2] )
          X_trans <- (1 / user_data[,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          Q_est_lm <- 1 / coef_lm[1]
          K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
          
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          # Compute mse
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[ ,1]
          res_lm <- user_data[ ,2] - 1/fitted_values_lm
          mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
          
          
          ### 2. MLE
          # 2-1 initial values from the previous linear model
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                                  sqrt(sigma2_initial))
          
          
          # 2-2 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle <- optim(fn = negLogLikelihood,
                             par = initial_value_mle,
                             hessian = T,
                             method = input$user_method,
                             data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle <- fit_mle$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle <- nlm(f = negLogLikelihood,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle <- fit_mle$estimate
            
          }
          
          ### MSE and AIC ###
          Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
          res_mle <- user_data[,2] - Y_original_pred_mle
          mse_lan <- (sum(res_mle^(2)) ) / (n_user_data  )
          
          AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x =  user_data[,1]  , data_y =  user_data[,2]  ) 
          ###
          
          # 2-3 get estimates
          Q_est_mle <- coef_mle[1]
          K_est_mle <- exp( coef_mle[2] )
          sigma_est_mle <- coef_mle[3]
          
          # 2-4 prepare for plotting
          
          Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
          res_mle <- user_data[,2] - Y_original_pred_mle
          mse_3 <- (sum(res_mle^(2)) ) / (n_user_data  )
          # compute mse
          
          
          fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_vals) / (1 + exp(coef_mle[2]) * x_vals)
          # fitted values based on mle for plotting 
          
          
          # 2-5 MLE Confidence band 
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if(input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          Cov_est <- solve(fit_mle$hessian)
          c_pred_vec <- x_vals
          
          mle_band_upper <- rep(NA, length(c_pred_vec))
          mle_band_lower <- rep(NA, length(c_pred_vec))
          
          for(i in 1:length(c_pred_vec)){
            c_now <- c_pred_vec[i]
            se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                    sigma = sigma_est_mle, cov_mat=Cov_est)
            
            term <- conf_val * se_est
            mle_band_upper[i] <- fitted_values_mle[i] + term
            mle_band_lower[i] <- fitted_values_mle[i] - term
          }
          
          
          fitted_values <- fitted_values_mle
          upper_bound <- mle_band_upper
          lower_bound <- mle_band_lower
          
          trajectory_result <- data.frame("Cw_values"=x_vals,
                                          "fitted_values"=fitted_values,
                                          "upper_bound"=upper_bound,
                                          "lower_bound"=lower_bound)
          
          
          ### MLE Table
          Q_est_mle <- coef_mle[1]
          K_est_mle <- exp( coef_mle[2] )
          sigma_est_mle <- coef_mle[3]
          
          
          Cov_est <- solve(fit_mle$hessian)
          std_error1 <- sqrt(Cov_est[1,1]  )
          std_error2 <- sqrt(Cov_est[2,2]  )
          std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
          std_error3 <- sqrt(Cov_est[3,3]  )
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
          SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
          upper <- c(coef_mle[1] + conf_val * std_error1,
                     coef_mle[2] + conf_val * std_error2, 
                     exp(coef_mle[2] + conf_val * std_error2),
                     # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                     coef_mle[3] + conf_val * std_error3 )
          lower <- c(coef_mle[1] - conf_val * std_error1,
                     coef_mle[2] - conf_val * std_error2, 
                     exp(coef_mle[2] - conf_val * std_error2),
                     # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                     coef_mle[3] - conf_val * std_error3 )
          
          
          
          table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                      Upper=format(upper, scientific = T), Lower=format(lower, scientific = T),
                                      AIC = c(format(AIC_lan, scientific = T), rep("",3) ),
                                      MSE = c(format(mse_lan, scientific = T), rep("",3 ) )
          )
        } # end of zeors == "No"
        
      } ## end of Langmuir
      else if( input$user_model == "Freundlich" ) {
        if( input$zeros == "Yes"){
          
          ### Freundlich model
          ## 1. initial values
          Y_trans <-  log(user_data[-1,2] )
          X_trans <- log(user_data[-1,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
          
          K_F_est_lm <- exp( coef_lm[1] )
          n_est_lm <- ( 1 / coef_lm[2] ) 
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          
          res_lm <- user_data[-1,2] - exp(fitted_values_lm)
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          
          
          ## 2. MLE
          initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
          
          
          # 2-1 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                                 par = initial_value_mle,
                                 hessian = T,
                                 method = input$user_method,
                                 data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle_fre <- fit_mle_fre$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                               p = initial_value_mle,
                               hessian = T,
                               data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle_fre <- fit_mle_fre$estimate
            
          }
          
          ### MSE and AIC ###
          Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          
          AIC_fre <- 2 * length(coef_mle_fre) + negLogLikelihood_fre_trans(coef_mle_fre, data_x = user_data[,1] , data_y = user_data[,2])
          ###
          
          # 2-3 get estimates
          K_F_prime_est_mle_fre <- coef_mle_fre[1]
          K_F_est_mle_fre <- exp(coef_mle_fre[1])
          
          n_est_mle_fre <- coef_mle_fre[2]
          sigma_est_mle_fre <- coef_mle_fre[3]
          
          # 2-4 prepare for plotting
          Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          # compute mse
          
          
          fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_vals^(1/coef_mle_fre[2])
          # fitted values based on mle for plotting
          
          
          # 2-5 MLE Confidence band
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if(input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          Cov_est_fre <- solve(fit_mle_fre$hessian)
          c_pred_vec <- x_vals
          
          mle_band_upper_fre <- rep(NA, length(c_pred_vec))
          mle_band_lower_fre <- rep(NA, length(c_pred_vec))
          
          for(i in 2:length(c_pred_vec)){
            c_now <- c_pred_vec[i]
            se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)
            
            term <- conf_val * se_est
            mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
            mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
          }
          
          mle_band_upper_fre[1] <- 0
          mle_band_lower_fre[1] <- 0
          # since derivative[2] does not exist when C= 0
          
          fitted_values <- fitted_values_mle_fre
          upper_bound <- mle_band_upper_fre
          lower_bound <- mle_band_lower_fre
          
          
          trajectory_result <- data.frame("Cw_values"=x_vals,
                                          "fitted_values"=fitted_values,
                                          "upper_bound"=upper_bound,
                                          "lower_bound"=lower_bound)
          
          #### MLE Table
          K_F_prime_est_mle_fre <- coef_mle_fre[1]
          K_F_est_mle_fre <- exp(coef_mle_fre[1])
          
          n_est_mle_fre <- coef_mle_fre[2]
          sigma_est_mle <- coef_mle_fre[3]
          
          
          Cov_est <- solve(fit_mle_fre$hessian)
          std_error1 <- sqrt(Cov_est[1,1])
          std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle_fre[1]))^(2))
          std_error2 <- sqrt(Cov_est[2,2])
          std_error3 <- sqrt(Cov_est[3,3])
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle_fre[1], exp(coef_mle_fre[1]), coef_mle_fre[2], coef_mle_fre[3])
          SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
          upper <- c(coef_mle_fre[1] + conf_val * std_error1, 
                     # exp(coef_mle_fre[1]) + conf_val * std_error1_delta_method,
                     exp(coef_mle_fre[1] + conf_val * std_error1),
                     coef_mle_fre[2] + conf_val * std_error2, 
                     coef_mle_fre[3] + conf_val * std_error3 )
          lower <- c(coef_mle_fre[1] - conf_val * std_error1, 
                     # exp(coef_mle_fre[1]) - conf_val * std_error1_delta_method,
                     exp(coef_mle_fre[1] - conf_val * std_error1),
                     coef_mle_fre[2] - conf_val * std_error2, 
                     coef_mle_fre[3] - conf_val * std_error3 )
          
          
          table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                     MLE = format(MLE, scientific = T),
                                     se = format( SE, scientific = T),
                                     Upper=format(upper, scientific = T),
                                     Lower=format(lower, scientific = T),
                                     AIC = c(format(AIC_fre, scientific = T), rep("",3) ),
                                     MSE = c(format(mse_fre, scientific = T), rep("",3 ) )
                                     
          )
        } else {
          
          ### Freundlich model
          ## 1. initial values
          Y_trans <-  log(user_data[,2] )
          X_trans <- log(user_data[,1] )
          
          data_trans <- data.frame(q.mol.kg = Y_trans, 
                                   Cw.mol.l = X_trans)
          
          fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
          coef_lm <- fit_lm$coefficients
          fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[,1]
          
          K_F_est_lm <- exp( coef_lm[1] )
          n_est_lm <- ( 1 / coef_lm[2] ) 
          x_vals <- (0:100)/100*max( user_data[,1] )
          
          
          res_lm <- user_data[,2] - exp(fitted_values_lm)
          sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
          
          
          ## 2. MLE
          initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
          
          
          # 2-1 fit a new model of reparameterization for Kd
          if( input$user_Rfunction == "optim"){
            fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                                 par = initial_value_mle,
                                 hessian = T,
                                 method = input$user_method,
                                 data_x = user_data[,1] , data_y = user_data[,2]
            )
            coef_mle_fre <- fit_mle_fre$par
            
          } else if( input$user_Rfunction == "nlm"){
            fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                               p = initial_value_mle,
                               hessian = T,
                               data_x = user_data[,1] , data_y = user_data[,2]
            )
            
            coef_mle_fre <- fit_mle_fre$estimate
            
          }
          
          ### MSE and AIC ###
          Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          
          AIC_fre <- 2 * length(coef_mle_fre) + negLogLikelihood_fre_trans(coef_mle_fre, data_x = user_data[,1] , data_y = user_data[,2])
          ###
          
          # 2-3 get estimates
          K_F_prime_est_mle_fre <- coef_mle_fre[1]
          K_F_est_mle_fre <- exp(coef_mle_fre[1])
          
          n_est_mle_fre <- coef_mle_fre[2]
          sigma_est_mle_fre <- coef_mle_fre[3]
          
          # 2-4 prepare for plotting
          Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
          res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
          mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
          # compute mse
          
          
          fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_vals^(1/coef_mle_fre[2])
          # fitted values based on mle for plotting
          
          
          # 2-5 MLE Confidence band
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if(input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          Cov_est_fre <- solve(fit_mle_fre$hessian)
          c_pred_vec <- x_vals
          
          mle_band_upper_fre <- rep(NA, length(c_pred_vec))
          mle_band_lower_fre <- rep(NA, length(c_pred_vec))
          
          for(i in 2:length(c_pred_vec)){
            c_now <- c_pred_vec[i]
            se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)
            
            term <- conf_val * se_est
            mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
            mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
          }
          
          mle_band_upper_fre[1] <- 0
          mle_band_lower_fre[1] <- 0
          # since derivative[2] does not exist when C= 0
          
          fitted_values <- fitted_values_mle_fre
          upper_bound <- mle_band_upper_fre
          lower_bound <- mle_band_lower_fre
          
          trajectory_result <- data.frame("Cw_values"=x_vals,
                                          "fitted_values"=fitted_values,
                                          "upper_bound"=upper_bound,
                                          "lower_bound"=lower_bound)
          
          #### MLE Table
          K_F_prime_est_mle_fre <- coef_mle_fre[1]
          K_F_est_mle_fre <- exp(coef_mle_fre[1])
          
          n_est_mle_fre <- coef_mle_fre[2]
          sigma_est_mle <- coef_mle_fre[3]
          
          
          Cov_est <- solve(fit_mle_fre$hessian)
          std_error1 <- sqrt(Cov_est[1,1])
          std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle_fre[1]))^(2))
          std_error2 <- sqrt(Cov_est[2,2])
          std_error3 <- sqrt(Cov_est[3,3])
          
          if( input$user_conf_level == "95%"){
            conf_val <- qnorm(0.975)
          } else if( input$user_conf_level == "90%"){
            conf_val <- qnorm(0.95)
          }  else if (input$user_conf_level == "80%"){
            conf_val<- qnorm(0.90)
          } else if (input$user_conf_level == "99%"){
            conf_val<- qnorm(0.995)
          }
          
          MLE <- c(coef_mle_fre[1], exp(coef_mle_fre[1]), coef_mle_fre[2], coef_mle_fre[3])
          SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
          upper <- c(coef_mle_fre[1] + conf_val * std_error1, 
                     # exp(coef_mle_fre[1]) + conf_val * std_error1_delta_method,
                     exp(coef_mle_fre[1] + conf_val * std_error1),
                     coef_mle_fre[2] + conf_val * std_error2, 
                     coef_mle_fre[3] + conf_val * std_error3 )
          lower <- c(coef_mle_fre[1] - conf_val * std_error1, 
                     # exp(coef_mle_fre[1]) - conf_val * std_error1_delta_method,
                     exp(coef_mle_fre[1] - conf_val * std_error1),
                     coef_mle_fre[2] - conf_val * std_error2, 
                     coef_mle_fre[3] - conf_val * std_error3 )
          
          
          table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                     MLE = format(MLE, scientific = T),
                                     se = format( SE, scientific = T),
                                     Upper=format(upper, scientific = T),
                                     Lower=format(lower, scientific = T),
                                     AIC = c(format(AIC_fre, scientific = T), rep("",3) ),
                                     MSE = c(format(mse_fre, scientific = T), rep("",3 ) )
                                     
          )
          
          
        } #end of zeros == "No"
      } ## end of Freundlich
      
      
      
      
      
      if( input$down_dataset == "Summary" ){
        data <- table_result
      } else {
        data <- trajectory_result
      }
      
      data
    }, bordered = T)

    
    
    
    # Downloadable csv of selected dataset ----

data_down_reactive <- reactive({
  value <- inputVar()
  
  file <- input$file1
  ext <- tools::file_ext(file$datapath)
  validate(need(ext == "csv", "Please upload a csv file"))
  
  req(input$file1)
  
  user_data <- read.csv(input$file1$datapath,
                        # header = input$header,
                        # sep = input$sep,
                        # quote = input$quote)
                        header = input$header == "Yes",
                        sep = input$sep,
                        quote = input$quote)
  
  colnames(user_data) <- c("Cw.mol.l", "q.mol.kg")
  n_user_data <- dim(user_data)[1]
  
  
  if(input$user_model == "Langmuir"){
    
    if( input$zeros == "Yes" ){
      #### Langmuir
      ### 1. initial values
      Y_trans <-  (1 / user_data[-1,2] )
      X_trans <- (1 / user_data[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vals <- (0:100)/100*max( user_data[,1] )
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
      res_lm <- user_data[-1,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
      
      
      ### 2. MLE
      # 2-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      # 2-2 fit a new model of reparameterization for Kd
      if(input$user_Rfunction == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = input$user_method,
                         data_x = user_data[,1] , data_y = user_data[,2]
        )
        coef_mle <- fit_mle$par
        
      } else if( input$user_Rfunction == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = user_data[,1] , data_y = user_data[,2]
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      ### MSE and AIC ###
      Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
      res_mle <- user_data[,2] - Y_original_pred_mle
      mse_lan <- (sum(res_mle^(2)) ) / (n_user_data  )
      
      AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x =  user_data[,1]  , data_y =  user_data[,2]  ) 
      ###
      
      # 2-3 get estimates
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]
      
      # 2-4 prepare for plotting
      
      Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
      res_mle <- user_data[,2] - Y_original_pred_mle
      mse_3 <- (sum(res_mle^(2)) ) / (n_user_data  )
      # compute mse
      
      
      fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_vals) / (1 + exp(coef_mle[2]) * x_vals)
      # fitted values based on mle for plotting 
      
      
      # 2-5 MLE Confidence band 
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est <- solve(fit_mle$hessian)
      c_pred_vec <- x_vals
      
      mle_band_upper <- rep(NA, length(c_pred_vec))
      mle_band_lower <- rep(NA, length(c_pred_vec))
      
      for(i in 1:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                sigma = sigma_est_mle, cov_mat=Cov_est)
        
        term <- conf_val * se_est
        mle_band_upper[i] <- fitted_values_mle[i] + term
        mle_band_lower[i] <- fitted_values_mle[i] - term
      }
      
      
      fitted_values <- fitted_values_mle
      upper_bound <- mle_band_upper
      lower_bound <- mle_band_lower
      
      
      trajectory_result <- data.frame("Cw_values"=x_vals,
                                      "fitted_values"=fitted_values,
                                      "upper_bound"=upper_bound,
                                      "lower_bound"=lower_bound)
      
      
      ### MLE Table
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]
      
      
      Cov_est <- solve(fit_mle$hessian)
      std_error1 <- sqrt(Cov_est[1,1]  )
      std_error2 <- sqrt(Cov_est[2,2]  )
      std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
      std_error3 <- sqrt(Cov_est[3,3]  )
      
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if (input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
      SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
      upper <- c(coef_mle[1] + conf_val * std_error1,
                 coef_mle[2] + conf_val * std_error2, 
                 exp(coef_mle[2] + conf_val * std_error2),
                 # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                 coef_mle[3] + conf_val * std_error3 )
      lower <- c(coef_mle[1] - conf_val * std_error1,
                 coef_mle[2] - conf_val * std_error2, 
                 exp(coef_mle[2] - conf_val * std_error2),
                 # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                 coef_mle[3] - conf_val * std_error3 )
      
      
      
      table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T), Lower=format(lower, scientific = T),
                                 AIC = c(format(AIC_lan, scientific = T), rep("",3) ),
                                 MSE = c(format(mse_lan, scientific = T), rep("",3 ) )
      )
      
      
      
    } else {
      
      ### 1. initial values
      Y_trans <-  (1 / user_data[,2] )
      X_trans <- (1 / user_data[,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      Q_est_lm <- 1 / coef_lm[1]
      K_est_lm <- ( coef_lm[1] / coef_lm[2] ) 
      
      x_vals <- (0:100)/100*max( user_data[,1] )
      
      # Compute mse
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[ ,1]
      res_lm <- user_data[ ,2] - 1/fitted_values_lm
      mse_2 <- (sum(res_lm^(2)) ) / (n_user_data  )
      
      
      ### 2. MLE
      # 2-1 initial values from the previous linear model
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
      initial_value_mle <- c( 1/coef_lm[1], log( coef_lm[1] / coef_lm[2]), 
                              sqrt(sigma2_initial))
      
      
      # 2-2 fit a new model of reparameterization for Kd
      if( input$user_Rfunction == "optim"){
        fit_mle <- optim(fn = negLogLikelihood,
                         par = initial_value_mle,
                         hessian = T,
                         method = input$user_method,
                         data_x = user_data[,1] , data_y = user_data[,2]
        )
        coef_mle <- fit_mle$par
        
      } else if( input$user_Rfunction == "nlm"){
        fit_mle <- nlm(f = negLogLikelihood,
                       p = initial_value_mle,
                       hessian = T,
                       data_x = user_data[,1] , data_y = user_data[,2]
        )
        
        coef_mle <- fit_mle$estimate
        
      }
      
      
      ### MSE and AIC ###
      Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
      res_mle <- user_data[,2] - Y_original_pred_mle
      mse_lan <- (sum(res_mle^(2)) ) / (n_user_data  )
      
      AIC_lan <- 2 * length(coef_mle) + negLogLikelihood(coef_mle, data_x =  user_data[,1]  , data_y =  user_data[,2]  ) 
      ###
      
      # 2-3 get estimates
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]
      
      # 2-4 prepare for plotting
      
      Y_original_pred_mle <- ( coef_mle[1] * exp(coef_mle[2]) * user_data[,1] ) / (1 + exp(coef_mle[2]) * user_data[,1])
      res_mle <- user_data[,2] - Y_original_pred_mle
      mse_3 <- (sum(res_mle^(2)) ) / (n_user_data  )
      # compute mse
      
      
      fitted_values_mle <-  ( coef_mle[1] * exp(coef_mle[2]) *  x_vals) / (1 + exp(coef_mle[2]) * x_vals)
      # fitted values based on mle for plotting 
      
      
      # 2-5 MLE Confidence band 
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est <- solve(fit_mle$hessian)
      c_pred_vec <- x_vals
      
      mle_band_upper <- rep(NA, length(c_pred_vec))
      mle_band_lower <- rep(NA, length(c_pred_vec))
      
      for(i in 1:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle(Q = Q_est_mle, K_new=coef_mle[2], C=c_now,
                                sigma = sigma_est_mle, cov_mat=Cov_est)
        
        term <- conf_val * se_est
        mle_band_upper[i] <- fitted_values_mle[i] + term
        mle_band_lower[i] <- fitted_values_mle[i] - term
      }
      
      
      fitted_values <- fitted_values_mle
      upper_bound <- mle_band_upper
      lower_bound <- mle_band_lower
      
      trajectory_result <- data.frame("Cw_values"=x_vals,
                                      "fitted_values"=fitted_values,
                                      "upper_bound"=upper_bound,
                                      "lower_bound"=lower_bound)
      
      
      ### MLE Table
      Q_est_mle <- coef_mle[1]
      K_est_mle <- exp( coef_mle[2] )
      sigma_est_mle <- coef_mle[3]
      
      
      Cov_est <- solve(fit_mle$hessian)
      std_error1 <- sqrt(Cov_est[1,1]  )
      std_error2 <- sqrt(Cov_est[2,2]  )
      std_error2_delta_method <- sqrt(Cov_est[2,2] * (exp(coef_mle[2]))^(2) )
      std_error3 <- sqrt(Cov_est[3,3]  )
      
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if (input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      MLE <- c(coef_mle[1], coef_mle[2], exp(coef_mle[2]), coef_mle[3])
      SE <- c(std_error1, std_error2, std_error2_delta_method, std_error3)
      upper <- c(coef_mle[1] + conf_val * std_error1,
                 coef_mle[2] + conf_val * std_error2, 
                 exp(coef_mle[2] + conf_val * std_error2),
                 # exp(coef_mle[2]) + conf_val * std_error2_delta_method,
                 coef_mle[3] + conf_val * std_error3 )
      lower <- c(coef_mle[1] - conf_val * std_error1,
                 coef_mle[2] - conf_val * std_error2, 
                 exp(coef_mle[2] - conf_val * std_error2),
                 # exp(coef_mle[2]) - conf_val * std_error2_delta_method, 
                 coef_mle[3] - conf_val * std_error3 )
      
      
      
      table_result <- data.frame(Par = c("Qmax", "K'd", "Kd", "sigma"),  MLE = format(MLE, scientific = T), se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T), Lower=format(lower, scientific = T),
                                 AIC = c(format(AIC_lan, scientific = T), rep("",3) ),
                                 MSE = c(format(mse_lan, scientific = T), rep("",3 ) )
      )
    } # end of zeors == "No"
    
  } ## end of Langmuir
  else if( input$user_model == "Freundlich" ) {
    if( input$zeros == "Yes"){
      
      ### Freundlich model
      ## 1. initial values
      Y_trans <-  log(user_data[-1,2] )
      X_trans <- log(user_data[-1,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[-1,1]
      
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      x_vals <- (0:100)/100*max( user_data[,1] )
      
      
      res_lm <- user_data[-1,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
      
      
      ## 2. MLE
      initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
      
      
      # 2-1 fit a new model of reparameterization for Kd
      if( input$user_Rfunction == "optim"){
        fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                             par = initial_value_mle,
                             hessian = T,
                             method = input$user_method,
                             data_x = user_data[,1] , data_y = user_data[,2]
        )
        coef_mle_fre <- fit_mle_fre$par
        
      } else if( input$user_Rfunction == "nlm"){
        fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = user_data[,1] , data_y = user_data[,2]
        )
        
        coef_mle_fre <- fit_mle_fre$estimate
        
      }
      
      ### MSE and AIC ###
      Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
      res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
      mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
      
      AIC_fre <- 2 * length(coef_mle_fre) + negLogLikelihood_fre_trans(coef_mle_fre, data_x = user_data[,1] , data_y = user_data[,2])
      ###
      
      # 2-3 get estimates
      K_F_prime_est_mle_fre <- coef_mle_fre[1]
      K_F_est_mle_fre <- exp(coef_mle_fre[1])
      
      n_est_mle_fre <- coef_mle_fre[2]
      sigma_est_mle_fre <- coef_mle_fre[3]
      
      # 2-4 prepare for plotting
      Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
      res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
      mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
      # compute mse
      
      
      fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_vals^(1/coef_mle_fre[2])
      # fitted values based on mle for plotting
      
      
      # 2-5 MLE Confidence band
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est_fre <- solve(fit_mle_fre$hessian)
      c_pred_vec <- x_vals
      
      mle_band_upper_fre <- rep(NA, length(c_pred_vec))
      mle_band_lower_fre <- rep(NA, length(c_pred_vec))
      
      for(i in 2:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)
        
        term <- conf_val * se_est
        mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
        mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
      }
      
      mle_band_upper_fre[1] <- 0
      mle_band_lower_fre[1] <- 0
      # since derivative[2] does not exist when C= 0
      
      fitted_values <- fitted_values_mle_fre
      upper_bound <- mle_band_upper_fre
      lower_bound <- mle_band_lower_fre
      
      
      trajectory_result <- data.frame("Cw_values"=x_vals,
                                      "fitted_values"=fitted_values,
                                      "upper_bound"=upper_bound,
                                      "lower_bound"=lower_bound)
      
      #### MLE Table
      K_F_prime_est_mle_fre <- coef_mle_fre[1]
      K_F_est_mle_fre <- exp(coef_mle_fre[1])
      
      n_est_mle_fre <- coef_mle_fre[2]
      sigma_est_mle <- coef_mle_fre[3]
      
      
      Cov_est <- solve(fit_mle_fre$hessian)
      std_error1 <- sqrt(Cov_est[1,1])
      std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle_fre[1]))^(2))
      std_error2 <- sqrt(Cov_est[2,2])
      std_error3 <- sqrt(Cov_est[3,3])
      
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if (input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      MLE <- c(coef_mle_fre[1], exp(coef_mle_fre[1]), coef_mle_fre[2], coef_mle_fre[3])
      SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
      upper <- c(coef_mle_fre[1] + conf_val * std_error1, 
                 # exp(coef_mle_fre[1]) + conf_val * std_error1_delta_method,
                 exp(coef_mle_fre[1] + conf_val * std_error1),
                 coef_mle_fre[2] + conf_val * std_error2, 
                 coef_mle_fre[3] + conf_val * std_error3 )
      lower <- c(coef_mle_fre[1] - conf_val * std_error1, 
                 # exp(coef_mle_fre[1]) - conf_val * std_error1_delta_method,
                 exp(coef_mle_fre[1] - conf_val * std_error1),
                 coef_mle_fre[2] - conf_val * std_error2, 
                 coef_mle_fre[3] - conf_val * std_error3 )
      
      
      table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                 MLE = format(MLE, scientific = T),
                                 se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T),
                                 Lower=format(lower, scientific = T),
                                 AIC = c(format(AIC_fre, scientific = T), rep("",3) ),
                                 MSE = c(format(mse_fre, scientific = T), rep("",3 ) )
                                 
      )
    } else {
      
      ### Freundlich model
      ## 1. initial values
      Y_trans <-  log(user_data[,2] )
      X_trans <- log(user_data[,1] )
      
      data_trans <- data.frame(q.mol.kg = Y_trans, 
                               Cw.mol.l = X_trans)
      
      fit_lm <- lm(q.mol.kg ~ . , data = data_trans)
      coef_lm <- fit_lm$coefficients
      fitted_values_lm <- coef_lm[1] + coef_lm[2] * user_data[,1]
      
      K_F_est_lm <- exp( coef_lm[1] )
      n_est_lm <- ( 1 / coef_lm[2] ) 
      x_vals <- (0:100)/100*max( user_data[,1] )
      
      
      res_lm <- user_data[,2] - exp(fitted_values_lm)
      sigma2_initial <- (sum(res_lm^(2)) ) / (n_user_data -1  )
      
      
      ## 2. MLE
      initial_value_mle <- c( log(K_F_est_lm), n_est_lm, sqrt(sigma2_initial))
      
      
      # 2-1 fit a new model of reparameterization for Kd
      if( input$user_Rfunction == "optim"){
        fit_mle_fre <- optim(fn = negLogLikelihood_fre_trans,
                             par = initial_value_mle,
                             hessian = T,
                             method = input$user_method,
                             data_x = user_data[,1] , data_y = user_data[,2]
        )
        coef_mle_fre <- fit_mle_fre$par
        
      } else if( input$user_Rfunction == "nlm"){
        fit_mle_fre <- nlm(f = negLogLikelihood_fre_trans,
                           p = initial_value_mle,
                           hessian = T,
                           data_x = user_data[,1] , data_y = user_data[,2]
        )
        
        coef_mle_fre <- fit_mle_fre$estimate
        
      }
      
      ### MSE and AIC ###
      Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
      res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
      mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
      
      AIC_fre <- 2 * length(coef_mle_fre) + negLogLikelihood_fre_trans(coef_mle_fre, data_x = user_data[,1] , data_y = user_data[,2])
      ###
      
      # 2-3 get estimates
      K_F_prime_est_mle_fre <- coef_mle_fre[1]
      K_F_est_mle_fre <- exp(coef_mle_fre[1])
      
      n_est_mle_fre <- coef_mle_fre[2]
      sigma_est_mle_fre <- coef_mle_fre[3]
      
      # 2-4 prepare for plotting
      Y_original_pred_mle_fre <-  exp(coef_mle_fre[1]) * user_data[,1]^(1/coef_mle_fre[2])
      res_mle_fre <- user_data[,2] - Y_original_pred_mle_fre
      mse_fre <- (sum(res_mle_fre^(2)) ) / (n_user_data  )
      # compute mse
      
      
      fitted_values_mle_fre <-   exp(coef_mle_fre[1]) *   x_vals^(1/coef_mle_fre[2])
      # fitted values based on mle for plotting
      
      
      # 2-5 MLE Confidence band
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if(input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      Cov_est_fre <- solve(fit_mle_fre$hessian)
      c_pred_vec <- x_vals
      
      mle_band_upper_fre <- rep(NA, length(c_pred_vec))
      mle_band_lower_fre <- rep(NA, length(c_pred_vec))
      
      for(i in 2:length(c_pred_vec)){
        c_now <- c_pred_vec[i]
        se_est <- se_for_CI_mle_fre_trans( K_F_prime = coef_mle_fre[1], n=coef_mle_fre[2], C=c_now, sigma = sigma_est_mle_fre, cov_mat=Cov_est_fre)
        
        term <- conf_val * se_est
        mle_band_upper_fre[i] <- fitted_values_mle_fre[i] + term
        mle_band_lower_fre[i] <- fitted_values_mle_fre[i] - term
      }
      
      mle_band_upper_fre[1] <- 0
      mle_band_lower_fre[1] <- 0
      # since derivative[2] does not exist when C= 0
      
      fitted_values <- fitted_values_mle_fre
      upper_bound <- mle_band_upper_fre
      lower_bound <- mle_band_lower_fre
      
      trajectory_result <- data.frame("Cw_values"=x_vals,
                                      "fitted_values"=fitted_values,
                                      "upper_bound"=upper_bound,
                                      "lower_bound"=lower_bound)
      
      #### MLE Table
      K_F_prime_est_mle_fre <- coef_mle_fre[1]
      K_F_est_mle_fre <- exp(coef_mle_fre[1])
      
      n_est_mle_fre <- coef_mle_fre[2]
      sigma_est_mle <- coef_mle_fre[3]
      
      
      Cov_est <- solve(fit_mle_fre$hessian)
      std_error1 <- sqrt(Cov_est[1,1])
      std_error1_delta_method <- sqrt(Cov_est[1,1] * (exp(coef_mle_fre[1]))^(2))
      std_error2 <- sqrt(Cov_est[2,2])
      std_error3 <- sqrt(Cov_est[3,3])
      
      if( input$user_conf_level == "95%"){
        conf_val <- qnorm(0.975)
      } else if( input$user_conf_level == "90%"){
        conf_val <- qnorm(0.95)
      }  else if (input$user_conf_level == "80%"){
        conf_val<- qnorm(0.90)
      } else if (input$user_conf_level == "99%"){
        conf_val<- qnorm(0.995)
      }
      
      MLE <- c(coef_mle_fre[1], exp(coef_mle_fre[1]), coef_mle_fre[2], coef_mle_fre[3])
      SE <- c(std_error1, std_error1_delta_method,std_error2, std_error3)
      upper <- c(coef_mle_fre[1] + conf_val * std_error1, 
                 # exp(coef_mle_fre[1]) + conf_val * std_error1_delta_method,
                 exp(coef_mle_fre[1] + conf_val * std_error1),
                 coef_mle_fre[2] + conf_val * std_error2, 
                 coef_mle_fre[3] + conf_val * std_error3 )
      lower <- c(coef_mle_fre[1] - conf_val * std_error1, 
                 # exp(coef_mle_fre[1]) - conf_val * std_error1_delta_method,
                 exp(coef_mle_fre[1] - conf_val * std_error1),
                 coef_mle_fre[2] - conf_val * std_error2, 
                 coef_mle_fre[3] - conf_val * std_error3 )
      
      
      table_result <- data.frame(Par = c("K'_F", "K_F", "n", "sigma"),
                                 MLE = format(MLE, scientific = T),
                                 se = format( SE, scientific = T),
                                 Upper=format(upper, scientific = T),
                                 Lower=format(lower, scientific = T),
                                 AIC = c(format(AIC_fre, scientific = T), rep("",3) ),
                                 MSE = c(format(mse_fre, scientific = T), rep("",3 ) )
                                 
      )
      
      
    } #end of zeros == "No"
  } ## end of Freundlich
  
  
  
  
  
  if( input$down_dataset == "Summary" ){
    data <- table_result
  } else {
    data <- trajectory_result
  }
  
  data
})
    
    output$downloadData <- downloadHandler(
      
      filename = function() {
        paste(input$user_model,"_" , input$down_dataset,".csv", sep = "")
      },
      content = function(file) {
        write.csv(data_down_reactive(), file, row.names = FALSE)
      }
    )
    
    
    output$downloadData_ex <- downloadHandler(
      
      filename = "example_data.csv",
      content = function(file) {
        file.copy('www/example_data.csv', file)
      }
    )


} # end of server

# Run the application 
shinyApp(ui = ui, server = server)
