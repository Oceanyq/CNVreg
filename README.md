
# title: "An Introduction to `CNVreg` Package to Perform Copy Number Variant Association Analysis with Penalied Regression"


## Introduction

The *CNVreg* package provides functions to perform copy number variants (CNV) association analysis with penalized regression model. 

This package convert CNVs over a genomic region as a piecewise constant curve to capture the dosage and length of CNVs. The association analysis is then evaluated by regressing outcome traits on all CNV fragments in the region while adjusting for covariates. The corresponding CNV effects are obtained at each genome position. The penalized regression model with Lasso and weighted fusion penalties would perform variable selection and encourage adjacent CNVs to share similar effect size.

This package has 3 main functions: 

* `prep()`: Data preprocessing and format conversion.

* `cvfit_WTSMTH()`: Model fitting and effect estimate with cross-validation(CV). The CV procedure is to tune an optimal model by selecting the best pair of candidate tuning parameters.

* `fit_WTSMTH()`: Model fitting and effect estimate with a given pair of tuning parameters. 


All functions use an example data included in the `CNVreg` package.


## Example data: CNVCOVY

The `CNVCOVY` dataset included in this package contains a small sample of data for demonstration purposes. It has 4 separate data files: copy number variants data in `CNV`, covariate data in `Cov`, and outcome traits `Y_QT` (quantitative) and `Y_BT`(Binary). 



 * `CNV`: A data frame describing CNV data in PLINK format with 5 variables `ID`, `CHR`, `BP1`, `BP2`, and `TYPE`.
 * `Cov`: A data frame contains 3 variables: `ID`, `Sex`, and `Age`.
 * `Y_QT` and `Y_BT`: each is a data frame for outcomes traits. `Y_QT` contains a quantitative trait. `Y_BT` contains a binary trait. Both have 2 variables: `ID` and `Y`.

Here is how you can load and view the summary of the datasets:


```{r}
# load the example dataset
data("CNVCOVY", package="CNVreg")
```


<details><summary>Here is the summary of each file in the dataset </summary>
```{r}
# view the dataset
summary(CNV)
summary(Cov)
summary(Y_QT)
summary(Y_BT)
```
</details>

Briefly, the dataset has the CNV (2680 records), covariates (`Sex` and `Age`), and outcome traits for 900 individuals. 


## Data preprocessing

The `prep()` function converts an individual's CNV events within a genomic region to fragments, and filter out rare events. It analyzes the adjacency relationship between CNV fragments and prepares different weight options for the penalized regression analysis. Here is how you can use the `prep()` function to preprocess `CNV`, `Cov` and a quantitative outcome `Y_QT`. You can replace `Y_QT` with a binary outcome `Y_BT` to prepare for a binary outcome regression.
     
```{r, warning=FALSE}
# data preprocessing 
frag_data_QT <- prep(CNV = CNV, Y = Y_QT, Z = Cov, rare.out = 0.05)
```

The result `frag_data_QT` is the output from the `prep()` function, which has a specially designed "WTsmth.data" format for easy application in the next step for CNV association analysis. It has 6 components.    
  
  * `design`: the CNV fragments in n by p dimensions, where n is the number of samples and p is the total number of CNV fragments. 
    
  * `Z`: a matrix of covariates with sample ID as rownames. The rownames are in the same order as in the outcomes. 
    
  * `Y`: a matrix of 1 column with sample ID as rownames. The rownames are in the same order as in the covariates.
    
  * `weight.structure`: a matrix that describes the adjacency structure of CNV fragments. The matrix is sparse and most values are zero, while adjacent non-zero values represent two adjacent CNV fragments that are overlapped by at least one CNV event in the population. 
    
  * `weight.options`: we provide 6 different options of weights that encourage differential information sharing based on the relationship between adjacent CNV fragments. Equal weight `eql`, Coscine-similarity based weight `wcs`, Inverse frequency weight `wif`, and these 3 weights further considering frequency of any CNV event (k) within each CNV-rich region ("keql", "kwcs", and "kwif").
             
  * `CNVR.info` summarizes the positions of all CNV fragments and their adjacency information. Each row represents a CNV fragment and the fragment names match the column names in `design`.


<details><summary> Here is a more detailed description of the output `frag_data_QT` </summary>
```{r}
# Format of `prep()` funtion output
 str(frag_data_QT)
 
```
</details>



## CNV association analysis with cross-validation

The `cvfit_WTSMTH()` function analyzes the association between a continuous/binary trait value and `CNV` while adjusting for the covariates `Cov`. Since we already have `Y_QT` prepared in the `prep()` step, we will fit a model to perform CNV association analysis for a continuous outcome first. 
  
 
```{r warning=FALSE}
  set.seed(12345)
  QT_TUNE <- cvfit_WTSMTH(data = frag_data_QT, 
                           lambda1=seq(-8, -3, 1), 
                           lambda2 = seq(12, 25, 2), 
                           weight="eql", 
                           family = "gaussian",
                           cv.control = list(n.fold = 5L, 
                                           n.core = 1L, 
                                           stratified = FALSE),
                          verbose = FALSE)
```


The `cvfit_WTSMTH()` function takes the output of `prep()` function as one of the major inputs.

`lambda1` and `lambda2` takes the candidate tuning parameters that control variable selection (`lambda1`) and effect smoothness (`lambda2`)

`weight` has six different options as described earlier. Since we only have a small dataset, varying the `weight` options will not have much influence on the model fitting results. With real CNV data with different similarity patterns and CNV frequencies, varying the `weight` option are expected to have different effects.

`family` has two options: `gaussian` for a continuous outcome, and `binomial` for a binary outcome.
                           
This function also supports parallel computing and change of n-folds in CV by adjust the `cv.control` list.

 * `n.fold` controls the number of folds used in CV.
 
 * `n.core` controls the the number of cores used in parallel computing.
 
 * `stratified` only has control for a binary outcome. We will skip it here and describe it later. 
 
If choose `verbose` = TRUE, it will print a message about what the program is currently working on. 
 
The output of the `cvfit_WTSMTH()` function is a list object containing 3 elements: `Loss`,  `lambda.selected`, and `coef`.


* `Loss`

  The `Loss` keeps track of the average validation loss in CV for each pair of candidate tuning parameters $\lambda_{1}$ and $\lambda_{2}$.      In the following table, the minimum loss is highlighted and the corresponding $\lambda_{1}$ and $\lambda_{2}$ values are selected to fit a     final model.
   
  In this simulated data, the variation of loss for different $\lambda_{2}$ with the same $\lambda_{1}$ is not very large. One reason is that    $\lambda_{2}$ controls the effect smoothness between adjacent CNVs, and the simulation data only has a small number of CNVs in adjacent that   share effects to other CNVs. The effect of changing $\lambda_{2}$ seems not prominence in this case. When we have more CNVs in adjacent and    share effects, it should have larger variance across $\lambda_{2}$.  

```{r echo=FALSE}
# loss matrix of candidate tuning parameters
 QT_TUNE$Loss %>% format( digits = 6) %>%
  mutate(
    across(2:ncol(QT_TUNE$Loss), 
               ~ cell_spec(.x, 
                           color = ifelse(.x > min(QT_TUNE$Loss[,2:ncol(QT_TUNE$Loss)])+0.0001, "black", "white"), 
                           background = ifelse(.x <= min(QT_TUNE$Loss[, 2:ncol(QT_TUNE$Loss)])+0.0001, "red", "white")
                           )
               )
)%>%
  kable(booktabs = FALSE, linesep = "",  align = "c", format = "html", escape = F,  caption = "Average loss for each pair of candidate tuning parameters") %>%add_header_above(c(" " = 1, "Lambda 1" = ncol(QT_TUNE$Loss)-1))
```

* `selected.lambda`

The `selected.lambda` are the optimal tuning parameters from the candidate lists that has the lowest loss, which can be confirmed with the `Loss` table.

```{r}
# selected optimal tuning parameters with minimum loss
 QT_TUNE$selected.lambda 
```
   
* `coef` 

The `coef` shows the estimated beta coefficients at the selected tuning parameters. It has `(intercept)`, CNV fragments (with detailed positions/type information), and covariate effects. In this small example, we can print all coefficient estimate, but you can modify the code to show only non-zero ones.  


Here lists the coefficients for `(Intercept)` and covariates. The characteristics for CNV (CHR, CNV.start, CNV.end, and deldup) are left as `NA`s intentionally in the original output. Here we only show the effect estimate. 

```{r}
##coefficients of intercept and covariates 
QT_TUNE$coef[c(1, 21:22), c("Vnames", "coef") ] 
```


Here lists the coefficients for CNVs and the corresponding plots.

We highlight the regions with adjacent CNVs. From the coefficient estimates, the model selectes several non-zero CNVs (data points) that are associated with the trait. Among the data points, the red ones have stronger effect than the black dots. The black dots are likely noise. 

We also zoom in to show the effect smoothness within the highlighted regions with adjacent CNVs. 

The results illustrate the variable selection and effect smoothness of the penalized regression method for CNV association analysis. 


<details><summary>  CNV coefficient estimate of a fine-tuned model for a continuous outcome.</summary>   
```{r}
# estimated coefficents for CNV
QT_TUNE$coef[2:20, ]
# non-zero coefficients 
# QT_TUNE$coef[which(abs(QT_TUNE$coef$coef)>0), ] 
```   
</details>


```{r warning=FALSE, echo=FALSE}
# plot the coefficients 
# keep CNV fragments and exclude intercept and covariates(Age, Sex) in ploting, row(2:20)

CNVR <- frag_data_QT$CNVR.info %>% group_by(CNV.id, deldup)%>%
  summarise(CNV.start = min(lower.boundary), 
            CNV.end = max(upper.boundary), 
            nfrag = length(CNV.id), 
            .groups = 'drop')
CNVR_adj <- CNVR[CNVR$nfrag > 1, ]
CNV_coef <- QT_TUNE$coef[2:20,]


P <- ggplot(CNV_coef, aes(x= CNV_coef$CNV.start , y=CNV_coef$coef))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef) > 0.001, "red", ifelse(abs(CNV_coef$coef) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("Genomic position", labels = as.character(CNVR_adj$CNV.start), breaks = CNVR_adj$CNV.start)+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start, xend=CNV_coef$CNV.end, y=CNV_coef$coef, yend=CNV_coef$coef)+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[1]-1000000, xmax=CNVR_adj$CNV.end[1]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x = CNVR_adj$CNV.start[1], y = max(CNV_coef$coef) +0.0002, label = "A", color="red")+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[2]-1000000, xmax=CNVR_adj$CNV.end[2]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[3]-1000000, xmax=CNVR_adj$CNV.end[3]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x = CNVR_adj$CNV.start[3], y = max(CNV_coef$coef)+0.0002, label = "B", color="red")+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[4]-1000000, xmax=CNVR_adj$CNV.end[4]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)

    

PA <- ggplot(CNV_coef[2:5,], aes(x= 1/2*(CNV_coef$CNV.start[2:5]+CNV_coef$CNV.end[2:5]) , y=CNV_coef$coef[2:5]))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef[2:5]) > 0.001, "red", ifelse(abs(CNV_coef$coef[2:5]) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("", labels = as.character(CNV_coef$CNV.start[2:5]), breaks = CNV_coef$CNV.start[2:5])+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start[2:5], xend=CNV_coef$CNV.end[2:5], y=CNV_coef$coef[2:5], yend=CNV_coef$coef[2:5])+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[1]-10, xmax=CNVR_adj$CNV.end[1]+10, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x =1/2*( CNVR_adj$CNV.start[1]+CNVR_adj$CNV.end[1]), y = max(CNV_coef$coef) +0.0002, label = "A", color="red")
  

PB <- ggplot(CNV_coef[11:14,], aes(x= 1/2*(CNV_coef$CNV.start[11:14]+CNV_coef$CNV.end[11:14]) , y=CNV_coef$coef[11:14]))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef[11:14]) > 0.001, "red", ifelse(abs(CNV_coef$coef[11:14]) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("Genomic position", labels = as.character(CNV_coef$CNV.start[11:14]), breaks = CNV_coef$CNV.start[11:14])+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start[11:14], xend=CNV_coef$CNV.end[11:14], y=CNV_coef$coef[11:14], yend=CNV_coef$coef[11:14])+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[3]-10, xmax=CNVR_adj$CNV.end[3]+10, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x =1/2*( CNVR_adj$CNV.start[3]+CNVR_adj$CNV.end[3]), y = max(CNV_coef$coef) +0.0002, label = "B", color="red")


P + (PA/PB) + 
  plot_annotation(title = "CNV coefficient estiamte across the genomic region - A fine-tuned model")+ 
  plot_layout(axes = "collect", widths = c(2,1)) + 
  plot_layout( guide = "collect") & theme(legend.position="Top",
                      legend.text = element_text(size=12), legend.title = element_text(size=12)) 
```
    




## CNV association analysis with a specific pair of tuning parameters

The user can choose function `fit_WTSMTH()` to show the regression result with some random combination of parameters without going through the CV process. 
Although, it is much faster to perform `fit_WTSMTH()`, we do recommend the user to stick with the parameter tuning procedure with `cvfit_WTSMTH()` and find the best parameters and the best fitted model. 
 
 
The `fit_WTSMTH()` function and the `cvfit_WTSMTH()` function uses the same analytical methods to perform CNV association analysis with penalized regression. Unlike the `cvfit_WTSMTH()` function that will fine-tune the parameters and select the optimal combination of $\lambda_{1}$ and $\lambda_{2}$ from a series of candidates, the `fit_WTSMTH()` function takes a user-specified value for $\lambda_{1}$ and $\lambda_{2}$ and estimate the coefficients for the given pair of parameters. 


Here we fit the model with some random parameters $\lambda_{1}$ and $\lambda_{2}$ . We can compare the regression results with the fine-tuned model above. This model has more noise signals (effect of `sex` is no longer 0) and missed some parts of the two consecutive CNV regions. The effect size between adjacent CNVs are quite different from each other. 

       
```{r}
# we know the optimal tuning parameters and directly apply it here.
QT_fit <- fit_WTSMTH(frag_data_QT, 
                      lambda1 = -5.5, 
                      lambda2 = 2, 
                      weight="eql",
                      family="gaussian")
```


Here lists the coefficients for `(Intercept)` and covariates. 

```{r}

QT_fit[c(1, 21:22), c("Vnames", "coef") ]

```


Here lists the coefficients for CNVs and the corresponding plots. 

<details><summary> CNV coefficient estimate with given parameters for a continuous outcome. <summary>
```{r}
QT_fit[2:20, ]
```
</details>


```{r warning=FALSE, echo=FALSE}
CNV_coef <- QT_fit[2:20,]


P <- ggplot(CNV_coef, aes(x= CNV_coef$CNV.start , y=CNV_coef$coef))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef) > 0.001, "red", ifelse(abs(CNV_coef$coef) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("Genomic position", labels = as.character(CNVR_adj$CNV.start), breaks = CNVR_adj$CNV.start)+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start, xend=CNV_coef$CNV.end, y=CNV_coef$coef, yend=CNV_coef$coef)+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[1]-1000000, xmax=CNVR_adj$CNV.end[1]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x = CNVR_adj$CNV.start[1], y = max(CNV_coef$coef) +0.001, label = "A", color="red")+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[2]-1000000, xmax=CNVR_adj$CNV.end[2]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[3]-1000000, xmax=CNVR_adj$CNV.end[3]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x = CNVR_adj$CNV.start[3], y = max(CNV_coef$coef)+0.001, label = "B", color="red")+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[4]-1000000, xmax=CNVR_adj$CNV.end[4]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)

    

PA <- ggplot(CNV_coef[2:5,], aes(x= 1/2*(CNV_coef$CNV.start[2:5]+CNV_coef$CNV.end[2:5]) , y=CNV_coef$coef[2:5]))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef[2:5]) > 0.001, "red", ifelse(abs(CNV_coef$coef[2:5]) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("", labels = as.character(CNV_coef$CNV.start[2:5]), breaks = CNV_coef$CNV.start[2:5])+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start[2:5], xend=CNV_coef$CNV.end[2:5], y=CNV_coef$coef[2:5], yend=CNV_coef$coef[2:5])+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[1]-10, xmax=CNVR_adj$CNV.end[1]+10, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x =1/2*( CNVR_adj$CNV.start[1]+CNVR_adj$CNV.end[1]), y = max(CNV_coef$coef) +0.001, label = "A", color="red")
  

PB <- ggplot(CNV_coef[11:14,], aes(x= 1/2*(CNV_coef$CNV.start[11:14]+CNV_coef$CNV.end[11:14]) , y=CNV_coef$coef[11:14]))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef[11:14]) > 0.001, "red", ifelse(abs(CNV_coef$coef[11:14]) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("Genomic position", labels = as.character(CNV_coef$CNV.start[11:14]), breaks = CNV_coef$CNV.start[11:14])+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start[11:14], xend=CNV_coef$CNV.end[11:14], y=CNV_coef$coef[11:14], yend=CNV_coef$coef[11:14])+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[3]-10, xmax=CNVR_adj$CNV.end[3]+10, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x =1/2*( CNVR_adj$CNV.start[3]+CNVR_adj$CNV.end[3]), y = max(CNV_coef$coef) +0.001, label = "B", color="red")


P + (PA/PB) +plot_annotation(title = "CNV coefficient estiamte across the genomic region - A random model")+
  plot_layout(axes = "collect", widths = c(2,1)) + plot_layout( guide = "collect") & theme(legend.position="Top",
                      legend.text = element_text(size=12), legend.title = element_text(size=12)) 
```


## Repeat the procedure with a binary outcome. 

#### Preprocess the data 

```{r, warning=FALSE}
# data preprocessing 
frag_data_BT <- prep(CNV = CNV, Y = Y_BT, Z = Cov, rare.out = 0.05)
```

The result `frag_data_QT` is the output from the `prep()` function, which has a specially designed "WTsmth.data" format for easy application in the next step for CNV association analysis. It contains the same list as mentioned earlier: `design`, `Z`, `Y`, `weight.structure`, `weight.options`, and `CNVR.info`.

Here is an alternative way to prepare data when performing CNV association analysis with the same set of CNVs for multiple outcome traits. Since we have the same `CNV` data, `Cov` data, and a different outcome trait `Y_BT`, we can manually format `Y_BT` to match the format in `frag_data_QT$Y`. 

<details><summary> Here we provide the command for reference. </summary>
```{r}
# It would be useful when we have large CNV data set and perform association analysis for multiple traits with the same set of CNV data. 

## copy frag_data_QT
#frag_data_BT <- frag_data_QT
#
### replace Y with Y_BT in the correct format: ordered named vector
### order the sample in Y_BT as in frag_data_QT$Y
#rownames(Y_BT) <- Y_BT$ID
#
#frag_data_BT$Y <- Y_BT[names(frag_data_QT$Y), "Y"] |> drop()
#names(frag_data_QT$Y) <- rownames(frag_data_QT$Y) 

## Directly replace frag_data_QT$Z is also possible, keep in mind to use the correct variable name and order.
```
</details>


#### CNV association analysis with cross-validation


For CNV association analysis with a binary outcome, choose `family` = "binomial". 

There are a few more options for the binary trait scenario compared to the continuous trait. 

1. `stratified` within the `cv.control` list: If one category of the binary outcome is considered "rare", `stratified` = TRUE is recommended to make sure the data splits are having the same proportion of outcomes. 

2. `iter.control`: For a binary outcome, we can also adjust the `iter.control` list with desired threshold that is deemed converged for coefficient estimate of a binary outcome.



```{r}
set.seed(12345)
BT_TUNE <- cvfit_WTSMTH(frag_data_BT, 
                         lambda1 = seq(-5.25, -4.75, 0.25), 
                         lambda2 = seq(2,  8, 2), 
                         weight="eql",
                         family="binomial", 
                         cv.control = list(n.fold = 5L, 
                                           n.core = 1L, 
                                           stratified = FALSE),
                         iter.control = list(max.iter = 8L, 
                                             tol.beta = 10^(-3), 
                                             tol.loss = 10^(-6)), 
                        verbose = FALSE)

```

The output of the `cvfit_WTSMTH()` function has the same list object containing 3 elements: `Loss`,  `lambda.selected`, and `coef`.

* `Loss`

`The `Loss` keeps track of the average validation loss in CV for each pair of candidate tuning parameters $\lambda_{1}$ and $\lambda_{2}$. In the following table, the minimum loss is highlighted and the corresponding $\lambda_{1}$ and $\lambda_{2}$ values are selected to fit a final model. 

Since the regression process for a binary trait takes longer time to converge, here we only use a short list of candidate tuning parameters for illustration purpose. 
 
```{r echo=FALSE}
# loss matrix of candidate tuning parameters
 
BT_TUNE$Loss %>% format( digits = 6) %>%
  mutate(
    across(2:ncol(BT_TUNE$Loss), 
               ~ cell_spec(.x, 
                           color = ifelse(.x > min(BT_TUNE$Loss[,2:ncol(BT_TUNE$Loss)])+0.000001, "black", "white"), 
                           background = ifelse(.x <= min(BT_TUNE$Loss[, 2:ncol(BT_TUNE$Loss)])+0.000001, "red", "white")
                           )
               )
)%>%
  kable(booktabs = FALSE, linesep = "",  align = "c", format = "html", escape = F,  caption = "Average loss for each pair of candidate tuning parameters") %>%add_header_above(c(" " = 1, "Lambda 1" = ncol(BT_TUNE$Loss)-1))
```

* `selected.lambda`

The `selected.lambda` are the optimal tuning parameters from the candidate lists that has the lowest loss, which can be confirmed with the `Loss` table.

```{r}
# selected optimal tuning parameters with minimum loss
 BT_TUNE$selected.lambda 
```
 
 * `coef`  
 
The estimated beta coefficients `coef` at the selected tuning parameters. It has `(intercept)`, CNV fragments (with detailed positions/type information), and covariate effects. In this small data example, we can print all coefficient estimate, but you can modify the code to show non-zero ones or the first few ones.  
    
    
    
Here lists the coefficients for `(Intercept)` and covariates. 

```{r}

BT_TUNE$coef[c(1, 21:22), c("Vnames", "coef") ]

```


Here lists the coefficients for CNVs and the corresponding plots. 

<details><summary> CNV coefficient estimate of a fine-tuned model for a binary outcome. </summary>
```{r}
BT_TUNE$coef[2:20, ]
```
</details>



```{r warning=FALSE, echo=FALSE}
CNV_coef <- BT_TUNE$coef[2:20,]


P <- ggplot(CNV_coef, aes(x= CNV_coef$CNV.start , y=CNV_coef$coef))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef) > 0.001, "red", ifelse(abs(CNV_coef$coef) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("Genomic position", labels = as.character(CNVR_adj$CNV.start), breaks = CNVR_adj$CNV.start)+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start, xend=CNV_coef$CNV.end, y=CNV_coef$coef, yend=CNV_coef$coef)+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[1]-1000000, xmax=CNVR_adj$CNV.end[1]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x = CNVR_adj$CNV.start[1], y = max(CNV_coef$coef) +0.001, label = "A", color="red")+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[2]-1000000, xmax=CNVR_adj$CNV.end[2]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[3]-1000000, xmax=CNVR_adj$CNV.end[3]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x = CNVR_adj$CNV.start[3], y = max(CNV_coef$coef)+0.001, label = "B", color="red")+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[4]-1000000, xmax=CNVR_adj$CNV.end[4]+1000000, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)

    

PA <- ggplot(CNV_coef[2:5,], aes(x= 1/2*(CNV_coef$CNV.start[2:5]+CNV_coef$CNV.end[2:5]) , y=CNV_coef$coef[2:5]))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef[2:5]) > 0.001, "red", ifelse(abs(CNV_coef$coef[2:5]) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("", labels = as.character(CNV_coef$CNV.start[2:5]), breaks = CNV_coef$CNV.start[2:5])+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start[2:5], xend=CNV_coef$CNV.end[2:5], y=CNV_coef$coef[2:5], yend=CNV_coef$coef[2:5])+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[1]-10, xmax=CNVR_adj$CNV.end[1]+10, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x =1/2*( CNVR_adj$CNV.start[1]+CNVR_adj$CNV.end[1]), y = max(CNV_coef$coef) +0.001, label = "A", color="red")
  

PB <- ggplot(CNV_coef[11:14,], aes(x= 1/2*(CNV_coef$CNV.start[11:14]+CNV_coef$CNV.end[11:14]) , y=CNV_coef$coef[11:14]))+theme_bw()+ 
  geom_point(aes(color = ifelse(abs(CNV_coef$coef[11:14]) > 0.001, "red", ifelse(abs(CNV_coef$coef[11:14]) < 10^(-8), "white", "black")))) + 
  scale_color_identity()+
  theme(axis.text.x = element_text(size = 10, angle = 25, vjust = 1, hjust = 1))+
  scale_x_continuous("Genomic position", labels = as.character(CNV_coef$CNV.start[11:14]), breaks = CNV_coef$CNV.start[11:14])+
  ylab("CNV coefficients")+
  geom_segment(x=CNV_coef$CNV.start[11:14], xend=CNV_coef$CNV.end[11:14], y=CNV_coef$coef[11:14], yend=CNV_coef$coef[11:14])+
  geom_rect(aes(xmin=CNVR_adj$CNV.start[3]-10, xmax=CNVR_adj$CNV.end[3]+10, ymin=min(CNV_coef$coef), ymax = max(CNV_coef$coef)), fill="yellow", alpha = 0.02)+
  annotate("text", x =1/2*( CNVR_adj$CNV.start[3]+CNVR_adj$CNV.end[3]), y = max(CNV_coef$coef) +0.001, label = "B", color="red")


P + (PA/PB) +plot_annotation(title = "CNV coefficient estiamte across the genomic region - A fine-tuned model")+
  plot_layout(axes = "collect", widths = c(2,1)) + plot_layout( guide = "collect") & theme(legend.position="Top",
                      legend.text = element_text(size=12), legend.title = element_text(size=12)) 
```
