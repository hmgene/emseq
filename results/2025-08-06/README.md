    knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE)
    library(ggplot2)

    ## Warning: package 'ggplot2' was built under R version 4.4.3

    library(dplyr)

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

    # Simulated example dataset
    df <- data.frame(
      Sample = paste0("S", 1:6),
      group = c("E", "E", "W", "W", "Y", "Y"),
      CpG = c(1200, 1350, 1100, 1080, 1150, 1170),
      uCpG = c(300, 250, 400, 420, 380, 360)
    )

    df <- df %>%
      mutate(
        total = CpG + uCpG,
        prop = CpG / total
      )

    knitr::kable(df, caption = "CpG and uCpG Counts with Methylation Proportions")

<table>
<caption>CpG and uCpG Counts with Methylation Proportions</caption>
<thead>
<tr class="header">
<th style="text-align: left;">Sample</th>
<th style="text-align: left;">group</th>
<th style="text-align: right;">CpG</th>
<th style="text-align: right;">uCpG</th>
<th style="text-align: right;">total</th>
<th style="text-align: right;">prop</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">S1</td>
<td style="text-align: left;">E</td>
<td style="text-align: right;">1200</td>
<td style="text-align: right;">300</td>
<td style="text-align: right;">1500</td>
<td style="text-align: right;">0.8000000</td>
</tr>
<tr class="even">
<td style="text-align: left;">S2</td>
<td style="text-align: left;">E</td>
<td style="text-align: right;">1350</td>
<td style="text-align: right;">250</td>
<td style="text-align: right;">1600</td>
<td style="text-align: right;">0.8437500</td>
</tr>
<tr class="odd">
<td style="text-align: left;">S3</td>
<td style="text-align: left;">W</td>
<td style="text-align: right;">1100</td>
<td style="text-align: right;">400</td>
<td style="text-align: right;">1500</td>
<td style="text-align: right;">0.7333333</td>
</tr>
<tr class="even">
<td style="text-align: left;">S4</td>
<td style="text-align: left;">W</td>
<td style="text-align: right;">1080</td>
<td style="text-align: right;">420</td>
<td style="text-align: right;">1500</td>
<td style="text-align: right;">0.7200000</td>
</tr>
<tr class="odd">
<td style="text-align: left;">S5</td>
<td style="text-align: left;">Y</td>
<td style="text-align: right;">1150</td>
<td style="text-align: right;">380</td>
<td style="text-align: right;">1530</td>
<td style="text-align: right;">0.7516340</td>
</tr>
<tr class="even">
<td style="text-align: left;">S6</td>
<td style="text-align: left;">Y</td>
<td style="text-align: right;">1170</td>
<td style="text-align: right;">360</td>
<td style="text-align: right;">1530</td>
<td style="text-align: right;">0.7647059</td>
</tr>
</tbody>
</table>

CpG and uCpG Counts with Methylation Proportions

    ggplot(df, aes(x = group, y = prop, fill = group)) +
      geom_boxplot(width = 0.4, alpha = 0.5, outlier.shape = NA) +
      geom_jitter(width = 0.1, size = 2) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(title = "Methylation Proportion by Group",
           y = "Methylation (%)", x = "Group") +
      theme_minimal()

![](README_files/figure-markdown_strict/setup-1.png)

    glm_fit <- glm(cbind(CpG, uCpG) ~ group, data = df, family = binomial)
    summary(glm_fit)

    ## 
    ## Call:
    ## glm(formula = cbind(CpG, uCpG) ~ group, family = binomial, data = df)
    ## 
    ## Coefficients:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  1.53393    0.04701  32.627  < 2e-16 ***
    ## groupW      -0.55615    0.06236  -8.919  < 2e-16 ***
    ## groupY      -0.39126    0.06319  -6.192 5.94e-10 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## (Dispersion parameter for binomial family taken to be 1)
    ## 
    ##     Null deviance: 95.830  on 5  degrees of freedom
    ## Residual deviance: 11.538  on 3  degrees of freedom
    ## AIC: 62.067
    ## 
    ## Number of Fisher Scoring iterations: 4

    cat("Null deviance:", round(glm_fit$null.deviance), "on", glm_fit$df.null, "degrees of freedom\n")

    ## Null deviance: 96 on 5 degrees of freedom

    cat("Residual deviance:", round(glm_fit$deviance), "on", glm_fit$df.residual, "degrees of freedom\n")

    ## Residual deviance: 12 on 3 degrees of freedom

    cat("AIC:", round(AIC(glm_fit)), "\n")

    ## AIC: 62

    invlogit <- function(x) exp(x) / (1 + exp(x))

    p_E <- invlogit(coef(glm_fit)["(Intercept)"])
    p_W <- invlogit(coef(glm_fit)["(Intercept)"] + coef(glm_fit)["groupW"])
    p_Y <- invlogit(coef(glm_fit)["(Intercept)"] + coef(glm_fit)["groupY"])

    rate_df <- data.frame(
      Group = c("E", "W", "Y"),
      Estimated_Proportion = round(c(p_E, p_W, p_Y), 4)
    )

    knitr::kable(rate_df, caption = "Estimated Methylation Proportions by Group (from GLM)")

<table>
<caption>Estimated Methylation Proportions by Group (from GLM)</caption>
<thead>
<tr class="header">
<th style="text-align: left;">Group</th>
<th style="text-align: right;">Estimated_Proportion</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">E</td>
<td style="text-align: right;">0.8226</td>
</tr>
<tr class="even">
<td style="text-align: left;">W</td>
<td style="text-align: right;">0.7267</td>
</tr>
<tr class="odd">
<td style="text-align: left;">Y</td>
<td style="text-align: right;">0.7582</td>
</tr>
</tbody>
</table>

Estimated Methylation Proportions by Group (from GLM)
