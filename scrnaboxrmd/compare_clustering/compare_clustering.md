---
title: "Compare Clustering Prediction "
author: "Saeid Amiri"
date: "2022-12-05"
output: 
  html_document:
    keep_md: true
---



# content
- Accuracy
   + Clustering accuracy
   + Association measure
- Computation

## Accuracy
### Clustering accurary
ACC: Computes one subtract of the error rate; it uses the hungarian method to macth to clusters.  
AMI: Adjusted Mutual Information
NMI: computes the normalized mutual information
RI: Rand Index
ARI: Adjusted Rand Index

###  Assocation measure
Cramér's V:  is a measure of association between two nominal variables, giving a value between 0 and +1 (inclusive)

## Computation


```r
X<-sample(c('a1', 'a2'), 1000, replace=TRUE)
Y<-sample(c('b1', 'b2', 'b3'), 1000, replace=TRUE)
compute_index(X,Y)
```

```
## Package 'mclust' version 6.0.0
## Type 'citation("mclust")' for citing this R package in publications.
```

```
## 
## Attaching package: 'mclust'
```

```
## The following object is masked from 'package:purrr':
## 
##     map
```

```
## $chisq_test
## 
## 	Pearson's Chi-squared test
## 
## data:  con1
## X-squared = 2.014, df = 2, p-value = 0.3653
## 
## 
## $contigency_table
##     Y
## X     b1  b2  b3
##   a1 156 176 162
##   a2 181 172 153
## 
## $prop_table
##     Y
## X       b1    b2    b3
##   a1 0.156 0.176 0.162
##   a2 0.181 0.172 0.153
## 
## $prop_table_row
##     Y
## X           b1        b2        b3
##   a1 0.3157895 0.3562753 0.3279352
##   a2 0.3577075 0.3399209 0.3023715
## 
## $prop_table_col
##     Y
## X           b1        b2        b3
##   a1 0.4629080 0.5057471 0.5142857
##   a2 0.5370920 0.4942529 0.4857143
## 
## $cramerv
## [1] 0.04487771
## 
## $accuracy
## [1] 0.357
## 
## $adjusted_mutual_information
## [1] 0.0009179856
## 
## $normalized_mutual_information
## [1] 0.0009179856
## 
## $rand_index
## [1] 0.5001502
## 
## $adjusted_rand_index
## [1] 1.451052e-05
```
