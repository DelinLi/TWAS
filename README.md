## Examples code for [my TWAS paper](https://academic.oup.com/plphys/advance-article/doi/10.1093/plphys/kiab161/6212071) on Plant Physiology 


#### Data Prepare
Convert the numeric expression into genotype-like (like numeric range from 0 to 2 for GAPIT) format required by GWAS tools. I could not find the best yet, but the Quantile one gave overall better results in cases of the [TWAS paper](https://academic.oup.com/plphys/advance-article/doi/10.1093/plphys/kiab161/6212071)
1. Linearly (used by the [eRD-GWAS](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1328-6)
2. Quantile (used in the paper and here to handle extreme expression values (outliers)): 
First, expression values smaller than quantile 5 were converted to 0 and values larger than quantile 95 were converted to 2. The remaining expression values were linearly transformed into values between 0 and 2.
<pre>
#R code for the Quantile solution (if you have a higher sample size, you can try Quantile10 Quantile90 or any value you want)
  Quantile<- apply(Expression.TPM.Matrix,1,
                         function(A){min_x=as.numeric(quantile(A,5));
                         max_x=as.numeric(quantile(A,95));
                         out<-(2*(A-min_x)/(max_x-min_x));
                         out[out>2]<-2;out[out< 0]<- 0;return(out)})
</pre>
3. QQNorm, a solution used by eQTL paper traiting expression as trait

#### Run the TWAS on two Traits


#### Visualization 

#### Citation:
Delin Li, Qiang Liu, Patrick S Schnable, TWAS results are complementary to and less affected by linkage disequilibrium than GWAS, Plant Physiology, 2021;, kiab161, https://doi.org/10.1093/plphys/kiab161
