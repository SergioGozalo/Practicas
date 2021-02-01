Functional tables analysis
================
Sergio Gozalo
20 de enero de 2021

## Loading necessary libraries

``` r
library("vegan")
```

    ## Loading required package: permute

    ## Loading required package: lattice

    ## This is vegan 2.5-7

``` r
library("ggplot2")
library("tidyr")
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

# EBI tables

## Table reading EBI

``` r
ebi_go_abundances <- read.table("functional.tables/ERP112966_GO_abundances_v4.1.tsv", header = TRUE, row.names = 1, sep ="\t")
ebi_go_abundances$description <- NULL #Need to remove this column to work
ebi_go_abundances$category <- NULL #Need to remove this column to work
ebi_go_abundances <- t(ebi_go_abundances)

ebi_go.slim_abundances <- read.table("functional.tables/ERP112966_GO-slim_abundances_v4.1.tsv", header = TRUE, row.names = 1, sep ="\t")
ebi_go.slim_abundances$description <- NULL
ebi_go.slim_abundances$category <- NULL
ebi_go.slim_abundances <- t(ebi_go.slim_abundances)

ebi_ipr_abundances <- read.table("functional.tables/ERP112966_IPR_abundances_v4.1.tsv", header = TRUE, row.names = 1, sep ="\t")
ebi_ipr_abundances$description <- NULL
ebi_ipr_abundances$category <- NULL
ebi_ipr_abundances <- t(ebi_ipr_abundances)
```

## Richness

``` r
# Subsampling
n1 <- min(rowSums(ebi_go_abundances))
ebi_go_abundances_ss <- rrarefy(ebi_go_abundances, n1)

n2 <- min(rowSums(ebi_go.slim_abundances))
ebi_go.slim_abundances_ss <- rrarefy(ebi_go.slim_abundances, n2)

n3 <- min(rowSums(ebi_ipr_abundances))
ebi_ipr_abundances_ss <- rrarefy(ebi_ipr_abundances, n3)

# Richness

ebi_go_abundances_richness <- estimateR(ebi_go_abundances_ss)
ebi_go.slim_abundances_richness <- estimateR(ebi_go.slim_abundances_ss)
ebi_ipr_abundances_richness <- estimateR(ebi_ipr_abundances_ss)
```

## NMDS

### GO

``` r
ebi_go_abundances_ss_bray <- vegdist(ebi_go_abundances_ss, method="bray")
ebi_go_abundances_ss_bray.nmds <- monoMDS(ebi_go_abundances_ss_bray)

plot(ebi_go_abundances_ss_bray.nmds, main = "GO")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-4-1.png)

Se observa un cluster, el resto son independientes.

``` r
stressplot(ebi_go_abundances_ss_bray.nmds, main = "GO")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-5-1.png)

### GO-SLIM

``` r
ebi_go.slim_abundances_ss_bray <- vegdist(ebi_go.slim_abundances_ss, method="bray")
ebi_go.slim_abundances_ss_bray.nmds <- monoMDS(ebi_go.slim_abundances_ss_bray)

plot(ebi_go.slim_abundances_ss_bray.nmds, main = "GO-SLIM")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-6-1.png)

Se observa un cluster, el resto son independientes.

``` r
stressplot(ebi_go.slim_abundances_ss_bray.nmds, main = "GO-SLIM")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-7-1.png)

### IPR

``` r
ebi_ipr_abundances_ss_bray <- vegdist(ebi_ipr_abundances_ss, method="bray")
ebi_ipr_abundances_ss_bray.nmds <- monoMDS(ebi_ipr_abundances_ss_bray)

plot(ebi_ipr_abundances_ss_bray.nmds, main = "IPR")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-8-1.png)

Se ve un cluster.

``` r
stressplot(ebi_ipr_abundances_ss_bray.nmds, main = "IPR")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-9-1.png)

## Dendograms

### GO

``` r
ebi_go_abundances_ss_bray_hclust <- hclust(ebi_go_abundances_ss_bray, method="average")
plot(ebi_go_abundances_ss_bray_hclust, main = "GO")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-10-1.png)

Se observan dos grupos principales, uno de ellos es mucho mas grande y tiene bastantes subdivisiones, siendo las dos primeras.

### GO-SLIM

``` r
ebi_go.slim_abundances_ss_bray_hclust <- hclust(ebi_go.slim_abundances_ss_bray, method="average")
plot(ebi_go.slim_abundances_ss_bray_hclust, main = "GO-SLIM")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-11-1.png)

Muy parecido al GO, pero con distinta longitud de ramas.

### IPR

``` r
ebi_ipr_abundances_ss_bray_hclust <- hclust(ebi_ipr_abundances_ss_bray, method="average")
plot(ebi_ipr_abundances_ss_bray_hclust, main = "IPR")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-12-1.png)

Se observan dos grupos principales, uno de ellos es mucho mas grande y tiene bastantes subdivisiones, distribuidas de distinta forma que en GO.

# ICM tables

## Table reading ICM

``` r
icm_cog_meta <- t(read.table("functional.tables/EMOSE-GC_ICM_250bp_COG.lengthNorm.metaGsizeNorm.counts.tbl", header = TRUE, sep = "\t", row.names = 1))
icm_cog_scg <- t(read.table("functional.tables/EMOSE-GC_ICM_250bp_COG.lengthNorm.SCGnorm.counts.tbl", header = TRUE, sep = "\t", row.names = 1))

icm_kegg_meta <- t(read.table("functional.tables/EMOSE-GC_ICM_250bp_KEGG.ko.lengthNorm.metaGsizeNorm.counts.tbl", header = TRUE, sep = "\t", row.names = 1))
icm_kegg_scg <- t(read.table("functional.tables/EMOSE-GC_ICM_250bp_KEGG.ko.lengthNorm.SCGnorm.counts.tbl", header = TRUE, sep = "\t", row.names = 1))

icm_pfam_meta <- t(read.table("functional.tables/EMOSE-GC_ICM_250bp_pfam.lengthNorm.metaGsizeNorm.counts.tbl", header = TRUE, sep = "\t", row.names = 1))
icm_pfam_scg <- t(read.table("functional.tables/EMOSE-GC_ICM_250bp_pfam.lengthNorm.SCGnorm.counts.tbl", header = TRUE, sep = "\t", row.names = 1))
```

## NMDs

### COG

``` r
icm_cog_meta_bray <- vegdist(icm_cog_meta, method="bray")
icm_cog_scg_bray <- vegdist(icm_cog_scg, method="bray")

icm_cog_meta_bray.nmds <- monoMDS(icm_cog_meta_bray)
icm_cog_scg_bray.nmds <- monoMDS(icm_cog_scg_bray)

# MetaGS
plot(icm_cog_meta_bray.nmds, main = "COG MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
# SC
plot(icm_cog_scg_bray.nmds, main = "COG SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-14-2.png)

COG MetaGS: Se pueden observar 4 clusters distintos, tres de ellos cercanos y uno mas distante y sin relación.

COG SCG: Se peuden observar 3 clusters, dos de ellos proximos y uno mas separado pero que aun parece mantener una relación.

``` r
stressplot(icm_cog_meta_bray.nmds, main = "COG MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-15-1.png)

``` r
stressplot(icm_cog_scg_bray.nmds, main = "COG SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-15-2.png)

### KEGG

``` r
icm_kegg_meta_bray <-vegdist(icm_kegg_meta, method="bray") 
icm_kegg_scg_bray <-vegdist(icm_kegg_scg, method="bray") 

icm_kegg_meta_bray.nmds <- monoMDS(icm_kegg_meta_bray)
icm_kegg_scg_bray.nmds <- monoMDS(icm_kegg_scg_bray)

# MetaGS
plot(icm_kegg_meta_bray.nmds, main = "KEGG MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-16-1.png)

``` r
# SC
plot(icm_kegg_scg_bray.nmds, main = "KEGG SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-16-2.png)

KEGG MetaGS: Se observan 4 clusters, en este caso parecen mejor separados, aunque dos de ellos son bastante proximos.

KEGG SCG: Se observan 4 clusters, dos de ellos muy proximos, cerca otro y el cuarto parece bastante alejado.

``` r
stressplot(icm_kegg_meta_bray.nmds, main = "KEGG MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-17-1.png)

``` r
stressplot(icm_kegg_scg_bray.nmds, main = "KEGG SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-17-2.png)

### PFAM

``` r
icm_pfam_meta_bray <-vegdist(icm_pfam_meta, method="bray") 
icm_pfam_scg_bray <-vegdist(icm_pfam_scg, method="bray")

icm_pfam_meta_bray.nmds <- monoMDS(icm_pfam_meta_bray)
icm_pfam_scg_bray.nmds <- monoMDS(icm_pfam_scg_bray)

# MetaGS
plot(icm_pfam_meta_bray.nmds, main = "PFAM MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-18-1.png)

``` r
# SC
plot(icm_pfam_scg_bray.nmds, main = "PFAM SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-18-2.png)

PFAM MetaGS: En este caso los clusters no estan tan bien definidos, pero se pueden llegar a distinguir 4, dos muy cercanos, un tercero que mantiene relacion y el cuarto totalmente aislado.

PFAM SCG: Se observan 4 clusters, dos de ellos casi solapados, el tercero bastante cerca (mas que en los casos anteriores) y el cuarto aislado.

``` r
stressplot(icm_pfam_meta_bray.nmds, main = "PFAM MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-19-1.png)

``` r
stressplot(icm_pfam_scg_bray.nmds, main = "PFAM SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-19-2.png)

## Dendograms

### COG

``` r
icm_cog_meta_bray_hclust <- hclust(icm_cog_meta_bray, method="average")
plot(icm_cog_meta_bray_hclust, main = "COG MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-20-1.png)

``` r
icm_cog_scg_bray_hclust <- hclust(icm_cog_scg_bray, method="average")
plot(icm_cog_scg_bray_hclust, main = "COG SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-20-2.png)

COG MetaGS: Se observan dos grupos principales que, a su vez, tienen dos subgrupos. Concuerda con los 4 clusters.

COG SCG: Se observan dos grupos principales, uno de ellos parece tener una subdivision importante, en el otro parece no ser tan importante. Concuerda con los clusters.

### KEGG

``` r
icm_kegg_meta_bray_hclust <-hclust(icm_kegg_meta_bray, method="average")
plot(icm_kegg_meta_bray_hclust, main = "KEGG MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-21-1.png)

``` r
icm_kegg_scg_bray_hclust <-hclust(icm_kegg_scg_bray, method="average")
plot(icm_kegg_scg_bray_hclust, main = "KEGG SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-21-2.png)

KEGG MetaGS: Se pueden observar dos grupos principales y cada uno de ellos dividido en otros dos grupos principales. Concuerda con los clusters.

KEGG SCG: Se pueden observar 2 grupos principales, divididos en otros dos grupos principales cada uno. Concuerda con los clusters.

Muy parecidos a los obtenidos usando COG.

### PFAM

``` r
icm_pfam_meta_bray_hclust <-hclust(icm_pfam_meta_bray, method="average")
plot(icm_pfam_meta_bray_hclust, main = "PFAM MetaGS")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-22-1.png)

``` r
icm_pfam_scg_bray_hclust <-hclust(icm_pfam_scg_bray, method="average")
plot(icm_pfam_scg_bray_hclust, main = "PFAM SCG")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-22-2.png)

PFAM MetaGS: Se pueden observar 2 grupos principales, divididos en otros dos grupos principales cada uno, pero a partir de ahi las divisiones son muy cortas, familias proximas.

PFAM SCG: En este caso, despues de la division principal, en uno de los grupos se puede observar que las relaciones son muy cercanas.

En ambos casos concuerda con los clusters poco definidos, que implican familias próximas.

## Conclusion EBI vs ICM

En ambos casos los NMDs coinciden con los dendogramas, son coherentes. Por otra parte, los datos obtenidos por el ICM tienen mucha mas capacidad para clasificar en diferentes familias o genes y manteniendo la coherencia entre NMDS y dendogramas.

Por tanto, considero que el tratamiento/pipeline utilizada en el ICM, con la finalidad de clasificación, es bastante más útil que la usada en el IBM, Mgnify.

# OTU tables MiTags

``` r
otu_mitags <- (read.delim("functional.tables/otu_table97.txt", header = TRUE, sep = "\t"))
otu_mitags$OTUId <- NULL
otu_mitags <- t(otu_mitags)
```

## Richness

``` r
otu_mitags_ss <- rrarefy(otu_mitags, min(rowSums(otu_mitags)))
richness <- estimateR(otu_mitags_ss)
richness
```

    ##          ERR2098365 ERR2098366 ERR2098367 ERR2098368 ERR2098369 ERR2098370
    ## S.obs    3417.00000 3327.00000 3320.00000  3152.0000 3214.00000 3035.00000
    ## S.chao1  5841.12429 5338.18068 5229.96625  5004.2059 4974.88136 4893.97655
    ## se.chao1  167.70276  141.51874  135.39392   136.2101  129.08239  140.29894
    ## S.ACE    5759.36504 5412.95297 5290.44311  4996.0848 4941.25272 4794.66249
    ## se.ACE     44.01762   42.50672   40.89097    39.6663   38.92132   38.91433
    ##          ERR2098371 ERR2098372 ERR2098373 ERR2098374 ERR2098375 ERR2098376
    ## S.obs    3106.00000 3118.00000 3008.00000 3038.00000 3173.00000 2907.00000
    ## S.chao1  4708.06157 4920.22332 4616.15789 4710.48155 4592.27778 4413.57092
    ## se.chao1  119.24140  133.53746  122.64906  124.98564  104.93939  112.18310
    ## S.ACE    4786.79213 4895.48293 4619.50972 4756.49735 4696.90633 4591.27497
    ## se.ACE     38.21722   39.00766   37.66176   38.57309   37.22395   38.56205
    ##          ERR2098377 ERR2098378 ERR2098379 ERR2098380 ERR2098381 ERR2098382
    ## S.obs    3109.00000  818.00000 2869.00000 2728.00000 3322.00000 3105.00000
    ## S.chao1  4777.89885 1979.21212 4459.66667 3896.33906 4886.67980 4499.15608
    ## se.chao1  124.10224  160.94758  125.41579   96.57242  112.58084  105.78399
    ## S.ACE    4867.90686 2067.47768 4403.15896 3923.60164 5021.74264 4524.32847
    ## se.ACE     39.60607   30.00777   36.36686   33.04346   38.82049   36.34368
    ##          ERR2098383 ERR2098384 ERR2098385 ERR2098386 ERR2098387 ERR2098388
    ## S.obs    3082.00000 2699.00000 1876.00000 1212.00000  553.00000 2918.00000
    ## S.chao1  4839.57143 3619.52601 3753.63323 2702.19780 1250.71429 4679.57456
    ## se.chao1  128.81029   77.58440  160.57726  160.29248  121.98546  135.45706
    ## S.ACE    4880.66091 3727.67804 4081.73675 3040.35679 1256.43690 4606.02075
    ## se.ACE     40.22799   31.81716   42.57281   37.46629   22.70129   37.36794
    ##          ERR2098389 ERR2098390 ERR2098391 ERR2098392 ERR2098393 ERR2098394
    ## S.obs    2832.00000 2786.00000 3624.00000 3329.00000 3689.00000 1281.00000
    ## S.chao1  4563.26941 4418.88732 5527.52761 4645.28788 5632.81138 3242.17714
    ## se.chao1  135.27068  130.08921  129.05317   98.85645  130.28101  205.89965
    ## S.ACE    4551.02440 4319.86434 5610.27291 4656.96154 5790.62260 3661.15290
    ## se.ACE     38.07818   36.05885   41.64987   35.25801   42.55709   42.02512
    ##          ERR2098395 ERR2098396 ERR2098397 ERR2098398 ERR2098399 ERR2098400
    ## S.obs    1151.00000 1069.00000 1893.00000 2683.00000 2808.00000 2800.00000
    ## S.chao1  2509.20690 2221.23494 3456.76033 4158.11084 4602.41606 4336.47477
    ## se.chao1  150.46554  133.02279  132.15932  120.72845  142.19370  122.88315
    ## S.ACE    2769.16244 2475.49592 3888.90618 4114.61731 4430.02096 4254.98139
    ## se.ACE     35.47151   32.96643   41.89195   35.37764   36.99592   35.55173
    ##          ERR2098401 ERR2098402 ERR2098403 ERR2098404 ERR2098405 ERR2098406
    ## S.obs     2786.0000 3542.00000 3584.00000 3735.00000 3655.00000  973.00000
    ## S.chao1   4489.0870 5336.05714 5483.66453 5743.15315 5852.03125 2107.32847
    ## se.chao1   135.8913  124.34151  130.64594  133.86566  148.11920  140.34932
    ## S.ACE     4336.1389 5402.75354 5525.55153 5800.60884 5820.33795 2250.72292
    ## se.ACE      35.9859   40.63569   40.81013   42.00796   42.49898   30.45351
    ##          ERR2098407 ERR2098408 ERR2098409 ERR2098410 ERR2098411 ERR2196983
    ## S.obs    1009.00000 2516.00000 2535.00000 3499.00000 3525.00000 1372.00000
    ## S.chao1  2289.21088 3843.65789 3930.51724 5302.11981 5330.81833 3236.08000
    ## se.chao1  151.81609  113.91522  116.19694  125.08768  126.16649  187.76793
    ## S.ACE    2638.43388 3762.80666 4001.26246 5338.00606 5317.64927 3481.99910
    ## se.ACE     35.05884   33.10677   35.47152   40.57472   39.37939   37.36963
    ##          ERR2196984 ERR2196985
    ## S.obs    2821.00000 3728.00000
    ## S.chao1  4406.43267 5634.08715
    ## se.chao1  124.55763  127.72516
    ## S.ACE    4399.67769 5718.95671
    ## se.ACE     36.54432   42.01198

## Rarefaction

``` r
#otu_rare <- rarefaction(otu_mitags_ss)
```

No he sido capaz de correr este plot.

## Rarecurve

``` r
#otu_rcurve <- rarecurve(otu_mitags_ss)
```

## Acumulation curve

``` r
#otu_acc <- specaccum(otu_mitags_ss, method = "exact", permutations = 1000)
#plot(otu_acc)
```

## Eveness

``` r
plot(radfit(colSums(otu_mitags_ss)))
```

    ## Error in glm.fit(x = structure(c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,  : 
    ##   NA/NaN/Inf in 'x'

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-28-1.png)

## Preston fit

``` r
otu_preston<-prestonfit(colSums(otu_mitags_ss))
otu_prestondistr<-prestondistr(colSums(otu_mitags_ss))
plot(otu_preston)
lines(otu_prestondistr, line.col="blue3")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-29-1.png)

## NMDS

``` r
otu_mitags_ss_bray <- vegdist(otu_mitags_ss, method="bray")
otu_mitags_ss_bray.nmds <- monoMDS(otu_mitags_ss_bray)

otu_mitags_ss_bray.nmds
```

    ## 
    ## Call:
    ## monoMDS(dist = otu_mitags_ss_bray) 
    ## 
    ## Non-metric Multidimensional Scaling
    ## 
    ## 50 points, dissimilarity 'bray', call 'vegdist(x = otu_mitags_ss, method = "bray")'
    ## 
    ## Dimensions: 2 
    ## Stress:     0.08454177 
    ## Stress type 1, weak ties
    ## Scores scaled to unit root mean square, rotated to principal components
    ## Stopped after 73 iterations: Stress nearly unchanged (ratio > sratmax)

``` r
plot(otu_mitags_ss_bray.nmds, main = "MiTags")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-30-1.png)

``` r
stressplot(otu_mitags_ss_bray.nmds, main = "MiTags")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-30-2.png)

## UPGMA

``` r
otu_mitags_ss_hclust <- hclust(otu_mitags_ss_bray, "average")
plot(otu_mitags_ss_hclust)
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-31-1.png)

``` r
otu_mitags_ss_hclust
```

    ## 
    ## Call:
    ## hclust(d = otu_mitags_ss_bray, method = "average")
    ## 
    ## Cluster method   : average 
    ## Distance         : bray 
    ## Number of objects: 50

## DCA

``` r
otu_dca <- decorana(otu_mitags_ss)
```

    ## Warning in decorana(otu_mitags_ss): some species were removed because they were
    ## missing in the data

``` r
plot(otu_dca, display = "sites")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-32-1.png)

``` r
plot(otu_dca, display = "species")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-32-2.png)

## Mantel test

``` r
otu_mitags_bray <- vegdist(otu_mitags, method="bray")
plot(otu_mitags_ss_bray,otu_mitags_bray)
abline(0,1, col = "red")
```

![](metagenomic_analysis_files/figure-markdown_github/unnamed-chunk-33-1.png)

``` r
mantel(otu_mitags_ss_bray,otu_mitags_bray)
```

    ## 
    ## Mantel statistic based on Pearson's product-moment correlation 
    ## 
    ## Call:
    ## mantel(xdis = otu_mitags_ss_bray, ydis = otu_mitags_bray) 
    ## 
    ## Mantel statistic r: 0.9757 
    ##       Significance: 0.001 
    ## 
    ## Upper quantiles of permutations (null model):
    ##    90%    95%  97.5%    99% 
    ## 0.0814 0.1173 0.1488 0.1687 
    ## Permutation: free
    ## Number of permutations: 999

Una buena correlación implica que el subsampling es correcto y que no afecta negativamente al analysis.

## Conclusión ITags vs MiTags

Comparando los plots obtenidos por MiTags e ITags (16S), puedo ver diferencias, pero no sabria interpretarlas.

Ahora bien, leyendo y utilizando el paper "Metagenomic 16S rDNA Illumina tags are a powerfulalternative to amplicon sequencing to explore diversityand structure of microbial communities", se puede deducir que usando MiTags la detección de OTUs es mas óptima.
