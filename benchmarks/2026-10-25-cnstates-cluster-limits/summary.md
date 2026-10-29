# cnstates benchmark: /data1/shahs3/users/myersm2/hatchet-cnstates-bench

## lowhet, mode 0

120 clusters (4 hot, 8 background), 12 samples, 3 true tumor clones, diploid, D = 2.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 3702 s | 2.21 GB | 814 | 2446.3 ms | 11895 / 1678 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| free_old | 0 | 3723 s | 2.24 GB | 814 | 2498.8 ms | n/a | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d | 0 | 4324 s | 2.16 GB | 752 | 3821.6 ms | 17250 / 2749 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d_exempt | 0 | 3837 s | 2.34 GB | 735 | 3447.6 ms | 17070 / 2713 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |
| d_exempt_partial | 0 | 3877 s | 2.17 GB | 740 | 3627.7 ms | 17160 / 2731 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| map_oracle | 0 | 4085 s | 2.31 GB | 484 | 371.4 ms | 15170 / 2333 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |
| d_exempt_pairwise | 0 | 5324 s | 2.19 GB | 738 | 13936.5 ms | 52665 / 7148 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |
| d1_exempt | 0 | 3252 s | 2.30 GB | 427 | 242.3 ms | 14770 / 2253 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d3_exempt | 0 | 4007 s | 2.22 GB | 777 | 4380.7 ms | 18795 / 3058 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | d_exempt | d_exempt_partial | map_oracle | d_exempt_pairwise | d1_exempt | d3_exempt |
|---|---|---|---|---|---|---|---|---|---|---|
| diploid | 2 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 |
| diploid | 3 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 48.4067 | 44.1722 |
| diploid | 4 | 37.5021 | 37.5021 | 41.0691 | 37.5021 | 38.6204 | 37.5021 | 37.5021 | 46.1666 | 37.5021 |
| diploid | 5 | 36.8188 | 36.8188 | 40.5592 | 37.3629 | 38.4820 | 37.3645 | 37.3629 | 46.0793 | 37.1543 |
| diploid | 6 | 36.1815 | 36.1815 | 39.4508 | 36.9480 | 38.0584 | 37.3602 | 36.9480 | 46.0850 | 36.0150 |

## lowhet, mode 2

120 clusters (4 hot, 8 background), 12 samples, 3 true tumor clones, diploid, D = 2.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 632 s | 0.73 GB | 814 | 2122.7 ms | 11895 / 1678 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| free_old | 0 | 12847 s | 0.74 GB | 814 | 2435.8 ms | n/a | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d | 0 | 912 s | 0.81 GB | 746 | 3758.2 ms | 17250 / 2749 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d_exempt | 0 | 840 s | 0.83 GB | 733 | 3490.1 ms | 17070 / 2713 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |
| d_exempt_partial | 0 | 861 s | 0.80 GB | 742 | 3659.0 ms | 17160 / 2731 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| map_oracle | 0 | 264 s | 0.47 GB | 482 | 364.0 ms | 15170 / 2333 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |
| d_exempt_pairwise | 0 | 2142 s | 1.65 GB | 736 | 13153.4 ms | 52665 / 7148 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |
| d1_exempt | 0 | 238 s | 0.46 GB | 427 | 244.4 ms | 14770 / 2253 | diploid.n2 | 1 | 0.0160 | 0.003 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d3_exempt | 0 | 962 s | 0.85 GB | 768 | 4528.0 ms | 18795 / 3058 | diploid.n4 | 3 | 0.0027 | 0.004 | 0.000 | 0.000 | 3.00 / 3.00 | 1.000 | 0 / 0 |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | d_exempt | d_exempt_partial | map_oracle | d_exempt_pairwise | d1_exempt | d3_exempt |
|---|---|---|---|---|---|---|---|---|---|---|
| diploid | 2 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 | 52.9956 |
| diploid | 3 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 44.1722 | 48.4067 | 44.1722 |
| diploid | 4 | 37.5021 | 37.5021 | 41.0691 | 37.5021 | 38.6204 | 37.5021 | 37.5021 | 46.1666 | 37.5021 |
| diploid | 5 | 36.8188 | 36.8188 | 40.5592 | 37.3629 | 38.4820 | 37.3646 | 37.3629 | 46.0793 | 37.1543 |
| diploid | 6 | 36.1815 | 36.1815 | 39.4508 | 36.9480 | 38.0584 | 37.3602 | 36.9480 | 46.0850 | 36.0150 |

## lowhet_noisy, mode 2

120 clusters (4 hot, 8 background), 12 samples, 3 true tumor clones, diploid, D = 2.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 981 s | 0.76 GB | 905 | 3656.8 ms | 11895 / 1678 | diploid.n2 | 1 | 0.0188 | 0.007 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| free_old | 0 | 981 s | 0.76 GB | 905 | 3677.3 ms | n/a | diploid.n2 | 1 | 0.0188 | 0.007 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d | 0 | 1164 s | 0.82 GB | 802 | 4488.3 ms | 17250 / 2749 | diploid.n2 | 1 | 0.0188 | 0.007 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d_exempt | 0 | 1097 s | 0.82 GB | 823 | 4418.1 ms | 17070 / 2713 | diploid.n2 | 1 | 0.0188 | 0.007 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| d_exempt_partial | 0 | 1103 s | 0.81 GB | 820 | 4431.9 ms | 17160 / 2731 | diploid.n2 | 1 | 0.0188 | 0.007 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |
| map_oracle | 0 | 248 s | 0.47 GB | 553 | 300.7 ms | 15170 / 2333 | diploid.n2 | 1 | 0.0188 | 0.007 | 0.000 | 0.056 | 1.00 / 3.00 | 1.000 | 0 / 2 |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | d_exempt | d_exempt_partial | map_oracle |
|---|---|---|---|---|---|---|---|
| diploid | 2 | 89.4630 | 89.4630 | 89.4630 | 89.4630 | 89.4630 | 89.4630 |
| diploid | 3 | 80.5649 | 80.5649 | 80.5649 | 80.5649 | 80.5649 | 80.5649 |
| diploid | 4 | 74.2573 | 74.2573 | 76.2107 | 74.9362 | 75.9358 | 74.9362 |
| diploid | 5 | 71.5870 | 71.5870 | 74.9712 | 72.3421 | 74.3733 | 74.6263 |
| diploid | 6 | 70.7931 | 70.7931 | 74.5857 | 71.3564 | 73.0612 | 74.4683 |

## moderate, mode 2

150 clusters (15 hot, 45 background), 10 samples, 4 true tumor clones, diploid, D = 2.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 784 s | 0.71 GB | 902 | 2504.8 ms | 12945 / 2098 | diploid.n3 | 2 | 0.0354 | 0.004 | 0.000 | 0.164 | 2.00 / 3.67 | 1.000 | 0 / 2 |
| free_old | 0 | 784 s | 0.70 GB | 902 | 2504.2 ms | n/a | diploid.n3 | 2 | 0.0354 | 0.004 | 0.000 | 0.164 | 2.00 / 3.67 | 1.000 | 0 / 2 |
| d | 0 | 1077 s | 0.91 GB | 819 | 3863.2 ms | 19650 / 3439 | diploid.n3 | 2 | 0.0354 | 0.004 | 0.000 | 0.164 | 2.00 / 3.67 | 1.000 | 0 / 2 |
| d_exempt | 0 | 925 s | 0.85 GB | 901 | 3103.5 ms | 18975 / 3304 | diploid.n5 | 4 | 0.0054 | 0.005 | 0.000 | 0.009 | 3.53 / 3.67 | 0.998 | 0 / 0 |
| d_exempt_partial | 0 | 997 s | 0.86 GB | 849 | 3517.4 ms | 19290 / 3367 | diploid.n3 | 2 | 0.0354 | 0.004 | 0.000 | 0.164 | 2.00 / 3.67 | 1.000 | 0 / 2 |
| map_oracle | 0 | 392 s | 0.58 GB | 674 | 1025.3 ms | 18195 / 3148 | diploid.n5 | 4 | 0.0054 | 0.005 | 0.000 | 0.009 | 3.53 / 3.67 | 0.998 | 0 / 0 |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | d_exempt | d_exempt_partial | map_oracle |
|---|---|---|---|---|---|---|---|
| diploid | 2 | 115.2699 | 115.2699 | 115.2699 | 115.2699 | 115.2699 | 115.2699 |
| diploid | 3 | 66.0165 | 66.0165 | 66.0165 | 66.0165 | 66.0165 | 66.0165 |
| diploid | 4 | 49.7781 | 49.7781 | 54.6024 | 49.7781 | 51.3389 | 49.7781 |
| diploid | 5 | 40.9749 | 40.9749 | 49.5780 | 38.4659 | 43.3371 | 38.4659 |
| diploid | 6 | 37.5302 | 37.5302 | 49.5780 | 37.7069 | 43.3371 | 38.1524 |

## nohot, mode 2

120 clusters (0 hot, 18 background), 10 samples, 3 true tumor clones, diploid, D = 2.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 665 s | 0.64 GB | 818 | 2225.8 ms | 10445 / 1678 | diploid.n4 | 3 | 0.0030 | 0.004 | 0.000 | 0.000 | n/a | 1.000 | 0 / 0 |
| free_old | 0 | 624 s | 0.63 GB | 818 | 2271.6 ms | n/a | diploid.n4 | 3 | 0.0030 | 0.004 | 0.000 | 0.000 | n/a | 1.000 | 0 / 0 |
| d | 0 | 761 s | 0.78 GB | 724 | 3056.2 ms | 15800 / 2749 | diploid.n4 | 3 | 0.0030 | 0.004 | 0.000 | 0.000 | n/a | 1.000 | 0 / 0 |
| map_oracle | 0 | 233 s | 0.44 GB | 524 | 357.4 ms | 13780 / 2345 | diploid.n4 | 3 | 0.0030 | 0.004 | 0.000 | 0.000 | n/a | 1.000 | 0 / 0 |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | map_oracle |
|---|---|---|---|---|---|
| diploid | 2 | 49.0722 | 49.0722 | 49.0722 | 49.0722 |
| diploid | 3 | 36.5581 | 36.5581 | 36.5581 | 36.5581 |
| diploid | 4 | 30.8132 | 30.8132 | 30.8132 | 30.8132 |
| diploid | 5 | 30.3090 | 30.3090 | 30.4351 | 30.6973 |
| diploid | 6 | 28.8774 | 28.8774 | 30.1794 | 30.5379 |

## scale, mode 0

300 clusters (15 hot, 45 background), 20 samples, 5 true tumor clones, diploid, D = 3.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |
| free_old | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |
| d | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |
| d_exempt | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |
| d_exempt_partial | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |
| map_oracle | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |
| d_exempt_pairwise | ? | n/a s | n/a GB | n/a | n/a ms | n/a | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |

## scale, mode 2

300 clusters (15 hot, 45 background), 20 samples, 5 true tumor clones, diploid, D = 3.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 7417 s | 2.05 GB | 965 | 13961.5 ms | 46598 / 5398 | diploid.n3 | 2 | 0.0237 | 0.003 | 0.000 | 0.110 | 2.00 / 4.07 | 0.969 | 0 / 3 |
| free_old | 0 | 6986 s | 2.02 GB | 965 | 12909.4 ms | n/a | diploid.n3 | 2 | 0.0237 | 0.003 | 0.000 | 0.110 | 2.00 / 4.07 | 0.969 | 0 / 3 |
| d | 0 | 11724 s | 2.60 GB | 960 | 29543.8 ms | 73508 / 10780 | diploid.n3 | 2 | 0.0237 | 0.003 | 0.000 | 0.110 | 2.00 / 4.07 | 0.969 | 0 / 3 |
| d_exempt | 0 | 11224 s | 2.52 GB | 959 | 28594.8 ms | 72158 / 10510 | diploid.n3 | 2 | 0.0237 | 0.003 | 0.000 | 0.110 | 2.00 / 4.07 | 0.969 | 0 / 3 |
| d_exempt_partial | 0 | 11586 s | 2.61 GB | 961 | 29321.9 ms | 72788 / 10636 | diploid.n3 | 2 | 0.0237 | 0.003 | 0.000 | 0.110 | 2.00 / 4.07 | 0.969 | 0 / 3 |
| map_oracle | 0 | 2509 s | 1.52 GB | 728 | 2370.6 ms | 60028 / 8084 | diploid.n6 | 5 | 0.0037 | 0.003 | 0.000 | 0.006 | 3.93 / 4.07 | 0.999 | 0 / 0 |
| d_exempt_pairwise | ? | n/a s | n/a GB | 199 | 35008.1 ms | 192294 / 19258 | none | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a / n/a |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | d_exempt | d_exempt_partial | map_oracle | d_exempt_pairwise |
|---|---|---|---|---|---|---|---|---|
| diploid | 2 | 144.2414 | 144.2414 | 144.2414 | 144.2414 | 144.2414 | 144.2414 | 144.2414 |
| diploid | 3 | 111.5411 | 111.5411 | 111.5411 | 111.5411 | 111.5411 | 111.5411 | 111.5411 |
| diploid | 4 | 102.6952 | 102.6952 | 102.6952 | 102.6952 | 102.6952 | 95.7760 | 102.6952 |
| diploid | 5 | 86.9896 | 86.9896 | 89.8343 | 88.9587 | 89.7671 | 84.1558 | 88.9587 |
| diploid | 6 | 78.4102 | 78.4102 | 83.9305 | 78.6559 | 80.6810 | 73.6883 | n/a |
| diploid | 7 | 75.7514 | 75.7514 | 82.0821 | 77.8402 | 78.3116 | 73.5632 | n/a |
| diploid | 8 | 75.0623 | 75.0623 | 79.3423 | 75.9918 | 78.1801 | 73.3556 | n/a |

## wgd, mode 2

120 clusters (4 hot, 8 background), 10 samples, 3 true tumor clones, WGD, D = 2.

| config | exit | wall | peak RSS | C steps | mean C step | mean C rows / discrete | chosen | tumor clones | frac MAE | purity MAE | over | under | hot states | exact CN | extra / missing |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| free | 0 | 577 s | 0.64 GB | 733 | 1402.3 ms | 10355 / 1672 | tetraploid.n2 | 1 | 0.0796 | 0.024 | 0.000 | 0.043 | 1.00 / 3.00 | 0.971 | 0 / 2 |
| free_old | 0 | 576 s | 0.64 GB | 733 | 1400.5 ms | n/a | tetraploid.n2 | 1 | 0.0796 | 0.024 | 0.000 | 0.043 | 1.00 / 3.00 | 0.971 | 0 / 2 |
| d | 0 | 652 s | 0.79 GB | 688 | 1901.7 ms | 15575 / 2716 | tetraploid.n2 | 1 | 0.0796 | 0.024 | 0.000 | 0.043 | 1.00 / 3.00 | 0.971 | 0 / 2 |
| d_exempt | 0 | 1169 s | 0.78 GB | 683 | 1876.7 ms | 15395 / 2680 | tetraploid.n2 | 1 | 0.0796 | 0.024 | 0.000 | 0.043 | 1.00 / 3.00 | 0.971 | 0 / 2 |
| d_exempt_partial | 0 | 644 s | 0.78 GB | 678 | 1842.0 ms | 15485 / 2698 | tetraploid.n2 | 1 | 0.0796 | 0.024 | 0.000 | 0.043 | 1.00 / 3.00 | 0.971 | 0 / 2 |
| map_oracle | 0 | 328 s | 0.52 GB | 634 | 356.3 ms | 13555 / 2312 | tetraploid.n2 | 1 | 0.0796 | 0.024 | 0.000 | 0.043 | 1.00 / 3.00 | 0.971 | 0 / 2 |

Objective per ploidy and clone count (compute-cn counts the normal clone):

| ploidy | clones | free | free_old | d | d_exempt | d_exempt_partial | map_oracle |
|---|---|---|---|---|---|---|---|
| tetraploid | 2 | 65.6162 | 65.6162 | 65.6162 | 65.6162 | 65.6162 | 65.6162 |
| tetraploid | 3 | 58.0146 | 58.0146 | 58.0146 | 58.0146 | 58.0146 | 58.0146 |
| tetraploid | 4 | 53.7875 | 53.7875 | 54.3070 | 54.0648 | 54.2711 | 54.9494 |
| tetraploid | 5 | 52.1525 | 52.1525 | 52.8427 | 52.5202 | 53.0775 | 54.7994 |
| tetraploid | 6 | 51.8529 | 51.8529 | 52.6668 | 51.6880 | 52.3795 | 54.6739 |
