# Round 1

## dev
```
> select.features.regularized(dev)
$lasso.features
 [1] "x1"   "x3"   "x7"   "x25"  "x26"  "x27"  "x28"  "x29"  "x30"  "x33" 
[11] "x34"  "x87"  "x88"  "x92"  "x97"  "x143" "x172" "x194" "x197"

$mcp.features
(Intercept)          x1          x2          x3          x7         x25 
          1           2           3           4           8          26 
        x26         x29         x30         x34 
         27          30          31          35 

$scad.features
(Intercept)          x1          x2          x3          x7         x25 
          1           2           3           4           8          26 
        x26         x27         x28         x29         x30         x31 
         27          28          29          30          31          32 
        x34         x88         x97        x143 
         35          89          98         144 


> check.importance(dev, xgb.params, xgb.nrounds)
[1] "SOIL important features"
          x1          x34           x3          x27           x7 
1.000000e+00 1.000000e+00 1.000000e+00 7.635998e-01 4.895164e-01 
         x25           x2          x26          x30          x29 
4.289833e-01 3.047652e-01 8.689867e-02 6.624677e-02 1.798096e-02 
         x28         x143          x97          x88          x92 
1.384403e-02 1.980348e-03 4.063207e-04 3.587983e-05 6.299630e-07 
        x175          x31          x87          x33         x194 
6.299619e-07 6.180155e-08 5.313724e-08 2.297487e-11 2.297487e-11 
        x197         x172          x65         x150         x124 
2.297487e-11 5.516044e-16 1.005230e-17 1.005230e-17 3.588327e-19 
         x84         x100          x16         x106         x193 
3.291517e-21 3.291517e-21 1.972826e-25 1.972826e-25 1.972826e-25 

```

# Round 2

## dev2

```
> select.features.abic(dev2)
$aic.mod

Call:
lm(formula = y ~ x1 + x2 + x3 + x4 + x7 + x25 + x26 + x28 + x29 + 
    x30 + x33 + x34 + x35 + x92 + x97 + x143 + x175 + x197, data = data)

$bic.mod

Call:
lm(formula = y ~ x1 + x2 + x3 + x7 + x25 + x26 + x28 + x29 + 
    x34 + x35 + x97 + x143 + x175, data = data)


> select.features.regularized(dev2)
$lasso.features
 [1] "x1"   "x2"   "x3"   "x4"   "x7"   "x22"  "x25"  "x26"  "x27"  "x28" 
[11] "x29"  "x30"  "x31"  "x33"  "x34"  "x35"  "x39"  "x44"  "x87"  "x88" 
[21] "x92"  "x97"  "x103" "x142" "x143" "x172" "x175" "x197"

$mcp.features
(Intercept)          x1          x2          x3          x4          x7 
          1           2           3           4           5           6 
        x22         x25         x26         x28         x29         x30 
          7           9          10          12          13          14 
        x33         x34         x35         x44         x88         x92 
         17          18          19          21          23          24 
        x97        x103        x143        x172        x175        x194 
         25          26          28          30          31          33 

$scad.features
(Intercept)          x1          x2          x3          x4          x7 
          1           2           3           4           5           6 
        x22         x25         x26         x27         x28         x29 
          7           9          10          11          12          13 
        x30         x33         x34         x35         x44         x87 
         14          17          18          19          21          22 
        x88         x92         x97        x103        x143        x172 
         23          24          25          26          28          30 
       x175        x197 
         31          34 


rf.top10: 1, 34, 3, 2, 4, 33, 171, 35, 88, 30, 44
xgb.top: 1, 34, 3, 4, 2, 171, 7, 87, 27, 33, 26, 32, 88, 172, 142, 193
```

## dev3

```
> select.features.abic(dev3)
$aic.mod

Call:
lm(formula = y ~ x1 + x2 + x3 + x4 + x7 + x25 + x26 + x28 + x29 + 
    x30 + x33 + x34 + x35 + x92 + x97 + x143 + x175 + x197, data = data)

$bic.mod

Call:
lm(formula = y ~ x1 + x2 + x3 + x7 + x25 + x26 + x28 + x29 + 
    x34 + x35 + x97 + x143 + x175, data = data)


> select.features.regularized(dev3)
$lasso.features
 [1] "x1"   "x2"   "x3"   "x4"   "x7"   "x22"  "x25"  "x26"  "x27"  "x28" 
[11] "x29"  "x30"  "x31"  "x33"  "x34"  "x35"  "x39"  "x44"  "x87"  "x88" 
[21] "x92"  "x97"  "x103" "x142" "x143" "x172" "x175" "x197"

$mcp.features
(Intercept)          x1          x2          x3          x4          x7 
          1           2           3           4           5           6 
        x22         x25         x26         x28         x29         x30 
          7           8           9          11          12          13 
        x33         x34         x35         x44         x88         x92 
         16          17          18          20          22          23 
        x97        x103        x143        x172        x175        x194 
         24          25          27          29          30          32 

$scad.features
(Intercept)          x1          x2          x3          x4          x7 
          1           2           3           4           5           6 
        x22         x25         x26         x27         x28         x29 
          7           8           9          10          11          12 
        x30         x31         x32         x33         x34         x35 
         13          14          15          16          17          18 
        x39         x44         x87         x88         x92         x97 
         19          20          21          22          23          24 
       x103        x142        x143        x171        x172        x175 
         25          26          27          28          29          30 
       x194        x197 
         32          33 


> check.importance(dev3, xgb.params, xgb.nrounds)
[1] "SOIL important features"
          x1          x34           x3          x25           x7 
1.000000e+00 1.000000e+00 1.000000e+00 9.275635e-01 8.975719e-01 
          x2          x27          x26          x29          x30 
7.438709e-01 5.648676e-01 5.411369e-01 3.930255e-01 3.028860e-01 
         x28          x97         x143          x88          x92 
1.268559e-01 1.072390e-01 4.409417e-02 3.879720e-02 2.333077e-02 
        x175          x31          x33         x194         x197 
2.083556e-02 1.183121e-02 6.379138e-03 4.190479e-03 1.126225e-03 
          x4          x35         x103          x87          x44 
7.692540e-04 5.581687e-04 4.367969e-04 3.156226e-04 2.444442e-04 
        x172          x22          x39         x142         x171 
1.538015e-04 1.401779e-04 6.586764e-06 1.366070e-06 1.491740e-07 



rf.top10: 1, 34, 3, 2, 171, 4, 35, 33, 27, 88, 7, 44, 30
xgb.top: 1, 34, 3, 4, 32, 27, 7, 33, 87, 26, 30, 175, 28, 92, 88
```

## dev4
```
> select.features.abic(dev4)
$aic.mod

Call:
lm(formula = y ~ x1 + x2 + x3 + x4 + x7 + x25 + x26 + x28 + x29 + 
    x30 + x33 + x34 + x35 + x92 + x97 + x143 + x175 + x197, data = data)

$bic.mod

Call:
lm(formula = y ~ x1 + x2 + x3 + x7 + x25 + x26 + x28 + x29 + 
    x34 + x35 + x97 + x143 + x175, data = data)

> select.features.regularized(dev4)
$lasso.features
 [1] "x1"   "x2"   "x3"   "x4"   "x7"   "x22"  "x25"  "x26"  "x27"  "x28" 
[11] "x29"  "x30"  "x31"  "x33"  "x34"  "x35"  "x39"  "x44"  "x87"  "x88" 
[21] "x92"  "x97"  "x103" "x142" "x143" "x172" "x175" "x197"

$mcp.features
(Intercept)          x1          x2          x3          x4          x7 
          1           2           3           4           5           6 
        x22         x25         x26         x27         x28         x29 
          7           8           9          10          11          12 
        x30         x31         x32         x33         x34         x35 
         13          14          15          16          17          18 
        x39         x44         x87         x88         x92         x97 
         19          20          21          22          23          24 
       x103        x142        x143        x171        x172        x175 
         25          26          27          28          29          30 
       x194        x197 
         31          32 

$scad.features
(Intercept)          x1          x2          x3          x4          x7 
          1           2           3           4           5           6 
        x22         x25         x26         x27         x28         x29 
          7           8           9          10          11          12 
        x30         x31         x32         x33         x34         x35 
         13          14          15          16          17          18 
        x39         x44         x87         x88         x92         x97 
         19          20          21          22          23          24 
       x103        x142        x143        x171        x172        x175 
         25          26          27          28          29          30 
       x194        x197 
         31          32 

```

