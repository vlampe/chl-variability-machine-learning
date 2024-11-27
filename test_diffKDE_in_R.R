# caluclate diffKDE in R

require(reticulate)

reticulate::py_version()


# use_python("/usr/local/homebrew/Cellar/python@3.9/3.9.17/bin/python3.9") # use the arm x86_64 version!! // not anymore, since i updated R

use_python("/opt/homebrew/Cellar/python@3.11/3.11.5/Frameworks/Python.framework/Versions/3.11/Resources/Python.app/Contents/MacOS/Python")

np      <- import("numpy",   convert=FALSE)
py_csv  <- import("csv",     convert=FALSE)
diffKDE <- import("diffKDE", convert=FALSE)
vals <- c(2,3,4,5,6,6,7,7,7,4,6,7,4,3,3)


diffKDE_rslt <- py_to_r(diffKDE$KDE(np$array(vals), 0, 10)) # diffusion KDE of sample data
diffKDE_rslt2 <- py_to_r(diffKDE$KDE(np$array(vals)))
diffKDE_rslt3 <- py_to_r(diffKDE$KDE(np$array(vals), -3, 60))

plot(diffKDE_rslt[[2]], diffKDE_rslt[[1]], 'l', col="red")
lines(diffKDE_rslt2[[2]], diffKDE_rslt2[[1]], col = "blue")
lines(diffKDE_rslt3[[2]], diffKDE_rslt3[[1]], col = "green")
