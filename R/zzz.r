utils::globalVariables(c("Individual", "X1", "X2", "X3", "X4","Var","Fractions"))

.onAttach <- function(libname, pkgname) {
  cite_info <- utils::citation(pkgname)
  cite_info <- "Bangken Ying, Yao Liu, Jiangshan Lai(2025). upset.hp: An R Package for Visualizing Hierarchical Partitioning and Variance Partitioning Results Using UpSet Plots. Data Intelligence,doi.org/10.3724/2096-7004.di.2025.0200"
  packageStartupMessage("Thank you for using this package! If you use this package in your research, please cite the following references \n")
  packageStartupMessage(cite_info)
}
