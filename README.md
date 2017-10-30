# CollabAssess
R package CollabAssess for modeling error-prone network data

How to use the CollabAssess package:

Open console, navigate to where the tar.gz file is saved

Type in
R CMD INSTALL CollabAssess_0.1.0.tar.gz

Package should install and you can load it in R using
library(CollabAssess)

There is one main function
CAgibbs()

It's documented and you can see examples using
?CAgibbs

Including producing the tables in the paper, and a form of the threshold plot
threshplot()

I've included data from Year 1 of the Q&C grant (2003-2004) as an example in the package. The adjacency matrix is called smYr1 and the hire list is called smYr1_hire.
