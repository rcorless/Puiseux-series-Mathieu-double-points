# Puiseux-series-Mathieu-double-points

The file BlanchClemmDoublePoints.maple is a Maple Workbook which contains code that will allow the user to locate double points to arbitrary precision and to compute Puiseux series expansions of the eigenvalue about those double points, using Newton's method in series.  This method is described in the paper Brimacombe, Corless, and Zamir: Computation of the Mathieu Functions, a historical perspective.

The file is intended to be used with Maple version 2020 or later, but might be usable with earlier versions.  If you do not have Maple, you can at least read the file by downloading the free Maple Player 

The BigDoublesTable.csv file contains the 32 digit precision coefficients for the Puiseux series about each of the Blanch and Clemm double points.

The columns are m, qstar, astar, alpha1, alpha2, ..., alpha6

The Puiseux expansion of a(q) about qstar is a(q) = astar + sum( alpha.k*(q-qstar)^(k/2), k=1..6) + ... 
