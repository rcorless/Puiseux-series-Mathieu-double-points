# Puiseux-series-Mathieu-double-points

This has turned into a repository for code related to solving the Mathieu equations, not just for double points.
-- ActiveDoublesc2p0.maple is a Maple Workbook that solves one particular case of hemodynamic flow with a double eigenvalue.  Code is intended for "supervised use"
-- HemodynamicsTalk.pdf are the slides for my talk on the subject, first given at U. Waterloo Tuesday January 17 2023
-- Activec1p0.maple is a Maple Workbook that solves several cases where the circumference is 1.0cm and plots various fluid mechanical quantities.  Has version 1.01
   of FlowQuantities.mpl tucked inside it.

The file BlanchClemmDoublePoints.maple is a Maple Workbook which contains code that will allow the user to locate double points to arbitrary precision and to compute Puiseux series expansions of the eigenvalue about those double points, using Newton's method in series.  This method is described in the paper Brimacombe, Corless, and Zamir: Computation of the Mathieu Functions, a historical perspective.https://arxiv.org/abs/2008.01812

The Workbook also contains my code to evaluate Mathieu eigenvalues/eigenfunctions, and to solve the Mathieu equation more generally by a Hermite-Obreschkoff method.  That code is not meant for automatic use but rather in a human-supervised way; currently that usage is not described.  I have now added the workbook MathieuCodeDemo to sketch out how to use it.

The file is intended to be used with Maple version 2020 or later, but might be usable with earlier versions.  If you do not have Maple, you can at least read the file by downloading the free Maple Player.

A similar file called MathieuTalk.maple contains the slides (and code etc) for my upcoming talk on Mathieu functions at Newcastle.

This repo also contains a copy of Robert Moir's translation to English of Mathieu's original 1868 Memoire.

The BigDoublesTable.csv file contains the 32 digit precision coefficients for the Puiseux series about each of the Blanch and Clemm double points.

The columns are m, qstar, astar, alpha1, alpha2, ..., alpha6

The Puiseux expansion of a(q) about qstar is a(q) = astar + sum( alpha.k*(q-qstar)^(k/2), k=1..6) + ... 
