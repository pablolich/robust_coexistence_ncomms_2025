The three files contain Mathematica Notebooks detailing how the equations in Sec. B1 were verified.

As shown in section B1, the expected variance depends on a combination of six terms, each depending in general on the initial guess $y_0$, the right-hand side $b$, and the matrix $M$ (the definition of M depends on the correlation).

For each of these six terms, first we compute the expectation symbolically, starting from a generic matrix. The calculation is performed for different sizes $n$. Then, the results of the calculations are contrasted with the formulas presented in Sec B1, showing perfect agreement.

The calculations are reported for small $n$, and have been verified up to $n = 10$. The maximum size to consider is stored in the first parameter $k$ in each notebook. Note that increasing $k$ increases the computing time (and memory use) exponentially.
