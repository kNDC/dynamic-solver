# Dynamic Solver
An R Shiny tool for solving systems of equations and inequalities, for use in macroeconomic planning.

## Workflow
The tool requires two input files stored in the `data` folder: `metadata.json` and `data.csv
* `metadata.json` lays out variables and equations (general inequality constraints are WIP) (see the next section for the syntax);
* `data.csv` contains data for initialising endogenous variables, exogenous variables and parameters;

## Semantic units and supported syntax
1. Equations (eqs) - each is determined by the formula field (eq1$formula).  Supported functions and operators (to be expanded):
   * Elementary algebraic operators: `+`, `-`, `*`, `/`, power operator `^` (`x^a`) + `sqrt(x) ≡ x^0.5`;
   * Logarithmic functions: general logarithm `log(x, a)`, natural logarithm `ln(x) ≡ log(x, exp(1))`, decenary logarithm `lg(x) ≡ log(x, 10)`;
   * Trigonometric functions: sine `sin(x)`, cosine `cos(x)`, tangent `tg(x)`, cotangent `ctg(x)`;
   * Inverse trigonometric functions: arctangent `arctg(x)`, `atan(x)`, arccotangent `arcctg(x)`;
   * Matrix multiplication: `%*%`;
   * Custom operators: first difference `dd(x)`, midpoint `mp(x)` (`mp(x)[i] ≡ (x[i] - x[i-1]) / 2`);
   * Indexing: numeric indices `x[1]`, `x[{1;2;3}]`, endpoint `x[T]`;
   * Named itemisation R style: `$` (`some.env$x`);

2. Endogenous variables (solved for in the model).  Supported properties:
   * Name: the name to show up in graphs and tables;
   * Non-negative: whether the variable is allowed to be negative `nonneg` = `TRUE`/`FALSE`;
   * Initialiser: whether initial values come from the dataset or are calculated straightaway `initialiser` = `TRUE`/`FALSE`;
   * Unit: the variable's unit;

4. Exogenous variables (taken as are from the data file):
   * Name: same as for endogenous;
   * Unit: same as for endogenous;

5. Parameters ((typically) uni-dimensional exogenous variables not set in the data file) (WIP):
   * Name: same as for endogenous;
   * Value: an integer/double `par1 = 10` or a vector `par1 = [1, 2, 3]`;
   * Unit: same as for endogenous;
