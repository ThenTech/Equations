#define _CRT_SECURE_NO_WARNINGS
//#define _MORE_PRECISION		// Uncomment for more precision in print functions.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

//	Epsilon, close to zero for comparison
#define EPS		1e-12

/**
*	Structure representing a root.
*/
typedef struct NUMBER {
	double real; double imag;
} ROOT;

/**
*	Return a new Equation given the degree.
*/
ROOT newRoot(double real, double imag) {
	return (ROOT) { real, imag };
}

/**
*	Print a root in form X+iY.
*/
void toStringRoot(ROOT r) {
#ifdef _MORE_PRECISION
	if		(r.imag == 0.0)	printf("% f", r.real);
	else if (r.real == 0.0)	printf("% fi", r.imag);
	else					printf("% fg %s %fi", r.real, r.imag < 0 ? "-" : "+", fabs(r.imag));
#else
	if		(r.imag == 0.0)	printf("% .8g", r.real);
	else if (r.real == 0.0)	printf("% .8gi", r.imag);
	else					printf("% .8g %s %.8gi", r.real, r.imag < 0 ? "-" : "+", fabs(r.imag));
#endif
}

/**
*	Structure representing an equation.
*/
typedef struct EQUATION {
	double *coef;	//	{lowest to highest degree} :: index==degree of coef
	size_t degree;	//	highest degree => (degree+1) coef
	ROOT *roots;	//	(degree) roots
} EQUATION;

/**
*	Structure representing the quotient after a division.
*/
typedef struct QUOTIENT {
	EQUATION quotient;
	EQUATION remainder;
} QUOTIENT;

/**
*	Return a new Equation given the degree.
*/
EQUATION newEquation(size_t degree) {
	return (EQUATION) { (double*)calloc((degree + 1), sizeof(double)),
						degree, 
						(ROOT*)calloc(degree, sizeof(ROOT))};
}

/**
*	Print an equation.
*/
void toStringEq(EQUATION eq) {
	for (int i = eq.degree; i >= 0; i--) {
#ifdef _MORE_PRECISION
		if (i == eq.degree)				printf("%f", eq.coef[i]);
		else if (eq.coef[i] != 0)		printf(" %s %f", eq.coef[i] < 0 ? "-" : "+", fabs(eq.coef[i]));
#else
		if (i == eq.degree)				printf("%.4g", eq.coef[i]);
		else if (eq.coef[i] != 0)		printf(" %s %.4g", eq.coef[i] < 0 ? "-" : "+", fabs(eq.coef[i]));
#endif
		if ((eq.coef[i] != 0) & (i != 0))	printf("x^%d", i);
	}
	putchar('\n');
}

/*
*	Prints alternate form of equation using roots.
*/
void toStringEqFromSol(EQUATION eq) {
	if (eq.coef[eq.degree] != 1) printf("%g * ", eq.coef[eq.degree]);
	for (int i = 0; i < eq.degree; i++) {
		printf("(x - (");
		toStringRoot(eq.roots[i]);
		printf("))%s", (i < eq.degree - 1 ? " * " : ""));
	}
	putchar('\n');
}

/**
*	Print all equation roots.
*/
void toStringRoots(EQUATION eq) {
	for (int i = 0; i < eq.degree; i++) {
		toStringRoot(eq.roots[i]);
		putchar('\n');
	}
}

/**
*	Returns true if difference is smaller than EPS or are equal.
*/
int fuzzyEquals(double x, double y) {
	return (fabs(x - y) <= EPS) || (x == y);
}

/*
*	Reads the coefficients and stores them in eqCoef.
*/
void readEquation(EQUATION *eq) {
	for (int i = eq->degree, cf = 97; i >= 0; i--, cf++) {
		if (i > 0)	printf("%cx^%d + ", cf, i);
		else		printf("%c = 0\n", cf, i);
	}

	for (int i = eq->degree, cf = 97; i >= 0; i--, cf++) {
		printf("%c = ", cf);
		scanf("%lf", &eq->coef[i]);
	}
	while (!eq->coef[eq->degree]) eq->degree--;		//	lower degree to compensate for leading zero coeff.
}

/*
*	Parses the coefficients from Cmd and stores them in eqCoef.
*/
void readEquationCmd(EQUATION *eq, char **in) {
	for (int i = eq->degree; i >= 0; i--)
		eq->coef[i] = atof(in[eq->degree - i + 1]);	//	Reverse store coefficients, convert string to float with 'atof(*char)'
	while (!eq->coef[eq->degree]) eq->degree--;		//	lower degree to compensate for leading zero coeff.
}

/*
*	Returns the derived equation.
*/
EQUATION deriveEq(EQUATION eq) {
	EQUATION deriv = newEquation(eq.degree - 1);
	for (int i = 0; i < eq.degree; i++)
		deriv.coef[i] = eq.coef[i + 1] * (i + 1);
	return deriv;
}

/*
*	Returns the integrated equation.
*/
EQUATION integrateEq(EQUATION eq) {
	EQUATION integ = newEquation(eq.degree + 1);
	integ.coef[0] = 0;
	for (int i = 0; i < integ.degree; i++)
		integ.coef[i + 1] = eq.coef[i] / (i + 1);
	return integ;
}

/*
*	Returns quotient of division of equation with another equation using:
*		- long division
*		- Horner's Scheme (div.degree == 0)
*/
QUOTIENT divideBy(EQUATION eq, EQUATION div) {
	if (div.degree > eq.degree) { printf("[ERROR] Divider degree to large!"); return (QUOTIENT) { NULL, NULL }; }
	int i, j;
	QUOTIENT q = (QUOTIENT) { newEquation(eq.degree - div.degree), newEquation(eq.degree) };

	if (div.degree == 1) {	// Horner
		div.coef[0] /= div.coef[1];	//	Make sure the divider is in form x+r
		q.quotient.coef[q.quotient.degree] = eq.coef[eq.degree];
		for (i = q.quotient.degree - 1; i >= 0; i--)
			q.quotient.coef[i] = div.coef[0] * q.quotient.coef[i + 1] + eq.coef[i + 1];
		q.remainder.degree = 0;
		q.remainder.coef[0] = div.coef[0] * q.quotient.coef[0] + eq.coef[0];
	} else {				// Long division
		memcpy(q.remainder.coef, eq.coef, sizeof(double) * (eq.degree + 1));

		for (i = q.quotient.degree; i >= 0; i--, q.remainder.degree--) {
			q.quotient.coef[i] = q.remainder.coef[i + div.degree] / div.coef[div.degree];
			for (j = 0; j < div.degree; j++)
				q.remainder.coef[i + j] -= div.coef[j] * q.quotient.coef[i];
		}
	}
	return q;
}

/*
*	Returns the value of the equation with given X using Horner's method.
*		== equation(X)
*/
double getValue(EQUATION eq, double x) {
	double y = 0.0;
	for (int i = eq.degree; i >= 0; i--)
		y = y * x + eq.coef[i];
	return y;
}

/*
*	Returns an array with all dividers of the given number. 
*/
long *getFactors(long n) {
	int l = log(abs(n)) * 9 + 10;
	long *facts = malloc(l * sizeof(long));
	for (int i = 0; i < l; facts[i++] = 0);
	facts[0] = 1; facts[1] = n;

	for (long i = 2; i * i <= n; i++)
		if (n % i == 0) {
			facts[i] = i;
			if (i != n / i)	facts[i++ + 1] = n / i;
		}
	return facts;
}

/*
*	Factorising
*		With integer coefficients, roots often are equal to:
*			+- (Any divider of highest degree coeff) / (Any divider of lowest degree coeff)
*/
double getXFactorising(EQUATION eq) {
	long *factA = getFactors(eq.coef[eq.degree]);
	long *factD = getFactors(eq.coef[0]);

	for (long i = 0, a; (a = factA[i]) != 0; i++)
		for (long j = 0, d; (d = factD[j]) != 0; j++) {
			double q = (double) d / a;
			if (fuzzyEquals(getValue(eq, q), 0.0))	return q;
			if (fuzzyEquals(getValue(eq, -q), 0.0))	return -q;
		}
	return 0.0;
}

/*
*	Newton-Raphson Method
*		Approximate a root by continuously subtracting the derivative value.
*/
double getXNewtonRaphson(EQUATION eq) {
	EQUATION deriv = deriveEq(eq);
	double x = 1, f, d;
	for (int i = 0; i < 1000; i++) {
		f = getValue(eq, x);
		d = getValue(deriv, x);
		if (fuzzyEquals(f, 0.0))	return x;
		if (!fuzzyEquals(d, 0.0))	x -= (f / d);
		else						x++;
	}
	return x;
}

/*
*	Range Halving method
*		Takes the average value in a range until range limits are equal thus giving a root.
*/
double getXRange(EQUATION eq, double a, double b) {
	double mid = (a + b) / 2.0;	//Midpoint
	if (getValue(eq, mid) == 0.0 || fuzzyEquals(a, b))
		return mid;
	if (getValue(eq, a) * getValue(eq, mid) < 0.0)
		return getXRange(eq, a, mid);
	return getXRange(eq, mid, b);
}

/*
*	Call function for Range Halving method.
*/
double getXRangeHalving(EQUATION eq) {
	return getXRange(eq, -100.0, 100.0);
}

/*
*	Set root if actual zero.
*/
int setSol(EQUATION *eq, size_t rootNmbr, double root) {
	if (fuzzyEquals(getValue(*eq, root), 0.0)) {
		eq->roots[rootNmbr] = newRoot(root, 0);
		return 1;
	}
	return 0;
}

/*
*	Solve and find the roots of a given equation.
*		eq		:	The start equations to add roots to.
*		div		:	The new equation with the previous root divided away.
*		rootNmbr:	The amount of roots found so far (as index for root array).
*	Returns true (1) if found a solution
*/
int solveEqRecur(EQUATION *eq, EQUATION div, size_t rootNmbr) {
	if (div.degree < 1)
		{ printf("No solutions, invalid equation.\n"); return 0; }
	else if (div.degree < 2)			//	ONE_S : 1 single solutions
		eq->roots[rootNmbr++] = newRoot(-div.coef[0] / div.coef[1], 0);
	else if (div.degree < 3) {
		double	a = div.coef[2], b = div.coef[1], c = div.coef[0],
				D = b * b - 4 * a * c;
		if (D >= 0) {				//	TWO_S : 2 separate solutions
			eq->roots[rootNmbr++] = newRoot((-b + sqrt(D)) / (2 * a), 0);
			eq->roots[rootNmbr++] = newRoot((-b - sqrt(D)) / (2 * a), 0);
		} else {					//	IMAGINARY_S : real part & imaginary part
			double real = (-b) / (2 * a), imag = sqrt(-D) / (2 * a);
			eq->roots[rootNmbr++] = newRoot(real, imag);
			eq->roots[rootNmbr++] = newRoot(real, -imag);
		}
	} else {
		int foundSol = 0;	 //bool
		EQUATION tmpSol = newEquation(1);
		tmpSol.coef[0] = div.coef[0];	// use -root for actual equation
		tmpSol.coef[1] = 1;

		//	No zero-coeff => 0 is root
		if (tmpSol.coef[0] == 0)	foundSol = setSol(eq, rootNmbr, 0);

		//	Try factorising zero-coeff and highest coeff
		if (!foundSol)		foundSol = setSol(eq, rootNmbr, (tmpSol.coef[0] = getXFactorising(div)));

		//	Try NewtonRaphson
		if (!foundSol)		foundSol = setSol(eq, rootNmbr, (tmpSol.coef[0] = getXNewtonRaphson(div)));

		//	Try rangeHalving
		if (!foundSol)		foundSol = setSol(eq, rootNmbr, (tmpSol.coef[0] = getXRangeHalving(div)));
	
		//	Solution found: divide away with horner
		if (foundSol)		return solveEqRecur(eq, divideBy(div, tmpSol).quotient, ++rootNmbr);
		else { printf("\nNo real solutions found (other complex solutions possible).\n"); return 0; }
	}
	return 1;
}

/*
*	Call solveEqRecur to start solving.
*/
int solveEq(EQUATION *eq) {
	return solveEqRecur(eq, *eq, 0);
}

/*
*	Print various info: equation, roots, derivative, integral
*/
void printInfo(EQUATION *eq) {
	printf("%s ", "Equation  :");
	toStringEq(*eq);

	printf("\n%s  %s\n", "Solving...", solveEq(eq) ? "Found roots!" : "No roots found (zeroed)!");

	printf("\n%s\n", "All roots :");
	toStringRoots(*eq);

	printf("\n%s ", "Alternate :");
	toStringEqFromSol(*eq);

	printf("\n%s ", "Derivative:");
	toStringEq(deriveEq(*eq));

	printf("\n%s ", "Integral  :");
	toStringEq(integrateEq(*eq));
}

/*
*	Makes the equation and print info.
*
*	degree > 0 ? parse cmd-line input
*	else asks for degree and reads the coefficients of an equation.
*
*	Returns the degree.
*/
int makeEquation(EQUATION *p, char **in, int degree) {
	if (degree > 1) {		//	Parse from cmd
		*p = newEquation(degree);
		readEquationCmd(p, in);
	} else {				//	Normal method
		printf("%s", "Degree? ");
		scanf("%d", &degree);
		if (degree < 1) return 0;
		
		*p = newEquation(degree);
		readEquation(p);

		putchar('\n');
	}

	printInfo(p);

	return p->degree;
}
