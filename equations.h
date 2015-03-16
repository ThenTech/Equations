#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//	Epsilon, close to zero for comparison
#define EPS		1e-12

// DEPRECATED :: Possible solutions for quadratic equation
typedef enum	{NO_S, ONE_S, TWO_S, DOUBLE_S, IMAGINARY_S} QuadEQsolution;

/**
*	Returns true if difference is smaller than EPS or are equal.
*/
int fuzzyEquals(double x, double y) {
	return (fabs(x - y) <= EPS) || (x == y);
}

/*
*	Reads the coefficients and stores them in eqCoef.
*/
void readEquation(double eqCoef[], size_t degree) {
	for (int i = degree, cf = 97; i >= 0; i--, cf++) {
		if (i > 0)	printf("%cx^%d + ", cf, i);
		else		printf("%c = 0\n", cf, i);
	}

	for (int i = degree, cf = 97; i >= 0; i--, cf++) {
		printf("%c = ", cf);
		scanf("%lf", &eqCoef[i]);
	}
}

/*
*	Parses the coefficients from Cmd and stores them in eqCoef.
*/
void readEquationCmd(double eqCoef[], char **in, size_t degree) {
	for (int i = degree; i >= 0; i--)
		eqCoef[i] = atof(in[degree - i + 1]);	//	Reverse store coefficients, convert string to float with 'atof(*char)'
}

/*
*	Prints the equation.
*/
void printEquation(double* eqCoef, size_t degree) {
	for (int i = degree; i >= 0; i--)
		if (i > 0)	printf("%.4gx^%d + ", eqCoef[i], i);
		else		printf("%.4g = 0\n", eqCoef[i], i);
}

/*	DEPRECATED
*	Prints alternate form of equation using roots.
*/
void printEqFromSol(double* eqCoef, double* solutions, size_t degree) {
	if (eqCoef[degree] != 1) printf("%g * ", eqCoef[degree]);
	for (int i = 1; i <= degree; i++) {
		switch ((int)solutions[0]) {
			case ONE_S: case TWO_S:	case DOUBLE_S: 
				printf("(x - (%g))", solutions[i]); break;
			case IMAGINARY_S:	
				printf("(x - (%g+%gi)) * (x - (%g+%gi))", solutions[i], solutions[i + 1], solutions[i], solutions[i + 1]); 
				i++; break;
			case NO_S: default:	break;
		}
		printf("%s", (i < degree ? " * " : " = 0"));
	}
}

/*
*	Returns the deriven equation.
*/
double* deriveEq(double* equation, size_t degree) {
	double *derivative = malloc(degree * sizeof(double));
	for (int i = 0; i < degree; i++)
		derivative[i] = equation[i + 1] * (i + 1);
	return derivative;
}

/*
*	Returns the integrated equation.
*/
double* integrateEq(double* equation, size_t degree) {
	double *integral = malloc((degree + 1) * sizeof(double));
	integral[0] = 0;
	for (int i = 0; i < degree + 1; i++)
		integral[i + 1] = equation[i] / (i + 1);
	return integral;
}

/*
*	Returns quotient of equation with its root using Horner's Scheme
*/
double* devideBy(double* equation, size_t degree, double root) {
	double *quotient = malloc(degree * sizeof(double));
	quotient[degree - 1] = equation[degree];
	for (int i = degree - 2; i >= 0; i--)
		quotient[i] = root * quotient[i + 1] + equation[i + 1];
	return quotient;
}

/*
*	Returns the value of the equation with given X using Horner's method.
*		== equation(X)
*/
double getValue(double *equation, size_t degree, double x) {
	double y = 0.0;
	for (int i = degree; i >= 0; i--)
		y = y * x + equation[i];
	return y;
}

/*
*	Returns an array with all deviders of the given number. 
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
*			+- (Any devider of highest degree coeff) / (Any devider of lowest degree coeff)
*/
double getXFactorising(double *equation, size_t degree) {
	long *factA = getFactors(equation[degree]);
	long *factD = getFactors(equation[0]);

	for (long i = 0, a; (a = factA[i]) != 0; i++)
		for (long j = 0, d; (d = factD[j]) != 0; j++) {
			double q = (double) d / a;
			if (fuzzyEquals(getValue(equation, degree, q), 0.0))	return q;
			if (fuzzyEquals(getValue(equation, degree, -q), 0.0))	return -q;
		}
	return 0.0;
}

/*
*	Newton-Raphson Method
*		Approximate a root by continiously subtracting the derivative value.
*/
double getXNewtonRaphson(double *equation, size_t degree) {
	double *deriv = deriveEq(equation, degree), x = 1, f, d;
	for (int i = 0; i < 1000; i++) {
		f = getValue(equation, degree, x);
		d = getValue(deriv, degree - 1, x);
		if (fuzzyEquals(f, 0.0))	return x;
		if (!fuzzyEquals(d, 0.0))	x -= (f / d);
		else						x++;
	}
	return x;
}

/*
*	Range Halving method
*		Takes the avarage value in a range until range limits are equal thus giving a root.
*/
double getXRange(double *equation, size_t degree, double a, double b) {
	double mid = (a + b) / 2.0;	//Midpoint
	if (getValue(equation, degree, mid) == 0.0 || fuzzyEquals(a, b))
		return mid;
	if (getValue(equation, degree, a) * getValue(equation, degree, mid) < 0.0)
		return getXRange(equation, degree, a, mid);
	return getXRange(equation, degree, mid, b);
}

/*
*	Call functin for Range Halving method.
*/
double getXRangeHalving(double *equation, size_t degree) {
	return getXRange(equation, degree, -100.0, 100.0);
}

/*
*	Print root if actual zero.
*/
int printSol(double *equation, size_t degree, double root) {
	if (fuzzyEquals(getValue(equation, degree, root), 0.0))
		return (printf("% .8g\n", root) > 0 ? 1 : 0);
	return 0;
}

/*
*	Solve and print out the roots of a given equation.
*/
void solvePrintEq(double* equation, size_t degree) {
	if (degree < 1)
		printf("No solutions, invalid equation.\n");
	else if (degree < 2)			//	ONE_S : 1 single solutions
		printf("% .8g\n", (-equation[0] / equation[1]));
	else if (degree < 3) {
		double	a = equation[2], b = equation[1], c = equation[0],
				D = b * b - 4 * a * c;
		if (D >= 0) {				//	TWO_S : 2 seperate solutions
			printf("% .8g\n", (-b + sqrt(D)) / (2 * a));
			printf("% .8g\n", (-b - sqrt(D)) / (2 * a));
		} else {					//	IMAGINARY_S : real part & imaginary part
			double real = (-b) / (2 * a), imag = sqrt(-D) / (2 * a);
			printf("% .8g + %.8gi\n", real, imag);
			printf("% .8g - %.8gi\n", real, imag);
		}
	} else {
		int foundSol = 0;	 //bool
		double tmpSol = equation[0];

		//	No zero-coeff => 0 is root
		if (tmpSol == 0)	foundSol = (printf("% d\n", 0) > 0 ? 1 : 0);

		//	Try factorising zero-coeff and highest coeff
		if (!foundSol)		foundSol = printSol(equation, degree, (tmpSol = getXFactorising(equation, degree)));

		//	Try NewtonRaphson
		if (!foundSol)		foundSol = printSol(equation, degree, (tmpSol = getXNewtonRaphson(equation, degree)));

		//	Try rangeHalving
		if (!foundSol)		foundSol = printSol(equation, degree, (tmpSol = getXRangeHalving(equation, degree)));
	
		//	Solution found: devide away with horner
		if (foundSol)		solvePrintEq(devideBy(equation, degree, tmpSol), degree - 1);
		else				printf("No real solutions found (other complex solutions possible).\n");
	}
}

/*
*	Print various info: equation, roots, derivative, integral
*/
void printInfo(double* equation, size_t degree) {
	printf("%s\n", "Equation:");
	printEquation(equation, degree);

	printf("\n%s\n", "Solutions:");
	solvePrintEq(equation, degree);

	printf("\n%s\n", "Derivative:");
	printEquation(deriveEq(equation, degree), degree - 1);

	printf("\n%s\n", "Integral:");
	printEquation(integrateEq(equation, degree), degree + 1);
}

/*
*	Makes the equation and print info.
*
*	degree > 0 ? parse cmd-line input
*	else asks for degree and reads the coefficients of an equation.
*
*	Returns the degree.
*/
int makeEquation(double* equation, char **in, int degree) {
	if (degree > 1) {		//	Parse from cmd
		equation = malloc((degree + 1) * sizeof(double));
		readEquationCmd(equation, in, degree);
	} else {				//	Normal method
		printf("%s", "Degree? ");
		scanf("%d", &degree);
		if (degree < 1) return 0;
		
		equation = malloc((degree + 1) * sizeof(double));
		readEquation(equation, degree);

		putchar('\n');
	}

	printInfo(equation, degree);

	return degree;
}
