# Equations
Program that solves polynomials read and parsed by user input or from command line.
Gives all solutions in real and complex plane if (possible) using various methods.
Base method: find a real solution and divide away with Horner, then solve quotient.

equations.h  
A header containing every method necessary for solving equations. With structs.

equations.old.h  
The initial header.

rpoly.h  
Uses Jenkins-Traub algorithm to solve equations in both real and complex planes.

_equationTester.c  
The  main method for testing.
