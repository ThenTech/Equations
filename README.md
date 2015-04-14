# Equations
Program that solves polynomials read and parsed by user input or from command line.
Gives all solutions in real and complex plane if (possible) using various methods.
Base method: find a real solution and divide away with Horner, then solve quotient.

<b>equations.h</b>  
A header containing every method necessary for solving equations. With structs.

<b>equations.old.h</b>  
The initial header.

<b>rpoly.h</b>  
Uses Jenkins-Traub algorithm to solve equations in both real and complex planes.

<b>_equationTester.c</b>  
The  main method for testing.
