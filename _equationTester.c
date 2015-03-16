#include "equations.h"

/**
*	Parse the input and give the equation output.
*		argc = lenght of cmd-line args
*		argv = array of args = array of strings
*
*	Cmd-line expects greatest degree coefficient first.
*	Parse cmd-line args only once, then user input.
*
*	No args? Read from user input.
*/
int main(int argc, char *argv[]) {
	double *equation;
	int readFromCmd = argc - 2;

	while (makeEquation(&equation, argv, readFromCmd) > 1) { readFromCmd = 0; putchar('\n'); };

	printf("\n");
	//system("PAUSE");
	return 0;
}
