/*
*	Toggle the comment on the header you want to use.
*		equations.h		Standard functions regarding equations;
*						e.g. reading, printing, basic solver and operations.
*
*		rpoly.h			Advanced eaquation solver with read/print functions. 
*						Finds all solutions in both real and complex planes.
*/
#include "equations.h"
//#include "rpoly.h"

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
	int readFromCmd = argc - 2;
	EQUATION e;

#ifdef _USE_RPOLY
	while (makeEquation_ak1(&e, argv, readFromCmd) > 0)
#else
	while (makeEquation(&e, argv, readFromCmd) > 0)		
#endif
		{ readFromCmd = 0; putchar('\n'); };

	printf("\n");
	system("PAUSE");
	return 0;
}
