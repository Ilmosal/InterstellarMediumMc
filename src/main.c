#include <stdio.h>
#include <string.h>
#include "Rtmc_core.h"

int main(int argc, char *argv[])
{
	if (argc < 2)
		return -1;

	int mode = strtol(argv[1], NULL, 10);

	if (mode < 0 || mode > 3)
		return -1;
	
	run_RTMC(mode);
	return 1;
}
