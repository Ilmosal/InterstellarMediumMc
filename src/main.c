#include <stdio.h>
#include <string.h>
#include "Rtmc_core.h"

int main(int argc, char *argv[])
{
	int mode = strtol(*argv, NULL, 10);

	if (mode < 0 || mode > 3)
		return -1;
	
	run_RTMC(0);
	return 1;
}
