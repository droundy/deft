#include <stdio.h>
#include <time.h>
#include "Monte-Carlo/monte-carlo.h"
#include <cassert>
#include <math.h>
#include <stdlib.h> 
#include <string.h>
#include <vector>
using std::vector;
using namespace std;
int main(){
	char outfilename[] = "test.dat";
	char pressureFile[] = strcat(outfilename,".prs");
	FILE *out = fopen((const char *)outfilename,"w");
        if (out == NULL) {
          printf("Error creating file %s\n", outfilename);
          return 1;
	}
        
        fprintf(out, "test outprint\n");
        fflush(stdout);
        fclose(out);

        FILE *pressuref = fopen((const char *)pressureFile,"a");
        if (pressuref == NULL) {
          printf("Error creating file %s\n", pressureFile);
          return 1;
        }
        fprintf(pressuref, "pressure test\n");
        fflush(stdout);
        fclose(pressuref);
	return 0;
}
