#include <stdio.h>

int main(void) 
{

//create file to plot data
FILE *newmeltplot;
newmeltplot = fopen("newmeltdataout.gnu", "w");
//fprintf(newmeltplot, "set title 'Graph of Free Energy vs Gaussian width'\n");
//fprintf(newmeltplot, "set title 'Graph of N_crystal vs Gaussian width'\n");
//fprintf(newmeltplot, "set xlabel 'width'\n");
//fprintf(newmeltplot, "set ylabel 'N_crystal'\n");
//fprintf(newmeltplot, "set ylabel 'Free Energy'\n");
//fprintf(newmeltplot, "plot 'newmeltdataout.dat' using ($1):($2) with lines title \"newmeltdataout.dat\" \n");
fprintf(newmeltplot, "plot 'newmeltdataout.dat' using ($1):($2) with lines title \"newmeltdataout.dat\" \n");
fprintf(newmeltplot, "pause -1\n");

return 0;
}
