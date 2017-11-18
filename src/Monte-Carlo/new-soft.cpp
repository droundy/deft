#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <chrono>
#include "vector3d.h"
#include <popt.h>
#include <ctime>
#include <cassert>

// ------------------------------------------------------------------------------
// Global Constants
// ------------------------------------------------------------------------------
    const int x = 0;
    const int y = 1;
    const int z = 2;
    const double sigma = 1;
    const double epsilon = 1;
    const double cutoff = sigma * pow(2,1./6.);

// ------------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------------
// Places a desired number of spheres into an FCC lattice
inline vector3d *FCCLattice(int numOfSpheres, double systemLength[3]);

// Applies periodic boundary conditions in x,y,z
static inline vector3d periodicBC(vector3d inputVector, double systemLength[3],bool wall[3]);

// Calculates total potential of system
double totalPotential(vector3d *sphereMatrix, int numOfSpheres,double systemLength[3],bool wall[3]);

// Checks whether the move of a sphere is valid
bool conditionCheck(vector3d *sphereMatrix, vector3d movedSpherePos, int numOfSpheres, int movedSphereNum,double Temperature);

// Counts the number of spheres within a certain radius. Does not take into consideration mirror image particles.
double *radialDistribution(vector3d *sphereMatrix, int numOfSpheres, double sizeOfSystem[3],bool wall[3]);

// Finds the nearest mirror image in adjacent cubes and specifies the vector displacement
vector3d nearestImage(vector3d R2,vector3d R1, double systemLength[3], bool wall[3]);

// Calculates the potential energy due to radial distances between spheres
double bondEnergy(vector3d R);

// Calculates product of radial forces and displacements of spheres for an ensemble
double forceTimesDist(vector3d *spheres, int numOfSpheres,double systemLength[3],bool wall[3]);


int main(int argc, const char *argv[])  {
	time_t mytime;
    clock_t start = clock();
    // Initialize Variables and Dummy Variables
    double systemLength[3] = {0,0,0};
    double xy = 0;
    double reducedDensity = 0;
    double reducedTemperature = 0;
    long totalIterations = 1000;
    double dr = 0.005;
    bool wall[3] = {false,false,false};
    int numOfSpheres = 0;
    
	char *data_dir = new char[1024];
	sprintf(data_dir,"none");
	char *filename = new char[1024];
	sprintf(filename, "none");
	char *filename_suffix = new char[1024];
	sprintf(filename_suffix, "none");
	
	// ----------------------------------------------------------------------------
	// Parse input options
	// ----------------------------------------------------------------------------
    // To assign values from command line ex.
    // in deft ./new-soft --lenx 69.0 --leny 42.1717 --lenz 99 ... ad nauseam
	poptContext optCon;
    poptOption optionsTable[] = {
        /*** System Dimensions ***/
        {"lenx", '\0', POPT_ARG_DOUBLE, &systemLength[x], 0, 
            "System Length in X. "
            "This command is used for NVT/DVT. If NDT is true, this value"
            " will be rewritten.", "DOUBLE"},
        {"leny", '\0', POPT_ARG_DOUBLE, &systemLength[y], 0, 
            "System Length in Y. "
            "This command is used for NVT/DVT. If NDT is true, this value"
            " will be rewritten.", "DOUBLE"},
        {"lenz", '\0', POPT_ARG_DOUBLE, &systemLength[z], 0, 
            "System Length in Z. "
            "This command is used for NVT/DVT. If NDT is true, this value"
            " will be rewritten.", "DOUBLE"},
		{"lenxy", '0', POPT_ARG_DOUBLE, &xy, 0,
			"Length of x and y. This command is used for NDT. "
			"If DVT/NVT is true, this value will"
			" be unused.", "DOUBLE"},
        
        /*** Thermodynamic Properties ***/ 
        {"sphereNum",'\0', POPT_ARG_INT, &numOfSpheres, 0,
			"Number of Spheres in System", "INT"},
        {"density", '\0', POPT_ARG_DOUBLE, &reducedDensity, 0, 
            "Reduced Density", "DOUBLE"},
        {"temp", '\0', POPT_ARG_DOUBLE, &reducedTemperature, 0,
            "Reduced Temperature", "DOUBLE"},
        {"wallx", '\0', POPT_ARG_NONE, &wall[x], 0,
            "If x dimension has a wall (int). ", "BOOL"},
        {"wally", '\0', POPT_ARG_NONE, &wall[y], 0,
            "If y dimension has a wall (int)", "BOOL"},
        {"wallz", '\0', POPT_ARG_NONE, &wall[z], 0,
            "If z dimension has a wall (int)", "BOOL"},
        
        /*** Simulation Characteristics ***/
        {"iters" ,'\0', POPT_ARG_LONG, &totalIterations, 0,
            "Total number of iterations", "LONG"},
        {"dr", '\0', POPT_ARG_DOUBLE, &dr, 0,
            "Maximum size of random move", "DOUBLE"},
            
		/*** PARAMETERS DETERMINING OUTPUT FILE DIRECTORY AND NAMES ***/

		{"dir", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &data_dir, 0,
		 "Directory in which to save data", "DIRNAME"},
		{"filename", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT, &filename, 0,
		 "Base of output file names", "STRING"},
		{"filename-suffix", '\0', POPT_ARG_STRING | POPT_ARGFLAG_SHOW_DEFAULT,
		 &filename_suffix, 0, "Output file name suffix", "STRING"},
        
        POPT_AUTOHELP
        POPT_TABLEEND
    };
    optCon = poptGetContext(NULL, argc, argv, optionsTable, 0);
    poptSetOtherOptionHelp(optCon, "[OPTION...]\n"
                         "Type of BOOL means that the command just needs"
                         " to be called. You don't need to associate any"
                         " value with that argument. :)");
    int c = 0;
    while((c = poptGetNextOpt(optCon)) >= 0);
    if (c < -1) {
        fprintf(stderr, "\n%s: %s\n", poptBadOption(optCon, 0), poptStrerror(c));
        return 1;
    }
    poptFreeContext(optCon);

    if ((systemLength[x] < 0) || (systemLength[y] < 0) || (systemLength[z] < 0)){
        printf("System lengths can't be less than or zero\n");
        printf("System Length: (%g %g %g)\n",systemLength[x],systemLength[y],systemLength[z]);
        return 1;
    }
    if (totalIterations < 0) {
        printf("Total iterations can't be less than zero: %ld\n",totalIterations);
        return 1;
    }
    if (reducedTemperature <= 0){
        printf("Despite the attempts and dreams of delusional humans," 
        "the temperature can't be less than or equal to zero: "
        "Temperature: %g\n", reducedTemperature);
        return 1;
    }
    if (reducedDensity <= 0){
        printf("Do you really want a system with zero or negative balls???\n"
        "Reduced Density: %g\n", reducedDensity);
        return 1;
    }
    // Changes to Input Conditions
    if (numOfSpheres == 0) {
		assert(reducedDensity > 0);
		printf("Performing DVT simulation.\n");
		numOfSpheres = int(systemLength[x] * systemLength[y] * systemLength[z] * reducedDensity);
		printf("N to: %d\n",numOfSpheres);
	} else if (systemLength[x] == 0 && systemLength[y] == 0 && systemLength[z] == 0) {
		assert(reducedDensity > 0);
		assert(numOfSpheres > 0);
		printf("Performing NDT simulation\n");
		double tempVol = (sigma*sigma*sigma*numOfSpheres)/(reducedDensity);
		systemLength[x] = systemLength[y] = systemLength[z] = pow(tempVol,1./3.);
	} else if (systemLength[x] == 0 || systemLength[y] == 0 || systemLength[z] == 0) {
		printf("You need to specify the all of the lengths or none!\n");
		exit(1);
	}
    // ----------------------------------------------------------------------------
    // Initial Conditions
    // ----------------------------------------------------------------------------
	double volume = systemLength[x] * systemLength[y] * systemLength[z];
    vector3d *spheres = FCCLattice(numOfSpheres,systemLength);
    vector3d *sphereMoveVec = new vector3d[numOfSpheres];
    // Energy
    double currentEnergy = totalPotential(spheres,numOfSpheres,systemLength,wall);
	double totalEnergy = 0.0;
    // Pressure
    double pressureIdeal = (numOfSpheres*reducedTemperature*epsilon)/volume;
    double virial = forceTimesDist(spheres, numOfSpheres, systemLength,wall);
    double exPressure = 0.0;
    // Radial Dist.
    long *runningRadial = new long[1000];
    long radialWrites = 0;
    // Diffusion
	double diffusion = 0;
    double sphereTotalMove[numOfSpheres] = {0};
    // Structure Factor
    double kmax = 24*M_PI;
    const double cellNumber = ceil(pow(numOfSpheres/4,1./3.));
    double cellLength = systemLength[x]/cellNumber;
    double dk = M_PI / (20*cellLength);
    int kpoints = ceil(((kmax-2*M_PI)/cellLength)/dk);
    long structureWrites = 0;
    double structureFactor[kpoints][kpoints] = {0};
    // File Saving and misc.
	long acceptedTrials = 0;
    int outputCount = 0;
    bool save = false;
    long saveIter = 1e3;
    int saveTimes[8] = {1,5,10,30,60,300,600,1800}; // seconds
	mytime = time(NULL);
	long drAdjust = numOfSpheres*numOfSpheres;
    int structureReset = 0;
    // Introduction
	printf("Filename: %s\n",filename);
	printf("Lattice Made\n");
    printf("Number of Spheres: %d\n", numOfSpheres);
    printf("Size of System: %g %g %g\n", systemLength[x],systemLength[y],systemLength[z]);
    printf("Wall existence on x,y,z %d %d %d\n",wall[x],wall[y],wall[z]);
    printf("dr: %.5g\n",dr);
    printf(ctime(&mytime));
    // ----------------------------------------------------------------------------
    // File Save Initiation
    // ----------------------------------------------------------------------------
    char *headerinfo = new char[4096];
    sprintf(headerinfo,
		"# System Lengths: %g %g %g\n"
		"# Density: %g\n"
		"# Temperature: %g\n"
		"# N: %d\n"
		"# walls: %d %d %d\n"
		"# volume: %g\n",
		systemLength[x],systemLength[y],systemLength[z],
		reducedDensity,reducedTemperature,numOfSpheres,
		wall[x],wall[y],wall[z], volume);
    char *countinfo = new char[4096];
	char *pos_fname = new char[1024];
	sprintf(pos_fname, "%s/%s-pos.dat", data_dir, filename);
	char *radial_fname = new char[1024];
	sprintf(radial_fname, "%s/%s-radial.dat", data_dir, filename);
	char *press_fname = new char[1024];
	sprintf(press_fname,"%s/%s-press.dat", data_dir, filename);
	char *energy_fname = new char[1024];
	sprintf(energy_fname, "%s/%s-energy.dat",data_dir,filename);
	char *dif_fname = new char[1024];
	sprintf(dif_fname, "%s/%s-dif.dat",data_dir,filename);
    char *struc_fname = new char[1024];
    sprintf(struc_fname, "%s/%s-struc.dat",data_dir,filename);
  // ----------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------
	printf("Main loop started\n\n");
	fflush(stdout);
    for (long currentIteration = 0; currentIteration < totalIterations; ++currentIteration) {
        exPressure += virial;
        totalEnergy += currentEnergy;
        bool trialAcceptance = false;
        // Picks a random sphere to move and performs a random move on that sphere.
        int movedSphereNum = random::ran64() % numOfSpheres;
        vector3d initialSpherePos = spheres[movedSphereNum];
        vector3d randMove = vector3d::ran(dr);
        vector3d movedSpherePos = periodicBC(spheres[movedSphereNum]+randMove,systemLength,wall);
        // Determines the difference in energy due to the random move
        double energyNew = 0.0, energyOld = 0.0;
        for (int currentSphere = 0; currentSphere < numOfSpheres; ++currentSphere) {
            if (currentSphere != movedSphereNum) {
                vector3d Rj = spheres[currentSphere];
                vector3d R2 = nearestImage(Rj,movedSpherePos,systemLength,wall);
                vector3d R1 = nearestImage(Rj,initialSpherePos,systemLength,wall);
                energyNew += bondEnergy(R2);
                energyOld += bondEnergy(R1);
            }
        }

        double energyChange = energyNew - energyOld;
        if (energyChange <= 0)  {
            trialAcceptance = true;
        }   else if (energyChange > 0)  {
                double p = exp(-energyChange / reducedTemperature);
                double r = random::ran(); // fixme
                if (p > r)    {
                    trialAcceptance = true;
                }   else    {
                    trialAcceptance = false;
                }
            }
        if (trialAcceptance){ 
            spheres[movedSphereNum] = movedSpherePos;
            currentEnergy = totalPotential(spheres,numOfSpheres,systemLength,wall); // update instead
            // try currentEnergy += energyChange; // rerun totalPotential when saving?
            virial = forceTimesDist(spheres,numOfSpheres,systemLength,wall); // update instead
            // TODO: try timing runs with different N values, and plotting time per move versus N
            sphereMoveVec[movedSphereNum] += randMove;
            // small speed improvement:  delay the following until saving
			sphereTotalMove[movedSphereNum] = sphereMoveVec[movedSphereNum].normsquared();
			diffusion = 0; // reset to zero
            for (int i = 0; i < numOfSpheres; i++){
				diffusion += sphereTotalMove[i]/numOfSpheres;
			}
            acceptedTrials += 1;
        }
        if ((currentIteration % (10*numOfSpheres)) == 0){ // try increasing 10 -> 100?
			radialWrites += 1;
            double *radialDistHist;
            radialDistHist = radialDistribution(spheres,numOfSpheres,systemLength,wall);
            for (int i = 0; i < 1000; ++i){
                runningRadial[i] += radialDistHist[i];
            }
            delete[] radialDistHist;
        if ((currentIteration % (10*numOfSpheres*kpoints*kpoints)) == 0){ // try increasing 10 -> 100?
            structureWrites += 1;
            for (int kx = 0; kx < kpoints -1; kx++) {
                    for (int ky = 0; ky < kpoints -1; ky++) {
                        double rho_k_real = 0, rho_k_imag = 0;
                        for (int n = 0; n < numOfSpheres; n++) {
                            double kr = dk*(kx*spheres[n].x + ky*spheres[n].y);
                            rho_k_real += cos(kr);
                            rho_k_imag += sin(kr);
                        }
                        structureFactor[kx][ky] += (rho_k_real*rho_k_real + rho_k_imag*rho_k_imag)/numOfSpheres;
                    }
                }
            }
            
        }
    // ---------------------------------------------------------------
    // Displacement Vector Adjustment
    // ---------------------------------------------------------------
    if ((currentIteration % drAdjust == 0) && (currentIteration != 0)) {
		double ratio = double(acceptedTrials) / double(currentIteration);
		if (ratio >= 0.6) {
            dr *= 1.25;
            drAdjust += drAdjust;
		} else if (ratio <= 0.4) {
            dr *= 0.8;
            drAdjust += drAdjust;
		} else {
            drAdjust *= 10;
		}
	}
	// ---------------------------------------------------------------
    // Save data to files
    // ---------------------------------------------------------------
		if ((currentIteration == saveIter) || (currentIteration + 1 == totalIterations)){
			clock_t end = clock();
			double clock = (end - start)/double(CLOCKS_PER_SEC);
			if (currentIteration + 1 == totalIterations){
				save = true;
			}	else if (clock < saveTimes[outputCount]){
				saveIter = ceil(saveTimes[outputCount]*double(saveIter)/clock);
				save = false;
			}	else if ((clock >= saveTimes[outputCount])&&(outputCount < 7)) {
				outputCount += 1;
				save = true;
				saveIter = ceil(saveTimes[outputCount]*double(saveIter)/clock);
			}	else if ((clock >= saveTimes[outputCount]) && (outputCount == 7)){
				saveIter = ceil((saveTimes[outputCount]+clock)*double(saveIter)/clock);
                structureReset += 1;
				save = true;
			}

			if (save) {
				long minutes = long(clock) / 60;
				long hours = minutes / 60;
				long days = hours / 24;
				double ratio = (double(currentIteration+1)/ 
								double(totalIterations));
				long estEndTime = (double(clock)/ratio)-clock;
				long fmin = estEndTime/60;
				long fhrs = fmin/60;
				long fdays = fhrs/24;
				printf("Filename: %s\n",filename);
				printf("Iteration: %ld of %ld\n",
					currentIteration+1,totalIterations);
				printf("Percent Complete: %.7f \n",100.0*ratio);
				printf("Acceptance Rate: %.5f \n", 
					100.0*double(acceptedTrials)/double(currentIteration+1));
				printf("dr: %.5f, drAdjust: %ld\n",dr,drAdjust);
				printf("Running for: %ld days, %ld hrs, %ld min, %ld sec \n"
						,days,hours-days*24,
						minutes-hours*60,long(clock)-minutes*60);
				printf("Estimated to Finish in: %ld days, %ld hrs, %ld min, %ld sec \n"
						,fdays,fhrs-fdays*24,
						fmin-fhrs*60,estEndTime-fmin*60);
				mytime = time(NULL);
				printf("Save time: %s\n\n",ctime(&mytime));
				fflush(stdout);
				sprintf(countinfo,
					"# iterations: %ld\n"
					"# accepted moves: %ld\n",
					currentIteration + 1, acceptedTrials);
				// Save Positions
				FILE *pos_out = fopen(pos_fname, "w");
				fprintf(pos_out, "%s%s", headerinfo,countinfo);
				for (int i=0; i<numOfSpheres; i++) {
					fprintf(pos_out, "%g\t%g\t%g\n",
					spheres[i].x, spheres[i].y, spheres[i].z);
				}
				// Save Radial
				FILE *radial_out = fopen(radial_fname,"w");
				fprintf(radial_out, "%s%s", headerinfo,countinfo);
				for (int i = 0; i < 1000; ++i){ 
                    const double dr = systemLength[x]/(2*1000.0);
                    double r_i = i*dr;
					fprintf(radial_out, "%g\t%g\n", r_i, 
						double(volume*runningRadial[i]/
							  (4*M_PI*r_i*r_i*dr*numOfSpheres*radialWrites)));
				}
				// Save Pressure
				FILE *press_out = fopen(press_fname,"w");
				fprintf(press_out, "%s%s", headerinfo,countinfo);
				fprintf(press_out, "%g\n",pressureIdeal + 
					(exPressure / (currentIteration*volume))); 
				// Save Energy
				FILE *energy_out = fopen(energy_fname,"w");
				fprintf(energy_out,"%s%s",headerinfo,countinfo);
				fprintf(energy_out,"%g\n",
						totalEnergy/double(currentIteration));
				// Save Diffusion
				FILE *dif_out = fopen(dif_fname,"w");
				fprintf(dif_out, "%s%s", headerinfo,countinfo);
                fprintf(dif_out,"%g\n",double(diffusion/(currentIteration+1)));
                // Save Structure Factor
                FILE *struc_out = fopen(struc_fname,"w");
                fprintf(struc_out, "%s%s", headerinfo,countinfo);
                fprintf(struc_out, "# Structure Writes %ld\n",structureWrites);
                fprintf(struc_out, "# kpoints %d\n", kpoints);
                for (int kx = 0; kx < kpoints -1; kx++) {
                    for (int ky = 0; ky < kpoints -1; ky++) {
                        fprintf(struc_out,"%g\t",structureFactor[kx][ky]/structureWrites);
                    }
                    fprintf(struc_out,"\n");
                }
                // Close files
				fclose(pos_out);
				fclose(radial_out);
				fclose(press_out);
				fclose(energy_out);
				fclose(dif_out);
                fclose(struc_out);
				save = false; 
			}
            if ((structureReset == 1) || (structureReset == 3)) {
                structureFactor[kpoints][kpoints] = {0};
                printf("Structure Factor Reset\n");
                printf("Reset Time: %s\n\n",ctime(&mytime));
                fflush(stdout);
            }
		}
    }
	// ----------------------------------------------------------------------------
	// END OF MAIN PROGRAM LOOP
	// ----------------------------------------------------------------------------
    printf("All Done :)\n");
    delete[] runningRadial;
    delete[] spheres;
    delete[] countinfo;
    delete[] pos_fname;
    delete[] radial_fname;
    delete[] press_fname;
    delete[] energy_fname;
    delete[] dif_fname;
    delete[] struc_fname;
}

inline vector3d *FCCLattice(int numOfSpheres, double systemLength[3])   {
  assert(systemLength[0] == systemLength[1]);
  assert(systemLength[0] == systemLength[2]);
  vector3d *sphereMatrix = new vector3d[numOfSpheres];
  double cellNumber = ceil(pow(numOfSpheres/4,1./3.));
  double cellLength = systemLength[x]/cellNumber;

  vector3d *offset = new vector3d[4];
  offset[0] = vector3d(0,0,0);
  offset[1] = vector3d(cellLength,0,cellLength)/2;
  offset[2] = vector3d(cellLength,cellLength,0)/2;
  offset[3] = vector3d(0,cellLength,cellLength)/2;

  int sphereNum = 0;
  for (int i=0; i<cellNumber; i++) {
    for (int j=0; j<cellNumber; j++) {
      for (int k=0; k<cellNumber; k++) {
        vector3d cornerSphere = { cellLength*i, cellLength*j, cellLength*k };
        for (int b=0; b<4; b++) {
          sphereMatrix[sphereNum++] = cornerSphere + offset[b];
          if (sphereNum == numOfSpheres) {
            delete[] offset;
            return sphereMatrix;
          }
        }
      }
    }
  }
  delete[] offset;
  delete[] sphereMatrix;
  return 0; // return null pointer if we could not fit the spheres
}

static inline vector3d periodicBC(vector3d inputVector, double systemLength[3], bool wall[3])   {
    for (int i = 0; i < 3; i++){
        while (inputVector[i] > systemLength[i]){
            if (wall[i] == true){
                inputVector[i] = 2*systemLength[i] - inputVector[i];
            } else {
                inputVector[i] -= systemLength[i];
            }
        }
        while (inputVector[i] < 0){
            if (wall[i] == true){
                inputVector[i] = fabs(inputVector[i]);
            } else{
                inputVector[i] += systemLength[i];
            }
        }
    }
     return inputVector;
}

double totalPotential(vector3d *sphereMatrix, int numOfSpheres, double systemLength[3],bool wall[3]) {
    double totalPotential = 0.0;

    for (int i = 0; i < numOfSpheres; ++i)  {
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j)  {
            vector3d Rj = sphereMatrix[j];
            vector3d R = nearestImage(Rj,Ri,systemLength,wall);
            totalPotential += bondEnergy(R);
        }
    }
    return totalPotential;
}

double *radialDistribution(vector3d *sphereMatrix, int numOfSpheres, double systemLength[3], bool wall[3]) {
    const int numOfBoxes = 1000;
    double *deltan = new double[numOfBoxes];
    
    for (int i = 0; i < numOfBoxes; i++){ 
		deltan[i] = 0;					  
	}
    for (int i = 0; i < numOfSpheres; ++i) {
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j) {
            vector3d Rj = sphereMatrix[j];
            vector3d R = nearestImage(Rj,Ri,systemLength,wall);
            double Rmag = R.norm();
            if (Rmag <= (systemLength[x]/2)){
                int box = floor(Rmag*numOfBoxes/(systemLength[x]/2));   // box = (R/dr) = (R / (sizeOfSystem/numOfBoxes))
                deltan[box] += 1;
            }
        }
    }
    return deltan;
}

vector3d nearestImage(vector3d R2,vector3d R1, double systemLength[3], bool wall[3]){
    vector3d R = R2 - R1;
    if (fabs(R.x) > systemLength[x]/2 && !wall[x]) {
        R.x -= copysign(systemLength[x],R.x);
    }
    if (fabs(R.y) > systemLength[y]/2 && !wall[y]) {
        R.y -= copysign(systemLength[y],R.y);
    }
    if (fabs(R.z) > systemLength[z]/2 && !wall[z]) {
        R.z -= copysign(systemLength[z],R.z);
    }
    return R;
}

double bondEnergy(vector3d R){
    if (R.norm() < cutoff){
        double Rsq = R.normsquared();
        double SR2 = (sigma*sigma)/Rsq;
        double SR6 = SR2*SR2*SR2;
        double SR12 = SR6*SR6;
        double bondEnergy = SR12 - SR6;
        return 4*epsilon*bondEnergy + epsilon;
    } else {
        return 0;
    }
}

double forceTimesDist(vector3d *spheres, int numOfSpheres, double systemLength[3], bool wall[3]) {
    double w = 0;
    for (int i = 0; i < numOfSpheres; ++i){
        vector3d Ri = spheres[i];
        for (int j = i + 1; j < numOfSpheres; ++j){
            vector3d Rj = spheres[j];
            vector3d R = nearestImage(Rj,Ri,systemLength,wall);
            double R2 = R.normsquared();
            if (R2 < cutoff){
                double SR2 = (sigma*sigma)/R2;
                double SR6 = SR2*SR2*SR2;
                double SR12 = SR6*SR6;
                w += 2*SR12 - SR6;
            }
        }
    }
    return 8*epsilon*w;
}
