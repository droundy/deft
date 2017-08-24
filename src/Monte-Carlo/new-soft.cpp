#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <chrono>
#include <cassert>
#include "vector3d.h"
#include <popt.h>
#include <ctime>

using namespace std;
using sec = chrono::seconds;
using ms = chrono::milliseconds;
using get_time = chrono::steady_clock;
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
    auto start = get_time::now();
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
	}	// NVT doesn't need any tweaks to initial conditions
  // ----------------------------------------------------------------------------
  // Initial Conditions
  // ----------------------------------------------------------------------------
	double volume = systemLength[x] * systemLength[y] * systemLength[z];
    vector3d *spheres = FCCLattice(numOfSpheres,systemLength);
    // rename this probably???
    vector3d *sphDif = new vector3d[numOfSpheres];
    double tempEnergy = totalPotential(spheres,numOfSpheres,systemLength,wall);
    double tempE2 = tempEnergy*tempEnergy;
    double pressureIdeal = (numOfSpheres*reducedTemperature*epsilon)/volume;
    double virial = forceTimesDist(spheres, numOfSpheres, systemLength,wall);
    double totalEnergy = 0.0;
    double energy2 = 0.0;
    double exPressure = 0.0;
    long acceptedTrials = 0;
    long *runningRadial = new long[1000];
    long radialWrites = 0;
    int outputCount = 0;
    bool save = false;
    long saveIter = 1e3;
    double *glass = new double[totalIterations];
    double Di[numOfSpheres] = {0};
    double tempD = 0;
    int saveTimes[8] = {1,5,10,30,60,300,600,1800}; // seconds
	mytime = time(NULL);
	
	printf("Filename: %s\n",filename);
	printf("Lattice Made\n");
    printf("Number of Spheres: %d\n", numOfSpheres);
    printf("Size of System: %g %g %g\n", systemLength[x],systemLength[y],systemLength[z]);
    printf("Wall existence on x,y,z %d %d %d\n",wall[x],wall[y],wall[z]);
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
	//~ FILE *dif_out = fopen((const char *)dif_fname,"w");//Rewrites
  // ----------------------------------------------------------------------------
  // MAIN PROGRAM LOOP
  // ----------------------------------------------------------------------------
	printf("Main loop started\n\n");
	fflush(stdout);
    for (long currentIteration = 0; currentIteration < totalIterations; ++currentIteration) {
        exPressure += virial;
        totalEnergy += tempEnergy;
        energy2 += tempE2;
        glass[currentIteration] = tempD;
        bool trialAcceptance = false;
        // Picks a random sphere to move and performs a random move on that sphere.
        int movedSphereNum = random::ran64() % numOfSpheres;
        vector3d initialSpherePos = spheres[movedSphereNum];
        vector3d randMove = vector3d::ran(dr);
        vector3d movedSpherePos = spheres[movedSphereNum] + randMove;
        movedSpherePos = periodicBC(movedSpherePos,systemLength,wall);
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
                double p = exp(-energyChange / reducedTemperature); // Need to get units correct in here.
                double r = random::ran();
                if ((p > r))    {
                    trialAcceptance = true;
                }   else    {
                    trialAcceptance = false;
                }
            }
        if (trialAcceptance == true){   // If accepted, update lattice, PE, and radial dist. func.
            spheres[movedSphereNum] = movedSpherePos;
            tempEnergy = totalPotential(spheres,numOfSpheres,systemLength,wall);
            tempE2 = tempEnergy*tempEnergy;
            virial = forceTimesDist(spheres, numOfSpheres,systemLength,wall);
            sphDif[movedSphereNum] += randMove;
            double paper = sphDif[movedSphereNum].normsquared();
            Di[movedSphereNum] += paper;
            for (int i = 0; i < numOfSpheres; i++){
				tempD += Di[i]/numOfSpheres;
			}
            acceptedTrials += 1;
        }
        if ((currentIteration % (10*numOfSpheres)) == 0){
			radialWrites += 1;
            double *radialDistHist;
            radialDistHist = radialDistribution(spheres,numOfSpheres,systemLength,wall);
            for (int i = 0; i < 1000; ++i){
                runningRadial[i] += radialDistHist[i];
            }
            delete[] radialDistHist;
        }
	// ---------------------------------------------------------------
    // Save data to files
    // ---------------------------------------------------------------
		if ((currentIteration == saveIter) || (currentIteration + 1 == totalIterations)){
			auto end = get_time::now();
			auto diff = end - start;
			double clockms = double(chrono::duration_cast<ms>(diff).count());
			long clock = long(chrono::duration_cast<sec>(diff).count());
			if (currentIteration + 1 == totalIterations){
				save = true;
			}	else if (clock < saveTimes[outputCount]){
				saveIter = ceil(1e3*saveTimes[outputCount]*double(saveIter)/clockms);
				save = false;
			}	else if ((clock >= saveTimes[outputCount])&&(outputCount < 7)) {
				outputCount += 1;
				save = true;
				saveIter = ceil(saveTimes[outputCount]*double(saveIter)/clock);
			}	else if ((clock >= saveTimes[outputCount]) && (outputCount == 7)){
				saveIter = ceil((saveTimes[outputCount]+clock)*double(saveIter)/clock);
				save = true;
			}

			if (save == true){
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
				printf("Percent Complete: %.7f %\n",100.0*ratio);
				printf("Acceptance Rate: %.5f %\n", 
					100.0*double(acceptedTrials)/double(currentIteration+1));
				printf("Running for: %ld days, %ld hrs, %ld min, %ld sec \n"
						,days,hours-days*24,
						minutes-hours*60,clock-minutes*60);
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
				FILE *pos_out = fopen((const char *)pos_fname, "w");
				fprintf(pos_out, "%s%s", headerinfo,countinfo);
				for (int i=0; i<numOfSpheres; i++) {
					fprintf(pos_out, "%g\t%g\t%g\n",
					spheres[i].x, spheres[i].y, spheres[i].z);
				}
				// Save Radial
				FILE *radial_out = fopen((const char *)radial_fname,"w");
				fprintf(radial_out, "%s%s", headerinfo,countinfo);
				for (int i = 0; i < 1000; ++i){ // This isn't being averaged properly
					fprintf(radial_out, "%g\t%g\n", 
						(i*systemLength[x])/(2*1000.0), 
						double(volume*runningRadial[i]/
							  (numOfSpheres*numOfSpheres*radialWrites)));
				}
				// Save Pressure
				FILE *press_out = fopen((const char *)press_fname,"w");
				fprintf(press_out, "%s%s", headerinfo,countinfo);
				fprintf(press_out, "%g\n",pressureIdeal + 
					(exPressure / (totalIterations*volume)));
				// Save Energy
				FILE *energy_out = fopen((const char *)energy_fname,"w");
				fprintf(energy_out,"%s%s",headerinfo,countinfo);
				fprintf(energy_out,"%g\t%g\n",
						totalEnergy/double(currentIteration),
						energy2/double(currentIteration));
				// Save Diffusion
				FILE *dif_out = fopen((const char *)dif_fname,"w");
				for (int i = 0; i < currentIteration; i++){
					fprintf(dif_out,"%g\n",double(glass[i]/currentIteration));
				}
				fclose(pos_out);
				fclose(radial_out);
				fclose(press_out);
				fclose(energy_out);
				fclose(dif_out);
				save = false;
			}
		}
    }
	// ----------------------------------------------------------------------------
	// END OF MAIN PROGRAM LOOP
	// ----------------------------------------------------------------------------
    printf("All Done :)\n");
    delete[] runningRadial;
    delete[] spheres;
}

inline vector3d *FCCLattice(int numOfSpheres, double systemLength[3])   {
    vector3d *sphereMatrix = new vector3d[numOfSpheres];
    double cellNumber = ceil(pow(numOfSpheres/4,1./3.));
    double cellLength = systemLength[x]/cellNumber;
    int ysteps = 0; int  zsteps = 0;
    
    vector3d cornerSphere = {sigma/2,sigma/2,sigma/2};   // Sphere's Center is radial distance from borders
    vector3d *offset = new vector3d[4];
    offset[0] = vector3d(0,0,0);
    offset[1] = vector3d(cellLength,0,cellLength)/2;
    offset[2] = vector3d(cellLength,cellLength,0)/2;
    offset[3] = vector3d(0,cellLength,cellLength)/2;
    
    for (int sphereNum = 0; sphereNum < numOfSpheres; sphereNum++) {
        sphereMatrix[sphereNum] = cornerSphere + offset[sphereNum%4];
        if (sphereNum%4 == 0) {
            if (cornerSphere.z + cellLength < systemLength[z]) {
                cornerSphere.z += cellLength;
                zsteps += 1;
            } else if ((cornerSphere.z + cellLength >= systemLength[z])
                    && (cornerSphere.y + cellLength  < systemLength[y])) {
                cornerSphere.z -= zsteps*cellLength;
                cornerSphere.y += cellLength;
                ysteps += 1;
                zsteps = 0;
            } else if ((cornerSphere.y + cellLength >= systemLength[y])
                    && (cornerSphere.z + cellLength >=systemLength[z])
                    && (cornerSphere.x + cellLength < systemLength[x])) {
                cornerSphere.z -= zsteps*cellLength;
                cornerSphere.y -= ysteps*cellLength;
                cornerSphere.x += cellLength;
                zsteps = ysteps = 0;
            } else if ((cornerSphere.x + cellLength >= systemLength[x])
                    && (cornerSphere.y + cellLength >= systemLength[y])
                    && (cornerSphere.z + cellLength >= systemLength[z])){
                sphereNum = numOfSpheres;
                }
        }
        if ((fabs(sphereMatrix[sphereNum].x) < 0.0001) & 
			(fabs(sphereMatrix[sphereNum].y) < 0.0001) & // Rescale Primitive Unit Cell
            (fabs(sphereMatrix[sphereNum].z) < 0.0001)){
            printf("Cell Dimensions have shrunk.\n");
            cornerSphere.x = cornerSphere.y = cornerSphere.z =  sigma/2;
            ysteps=zsteps = 0;
            cellNumber += 1;
            cellLength = systemLength[x]/cellNumber;
            offset[1] = vector3d(cellLength,0,cellLength)/2;
            offset[2] = vector3d(cellLength,cellLength,0)/2;
            offset[3] = vector3d(0,cellLength,cellLength)/2;
            sphereNum = -1;
        }
    }
    return sphereMatrix;
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

    for (int i = 0; i < numOfSpheres; ++i) {
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j) {
            vector3d Rj = sphereMatrix[j];
            vector3d R = nearestImage(Rj,Ri,systemLength,wall);
            double Rmag = R.norm();
            int box = int(Rmag*numOfBoxes/(systemLength[x]/2));   // box = (R/dr) = (R / (sizeOfSystem/numOfBoxes))
            if (Rmag <= systemLength[x]/2){
                deltan[box] += 1;
            }
        }
    }
    return deltan;
}

vector3d nearestImage(vector3d R2,vector3d R1, double systemLength[3], bool wall[3]){
    vector3d R = R2 - R1;
    if ((fabs(R.x) > (systemLength[x]/2))
        && (wall[x] == false)){
        R.x -= copysign(systemLength[x],R.x);
    }
    if ((fabs(R.y) > (systemLength[y]/2))
        && (wall[y] == false)){
        R.y -= copysign(systemLength[y],R.y);
    }
    if ((fabs(R.z) > (systemLength[z]/2))
        && (wall[z] == false)){
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
