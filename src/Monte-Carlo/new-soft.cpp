
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <chrono>
#include "vector3d.h"

using namespace std;
using sec = chrono::seconds;
using get_time = chrono::steady_clock;

// ------------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------------

// Places a desired number of spheres into an FCC lattice
inline vector3d *FCCLattice(int numOfSpheres, double sphereRadius);

// Applies periodic boundary conditions in x,y,z
static inline vector3d periodicBC(vector3d inputVector, double sizeOfSystem);

// Calculates total potential of system
double totalPotential(vector3d *sphereMatrix, int numOfSpheres,int n,double dVr);

// Checks whether the move of a sphere is valid
bool conditionCheck(vector3d *sphereMatrix, vector3d movedSpherePos, int numOfSpheres, int movedSphereNum,double Temperature);

// Counts the number of spheres within a certain radius. Does not take into consideration mirror image particles.
int *radialDistribution(vector3d *sphereMatrix, int numOfSpheres, double sizeOfSystem);

// Finds the nearest mirror image in adjacent cubes and specifies the vector displacement
vector3d nearestImage(vector3d R2,vector3d R1, double sizeOfSystem);

// Calculates the potential energy due to radial distances between spheres
double bondEnergy(vector3d R);

// Calculates product of radial forces and displacements of spheres for an ensemble
double pairVirialFunction(vector3d *spheres);

// ------------------------------------------------------------------------------
// Test Variables
// ------------------------------------------------------------------------------

int numOfSpheres = 32;
long totalIterations = 100000;
double dr = 0.005;
double sigma = 3.405;
double epsilon = 1;
double reducedDensity = 1;
double reducedTemperature = 1 ;  // (epsilon * T / k_B) unit of (kelvin * joules / kelvin) = joule.
// epsilon / k_B = 119.8K from Verlet and Hansen
// I believe that without reducing the Temperature I was not making a dimensionless value and was creating
// artificially high acceptance rates.

// ------------------------------------------------------------------------------
// Global Constants
// ------------------------------------------------------------------------------
double volume = (numOfSpheres*sigma*sigma*sigma) / reducedDensity;
double sizeOfSystem = pow(volume, 1./3.);

int main()
{
    auto start = get_time::now();
    vector3d *spheres = FCCLattice(numOfSpheres,sigma/2);
    cout << "Number of Spheres in System: " << numOfSpheres << ", Size of System: " << sizeOfSystem << endl;

    double dVr = 0.01*sizeOfSystem;
    int dVsteps = int((1.2*sizeOfSystem - sizeOfSystem) / (dVr));

    double energyChange[dVsteps];
    double tempEnergy[dVsteps];

    for (int i = 0; i < dVsteps; ++i){
        energyChange[i] = 0;
        tempEnergy[i] = totalPotential(spheres,numOfSpheres,i,dVr)/totalIterations;
    }

    double totalEnergy = totalPotential(spheres,numOfSpheres,0,dVr);

    double pressureIdeal = (numOfSpheres*reducedTemperature*epsilon) / volume;
    double virial = pairVirialFunction(spheres);
    double exPressure = 0.0;

    long acceptedTrials = 0;
    int runningRadial[1000];

    FILE *positions_file = fopen("MonteCarloSS.positions", "w");
    FILE *radial_file = fopen("MonteCarloSS.radial", "w");
    FILE *energy_file = fopen("MonteCarloSS.energies", "w");
    FILE *pressure_file = fopen("MonteCarloSS.pressure", "w");
    FILE *energyArray_file = fopen("MonteCarloSS.energyArray","w");
    string filename[5] = {"MonteCarloSS.positions","MonteCarloSS.radial","MonteCarloSS.energies","MonteCarloSS.pressure","MonteCarloSS.energyArray"};
    ofstream outputFile[5];
    outputFile[2].open(filename[2].c_str());
    outputFile[1].open(filename[1].c_str());
    outputFile[3].open(filename[3].c_str());
    outputFile[4].open(filename[4].c_str());

    // Performs the random move and checking
    for (long currentIteration = 0; currentIteration < totalIterations; ++currentIteration) {
        if (currentIteration >= 1000000){
            for (int n = 0; n < dVsteps; ++n){
                energyChange[n] += tempEnergy[n];
            }
        }
        
        exPressure += virial;
        
        bool trialAcceptance = false;
        // Picks a random sphere to move and performs a random move on that sphere.
        int movedSphereNum = random::ran64() % numOfSpheres;
        vector3d movedSpherePos = spheres[movedSphereNum] + vector3d::ran(dr);

        movedSpherePos = periodicBC(movedSpherePos,sizeOfSystem);
        vector3d initialSpherePos = spheres[movedSphereNum];

        double energyNew = 0.0;
        double energyOld = 0.0;
        // Determines the difference in energy due to the random move
        for (int currentSphere = 0; currentSphere < numOfSpheres; ++currentSphere) {
            if (currentSphere != movedSphereNum) {
                vector3d Rj = spheres[currentSphere];
                vector3d R2 = nearestImage(Rj,movedSpherePos,sizeOfSystem);
                vector3d R1 = nearestImage(Rj,initialSpherePos,sizeOfSystem);
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
            totalEnergy = totalPotential(spheres,numOfSpheres,0,dVr);
            virial = pairVirialFunction(spheres);
            
            if (currentIteration >= 1000000){
                for (int n = 0; n < dVsteps; ++n){
                    tempEnergy[n] = totalPotential(spheres,numOfSpheres,n,dVr);
                }
            }
            acceptedTrials += 1;
            
            if ((acceptedTrials % (totalIterations/10)) == 0){
                int *radialDistHist;
                radialDistHist = radialDistribution(spheres,numOfSpheres,sizeOfSystem);
                for (int i = 0; i < 1000; ++i){
                    runningRadial[i] = radialDistHist[i] ;
                }
            }
        }

        fprintf(energy_file, "%g\n", totalEnergy);
        if ((currentIteration == 9*totalIterations/10) ||(currentIteration == totalIterations/10) ||(currentIteration == totalIterations/4) || (currentIteration == totalIterations/2) || (currentIteration == 3*totalIterations/4))
        { // Writes out status of simulation.
            auto end = get_time::now();
            auto diff = end - start;
            cout << "Iteration: " << currentIteration << " of " << totalIterations << endl;
            cout << "Time Elapsed: " << chrono::duration_cast<sec>(diff).count()<<" sec " << endl;
            cout << "Ratio of Acceptance: " << acceptedTrials << "/" << totalIterations << endl;
            cout << endl;
        }
    }

    cout << "Ratio of Accepted to Total: " << acceptedTrials << "/" << totalIterations << endl;

    // Write Positions to File
    outputFile[0].open(filename[0].c_str());
    outputFile[0] << spheres << endl;
    for (int i=0; i<numOfSpheres; i++) {
        fprintf(positions_file, "%g\t%g\t%g\n",
                spheres[i].x, spheres[i].y, spheres[i].z);
    }

    for (int i = 0; i < 1000; ++i){
        fprintf(radial_file, "%d\t",runningRadial[i]);
    }
    fprintf(radial_file,"\n");

    for (int i = 0; i < dVsteps; ++i){
        fprintf(energyArray_file,"%g\t",energyChange[i]);
    }
    
    fprintf(pressure_file, "%g\t",pressureIdeal + (exPressure / (totalIterations*volume)));
    fprintf(energyArray_file,"\n");

    fclose(energy_file);
    fclose(positions_file);
    fclose(radial_file);
    fclose(pressure_file);
    fclose(energyArray_file);
    auto end = get_time::now();
    auto diff = end - start;
    cout << "Total Time to Completion: " << chrono::duration_cast<sec>(diff).count() << " sec " <<endl;
}

inline vector3d *FCCLattice(int totalNumOfSpheres,double sphereRadius)   {
    vector3d *sphereMatrix = new vector3d[totalNumOfSpheres];
    double xRef,yRef,zRef,xNeighbor,yNeighbor,zNeighbor;
    int xsteps,ysteps,zsteps,smallCell, breaker;
    double cubeSideLengthFactor = ceil(pow(numOfSpheres,1./3.));

    xRef = yRef = zRef = (-sizeOfSystem/2) + sphereRadius;
    xNeighbor = yNeighbor = zNeighbor = 0;
    xsteps=ysteps=zsteps = 0;

    for (int sphereNum = 0; sphereNum < totalNumOfSpheres; sphereNum++) {
        double cubeSideLength = sizeOfSystem / cubeSideLengthFactor;
        sphereMatrix[sphereNum].x = xRef + xNeighbor;
        sphereMatrix[sphereNum].y = yRef + yNeighbor;
        sphereMatrix[sphereNum].z = zRef + zNeighbor;
        breaker = sphereNum;

        xNeighbor = yNeighbor = zNeighbor = 0;
        smallCell = sphereNum%4;

        if (smallCell == 0) {
            xNeighbor = cubeSideLength / 2;
            zNeighbor = cubeSideLength / 2;
        } else if (smallCell == 1) {
            xNeighbor = cubeSideLength / 2;
            yNeighbor = cubeSideLength / 2;
        } else if (smallCell == 2) {
            yNeighbor = cubeSideLength / 2;
            zNeighbor = cubeSideLength / 2;
        } else if (smallCell == 3) {
            if (xRef + cubeSideLength < sizeOfSystem/2) {
                xRef += cubeSideLength;
                xsteps += 1;
            } else if ((xRef + cubeSideLength >= sizeOfSystem/2)
                       && (yRef + cubeSideLength  < sizeOfSystem/2)) {
                xRef -= xsteps*cubeSideLength;
                yRef += cubeSideLength;
                ysteps += 1;
                xsteps = 0;
            } else if ((yRef + cubeSideLength >= sizeOfSystem/2)
                       && (xRef + cubeSideLength >=sizeOfSystem/2)
                       && (zRef + cubeSideLength < sizeOfSystem/2)) {
                xRef -= xsteps*cubeSideLength;
                yRef -= ysteps*cubeSideLength;
                zRef += cubeSideLength;
                xsteps = ysteps = 0;
            } else if ((yRef + cubeSideLength >= sizeOfSystem/2)
                        && (xRef + cubeSideLength >=sizeOfSystem/2)
                        && (zRef + cubeSideLength >= sizeOfSystem/2)) {
                sphereNum = totalNumOfSpheres;
            }
        }
        if ((fabs(sphereMatrix[sphereNum].x) <= 0.0001) & (fabs(sphereMatrix[sphereNum].y) <= 0.0001) &
            (fabs(sphereMatrix[sphereNum].z) <=0.0001)){
            cubeSideLengthFactor += 1;
            xRef = yRef = zRef = (-sizeOfSystem/2) + sphereRadius;
            xNeighbor = yNeighbor = zNeighbor = 0;
            xsteps=ysteps=zsteps = 0;
            sphereNum = -1;
        }
        else {
            sphereNum = sphereNum;
        }
    }
    return sphereMatrix;
}

static inline vector3d periodicBC(vector3d inputVector, double sizeOfSystem)   {
    for (int i = 0; i < 3; i++){
        if (inputVector[i] > sizeOfSystem/2){
            inputVector[i] -= sizeOfSystem;
        }
        else if (inputVector[i] < -sizeOfSystem/2){
            inputVector[i] += sizeOfSystem;
        }
    }
     return inputVector;
}

double totalPotential(vector3d *sphereMatrix, int numOfSpheres, int n,double dVr) {
    double totalPotential = 0.0;

    for (int i = 0; i < numOfSpheres; ++i)  {   // Sphere 1
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j)  {   // Vector from Sphere 2
            vector3d Rj = sphereMatrix[j];
            vector3d R = (1+n*dVr)*nearestImage(Rj,Ri,sizeOfSystem);
            totalPotential += bondEnergy(R);
        }
    }
    return totalPotential;
}

int *radialDistribution(vector3d *sphereMatrix, int numOfSpheres, double sizeOfSystem) {
    // I am having an issue making the size of the dr for the histogram a variable. So I hard set it. Don't judge.
    const int numOfBoxes = 1000;
    static int deltan[numOfBoxes];

    for (int i = 0; i < numOfSpheres; ++i)
    {
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j)
        {
            vector3d Rj = sphereMatrix[j];
            vector3d R = nearestImage(Rj,Ri,sizeOfSystem);
            double Rmag = R.norm();
            int box = int(Rmag*numOfBoxes/sizeOfSystem);   // box = (R/dr) = (R / (sizeOfSystem/numOfBoxes))
            if (Rmag < sizeOfSystem){
                deltan[box] += 1;
            }
        }
    }
    return deltan;
}

vector3d nearestImage(vector3d R2,vector3d R1, double sizeOfSystem){
    vector3d R = R2 - R1;
    if (fabs(R.x) > (sizeOfSystem/2)){
        R.x -= copysign(sizeOfSystem,R.x);
    }
    if (fabs(R.y) > (sizeOfSystem/2)){
        R.y -= copysign(sizeOfSystem,R.y);
    }
    if (fabs(R.z) > (sizeOfSystem/2)){
        R.z -= copysign(sizeOfSystem,R.z);
    }
    return R;
}

double bondEnergy(vector3d R){
    double Rsq = R.normsquared();
    double SR2 = (sigma*sigma)/Rsq;
    double SR6 = SR2*SR2*SR2;
    double SR12 = SR6*SR6;
    double bondEnergy = SR12 - SR6;
    return 4*epsilon*bondEnergy + epsilon;    // Carries units of energy.
}

double pairVirialFunction(vector3d *spheres) {
    double w = 0;
    for (int i = 0; i < numOfSpheres; ++i){
        vector3d Ri = spheres[i];
        for (int j = i + 1; j < numOfSpheres; ++j){
            vector3d Rj = spheres[j];
            vector3d R = nearestImage(Rj,Ri,sizeOfSystem);
            double R2 = R.normsquared();
            if (R2 < 2*sizeOfSystem){
                double SR2 = (sigma*sigma)/R2;
                double SR6 = SR2*SR2*SR2;
                double SR12 = SR6*SR6;
                w += 2*SR12 - SR6;
            } else {
                w += 0;
            }
        }
    }
    return 8*epsilon*w;
}


