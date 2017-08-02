
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <chrono>
#include "vector3d.h"

#include "MersenneTwister.h"

using namespace std;
using sec = chrono::seconds;
using get_time = chrono::steady_clock;
//random_device rd;
//mt19937 gen(rd());
// ------------------------------------------------------------------------------
// Global Constants
// ------------------------------------------------------------------------------
    const int x = 0;
    const int y = 1;
    const int z = 2;

// ------------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------------

// Places a desired number of spheres into an FCC lattice
inline vector3d *FCCLattice();

// Applies periodic boundary conditions in x,y,z
static inline vector3d periodicBC(vector3d inputVector, double sizeOfSystem);

// Calculates total potential of system
double totalPotential(vector3d *sphereMatrix, int numOfSpheres);

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

//int numOfSpheres = 105;
long totalIterations = 0;
double dr = 0.005;
double sigma = 1;
double epsilon = 1;
double reducedDensity = 4;
double reducedTemperature = 1 ;  // (epsilon * T / k_B) unit of (kelvin * joules / kelvin) = joule.

double cutoff = sigma * pow(2,5./6.);
double systemLength[3] = {5,5,5};
double volume = systemLength[x] * systemLength[y] * systemLength[z];
int numOfSpheres = int(volume * reducedDensity);
double sizeOfSystem = pow(volume, 1./3.);

// Stuff to Remove before pushing
//vector3d ran3();
//double ran();

int main()
{
    // Stuff to Remove before pushing
//    uniform_int_distribution<int> randSphere(0,numOfSpheres-1);
//    normal_distribution<double> randMove(0,dr);

    auto start = get_time::now();
    vector3d *spheres = FCCLattice();
    cout << "Number of Spheres in System: " << numOfSpheres << ", Size of System: "
        << systemLength[x] << ", " << systemLength[y] << ", " << systemLength[z] << endl;


    double totalEnergy = totalPotential(spheres,numOfSpheres);

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

    // Performs the random move and checking
    for (long currentIteration = 0; currentIteration < totalIterations; ++currentIteration) {

        exPressure += virial;

        bool trialAcceptance = false;
        // Picks a random sphere to move and performs a random move on that sphere.
        int movedSphereNum = random::ran64() % numOfSpheres;//randSphere(gen);
        vector3d movedSpherePos = spheres[movedSphereNum] + vector3d::ran(dr);//ran3()*dr;

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
                double r = random::ran();//randMove(gen);
                if ((p > r))    {
                    trialAcceptance = true;
                }   else    {
                    trialAcceptance = false;
                }
            }
        if (trialAcceptance == true){   // If accepted, update lattice, PE, and radial dist. func.
            spheres[movedSphereNum] = movedSpherePos;
            totalEnergy = totalPotential(spheres,numOfSpheres);
            virial = pairVirialFunction(spheres);


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
    for (int i=0; i<numOfSpheres; i++) {
        fprintf(positions_file, "%g\t%g\t%g\n",
                spheres[i].x, spheres[i].y, spheres[i].z);
    }

    for (int i = 0; i < 1000; ++i){
        fprintf(radial_file, "%d\t",runningRadial[i]);
    }
    fprintf(radial_file,"\n");



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

inline vector3d *FCCLattice()   {
    vector3d *sphereMatrix = new vector3d[numOfSpheres];

    double cellNumber = ceil(pow(numOfSpheres/4,1./3.));
    double cellLength = systemLength[x]/cellNumber;

    double sphereRadius = sigma/2;
    double sizeOfSystem = 2;
    double cubeSideLength = 1;//sizeOfSystem / cubeSideLengthFactor;

    vector3d cornerSphere = {sphereRadius,sphereRadius,sphereRadius};   // Sphere's Center is radial distance from borders
    vector3d *offset = new vector3d[4];
    offset[0] = vector3d(0,0,0);
    offset[1] = vector3d(cellLength,0,cellLength)/2;
    offset[2] = vector3d(cellLength,cellLength,0)/2;
    offset[3] = vector3d(0,cellLength,cellLength)/2;

    int xsteps = 0; int ysteps = 0;int  zsteps = 0;

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
                }
                else if ((cornerSphere.x + cellLength >= systemLength[x])
                         && (cornerSphere.y + cellLength >= systemLength[y])
                         && (cornerSphere.z + cellLength >= systemLength[z])){
                            sphereNum = numOfSpheres;
                         }
        }
        if ((fabs(sphereMatrix[sphereNum].x) <= 0.0001) & (fabs(sphereMatrix[sphereNum].y) <= 0.0001) &
            (fabs(sphereMatrix[sphereNum].z) <=0.0001)){
            cout << "Cell Dimensions have shrunk " << endl;
            cornerSphere.x = cornerSphere.y = cornerSphere.z =  sphereRadius;
            xsteps=ysteps=zsteps = 0;
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

// remove n and dVr?
double totalPotential(vector3d *sphereMatrix, int numOfSpheres) {
    double totalPotential = 0.0;

    for (int i = 0; i < numOfSpheres; ++i)  {   // Sphere 1
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j)  {   // Vector from Sphere 2
            vector3d Rj = sphereMatrix[j];
            vector3d R = nearestImage(Rj,Ri,sizeOfSystem);
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
    if (R.norm() < cutoff){
        double Rsq = R.normsquared();
        double SR2 = (sigma*sigma)/Rsq;
        double SR6 = SR2*SR2*SR2;
        double SR12 = SR6*SR6;
        double bondEnergy = SR12 - SR6;
        return 4*epsilon*bondEnergy + epsilon;
    }   else {
        return 0;
        }
}

double pairVirialFunction(vector3d *spheres) {
    double w = 0;
    for (int i = 0; i < numOfSpheres; ++i){
        vector3d Ri = spheres[i];
        for (int j = i + 1; j < numOfSpheres; ++j){
            vector3d Rj = spheres[j];
            vector3d R = nearestImage(Rj,Ri,sizeOfSystem);
            double R2 = R.normsquared();
            // fixme add cutoff for WCA
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

//double ran(){
//  const long unsigned int x =0;
//  static MTRand my_mtrand(x); // always use the same random number generator (for debugging)!
//  return my_mtrand.randExc(); // which is the range of [0,1)
//}
//
//vector3d ran3(){
//  double x, y, r2,z;
//  do{
//    x = 2 * ran() - 1;
//    y = 2 * ran() - 1;
//    z = 2 * ran() - 1;
//    r2 = x * x + y * y + z * z;
//  } while(r2 >= 1 || r2 == 0);
//  double fac = sqrt(-2*log(r2)/r2);
//  vector3d out(x*fac,y*fac,z*fac);
//}

