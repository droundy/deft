#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <chrono>
#include <vector3d.h>
#include <MersenneTwister.h>

using namespace std;
using sec = chrono::seconds;
using get_time = chrono::steady_clock;

// ------------------------------------------------------------------------------
// Functions
// ------------------------------------------------------------------------------

// Places a desired number of spheres into an FCC lattice
inline vector3d *FCCLattice(int numOfSpheres,double cubeSideLength, double sizeOfSystem, double sphereRadius);

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

// Random Number Stuff which can be changed
random_device rd;
mt19937 gen(rd());

// ------------------------------------------------------------------------------
// Test Variables
// ------------------------------------------------------------------------------

int numOfSpheres = 32;   // Subject to change depending on cell size parameters
long totalIterations = 1000000;

double sigma = 1;//1;
double cellLength = 1.76*sigma;//1.8*sigma;
double sizeOfSystem = 2*cellLength;
double epsilon = 1; //
double Temperature = 85 * epsilon  /  119.8;  // (epsilon * T / k_B) unit of (kelvin * joules / kelvin) = joule.
// epsilon / k_B = 119.8K from Verlet and Hansen

double dr = 0.005;

int main()
{
    auto start = get_time::now();
    vector3d *spheres = FCCLattice(numOfSpheres,cellLength,sizeOfSystem,sigma/2);
    cout << "Number of Spheres in System: " << numOfSpheres << ". The system may have been revised due to an unrealistic amount of spheres in the desired space." << endl;

    uniform_int_distribution<int> randSphere(0,numOfSpheres-1);
    normal_distribution<double> randMove(0,dr);
    uniform_real_distribution<double> randProb(0.0,1.0);

    double totalEnergy = totalPotential(spheres,numOfSpheres);
    long acceptedTrials = 0;
    int runningRadial[1000];

    FILE *positions_file = fopen("MonteCarloSS.positions", "w");
    FILE *radial_file = fopen("MonteCarloSS.radial", "w");
    FILE *energy_file = fopen("MonteCarloSS.energies", "w");
    string filename[3] = {"MonteCarloSS.positions","MonteCarloSS.radial","MonteCarloSS.energies"};
    ofstream outputFile[3];
    outputFile[2].open(filename[2].c_str());
    outputFile[1].open(filename[1].c_str());


    // Performs the random move and checking
    for (long currentIteration = 0; currentIteration < totalIterations; ++currentIteration) {
        bool overlap = false;
        bool trialAcceptance = false;
        int movedSphereNum = randSphere(gen);
        vector3d movedSpherePos = spheres[movedSphereNum];

        for (int i = 0; i < 3; ++i){   // Store to temporary array for comparison
            movedSpherePos[i] += randMove(gen);//dr*randGauss();
        }
        movedSpherePos = periodicBC(movedSpherePos,sizeOfSystem);
        vector3d initialSpherePos = spheres[movedSphereNum];

        double energyNew = 0.0;
        double energyOld = 0.0;
        for (int currentSphere = 0; currentSphere < numOfSpheres; ++currentSphere) {   // Determines the difference in energy due to the random move
            if (currentSphere != movedSphereNum) {
                vector3d Rj = spheres[currentSphere];
                vector3d R2 = nearestImage(Rj,movedSpherePos,sizeOfSystem);
                vector3d R1 = nearestImage(Rj,initialSpherePos,sizeOfSystem);
                energyNew += bondEnergy(R2);
                energyOld += bondEnergy(R1);
            }
        }

        double energyChange = energyNew - energyOld;

        if ((currentIteration % (totalIterations / 10)) ==0){
            cout << "Current Iteration: " << currentIteration << endl;
            cout << "Energy Change: " << energyChange << endl;
            double p = exp(-energyChange / Temperature);
            cout << "If energy change > 0, p: " << p << endl;
            cout << endl;
        }


        if (energyChange <= 0)  {
            trialAcceptance = true;
        }
        else if (energyChange > 0)  {
            double p = exp(-energyChange / Temperature); // Need to get units correct in here.
            double r = randProb(gen);
            if ((p > r))    {
                trialAcceptance = true;
            }
            else    {
                trialAcceptance = false;
            }
        }
        if (trialAcceptance == true){   // Rewrites position and energy files
            spheres[movedSphereNum] = movedSpherePos;
            totalEnergy = totalPotential(spheres,numOfSpheres);
            acceptedTrials += 1;

            if ((acceptedTrials%1000) == 0){
                int *radialDistHist;
                radialDistHist = radialDistribution(spheres,numOfSpheres,sizeOfSystem);
                for (int i = 0; i < 1000; ++i){
                    runningRadial[i] = radialDistHist[i];
                    fprintf(radial_file, "%d\t",runningRadial[i]);
                }
                fprintf(radial_file,"\n");
            }
        }

        fprintf(energy_file, "%g\n", totalEnergy);
        outputFile[2] << totalEnergy << endl;
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

    fclose(energy_file);
    fclose(positions_file);
    fclose(radial_file);
}

inline vector3d *FCCLattice(int totalNumOfSpheres,double cubeSideLength, double sizeOfSystem,double sphereRadius)   {
    vector3d *sphereMatrix = new vector3d[totalNumOfSpheres];
    double xRef,yRef,zRef,xNeighbor,yNeighbor,zNeighbor;
    int xsteps,ysteps,zsteps,smallCell, breaker;

    xRef = yRef = zRef = sphereRadius;
    xNeighbor = yNeighbor = zNeighbor = 0;
    xsteps=ysteps=zsteps = 0;

    for (int sphereNum = 0; sphereNum < totalNumOfSpheres; sphereNum++) {
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
            if (xRef + cubeSideLength < sizeOfSystem) {
                xRef += cubeSideLength;
                xsteps += 1;
            } else if ((xRef + cubeSideLength >= sizeOfSystem)
                       && (yRef + cubeSideLength  < sizeOfSystem)) {
                xRef -= xsteps*cubeSideLength;
                yRef += cubeSideLength;
                ysteps += 1;
                xsteps = 0;
            } else if ((yRef + cubeSideLength >= sizeOfSystem)
                       && (xRef + cubeSideLength >=sizeOfSystem)
                       && (zRef + cubeSideLength < sizeOfSystem)) {
                xRef -= xsteps*cubeSideLength;
                yRef -= ysteps*cubeSideLength;
                zRef += cubeSideLength;
                xsteps = ysteps = 0;
            } else if ((yRef + cubeSideLength >= sizeOfSystem)
                        && (xRef + cubeSideLength >=sizeOfSystem)
                        && (zRef + cubeSideLength >= sizeOfSystem)) {
                sphereNum = totalNumOfSpheres;
            } else {
                printf(" Help!!!\n");
                cout << " Help!!!" << endl;
            }
        }
    }
    return sphereMatrix;
}

static inline vector3d periodicBC(vector3d inputVector, double sizeOfSystem)   {
    for (int i = 0; i < 3; i++){
        if (inputVector[i] > sizeOfSystem){
            inputVector[i] -= sizeOfSystem;
        }
        else if (inputVector[i] < 0){
            inputVector[i] += sizeOfSystem;
        }
    }
     return inputVector;
}

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

    for (int currentSphere = 0; currentSphere < numOfSpheres; ++currentSphere)
    {
        vector3d R1 = sphereMatrix[currentSphere];
        for (int testSphere = currentSphere + 1; testSphere < numOfSpheres; ++testSphere)
        {
            vector3d R2 = sphereMatrix[testSphere];
            vector3d R = nearestImage(R2,R1,sizeOfSystem);
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
