#define _USE_MATH_DEFINES
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <random>
#include <math.h>
#include <chrono>
#include <vector3d.h>

using namespace std;
using sec = chrono::seconds;
using get_time = chrono::steady_clock;

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
inline vector3d *FCCLattice(int numOfSpheres,double cubeSideLength, double sizeOfSystem, double sphereRadius);

// Applies periodic boundary conditions in x,y,z
inline vector3d periodicBC(vector3d inputVector, double sizeOfSystem);

// Calculates total potential of system
double totalPotential(vector3d *sphereMatrix, int numOfSpheres);

// Checks whether the move of a sphere is valid
bool conditionCheck(vector3d *sphereMatrix, vector3d movedSpherePos, int numOfSpheres, int movedSphereNum,double Temperature);

// Counts the number of spheres within a certain radius. Does not take into consideration mirror image particles.
int *radialDistribution(vector3d *sphereMatrix, int numOfSpheres, double sizeOfSystem);

// Random Number Stuff which can be changed
random_device rd;
mt19937 gen(rd());

// ------------------------------------------------------------------------------
// Test Variables
// ------------------------------------------------------------------------------

int numOfSpheres = 32;   // Subject to change depending on cell size parameters
int totalIterations = 10000000;

double sigma = 3.45;
double cellLength = 6*sigma;
double sizeOfSystem = 2*cellLength;
double Temperature = 84;

double epsilon = 1;
double WCA = 0;

double dr = 0.01;

int main()
{
    auto start = get_time::now();
    vector3d *sphereMatrix = FCCLattice(numOfSpheres,cellLength,sizeOfSystem,sigma/2);
    cout << "Number of Spheres in System: " << numOfSpheres << ". The system may have been revised due to an unrealistic amount of spheres in the desired space." << endl;

    uniform_int_distribution<int> randSphere(0,numOfSpheres-1);
    normal_distribution<double> randMove(0,dr);
    uniform_real_distribution<double> randProb(0.0,1.0);

    int acceptedTrials = 0;

    double totalEnergy = totalPotential(sphereMatrix,numOfSpheres);
    FILE *positions_file = fopen("MonteCarloSS.positions", "w");
    FILE *radial_file = fopen("MonteCarloSS.radial", "w");
    FILE *energy_file = fopen("MonteCarloSS.energies", "w");
    string filename[3] = {"MonteCarloSS.positions","MonteCarloSS.radial","MonteCarloSS.energies"};
    ofstream outputFile[3];
    outputFile[2].open(filename[2].c_str());
    outputFile[1].open(filename[1].c_str());

    double energy1 = 0.0;
    double energy2 = 0.0;
    double R1sq, SR12, SR16, SR112;
    double R2sq, SR22, SR26, SR212;
    vector3d R1,R2,Rj,initialRow;
    bool trialAcceptance;

    // Performs the random move and checking
    for (int currentIteration = 0; currentIteration < totalIterations; ++currentIteration) {
        bool overlap = false;
        int movedSphereNum = randSphere(gen);
        vector3d movedSpherePos = sphereMatrix[movedSphereNum];

        for (int i = 0; i < 3; ++i){   // Store to temporary array for comparison
            movedSpherePos[i] += randMove(gen);
        }
        movedSpherePos = periodicBC(movedSpherePos,sizeOfSystem);
        initialRow = sphereMatrix[movedSphereNum];

        for (int currentSphere = 0; currentSphere < numOfSpheres; ++currentSphere) {   // Determines the difference in energy due to the random move
            if (currentSphere != movedSphereNum) {
                Rj = sphereMatrix[currentSphere];
                // Calculates Energy after Sphere move
                R2 = Rj - movedSpherePos;
                R2sq = R2.norm();
                SR22 = (sigma*sigma)/R2sq;
                SR26 = SR22*SR22*SR22;
                SR212 = SR26*SR26;
                energy2 += SR212 - SR26;
                // Calculates Energy before Sphere move.
                R1 = Rj - initialRow;
                R1sq = R1.norm();
                SR12 = (sigma*sigma)/R1sq;
                SR16 = SR12*SR12*SR12;
                SR112 = SR16*SR16;
                energy1 += SR112 - SR16;

                if (sqrt(R2sq) < sigma/2){
                    overlap = true;
                }
            }
        }

        double energyChange = energy2 - energy1;
        if (energyChange <= 0)  {
            trialAcceptance = true;
        }
        else if (energyChange > 0)  {
            double p = exp(-energyChange / Temperature);
            double r = randProb(gen);
            if ((p > r)&&(overlap == false))    {
                trialAcceptance = true;
            }
            else    {
                trialAcceptance = false;
            }
        }


        if (trialAcceptance == true){   // Rewrites position and energy files
            sphereMatrix[movedSphereNum] = movedSpherePos;
            totalEnergy = totalPotential(sphereMatrix,numOfSpheres);
            acceptedTrials += 1;
        }

        if ((currentIteration % 1000) == true){
           int *radialDistHist;

           radialDistHist = radialDistribution(sphereMatrix,numOfSpheres,sizeOfSystem);

           for (int i = 0; i < 1000; i++ ) {
                fprintf(radial_file, "%d\t",radialDistHist[i]);
           }
           fprintf(radial_file,"\n");
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
    outputFile[0] << sphereMatrix << endl;
    for (int i=0; i<numOfSpheres; i++) {
        fprintf(positions_file, "%g\t%g\t%g\n",
                sphereMatrix[i].x, sphereMatrix[i].y, sphereMatrix[i].z);
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

inline vector3d periodicBC(vector3d inputVector, double sizeOfSystem)   {
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
    double Rsq, SR2, SR6, SR12;
    vector3d R,Rj;

    for (int i = 0; i < numOfSpheres; ++i)  {   // Sphere 1
        vector3d Ri = sphereMatrix[i];
        for (int j = i + 1; j < numOfSpheres; ++j)  {   // Vector from Sphere 2
            Rj = sphereMatrix[j];
            R = Rj - Ri;
            Rsq = R.norm();
            SR2 = (sigma*sigma)/Rsq;
            SR6 = SR2*SR2*SR2;
            SR12 = SR6*SR6;
            totalPotential += SR12 - SR6;
        }
    }
    return totalPotential = 4*epsilon*totalPotential + WCA;
}

int *radialDistribution(vector3d *sphereMatrix, int numOfSpheres, double sizeOfSystem)   {
    vector3d R1,R2;
    double R;
    // I am having an issue making the size of the dr for the histogram a variable.
    const int numOfBoxes = 1000;
    static int deltan[numOfBoxes];

    for (int currentSphere = 0; currentSphere < numOfSpheres; ++currentSphere)
    {
        R1 = sphereMatrix[currentSphere];
        for (int testSphere = currentSphere + 1; testSphere < numOfSpheres; ++testSphere)
        {
            R2 = sphereMatrix[testSphere];
            R = (R2-R1).norm();
            int box = int(R*numOfBoxes/sizeOfSystem);   // box = (R/dr) = (R / (sizeOfSystem/numOfBoxes))
            if (R < sizeOfSystem){
                deltan[box] += 1;
            }
        }
    }
    return deltan;
}



