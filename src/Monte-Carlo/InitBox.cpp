//#include "stdafx.h"
#include <random>
#include <iostream>
#include <math.h>
#include <string.h>
#include <fstream>
#include <ctime>
#include "InitBox.h"

using namespace std;
std::random_device dummy;
std::seed_seq seed{ dummy(),dummy(),dummy(),dummy(),dummy(),dummy(),dummy(),dummy() };
std::mt19937 generator(seed);
std::uniform_real_distribution<double> myRand(0.0, 1.0);
std::uniform_int_distribution<int> myIntRand(0, 9);
//###########################

ATOM::ATOM(double a, double b, double c, double sig, double lam) { x = a; y = b; z = c; sigma = sig; lambda = lam; };
//ATOM::ATOM() { ; };
double ATOM::distanceToAtom(ATOM &atom) {
  double distance = (atom.x - x)*(atom.x - x) + (atom.y - y)*(atom.y - y) + (atom.z - z)*(atom.z - z);
  return distance;//must work with r^2 now...
}
//############################




INITBOX::INITBOX() {//not used... not tested...
  lx = 10.154912976;
  ly = 10.154912976;
  lz = 10.154912976;
  temperature = 1.5;
  sigma = 2.0;//hard coded for now  sigma=2*r
  lambda = 1.5;
  dStep = sigma / 5.0*3.0;
  maxStep = lx*2;//have to double check lx? or use int conversion and multiply? hmm
  binSizeX = lx/sigma/lambda;
  binSizeY = ly/sigma/lambda;
  binSizeZ = lz/sigma/lambda;
  resetNumberOfBonds();
  numAtoms=0;
  numAvailableInAtomList=500;
  list=new ATOM[numAvailableInAtomList];
  for(int n=0;n<numAvailableInAtomList;n++) list[n].id=n;
  return;
}
INITBOX::INITBOX(double L,int N){
  numAtoms=0;
  lx=L;
  ly=L;
  lz=L;
  sigma = 2.0;//hard coded for now  sigma=2*r
  lambda = 1.5;
  dStep = sigma / 5.0*3.0;
  maxStep = lx*2;//have to double check lx? or use int conversion and multiply? hmm
  binSizeX = lx/sigma/lambda;
  binSizeY = ly/sigma/lambda;
  binSizeZ = lz/sigma/lambda;
  resetNumberOfBonds();
  list=new ATOM[2];
  for(int count=0;count<10 && numAtoms!=N;count++){
    //INITBOX(L);
    delete [] list;
    numAtoms=0;
    temperature=1.1;
    fillBoxWithRain(N);
    //delete box;
  }
  temperature=1000.0;
  if(numAtoms!=0) for(int count=0;count<numAtoms+30;count++) simulate(numAtoms*1000);
  numAvailableInAtomList=0;
}
INITBOX::INITBOX(double L) {
  lx = L;
  ly = L;
  lz = L;
  temperature = 1.5;
  sigma = 2.0;//hard coded for now  sigma=2*r
  lambda = 1.5;
  dStep = sigma / 5.0*3.0;
  maxStep = lx*2;//have to double check lx? or use int conversion and multiply? hmm
  binSizeX = lx/sigma/lambda;
  binSizeY = ly/sigma/lambda;
  binSizeZ = lz/sigma/lambda;
  resetNumberOfBonds();
  numAtoms=0;
  numAvailableInAtomList=500;
  list=new ATOM[numAvailableInAtomList];
  for(int n=0;n<numAvailableInAtomList;n++) list[n].id=n;
  for (int nx = 0; nx < 30; nx++)
    for (int ny = 0; ny < 30; ny++)
      for (int nz = 0; nz < 30; nz++) {
        bins[nx][ny][nz] = -1;
        numInBins[nx][ny][nz] = 0;
      }
}
void INITBOX::resetNumberOfBonds(void) {
  for (unsigned int i = 0; i < sizeof(numberOfBonds) / sizeof(*numberOfBonds); i++) {
    numberOfBonds[i] = 0.0;
    for (int j = 0; j < 1024; j++) bondRadiusBins[i][j] = 0;
  }
}
void INITBOX::addAtom(double x, double y, double z){
  if(numAvailableInAtomList<=0){
    class ATOM *list2;
    list2=new ATOM[numAtoms+501];
    for(int n=0;n<numAtoms;n++) list2[n]=list[n];//save old list
    for(int n=0;n<numAtoms+501;n++) list2[n].id=n;//make sure new list has id's initialized
    for(int n=0;n<numAtoms;n++) removeAtomFromBins(list[n]);//take old list out of the bins
    numAvailableInAtomList=501;
    delete[] list;
    list=list2;
    for(int n=0;n<numAtoms;n++) addAtomToBins(list[n]);//put new list back into the bins
  }
  numAtoms++;
  numAvailableInAtomList--;
  list[numAtoms-1].x=x;
  list[numAtoms-1].y=y;
  list[numAtoms-1].z=z;
  list[numAtoms-1].id=numAtoms-1;
  addAtomToBins(list[numAtoms-1]);
  return;
}
bool INITBOX::randStep() {
  std::uniform_int_distribution<int> randAtom(0, numAtoms-1);
  int n=randAtom(generator);//pick random atom
  double prevEnergy = atomEnergy(n);//energy to take the atom out of the box
  std::normal_distribution<double>randStep(0.0, dStep);//bad overloaded name, sorry
  double dx = maxStep * 3;
  while (abs(dx)>maxStep) dx=randStep(generator);
  double dy = 3 * maxStep;
  while (abs(dy)>maxStep) dy=randStep(generator);
  double dz = 3 * maxStep;
  while(abs(dz)>maxStep) dz=randStep(generator);
  /*double r = sqrt(dx*dx + dy*dy + dz*dz);
    dx = dStep*dx / r;
    dy = dStep*dy / r;
    dz = dStep*dz / r;*/

  removeAtomFromBins(list[n]);
  list[n].x += dx;
  list[n].y += dy;
  list[n].z += dz;
  addAtomToBins(list[n]);
  double nextEnergy = atomEnergy(n);
  /*if (nextEnergy < 0.0) {
    removeAtomFromBins(list[n]);
    list[n].x -= dx;
    list[n].y -= dy;
    list[n].z -= dz;
    addAtomToBins(list[n]);
    if (atomEnergy(n) < 0.0) cout << "critical FAILURE WITH RANDSTEP" << endl;
    return false;
    }*/
  double acceptanceRate = exp((nextEnergy - prevEnergy) / temperature)*(1-2*(nextEnergy<0));
  std::uniform_real_distribution<double> randRate(0.0, 1.0);
  int reality = (randRate(generator) <= acceptanceRate);
  dx = (1 - reality)*dx;
  dy = (1 - reality)*dy;
  dz = (1 - reality)*dz;
  removeAtomFromBins(list[n]);
  list[n].x -= dx;
  list[n].y -= dy;
  list[n].z -= dz;
  addAtomToBins(list[n]);

  //cout << "modified ds: " << sqrt(dx*dx+dy*dy+dz*dz) << endl;
  /*if (nextEnergy < prevEnergy && reality==1) {
    cout << "reality: " << reality << endl;
    cout << "prevE: " << prevEnergy << " nextE: " << nextEnergy << endl;
    }*/
  return reality;//1=success 0=fail
}
double INITBOX::totalWallEnergyX(double x0) {
  double energy = 0.0;
  for (int n = 0; n < numAtoms; n++) {
    energy += atomWallEnergyX(n, x0);
  }
  return energy;
}
double INITBOX::totalWallEnergyYplane(double y0) {
  double energy = 0.0;
  for (int n = 0; n < numAtoms; n++) {
    energy += atomWallEnergyYplaneY(n, y0);
  }
  return energy;
}
double INITBOX::totalWallEnergy(void) {
  int nx = binSizeX - 1;
  double energy = 0.0;
  for (int ny = 0; ny < binSizeY; ny++) {
    for (int nz = 0; nz < binSizeZ; nz++) {
      if (numInBins[nx][ny][nz] <= 0) continue;
      int n = bins[nx][ny][nz];
      energy += atomWallEnergy(n);
      ATOM *current = &list[n];
      while (current->next != NULL) {
        current = current->next;
        energy += atomWallEnergy(current->id);
      }
    }
  }
  return energy;
}
double INITBOX::totalWallEnergyYplane(void) {
  int ny = binSizeY - 1;
  double energy = 0.0;
  for (int nx = 0; nx < binSizeX; nx++) {
    for (int nz = 0; nz < binSizeZ; nz++) {
      if (numInBins[nx][ny][nz] <= 0) continue;
      int n = bins[nx][ny][nz];
      energy += atomWallEnergyYplane(n);
      ATOM *current = &list[n];
      while (current->next != NULL) {
        current = current->next;
        energy += atomWallEnergyYplane(current->id);
      }
    }
  }
  return energy;
}

double INITBOX::atomWallEnergyYplaneY(int n, double y0) {
  //special case of energy method, checks if against the right wall
  //shifts normally in the y-z direction, but only shifts nx to the right
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];

  double modLy = -(atom->y > y0)*ly;//I think I can just shift the boundary to the right by ly, try later
  atom->y += modLy;//make sure target atom is below the boundary
  if (atom->y + sigma*lambda < y0) {//can't reach the boundary
    atom->y -= modLy;
    return 0.0;
  }
  //shift the atom back into the box, but also shift the boundary to compensate
  double modBoundaryLy = (atom->y < 0)*ly;
  atom->y += modBoundaryLy;//shift the atom back into the box
  y0 += modBoundaryLy;//maintain boundary relative to atom

  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;

  for (int di = -1; di <= 1; di++) {
    int modNx = -(nx + di >= binSizeX);// look at bin[x=nx+di] but ensure nx+di is in box
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = 0; dj <= 1; dj++) {//boundary is either in current box or up by one

      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      y0 += modNy*ly;//keep boundary relative to atom
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = -1; dk <= 1; dk++) {
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          if (current->y <= y0) {
            current = current->next;
            continue;
          }
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          if (current->id == atom->id) r = (sigma*sigma + 1)*lambda*lambda;//no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;
            atom->y -= modNy*ly+modBoundaryLy+modLy;//unshift atom
            atom->z -= modNz*lz;
            y0 -= modNy*ly + modBoundaryLy;//unshift boundary (why? I'm returning -1.0)
            cout << "found hard shell impact, BAAAAAAD" << endl;
            return -1.0;//bad hard shell impact
          }
          if (r < sigma*sigma*lambda*lambda) energy += wellDepth;
          current = current->next;//###########################this needs to be in loop iteration
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
      y0 -= modNy*ly;
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  atom->y -= modLy + modBoundaryLy;//undo shift of atom to left side, and shift to get back into bin
  y0 -= modBoundaryLy;
  return energy;
}
double INITBOX::atomWallEnergyYplane(int n) {
  //special case of energy method, checks if against the top wall
  //shifts normally in the x-z direction, but only shifts ny up
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];
  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;
  if (ny < binSizeY - 1) return 0;//want it agains the top wall
  for (int di = -1; di <= 1; di++) {
    int modNx = -(nx + di >= binSizeX);//shift in X dir
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = 1; dj <= 1; dj++) {
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = -1; dk <= 1; dk++) {
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          if (current->id == atom->id) r = (sigma*sigma + 1)*lambda*lambda;//no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz;
            return -1.0;//bad hard shell impact
          }
          if (r < sigma*sigma*lambda*lambda) energy += wellDepth;
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  return energy;
}
double INITBOX::totalWallEnergyZplane(double y0) {
  double energy = 0.0;
  for (int n = 0; n < numAtoms; n++) {
    energy += atomWallEnergyZplaneZ(n, y0);
  }
  return energy;
}
double INITBOX::totalWallEnergyZplane(void) {
  int nz = binSizeZ - 1;
  double energy = 0.0;
  for (int nx = 0; nx < binSizeX; nx++) {
    for (int ny = 0; ny < binSizeY; ny++) {
      if (numInBins[nx][ny][nz] <= 0) continue;
      int n = bins[nx][ny][nz];
      energy += atomWallEnergyZplane(n);
      ATOM *current = &list[n];
      while (current->next != NULL) {
        current = current->next;
        energy += atomWallEnergyZplane(current->id);
      }
    }
  }
  return energy;
}
double INITBOX::atomWallEnergyZplaneZ(int n, double z0) {
  //special case of energy method, checks if against the right wall
  //shifts normally in the y-z direction, but only shifts nx to the right
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];

  double modLz = -(atom->z > z0)*lz;//I think I can just shift the boundary to the back by lz, try later
  atom->z += modLz;//make sure target atom is t of the boundary
  if (atom->z + sigma*lambda < z0) {//can't reach the boundary
    atom->z -= modLz;
    return 0.0;
  }
  //shift the atom back into the box, but also shift the boundary to compensate
  double modBoundaryLz = (atom->z < 0)*lz;
  atom->z += modBoundaryLz;//shift the atom back into the box
  z0 += modBoundaryLz;//maintain boundary relative to atom

  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;

  for (int di = -1; di <= 1; di++) {
    int modNx = -(nx + di >= binSizeX);// look at bin[x=nx+di] but ensure nx+di is in box
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box

    //now boundary is at x0=x0+modBoundary*lx+modNx*lx
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = -1; dj <= 1; dj++) {
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = 0; dk <= 1; dk++) {//boundary is either in current box or to the back

        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        z0 += modNz*lz;//keep boundary relative to atom
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          if (current->z <= z0) {
            current = current->next;
            continue;
          }

          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          if (current->id == atom->id) r = (sigma*sigma + 1)*lambda*lambda;//no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz + modBoundaryLz + modLz;
            z0 -= modNz*lz + modBoundaryLz;//unshift boundary (why? returning -1.0)
            cout << "found hard shell impact, BAAAAAAD" << endl;
            return -1.0;//bad hard shell impact
          }
          if (r < sigma*sigma*lambda*lambda) energy += wellDepth;
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
        z0 -= modNz*lz;
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  atom->z -= modLz + modBoundaryLz;//undo shift of atom to left side, and shift to get back into bin
  z0 -= modBoundaryLz;  //why?
  return energy;
}
double INITBOX::atomWallEnergyZplane(int n) {
  //special case of energy method, checks if against back z-wall
  //shifts normally in the x-y direction, but only shifts nz to the back
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];
  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;
  if (nz < binSizeZ - 1) return 0;//want it agains the back wall
  for (int di = -1; di <= 1; di++) {
    int modNx = -(nx + di >= binSizeX);//shift in X dir
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = -1; dj <= 1; dj++) {
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = 1; dk <= 1; dk++) {
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          if (current->id == atom->id) r = (sigma*sigma + 1)*lambda*lambda;//no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz;
            return -1.0;//bad hard shell impact
          }
          if (r < sigma*sigma*lambda*lambda) energy += wellDepth;
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  return energy;
}
double INITBOX::atomWallEnergyX(int n, double x0) {
  //special case of energy method, checks if against the right wall
  //shifts normally in the y-z direction, but only shifts nx to the right
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];
  
  double modLx = -(atom->x > x0)*lx;//I think I can just shift the boundary to the right by lx, try later
  atom->x += modLx;//make sure target atom is to the left of the boundary
  if (atom->x + sigma*lambda < x0) {//can't reach the boundary
    atom->x -= modLx;
    return 0.0;
  }
  //shift the atom back into the box, but also shift the boundary to compensate
  double modBoundaryLx = (atom->x < 0)*lx;
  atom->x += modBoundaryLx;//shift the atom back into the box
  x0 += modBoundaryLx;//maintain boundary relative to atom
  
  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;
  
  for (int di = 0; di <= 1; di++) {//boundary is either in the current box or to the right
    int modNx = -(nx + di >= binSizeX);// look at bin[x=nx+di] but ensure nx+di is in box
    //modNx += (nx + di < 0);//this won't happen (only when di=-1..1
    atom->x += modNx*lx;//shift atom rather than shift box
    x0 += modNx*lx;//keep boundary relative to atom
    //now boundary is at x0=x0+modBoundary*lx+modNx*lx
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = -1; dj <= 1; dj++) {
      
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = -1; dk <= 1; dk++) {
        
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          if (current->x <= x0) {
            current = current->next;
            continue;
          }
          
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          if (current->id == atom->id) r = (sigma*sigma + 1)*lambda*lambda;//no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx + modBoundaryLx + modLx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz;
            x0 -= modNx*lx + modBoundaryLx;//unshift boundary
            cout << "found hard shell impact, BAAAAAAD" << endl;
            return -1.0;//bad hard shell impact
          }
          if (r < sigma*sigma*lambda*lambda) energy += wellDepth;
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
    x0 -= modNx*lx;
  }
  atom->x -= modLx + modBoundaryLx;//undo shift of atom to left side, and shift to get back into bin
  x0 -= modBoundaryLx;
  return energy;
}
double INITBOX::atomWallEnergy(int n) {
  //special case of energy method, checks if against the right wall
  //shifts normally in the y-z direction, but only shifts nx to the right
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];
  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;
  if (nx < binSizeX - 1) return 0;//want it agains the right wall##############changed#########
  for (int di = 1; di <= 1; di++) {
    int modNx = -(nx + di >= binSizeX);//shift in X dir
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = -1; dj <= 1; dj++) {
      
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = -1; dk <= 1; dk++) {
        
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          if (current->id == atom->id) r = (sigma*sigma + 1)*lambda*lambda;//no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz;
            return -1.0;//bad hard shell impact
          }
          if (r < sigma*sigma*lambda*lambda) energy += wellDepth;
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  return energy;
}
void INITBOX::testEnergy(void) {
  for (int n = 0; n < numAtoms; n++) {
    if(atomEnergy(n)<0) cout<<"found negative energy, FAIL"<<endl;
  }
}
void INITBOX::testAtomEnergy() {
  double energy = 0.0;
  for (int n = 0; n < numAtoms; n++) {
    ATOM *atom = &list[n];
    double prevX = atom->x, prevY = atom->y, prevZ = atom->z;
    double energyTemp = atomEnergy(n);
    energy += energyTemp;
    if (abs(prevX - atom->x)>1e-15) { cout << "bad energy alg." << endl; break; };
    if (abs(prevY - atom->y)>1e-15) { cout << "bad energy alg." << endl; break; };
    if (abs(prevZ - atom->z)>1e-15) { cout << "bad energy alg." << endl; break; };
    cout << energyTemp << endl;
  }
  cout << "total Energy: " << energy << endl;
}
bool INITBOX::testAllInBox(void) {
  //tests if all atoms are in the box
  for (int i = 0; i < numAtoms; i++) {
    if (list[i].x < 0) return false;
    if (list[i].x > lx) return false;
    if (list[i].y < 0) return false;
    if (list[i].y > ly) return false;
    if (list[i].z < 0) return false;
    if (list[i].z > lz) return false;
  }
  return true;
}
double INITBOX::totalEnergy(void) {
  double energy = 0.0;
  for (int n = 0; n < numAtoms; n++) {
    double atomE = atomEnergyWithRadiusTracking(n);
    energy += atomE;
    int nBond = (atomE / wellDepth + 0.5);
    if (nBond > 30) cout << "found bondNumber>30, BAAD        "<<endl;
    else numberOfBonds[nBond] += 1;
  }
  return energy / 2.0;
}
void INITBOX::testWallEnergy(void) {
  std::uniform_real_distribution<double> lxRand(0.0, lx);
  std::uniform_real_distribution<double> lyRand(0.0, ly);
  double averageEnergy = 0.0;
  double averageEnergyY = 0.0;
  double averageEnergyZ = 0.0;
  for (int i = 0; i < 1000; i++) {
    averageEnergy += randWallEnergyX();
    averageEnergyY += randWallEnergyY();
    averageEnergyZ += randWallEnergyZ();
  }
  averageEnergy = averageEnergy / 1000.0;
  averageEnergyY = averageEnergyY / 1000.0;
  averageEnergyZ = averageEnergyZ / 1000.0;
  cout << "regular: " << totalWallEnergy() << " right wall: " << totalWallEnergyX(lx)
       << " left wall: " << totalWallEnergyX(0.0)
       << " averageOf1000: " << averageEnergy << endl;
  cout << "yPlaneTop: " << totalWallEnergyYplane();
  cout << " TopY(ly): " << totalWallEnergyYplane(ly);
  cout << " bottomY(ly): " << totalWallEnergyYplane(0.0);
  cout << " averageOf1000Y: " << averageEnergyY << endl;
  cout << "zPlaneBack: " << totalWallEnergyZplane();
  cout << " Back(lz): " << totalWallEnergyZplane(lz);
  cout << " front(lz): " << totalWallEnergyZplane(0.0);
  cout << " averageOf1000Z: " << averageEnergyZ << endl << endl;
  /*cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
    cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
    cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
    cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
    cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
    cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
    cout << "some random x0: " << totalWallEnergyX(lxRand(generator)) << endl;
  */
}
bool INITBOX::testBinnedAtoms(void) {
  //tests if all atoms are in the bins
  int *atomNumber = new int[numAtoms];
  for (int n = 0; n < numAtoms; n++) atomNumber[n] = 0;
  int totalAtoms = 0;
  //first check if there is exactly one entry per atom in the bins
  for (int nx = 0; nx < binSizeX; nx++) {
    for (int ny = 0; ny < binSizeY; ny++) {
      for (int nz = 0; nz < binSizeZ; nz++) {
        totalAtoms += numInBins[nx][ny][nz];
        if (bins[nx][ny][nz] < 0) continue;
        //cout << nx << " " << ny << " " << nz <<" numInBins="<<numInBins[nx][ny][nz]<< " n=" << bins[nx][ny][nz] << endl;
        ATOM *current = &list[bins[nx][ny][nz]];
        for (int n = 1; n <= numInBins[nx][ny][nz]; n++) {
          atomNumber[current->id]++;
          current = current->next;
        }
      }
    }
  }
  //cout << "totalAtoms: " << totalAtoms << endl;
  for (int n = 0; n < numAtoms; n++) if (atomNumber[n] != 1) { cout << n << " num:" << atomNumber[n] << endl; delete[] atomNumber; return false; }
  delete[] atomNumber;
  //now that I know there is exactly one entry per atom
  //test if each atom is actually in the proper bin
  //I assume the coordinates will work without checking
  //will crash for bad coordinates
  for (int n = 0; n < numAtoms; n++) {
    int nx, ny, nz;
    ATOM *atom = &list[n];
    nx = (atom->x / lx)*binSizeX;
    ny = (atom->y / ly)*binSizeY;
    nz = (atom->z / lz)*binSizeZ;
    //the atom should be at bin [nx,ny,nz]...
    bool testIfInBox = false;//starts off as not verified
    if (bins[nx][ny][nz] < 0) return false;//none in the box.... failed
    ATOM *current = &list[bins[nx][ny][nz]];
    for (int i = 0; i < numInBins[nx][ny][nz]; i++) {
      if (current->id == n) {
        testIfInBox = true; break;
      }
      current = current->next;
    }
    if (testIfInBox == false) return false;
    //cout << "verified: " << n << endl;
  }
  return true;
}
double INITBOX::atomEnergyWithRadiusTracking(int n) {
  double energy = 0.0;
  double radiusBonds[30] = { 0 };//stores the radius of the nth bond, maximum of 30 bonds assumed
  int numBonds = 0;
  int nx, ny, nz;
  ATOM *atom = &list[n];
  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;
  for (int di = -1; di <= 1; di++) {
    
    int modNx = -(nx + di >= binSizeX);//shift in X dir
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = -1; dj <= 1; dj++) {
      
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = -1; dk <= 1; dk++) {
        
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          //if (current->id == atom->id) r=(sigma+1)*lambda;//no self interactions
          r += (current->id == atom->id)*(sigma*sigma + 1)*lambda*lambda;//ensures no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz;
            //this is official totalEnergy method
            //informing end user of hard shell impact...
            cout << "found hard shell IMPACT!!!!!!!!!!!!!" << endl;
            return -1.0;//bad hard shell impact
          }
          //if (r < sigma*lambda) energy += wellDepth;
          if (r < sigma*sigma*lambda*lambda) {
            energy += wellDepth;
            radiusBonds[numBonds] = sqrt(r);
            numBonds++;
          }
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  for (int n = 0; n < numBonds; n++) {
    int binNumber = abs(radiusBonds[n] - sigma) / (sigma*lambda - sigma) * 1024;
    binNumber -= (binNumber >= 1024);//odd cases where round off doesn't work for us
    binNumber += (binNumber < 0);//don't need since abs()
    bondRadiusBins[numBonds][binNumber]++;
  }
  return energy;
}
double INITBOX::atomEnergy(int n) {
  double energy = 0.0;
  int nx, ny, nz;
  ATOM *atom = &list[n];
  nx = (atom->x / lx)*binSizeX;
  ny = (atom->y / ly)*binSizeY;
  nz = (atom->z / lz)*binSizeZ;
  for (int di = -1; di <= 1; di++) {
    
    int modNx = -(nx + di >= binSizeX);//shift in X dir
    modNx += (nx + di < 0);
    atom->x += modNx*lx;//shift atom rather than shift box
    int NX = nx + di + binSizeX*modNx;//shift to the proper bin
    for (int dj = -1; dj <= 1; dj++) {
      
      int modNy = -(ny + dj >= binSizeY);//shift in Y dir
      modNy += (ny + dj < 0);
      atom->y += modNy*ly;
      int NY = ny + dj + binSizeY*modNy;
      for (int dk = -1; dk <= 1; dk++) {
        
        int modNz = -(nz + dk >= binSizeZ);//shift in Z dir
        modNz += (nz + dk < 0);
        int NZ = nz + dk + binSizeZ*modNz;
        if (numInBins[NX][NY][NZ] <= 0) continue;//verify there is work to be done
        atom->z += modNz*lz;
        ATOM *current = &list[bins[NX][NY][NZ]];//grab first on bin stack
        for (int i = 1; i <= numInBins[NX][NY][NZ]; i++) {//start summing up energies
          
          double r = atom->distanceToAtom(*current);//code changed to work with r^2
          //if (current->id == atom->id) r=(sigma+1)*lambda;//no self interactions
          r += (current->id == atom->id)*(sigma*sigma + 1)*lambda*lambda;//ensures no self interactions
          //if (r<3) cout << "r: " << r << endl;
          if (r < sigma*sigma) {
            atom->x -= modNx*lx;//unshift atom
            atom->y -= modNy*ly;
            atom->z -= modNz*lz;
            return -1.0;//bad hard shell impact
          }
          //if (r < sigma*lambda) energy += wellDepth;
          energy += wellDepth*(r < sigma*sigma*lambda*lambda);
          current = current->next;
        }
        atom->z -= modNz*lz;//unshift atom z
      }
      atom->y -= modNy*ly;//unshift atom y
    }
    atom->x -= modNx*lx;//unshift atom x
  }
  
  return energy;
}
bool INITBOX::randStepRain(double ddStep) {
  std::uniform_int_distribution<int> randAtom(0, numAtoms-1);
  int n=randAtom(generator);//pick random atom
  double radius2=(list[n].x-lx/2.0)*(list[n].x-lx/2.0);
  radius2+=(list[n].y-ly/2.0)*(list[n].y-ly/2.0);
  radius2+=(list[n].z-lz/2.0)*(list[n].z-lz/2.0);
  radius2=sqrt(radius2);
  double prevEnergy = atomEnergy(n)+30*radius2;//energy to take the atom out of the box
  std::normal_distribution<double>randStep(0.0, ddStep);
  double dx = maxStep * 3;
  while (abs(dx)>maxStep) dx=randStep(generator);
  double dy = 3 * maxStep;
  while (abs(dy)>maxStep) dy=randStep(generator);
  double dz = 3 * maxStep;
  while(abs(dz)>maxStep) dz=randStep(generator);
  /*double r = sqrt(dx*dx + dy*dy + dz*dz);
    dx = dStep*dx / r;
    dy = dStep*dy / r;
    dz = dStep*dz / r;*/
  
  removeAtomFromBins(list[n]);
  list[n].x += dx;
  list[n].y += dy;
  list[n].z += dz;
  addAtomToBins(list[n]);
  double nextEnergy = atomEnergy(n);
  if(nextEnergy<0.0){
    removeAtomFromBins(list[n]);
    list[n].x-=dx;
    list[n].y-=dy;
    list[n].z-=dz;
    addAtomToBins(list[n]);
    return false;
  }
  radius2=(list[n].x-lx/2.0)*(list[n].x-lx/2.0);
  radius2+=(list[n].y-ly/2.0)*(list[n].y-ly/2.0);
  radius2+=(list[n].z-lz/2.0)*(list[n].z-lz/2.0);
  radius2=sqrt(radius2);
  nextEnergy+=30*radius2;
  /*if (nextEnergy < 0.0) {
    removeAtomFromBins(list[n]);
    list[n].x -= dx;
    list[n].y -= dy;
    list[n].z -= dz;
    addAtomToBins(list[n]);
    if (atomEnergy(n) < 0.0) cout << "critical FAILURE WITH RANDSTEP" << endl;
    return false;
    }*/
  double acceptanceRate = exp((nextEnergy - prevEnergy) / temperature)*(1-2*(nextEnergy<0));
  std::uniform_real_distribution<double> randRate(0.0, 1.0);
  int reality = (randRate(generator) <= acceptanceRate);
  dx = (1 - reality)*dx;
  dy = (1 - reality)*dy;
  dz = (1 - reality)*dz;
  removeAtomFromBins(list[n]);
  list[n].x -= dx;
  list[n].y -= dy;
  list[n].z -= dz;
  addAtomToBins(list[n]);
  
  //cout << "modified ds: " << sqrt(dx*dx+dy*dy+dz*dz) << endl;
  /*if (nextEnergy < prevEnergy && reality==1) {
    
    cout << "reality: " << reality << endl;
    cout << "prevE: " << prevEnergy << " nextE: " << nextEnergy << endl;
    }*/
  return reality;//1=success 0=fail
  
}
int INITBOX::fillBoxWithRain(int N){
  list=new ATOM[N];
  numAvailableInAtomList=0;
  for(int n=0;n<N;n++){
    list[n].id=n;
    //cout<<list[n].x<<endl;
    //cout<<lx<<endl;
  }
  for (int nx = 0; nx < 30; nx++)
    for (int ny = 0; ny < 30; ny++)
      for (int nz = 0; nz < 30; nz++) {
        bins[nx][ny][nz] = -1;
        numInBins[nx][ny][nz] = 0;
      }
  std::normal_distribution<double>randStepX(0.0, lx/2.0/3.0);
  std::normal_distribution<double>randStepY(0.0, ly/2.0/3.0);
  std::normal_distribution<double>randStepZ(0.0, lz/2.0/3.0);
  numAtoms=0;
  //cout<<"putting in first atom"<<endl;
  //cout <<"lx="<<lx<<" x="<<list[numAtoms].x<<endl;
  list[numAtoms].x=lx*3.0;
  //cout <<"lx="<<lx<<" x="<<list[numAtoms].x<<endl;
  while(list[numAtoms].x<0.0 || list[numAtoms].x>lx) {
    list[numAtoms].x=lx/2.0+randStepX(generator);
    //cout <<"lx="<<lx<<" x="<<list[numAtoms].x<<endl;
  }
  list[numAtoms].y=ly*3.0;
  while(list[numAtoms].y<0.0 || list[numAtoms].y>ly) list[numAtoms].y=ly/2.0+randStepY(generator);
  list[numAtoms].z=lz*3.0;
  while(list[numAtoms].z<0.0 || list[numAtoms].z>lz) list[numAtoms].z=lz/2.0+randStepZ(generator);
  numAtoms=1;
  
  double ddstep=maxStep/6.0;
  double T=temperature;
  for(temperature=T;temperature>0.1;temperature-=0.1){
    int attemptsLeft=100;
    while(numAtoms<N && attemptsLeft>0){
      //cout<<"starting loop"<<endl;
      int acceptance=0;
      for(int count=0;count<numAtoms*1000;count++) acceptance+=randStepRain(ddstep);
      if( (acceptance/1000.0)/numAtoms>0.3 && ddstep<maxStep/3.0) ddstep*=1.1;
      if( (acceptance/1000.0)/numAtoms<0.2 ) ddstep*=0.9;
      bool failFlag=true;
      for(int count=0; count<100; count++){
        list[numAtoms].x=lx*3.0;
        while(list[numAtoms].x<0.0 || list[numAtoms].x>lx) list[numAtoms].x=lx/2.0+randStepX(generator);
        list[numAtoms].y=ly*3.0;
        while(list[numAtoms].y<0.0 || list[numAtoms].y>ly) list[numAtoms].y=ly/2.0+randStepY(generator);
        list[numAtoms].z=lz*3.0;
        while(list[numAtoms].z<0.0 || list[numAtoms].z>lz) list[numAtoms].z=lz/2.0+randStepZ(generator);
        if(atomEnergy(numAtoms)>=0.0){
          addAtomToBins(list[numAtoms]);
          numAtoms++;
          //cout<<"numAtoms="<<numAtoms<<endl;
          failFlag=false;
          break;
        }
      }
      if(failFlag==true) attemptsLeft--;
    }
  }
  temperature=T;
  return numAtoms;
}
int INITBOX::fillBoxFccPlus(int N){
  int nLeft=fillBoxFcc(N);
  if (nLeft<=0) return 0;
  int current=N-nLeft;
  double lScale=sqrt(2.0)*2+1e-14;
  double L=lx;
  int nMax=(int)(L/lScale);
  //zplane
  {
    int nz=nMax;
    double z=lScale*nz;
    for(int nx=0;nx<=nMax && current<N;nx++){
      double x=lScale*nx;
      for(int ny=0;ny<=nMax && current<N;ny++){
        double y=lScale*ny;
        list[current].x=x;
        list[current].y=y;
        list[current].z=z;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x+lScale/2.0;
        list[current].y=y+lScale/2.0;
        list[current].z=z;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x+lScale/2.0;
        list[current].y=y;
        list[current].z=z+lScale/2.0;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x;
        list[current].y=y+lScale/2.0;
        list[current].z=z+lScale/2.0;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        
      }
    }
  }
  //yPlane
  {
    int ny=nMax;
    double y=lScale*ny;
    for(int nx=0;nx<=nMax && current<N;nx++){
      double x=lScale*nx;
      for(int nz=0;nz<=nMax && current<N;nz++){
        double z=lScale*nz;
        
        list[current].x=x;
        list[current].y=y;
        list[current].z=z;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x+lScale/2.0;
        list[current].y=y+lScale/2.0;
        list[current].z=z;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x+lScale/2.0;
        list[current].y=y;
        list[current].z=z+lScale/2.0;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x;
        list[current].y=y+lScale/2.0;
        list[current].z=z+lScale/2.0;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
	
      }
    }
  }
  //xPlane
  {
    int nx=nMax;
    double x=lScale*nx;
    
    for(int nz=0;nz<=nMax && current<N;nz++){
      double z=lScale*nz;
      for(int ny=0;ny<=nMax && current<N;ny++){
        double y=lScale*ny;
        list[current].x=x;
        list[current].y=y;
        list[current].z=z;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x+lScale/2.0;
        list[current].y=y+lScale/2.0;
        list[current].z=z;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x+lScale/2.0;
        list[current].y=y;
        list[current].z=z+lScale/2.0;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
        list[current].x=x;
        list[current].y=y+lScale/2.0;
        list[current].z=z+lScale/2.0;
        while(list[current].x>L) list[current].x-=L;
        while(list[current].y>L) list[current].y-=L;
        while(list[current].z>L) list[current].z-=L;
        list[current].id=current;
        if(atomEnergy(current)>=0.0){
          addAtomToBins(list[current]);
          current++;
          if(current>=N) break;
        }
      }
    }
  }
  return N-current;
}
int INITBOX::fillBoxFcc(int N){
  //assumes a square box
  numAtoms=N;
  list=new ATOM[N];
  numAvailableInAtomList=0;
  for (int nx = 0; nx < 30; nx++)
    for (int ny = 0; ny < 30; ny++)
      for (int nz = 0; nz < 30; nz++) {
        bins[nx][ny][nz] = -1;
        numInBins[nx][ny][nz] = 0;
      }
  
  double lScale=sqrt(2.0)*2+1e-14;
  double L=lx;
  if(L<lScale)
    return N;
  int nMax=(int)(L/lScale);
  
  bool ****fccAtoms;//wooooo, should have used one page of pixels
  fccAtoms=new  bool***[nMax];
  for(int i=0;i<nMax;i++){//allocate bool for fcc array
    //not all lattice sites will be poplated
    //deterimine which ones will be populated first
    //then create the atoms at those lattice sites
    fccAtoms[i]=new bool**[nMax];
    for(int j=0;j<nMax;j++){
      fccAtoms[i][j]=new bool*[nMax];
      for(int k=0;k<nMax;k++){
        fccAtoms[i][j][k]=new bool[4];
        fccAtoms[i][j][k][0]=false;
        fccAtoms[i][j][k][1]=false;
        fccAtoms[i][j][k][2]=false;
        fccAtoms[i][j][k][3]=false;
      }
    }
  }
  int count=0;//now randomly choose lattice sites to populate
  while(count<N && count<nMax*nMax*nMax*4){
    count++;
    std::uniform_int_distribution<int> randN(1,nMax*nMax*nMax*4-count+1);
    bool foundIndex=false;
    int baseCount=0;
    int countId=randN(generator);
    for(int nx=0;nx<nMax && foundIndex==false;nx++){
      for(int ny=0;ny<nMax && foundIndex==false;ny++){
        for(int nz=0;nz<nMax && foundIndex==false;nz++){
          for(int offset1234=0;offset1234<4 ;offset1234++){
            baseCount+=(int)(fccAtoms[nx][ny][nz][offset1234]==false);
            if(baseCount==countId) {
              fccAtoms[nx][ny][nz][offset1234]=true;
              foundIndex=true; 
              break;}
            
          }
        }
      }
    }
  }
  
  //lattice sites chosen, now populate
  int baseCount=0;
  for(int nx=0;nx<nMax;nx++){
    double x=nx*lScale;
    for(int ny=0;ny<nMax;ny++){
      double y=ny*lScale;
      for(int nz=0;nz<nMax;nz++){
        double z=nz*lScale;
        if(fccAtoms[nx][ny][nz][0]==true){
          list[baseCount].x=x;
          list[baseCount].y=y;
          list[baseCount].z=z;
          list[baseCount].id=baseCount;
          addAtomToBins(list[baseCount]);
          baseCount++;
        }
        if(fccAtoms[nx][ny][nz][1]==true){
          list[baseCount].x=x+lScale/2.0;
          list[baseCount].y=y+lScale/2.0;
          list[baseCount].z=z;
          list[baseCount].id=baseCount;
          addAtomToBins(list[baseCount]);
          baseCount++;
        }
        if(fccAtoms[nx][ny][nz][2]==true){
          list[baseCount].x=x+lScale/2.0;
          list[baseCount].y=y;
          list[baseCount].z=z+lScale/2.0;
          list[baseCount].id=baseCount;
          addAtomToBins(list[baseCount]);
          baseCount++;
        }
        if(fccAtoms[nx][ny][nz][3]==true){
          list[baseCount].x=x;
          list[baseCount].y=y+lScale/2.0;
          list[baseCount].z=z+lScale/2.0;
          list[baseCount].id=baseCount;
          addAtomToBins(list[baseCount]);
          baseCount++;
        }
      }
    }
  }
  
  for(int i=0;i<nMax;i++){
    for(int j=0;j<nMax;j++){
      
      for(int k=0;k<nMax;k++){
        delete [] fccAtoms[i][j][k];
      }
      delete [] fccAtoms[i][j];
    }
    delete [] fccAtoms[i];
  }
  delete [] fccAtoms;
  return N-baseCount;
}
double INITBOX::fillBox(int N) {
  numAtoms = N;
  list = new ATOM[N];
  numAvailableInAtomList=0;
  std::uniform_real_distribution<double> lxRand(0.0, lx);
  std::uniform_real_distribution<double> lyRand(0.0, ly);
  std::uniform_real_distribution<double> lzRand(0.0, lz);
  
  //memset(bins, -1, sizeof(bins));
  //memset(numInBins, 0, sizeof(numInBins));
  for (int nx = 0; nx < 30; nx++)
    for (int ny = 0; ny < 30; ny++)
      for (int nz = 0; nz < 30; nz++) {
        bins[nx][ny][nz] = -1;
        numInBins[nx][ny][nz] = 0;
      }
  for (int n = 0; n < N; n++) {
    list[n].x = lxRand(generator);
    list[n].y = lyRand(generator);
    list[n].z = lzRand(generator);
    list[n].id = n;
    addAtomToBins(list[n]);
  }
  return attemptRelaxation();
  
}
double INITBOX::attemptRelaxation(void) {
  //trys to lower the contact area
  //TO DO: add stochastic noise like Brownian motion
  int failLower = 0;
  double prevError = 1.1;
  int stepDown = 10;
  for (int i = 0; i < 100; i += 1) {
    //cout << "testing for fail TRY: " << i << " stepDown:" << stepDown << endl;
    double error = percentFail();
    if (error > prevError) {
      failLower++;
      if (failLower >= 1) {
        failLower = 0;
        stepDown *= 2;
      }
    }
    prevError = error;
    cout <<"attempting relaxation, step: "<<i<<" of 100, error Fraction: "<< error <<flush<<'\r';
    if (error == 0.0) return 0.0;
    //cout << "moving boxes to minimize fail" << endl;
    for (int j = 0; j < 800; j += 1) {
      stepAway(stepDown);
    }
  }
  if (percentFail() > 0.0) cout << "detected failure in the relaxation, sorry..." << endl;
  return percentFail();
}

double INITBOX::randWallEnergy(void) {
  std::uniform_int_distribution<int> randAlgorithm(0, 2);
  switch (randAlgorithm(generator)) {
  case 0: return randWallEnergyX();
  case 1: return randWallEnergyY();
  case 2: return randWallEnergyZ();
  }
  return 0.0;//?shouldn't go here...
}
double INITBOX::randWallEnergyX(void) {
  std::uniform_real_distribution<double> lxRand(0.0, lx);
  return totalWallEnergyX(lxRand(generator));
}
double INITBOX::randWallEnergyY(void) {
  std::uniform_real_distribution<double> lyRand(0.0, ly);
  return totalWallEnergyYplane(lyRand(generator));
}
double INITBOX::randWallEnergyZ(void) {
  std::uniform_real_distribution<double> lzRand(0.0, lz);
  return totalWallEnergyZplane(lzRand(generator));
}


double INITBOX::percentFail(void) {
  int fail = 0;
  for (int n = 0; n < numAtoms; n++) {
    double minR = 10.0;
    for (int i = 0; i < numAtoms; i++) {
      if (i == n) continue; //no self interactions
      bool testFail = false;
      for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
          for (int dk = -1; dk <= 1; dk++) {
            list[i].x += di*lx;
            list[i].y += dj*ly;
            list[i].z += dk*lz;
            double r = list[n].distanceToAtom(list[i]);//code changed to work with r^2
            if (r < minR) minR = r;//used to display minimum radius per atom, not used atm
            //if(r<1.0) cout << r << endl;
            list[i].x -= di*lx;
            list[i].y -= dj*ly;
            list[i].z -= dk*lz;
            if (r < sigma*sigma) testFail = true;
          }
        }
      }
      
      if (testFail==true) {
        fail++;
        break;
      }
    }
    //cout << "minR=" << minR << endl;
  }
  return ((double)fail)/numAtoms;
}
INITBOX::~INITBOX() {
  delete [] list;
}
void INITBOX::stepAway(int downStep) {
  double ds = lx / binSizeX/downStep;
  for (int n = 0; n < numAtoms; n++) {
    double fx = 0.0;
    double fy = 0.0;
    double fz = 0.0;
    for (int i = 0; i < numAtoms; i++) {
      if (i == n) continue;
      
      
      //scan over all possible shifts (periodic boundary)
      for (int di = -1; di <= 1; di++) {
        for (int dj = -1; dj <= 1; dj++) {
          for (int dk = -1; dk <= 1; dk++) {
            list[n].x += di*lx;
            list[n].y += dj*ly;
            list[n].z += dk*lz;
            double r = list[n].distanceToAtom(list[i]);//code changed to work with r^2
            if (r < 2.0*2.0* sigma*sigma) {
              r = r*r*r*r;
              fx += (list[n].x - list[i].x) / (r + 1);
              fy += (list[n].y - list[i].y) / (r + 1);
              fz += (list[n].z - list[i].z) / (r + 1);
            }
            list[n].x -= di*lx;
            list[n].y -= dj*ly;
            list[n].z -= dk*lz;
          }
        }
      }
      
    }
    removeAtomFromBins(list[n]);
    list[n].x += fx*ds;
    list[n].y += fy*ds;
    list[n].z += fz*ds;
    if (list[n].x > lx) list[n].x -= lx;
    if (list[n].x < 0) list[n].x += lx;
    if (list[n].y > ly) list[n].y -= ly;
    if (list[n].y < 0) list[n].y += ly;
    if (list[n].z > lz) list[n].z -= lz;
    if (list[n].z < 0) list[n].z += lz;
    addAtomToBins(list[n]);
  }
}
void INITBOX::simulate(int sample, int iterations,double expectedFinish) {
  int acceptanceCount = 1;
  int count = 4;//starts at a ratio of 25% which is good
  for (int n = 0; n < iterations / 8; n++) {
    //if ((n+1) % (100000/8) == 0) {
    /*cout << "sample: " << sample //<< " donePercent: " << ((1.0*n) / iterations)
      << " acceptanceRate: " << (1.0*acceptanceCount) / (count) << " dstep: " << dStep
      << flush << '\r';*/
    //box->testWallEnergy();
    //fflush(stdout);
    //}
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    //box->testEnergy();
    if ((n + 1) % (100000 / 8) == 0) {
      if ((1.0*acceptanceCount) / count>0.3 && dStep<maxStep / 3.0) dStep *= 1.1;
      if ((1.0*acceptanceCount) / count < 0.2) dStep *= 0.9;
    }
  }
  if ((1.0*acceptanceCount) / count>0.3 && dStep<maxStep / 3.0) dStep *= 1.1;
  if ((1.0*acceptanceCount) / count < 0.2) dStep *= 0.9;
  cout << "sample: " << sample << " ETA: " << expectedFinish
       << " acceptanceRate: " << (1.0*acceptanceCount) / (count) << " dstep: " << dStep << "     "
       << flush << '\r';
  //residual iterations
  for (int n = 8 * iterations; n < iterations; n++) randStep();
}
void INITBOX::simulate(int sample,int iterations) {
  int acceptanceCount = 1;
  int count = 4;//starts at a ratio of 25% which is good
  for (int n = 0; n < iterations/8; n++) {
    //if ((n+1) % (100000/8) == 0) {
    /*cout << "sample: " << sample //<< " donePercent: " << ((1.0*n) / iterations)
      << " acceptanceRate: " << (1.0*acceptanceCount) / (count) << " dstep: " << dStep
      << flush << '\r';*/
    //box->testWallEnergy();
    //fflush(stdout);
    //}
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    //box->testEnergy();
    if ((n + 1) % (100000/8) == 0) {
      if ((1.0*acceptanceCount) / count>0.3 && dStep<maxStep / 3.0) dStep *= 1.1;
      if ((1.0*acceptanceCount) / count < 0.2) dStep *= 0.9;
    }
  }
  if ((1.0*acceptanceCount) / count>0.3 && dStep<maxStep / 3.0) dStep *= 1.1;
  if ((1.0*acceptanceCount) / count < 0.2) dStep *= 0.9;
  cout << "sample: " << sample //<< " donePercent: " << ((1.0))
       << " acceptanceRate: " << (1.0*acceptanceCount) / (count) << " dstep: " << dStep<<"               "
       << flush << '\r';
  //residual iterations
  for (int n = 8 * iterations; n < iterations; n++) randStep();
}
void INITBOX::simulate(int iterations) {
  int acceptanceCount = 1;
  int count = 4;//starts at a ratio of 25% which is good
  
  for (int n = 0; n < iterations/8; n++) {
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
    count++;
    if (randStep() == 1) acceptanceCount++;
  }
  for (int n = 8 * iterations; n < iterations; n++) randStep();//don't care about acceptance on residual
  if ((1.0*acceptanceCount) / count>0.3 && dStep<maxStep / 3.0) dStep *= 1.1;
  if ((1.0*acceptanceCount) / count < 0.2) dStep *= 0.9;
}

void INITBOX::removeAtomFromBins(ATOM &atom) {
  //assumes the atom was already in the bins
  int nx, ny, nz;
  nx = (atom.x / lx)*binSizeX;
  ny = (atom.y / ly)*binSizeY;
  nz = (atom.z / lz)*binSizeZ;
  numInBins[nx][ny][nz]--;
  //if (numInBins[nx][ny][nz] < 0) cout << "detected critical error in removeAtoms" << endl;
  if (numInBins[nx][ny][nz] <=0) {
    bins[nx][ny][nz] = -1;
    return;
  }
  if (bins[nx][ny][nz] == atom.id) {
    atom.next->prev = NULL;
    bins[nx][ny][nz] = atom.next->id;
    return;
  }
  atom.prev->next = atom.next;
  if (atom.next == NULL) return;
  atom.next->prev = atom.prev;
  return;
}
void INITBOX::addAtomToBins(ATOM &atom) {
  int nx, ny, nz;
  atom.x += (atom.x < 0)*lx;
  atom.y += (atom.y < 0)*ly;
  atom.z += (atom.z < 0)*lz;
  atom.x -= (atom.x >= lx)*lx;
  atom.y -= (atom.y >= ly)*ly;
  atom.z -= (atom.z >= lz)*lz;
  //added twice to support maxStep=2*lx , should really use integer arthimatic to be valid for maxStep=1000000*lx
  atom.x += (atom.x < 0)*lx;
  atom.y += (atom.y < 0)*ly;
  atom.z += (atom.z < 0)*lz;
  atom.x -= (atom.x >= lx)*lx;
  atom.y -= (atom.y >= ly)*ly;
  atom.z -= (atom.z >= lz)*lz;
  
  nx = (atom.x / lx)*binSizeX;
  ny = (atom.y / ly)*binSizeY;
  nz = (atom.z / lz)*binSizeZ;
  if (bins[nx][ny][nz] == -1) {
    bins[nx][ny][nz] = atom.id;
    atom.prev = NULL;
    atom.next = NULL;
    numInBins[nx][ny][nz]=1;
    return;
  }
  ATOM *current = &list[bins[nx][ny][nz]];
  for (int i = 1; i < numInBins[nx][ny][nz]; i++) {
    current = current->next;
  }
  current->next = &atom;
  atom.prev = current;
  atom.next = NULL;
  numInBins[nx][ny][nz]++;
  return;
}
