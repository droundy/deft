//#include "stdafx.h"



//###########################
class ATOM { 
public:
	class ATOM *prev, *next;
	int id;
	double x, y, z,sigma,lambda;
	ATOM(double a, double b, double c, double sig, double lam) ;
	ATOM() { ; };
	double distanceToAtom(ATOM &atom) ;
};
//############################

class INITBOX {
public:

	class ATOM *list;
	int numAtoms,numInBins[30][30][30];
	int binSizeX, binSizeY, binSizeZ;
	int bins[30][30][30];//hard coded for now, contains the atom id of the first atom
	double lx, ly, lz;
	double sigma, lambda,dStep,maxStep;
	double wellDepth;
	double temperature;
	int numberOfBonds[30];
	int bondRadiusBins[30][1024];
	INITBOX();
	INITBOX(double L,int N);
	INITBOX(double L);
	void resetNumberOfBonds(void);
	void addAtom(double x, double y, double z);
	bool randStep();
	double totalWallEnergyX(double x0);
	double totalWallEnergyYplane(double y0);
	double totalWallEnergy(void);
	double totalWallEnergyYplane(void);
	
	double atomWallEnergyYplaneY(int n, double y0);
	double atomWallEnergyYplane(int n);
	double totalWallEnergyZplane(double y0);
	double totalWallEnergyZplane(void);
	double atomWallEnergyZplaneZ(int n, double z0);
	double atomWallEnergyZplane(int n);
	double atomWallEnergyX(int n, double x0);
	double atomWallEnergy(int n);
	void testEnergy(void);
	void testAtomEnergy();
	bool testAllInBox(void);
	double totalEnergy(void);
	void testWallEnergy(void);
	bool testBinnedAtoms(void);
	double atomEnergyWithRadiusTracking(int n);
	double atomEnergy(int n);
	bool randStepRain(double ddStep);
	int fillBoxWithRain(int N);
	int fillBoxFccPlus(int N);
	int fillBoxFcc(int N);
	double fillBox(int N);
	double attemptRelaxation(void);
	
	double randWallEnergy(void);
	double randWallEnergyX(void);
	double randWallEnergyY(void);
	double randWallEnergyZ(void);
	
	
	double percentFail(void);
	~INITBOX();
	void stepAway(int downStep);
	void simulate(int sample, int iterations,double expectedFinish);
	void simulate(int sample,int iterations);
	void simulate(int iterations);

private:
	void removeAtomFromBins(ATOM &atom);
	void addAtomToBins(ATOM &atom);
	int numAvailableInAtomList;
};
