/*
 ============================================================================
 Name        : PDIRadTianSCC3.cpp
 Author      : Abdullah-Al Mamun
 Date		 : April 9, 2013
 Version     :
 Copyright   : BECAT
 Description : MPI C++ project

 sort (remove duplicacy)
 blocking
 sorting of blocks
 shuffle of blocks
 square prefix sum of blocks
 generating edgelist
 find connected components

 last editted on Aug 13, 2013 (final)
 ============================================================================
 */



#include <iostream>
#include <mpi.h>
#include <fstream>
#include <libxml/parser.h>
#include <libxml/xpath.h>
#include <boost/algorithm/string.hpp>
#include <algorithm>
#include <iterator>
#include <vector>
#include <boost/regex.hpp>
#include <math.h>
#include <string>

using namespace std;
using namespace MPI;
using namespace boost;


// MACRO
#define FILE_CONFIG			"sample.xml"
#define	RECORD_TOTAL		6000000
#define	FILE_SAME			0
#define L_MER				4
#define ALPHABET_SIZE		256
#define ROOT_ID				0
#define LINE_LENGTH_MIN		10
#define CLUSTER_FOR_MASTER	100
#define DATA_INVALID		-1


// structure
typedef struct
{
	string str;
	int ind;
}StrPacket;

typedef struct
{
	int ind, clusterInd;
}recordPacket;


// function declaration
vector<int> arrangeBlockListRand(int dataTotal);
vector<int> arrangeBlockListSorted(vector<int> dataArr);
vector<int> arrangeBlockList(vector<int> dataArr);
vector<int> generateRandSeq();

void outputCluster(vector<vector<int> > clusterExactIndMasterArr, vector<int> rootArr);
vector<int> findConnComp(vector<int> indArr, int pointTotal);
int findRoot(int pointID, vector<int> &parentArr);
void makeEquivalent(int rootU, int rootV, vector<int> &parentArr, vector<int> &weightArr);

string convertToLower(string str);
int convertStrToInd(string str);
void getInputFileNameList(vector<string> &fileNameArr);
void getInputComparisonPara(vector<vector<int> > &attrEditArr);
int getInputThreshold();
int getLNInd();
void readDataFromFile(vector<string> fileNameArr);

void sortData();
void findExactCluster(vector<int> &clusterExactIndArr);


vector<vector<int> > createBlock();
vector<vector<int> > splitIntoBlock(vector<int> blockMasterArr, vector<int> blockIndArr, vector<int> indArr, vector<int> fracArr, vector<int> blockCountArr);
void createClusterEdgeList(vector<vector<int> > blockArr);
void generateEdgilist(vector<int> blockRowArr, int blockInd, int blockTotal);
int calculateDistTotal(vector<string> a, vector<string> b);
int calculateEditDist(vector<string> a, vector<string> b, vector<vector<int> > attrEditArr);
int calculateBasicED(string str1, string str2, int threshRem);
int calculateBasicED2(string str1, string str2, int threshRem);


// global variable
unsigned long long checkTemp, checkTemp2;
double checkTemp3, checkTemp4, checkTemp5;
int clusterExactTotalThis, clusterTotal, procID, procTotal, threshold, fracStart, fracLast;
vector<vector<int> > attrEditArr;
vector<vector<string> > recordArr;
vector<recordPacket> recordExactIndArr;
vector<int> edgeArr, rootArr;
vector<StrPacket> recordStrArr;

double loopCount = 0, editCount = 0, loopExp = 0;
int blockcountTest = 0, blockcountMTest = 0;

// main function
int main(int argc, char* argv[])
{
	// variable declaration
	clock_t startT, startExactT, startApproxT, startSortT, startBlockT, startSCCT, startDataT, startBCastT, startCommT, startGenListT, startMergeT, startDistT, startMasterT;
	double totalT, totalSortT, totalBlockT, totalDataT, totalExactT, totalApproxT, totalSCCT, totalBCastT, totalCommT, totalGenListT, totalMergeT, totalDistT, totalMasterT;
	int recordTotal, blockTotal;
	int i, j, k, procIdSource, dataTemp, dataTemp2, clusterMasterSize, dataSum, recordExactTotal, count, ind;
	vector<string> fileNameArr;
	vector<int> blockArr, blockMasterArr, tempDataMasterArr, fracArr;
	vector<vector<int> > clusterApproxIndArr, clusterMasterArr, clusterExactIndMasterArr, blockTempArr;
	vector<int> blockStartIndArr, blockIndArr, blockIndTempArr, clusterExactIndArr, tempIndArr, tempDataArr, clusterSerialIndArr, procSerialLenArr, clusterPointMasterArr, clusterPointNextArr, strideArr, clusterProcArr;

	Status status;
	Datatype typeRecordPacket, typeOldArr[1];
	Aint offsetArr[1];
	int blockCountArr[1];

	// initialize MPI
	Init(argc, argv);

	procTotal	= COMM_WORLD.Get_size();
	procID		= COMM_WORLD.Get_rank();


	// custom datatype for record
	offsetArr[0]		= 0;
	typeOldArr[0]		= INT;
	blockCountArr[0]	= 2;
	typeRecordPacket	= Datatype::Create_struct(1, blockCountArr, offsetArr, typeOldArr);
	typeRecordPacket.Commit();

	// read data and config files
	if(procID == ROOT_ID)
		startDataT		= clock();
	getInputFileNameList(fileNameArr);
	threshold	= getInputThreshold();
	getInputComparisonPara(attrEditArr);

	readDataFromFile(fileNameArr);
	recordTotal		= recordArr.size();

	// initialize time
	totalSortT		= 0.0;
	totalBlockT		= 0.0;
	totalDataT		= 0.0;
	totalExactT		= 0.0;
	totalApproxT	= 0.0;
	totalSCCT		= 0.0;

	totalBCastT		= 0.0;
	totalCommT		= 0.0;
	totalGenListT	= 0.0;
	totalMergeT		= 0.0;
	totalDistT		= 0.0;
	totalMasterT	= 0.0;

	if(procID == ROOT_ID)
	{
		totalT			= 0.0;


		totalDataT		= (double)(clock() - startDataT) / CLOCKS_PER_SEC;
		startT	= clock();
	}



	startExactT		= clock();
	startSortT		= clock();

	// sort data using radix sort of whole string considering first 2-characters
	sortData();
	totalSortT		= (double)(clock() - startSortT) / CLOCKS_PER_SEC;

	// find exact clusters using duplicates
	findExactCluster(clusterExactIndArr);
	totalExactT		= (double)(clock() - startExactT) / CLOCKS_PER_SEC;

	//cout << procID << " exacct " << clusterExactIndArr.size() << endl;
	if(procID == ROOT_ID)
	{

		startApproxT	= clock();

	}

	if(procID == ROOT_ID)
	{
		procSerialLenArr.resize(procTotal, 0);
		//clusterSerialIndArr.resize(procTotal);
		//clusterProcArr.resize(procTotal);

	}

	dataTemp	= clusterExactIndArr.size();


	COMM_WORLD.Barrier();
	startCommT	= clock();
	// gather exact cluster length
	COMM_WORLD.Gather(&dataTemp, 1, INT, &procSerialLenArr.front(), 1, INT, ROOT_ID);
	totalCommT	= (double)(clock() - startCommT) / CLOCKS_PER_SEC;

	if(procID == ROOT_ID)
	{
		startMasterT	= clock();

		dataSum		= 0;
		for(i = 0; i < procTotal; ++i)
		{
			strideArr.push_back(dataSum);
			dataSum		+= procSerialLenArr.at(i);
			//cout << dataSum << " total " << procSerialLenArr.at(i) << endl;
		}
		tempIndArr.resize(dataSum);

		totalMasterT	+= (double)(clock() - startMasterT) / CLOCKS_PER_SEC;
	}


	COMM_WORLD.Barrier();
	startCommT	= clock();
	// gather exact cluster indices
	COMM_WORLD.Gatherv(&clusterExactIndArr.front(), dataTemp, INT, &tempIndArr.front(), &procSerialLenArr.front(), &strideArr.front(), INT, ROOT_ID);
	totalCommT	+= (double)(clock() - startCommT) / CLOCKS_PER_SEC;
	//COMM_WORLD.Send(&clusterSerialIndArr.front(), clusterSerialIndArr.size(), INT, ROOT_ID, 0);

	if(procID == ROOT_ID)
	{
		startMasterT	= clock();

//		clusterProcArr.resize(procTotal);
//		for(i = 0; i < procTotal; ++i)
//			COMM_WORLD.Recv(&clusterProcArr.front(), procSerialLenArr.at(i), INT, i, ANY_TAG, status);
		//cout << tempIndArr.size() << " clusterExactIndMasterArr " << endl;

		dataTemp2	= tempIndArr.size();
		dataTemp	= 0;
		for(i = 0; i < dataTemp2; )
		{
			dataTemp	= tempIndArr.at(i++);

			while(dataTemp > 0)
			{
				tempDataArr.push_back(tempIndArr.at(i++));
				--dataTemp;
			}
			clusterExactIndMasterArr.push_back(tempDataArr);
			tempDataArr.clear();
		}

		recordExactTotal	= clusterExactIndMasterArr.size();
//		for(i = 0; i < recordExactTotal; ++i)
//		{
//			cout << "cluster " << i << ": ";
//			for(j = 0; j < clusterExactIndMasterArr.at(i).size(); ++j)
//			{
//				cout << clusterExactIndMasterArr.at(i).at(j) << "\t";
//			}
//			cout << endl;
//		}
//		cout  << " clusterProcArr " << endl;

		for(i = 0; i < recordExactTotal; ++i)
			recordExactIndArr.push_back((recordPacket){clusterExactIndMasterArr.at(i).at(0), i});

		tempIndArr.clear();
		procSerialLenArr.clear();
		procSerialLenArr.resize(procTotal, 0);

		totalMasterT	+= (double)(clock() - startMasterT) / CLOCKS_PER_SEC;

	}

	startBCastT	= clock();
	COMM_WORLD.Bcast(&recordExactTotal, 1, INT, ROOT_ID);
	recordExactIndArr.resize(recordExactTotal);
	// broadcast representative record index with cluster index of each exact cluster
	COMM_WORLD.Bcast(&recordExactIndArr.front(), recordExactTotal, typeRecordPacket, ROOT_ID);
	totalBCastT	+= (double)(clock() - startBCastT) / CLOCKS_PER_SEC;
	//cout << procID << " ex " << recordExactIndArr.size() << endl;

	if(procID == ROOT_ID)
	{

		procSerialLenArr.clear();
		procSerialLenArr.resize(procTotal);
	}

	startBlockT		= clock();
	// create 26^l blocks
	blockTempArr	= createBlock();
	blockTotal		= blockTempArr.size();
	count			= 0;
	blockArr.clear();
	tempDataArr.clear();
	totalBlockT		= (double)(clock() - startBlockT) / CLOCKS_PER_SEC;

	startMergeT	= clock();
	for(i = 0; i < blockTotal; ++i)
	{
		dataTemp2	= blockTempArr[i].size();
		tempDataArr.push_back(dataTemp2);
		blockArr.push_back(dataTemp2);
		for(j = 0; j < dataTemp2; ++j)
		{
			blockArr.push_back(blockTempArr[i][j]);
		}
		count		= count + 1 + dataTemp2;
	}
	totalMergeT	+= (double)(clock() - startMergeT) / CLOCKS_PER_SEC;

	COMM_WORLD.Barrier();
	startCommT	= clock();
	COMM_WORLD.Gather(&count, 1, INT, &procSerialLenArr.front(), 1, INT, ROOT_ID);
	totalCommT	+= (double)(clock() - startCommT) / CLOCKS_PER_SEC;

	if(procID == ROOT_ID)
	{
		startMergeT	= clock();

		tempDataMasterArr.resize(blockTotal);
		strideArr.clear();
		dataSum		= 0;
		for(i = 0; i < procTotal; ++i)
		{
			strideArr.push_back(dataSum);
			dataSum		+= procSerialLenArr[i];
			//cout << dataSum << " total " << procSerialLenArr.at(i) << endl;
		}
		blockMasterArr.resize(dataSum);

		totalMergeT	+= (double)(clock() - startMergeT) / CLOCKS_PER_SEC;
	}

	COMM_WORLD.Barrier();
	startCommT	= clock();
	COMM_WORLD.Gatherv(&blockArr.front(), count, INT, &blockMasterArr.front(), &procSerialLenArr.front(), &strideArr.front(), INT, ROOT_ID);
	COMM_WORLD.Reduce(&tempDataArr.front(), &tempDataMasterArr.front(), blockTotal, INT, SUM, ROOT_ID);
	totalCommT	+= (double)(clock() - startCommT) / CLOCKS_PER_SEC;

	blockArr.clear();


	if(procID == ROOT_ID)
	{
		startDistT	= clock();

		double tempDataSqArr[blockTotal];
		double countData, tempSum, dataPerProc;
		int blockIndNew;
		vector<int> blockTestArr;

		//blockIndArr	= arrangeBlockList(tempDataMasterArr);
		blockIndTempArr	= arrangeBlockListSorted(tempDataMasterArr);
		//blockIndArr		= generateRandSeq();

		tempSum	= 0;
		for(i = 0; i < blockTotal; ++i)
		{
			blockTestArr.push_back(i);
			blockIndNew					= blockIndTempArr[i];
			tempDataSqArr[blockIndNew]	= (double)(tempDataMasterArr[blockIndNew] * tempDataMasterArr[blockIndNew]);
			tempSum						+= tempDataSqArr[blockIndNew];

			//if(i <= 7 || i > (blockTotal - 1))
				//cout << i << " ind :" << blockIndNew << " data " << tempDataMasterArr[blockIndNew] << " sq " << tempDataSqArr[blockIndNew] << endl;
		}


		dataPerProc	= ceil(tempSum / procTotal);

		cout << tempSum << " perpc " << dataPerProc << endl;
		tempDataArr.clear();
		fracArr.clear();
		countData	= 0;
		ind			= 1;
		int k = 0, countT1, countT2;
		double ratioT1;
		for(i = 0; i < blockTotal; ++i)
		{
			countData	+= tempDataSqArr[blockIndTempArr[blockTestArr[0]]];

			if(countData >= (dataPerProc * ind))
			{
				countT1	= countData - dataPerProc * ind;
				ratioT1	= (double)countT1 / dataPerProc;

				cout << "#proc: " << (ind - 1) << " rem: " << countT1 << " rat: " << ratioT1 << " val: " << tempDataSqArr[blockIndTempArr[blockTestArr[0]]] << endl;
				if(ratioT1 > .01)
				{
					countData	-= tempDataSqArr[blockIndTempArr[blockTestArr[0]]];
					for(k = 0; k < blockTestArr.size();)
					{
						//if(ind > 22)
						{
							//cout << "i: " << i << " s: " << blockTestArr.size() << " " << blockTestArr[blockTestArr.size() - 1] << endl;
							//cout << blockIndTempArr[blockTestArr[blockTestArr.size() - 1]] << " t " << tempDataSqArr[blockIndTempArr[blockTestArr[blockTestArr.size() - 1]]] << " c: " << countData << endl;
						}

						countData	+= tempDataSqArr[blockIndTempArr[blockTestArr[blockTestArr.size() - 1]]];
						blockIndArr.push_back(blockIndTempArr[blockTestArr[blockTestArr.size() - 1]]);
						//if(ind > 22)
						{
							//cout << "ii: " << i << " s: " << blockTestArr.size() << " " << blockTestArr[blockTestArr.size() - 1] << endl;
							//cout << blockIndTempArr[blockTestArr[blockTestArr.size() - 1]] << endl;
						}


						if(countData >= (dataPerProc * ind))
						{
							//if(ind > 20)
							cout << "proc: " << (ind - 1) << " frac: " << (countData - dataPerProc * ind) << endl;
							//fracArr.push_back(countData - dataPerProc * ind);
							fracArr.push_back(dataPerProc * ind - countData + tempDataSqArr[blockIndTempArr[blockTestArr[blockTestArr.size() - 1]]]);
							blockTestArr.erase(blockTestArr.begin() + blockTestArr.size() - 1);
							tempDataArr.push_back(blockIndArr.size() - 1);
							break;
						}
						blockTestArr.erase(blockTestArr.begin() + blockTestArr.size() - 1);
						++i;
					}
				}
				else
				{

					//fracArr.push_back(countData - dataPerProc * ind);
					fracArr.push_back(dataPerProc * ind - countData + tempDataSqArr[blockIndTempArr[blockTestArr[0]]]);
					tempDataArr.push_back(blockIndArr.size() - 1);
					blockIndArr.push_back(blockIndTempArr[blockTestArr[0]]);
					blockTestArr.erase(blockTestArr.begin());
				}
				//fracArr.push_back(tempDataSqArr[blockIndArr[i]] - (countData - dataPerProc * ind));


				//cout << "proc: " << (ind - 1) << " frac: " << (countData - dataPerProc * ind) << " last blk: " <<  i << " count " << countData << " last " << tempDataMasterArr[blockIndArr[i]] << " lastsq " << tempDataSqArr[blockIndArr[i]] << " fr " << fracArr[ind - 1] <<  endl;
				cout << "proc: " << (ind - 1) << " frac: " << (countData - dataPerProc * ind) << " last blk: " <<  i << " count " << countData << " last " << blockTestArr.size() << " fr " << fracArr[ind - 1] << " size: " << tempDataSqArr[blockIndArr[blockIndArr.size() - 1]] <<  endl;

				ind++;

				//countData	= 0;
			}
			else
			{
				blockIndArr.push_back(blockIndTempArr[blockTestArr[0]]);
				blockTestArr.erase(blockTestArr.begin());
			}
		}

		if((dataPerProc * procTotal) != tempSum)
		{
			fracArr.push_back(tempDataSqArr[blockIndTempArr[blockTestArr[0]]]);
			tempDataArr.push_back(blockTotal - 1);
		}


		count		= blockMasterArr.size();


		cout << endl << tempDataArr.size() << " co " << blockIndArr.size() << " t: " << blockTestArr.size() << " blk: " << blockTotal << endl;

		totalDistT	+= (double)(clock() - startDistT) / CLOCKS_PER_SEC;
	}
	else
	{
		blockIndArr.resize(blockTotal);
	}

	startBCastT	= clock();
	COMM_WORLD.Bcast(&count, 1, INT, ROOT_ID);
	totalBCastT	+= (double)(clock() - startBCastT) / CLOCKS_PER_SEC;

	if(procID != ROOT_ID)
	{
		fracArr.resize(procTotal);
		tempDataArr.resize(procTotal);
		blockMasterArr.resize(count);
		tempDataMasterArr.resize(blockTotal);
		//cout << procID << " count " << count << endl;
	}

	startBCastT	= clock();
	COMM_WORLD.Bcast(&fracArr.front(), procTotal, INT, ROOT_ID);
	COMM_WORLD.Bcast(&tempDataArr.front(), procTotal, INT, ROOT_ID);
	COMM_WORLD.Bcast(&blockIndArr.front(), blockTotal, INT, ROOT_ID);
	COMM_WORLD.Bcast(&blockMasterArr.front(), count, INT, ROOT_ID);
	COMM_WORLD.Bcast(&tempDataMasterArr.front(), blockTotal, INT, ROOT_ID);
	totalBCastT	+= (double)(clock() - startBCastT) / CLOCKS_PER_SEC;

	blockTempArr.clear();
	blockTempArr	= splitIntoBlock(blockMasterArr, blockIndArr, tempDataArr, fracArr, tempDataMasterArr);

	//cout<< procID << " bt " << blockTempArr.size() << endl;

	if(procID == ROOT_ID)
	{
		//cout << "proc total " << procTotal << "\t record total " << recordArr.size() << "\t total " << totalT << "\t data " << totalDataT << "\t exact" << totalExactT << "\t sorting " << totalSortT << "\t approx " << totalApproxT << "\t block " << totalBlockT << "\t scc " << totalSCCT << endl;

	}

	startGenListT	= clock();
	createClusterEdgeList(blockTempArr);
	totalGenListT	= (double)(clock() - startGenListT) / CLOCKS_PER_SEC;

	dataTemp	= edgeArr.size();

	COMM_WORLD.Barrier();
	startCommT	= clock();
	COMM_WORLD.Gather(&dataTemp, 1, INT, &procSerialLenArr.front(), 1, INT, ROOT_ID);
	totalCommT	+= (double)(clock() - startCommT) / CLOCKS_PER_SEC;

	if(procID == ROOT_ID)
	{
		startMasterT	= clock();

		strideArr.clear();
		dataSum		= 0;
		for(i = 0; i < procTotal; ++i)
		{
			strideArr.push_back(dataSum);
			dataSum		+= procSerialLenArr.at(i);
		//	cout << dataSum << " total " << procSerialLenArr.at(i) << endl;
		}
		tempIndArr.clear();
		tempIndArr.resize(dataSum);

		totalMasterT	+= (double)(clock() - startMasterT) / CLOCKS_PER_SEC;
	}

	COMM_WORLD.Barrier();
	startCommT	= clock();
	COMM_WORLD.Gatherv(&edgeArr.front(), dataTemp, INT, &tempIndArr.front(), &procSerialLenArr.front(), &strideArr.front(), INT, ROOT_ID);
	totalCommT	+= (double)(clock() - startCommT) / CLOCKS_PER_SEC;

	COMM_WORLD.Reduce(&blockcountTest, &blockcountMTest, 1, INT, SUM, ROOT_ID);

	if(procID == ROOT_ID)
	{
		//cout << "master edgelist: " << endl;
//		for(i = 0; i < dataSum; i = i + 2)
//		{
//			cout << recordExactIndArr.at(tempIndArr.at(i)).ind << " : " << recordExactIndArr.at(tempIndArr.at(i + 1)).ind << "\t";
//		}
//		cout << endl;


		startSCCT		= clock();

		dataTemp	= clusterExactIndMasterArr.size();
		rootArr		= findConnComp(tempIndArr, dataTemp);

		totalSCCT		= (double)(clock() - startSCCT) / CLOCKS_PER_SEC;
		totalApproxT	= (double)(clock() - startApproxT) / CLOCKS_PER_SEC;
		totalT			= (double)(clock() - startT) / CLOCKS_PER_SEC;
		startMasterT	= clock();
		outputCluster(clusterExactIndMasterArr, rootArr);

		totalMasterT	+= (double)(clock() - startMasterT) / CLOCKS_PER_SEC;


		cout << "Master: " << procID << ":" << blockcountMTest << " exp: " << loopExp << " count: " << loopCount << " edit: " << editCount << " bcast: " << totalBCastT << " comm: " << totalCommT << " master: " << totalMasterT << " exact: " << totalExactT << " blk: " << totalBlockT << " merge: " << totalMergeT << " dist: " << totalDistT << " gen: " << totalGenListT << " conn: " << totalSCCT << " data: " << totalDataT << " approx: " << totalApproxT << " total: " << totalT << endl;
		//cout << "proc total " << procTotal << "\t record total " << recordArr.size() << "\t total " << totalT << "\t data " << totalDataT << "\t exact" << totalExactT << "\t sorting " << totalSortT << "\t approx " << totalApproxT << "\t block " << totalBlockT << "\t scc " << totalSCCT << endl;
	}
	else
	//cout << "proc total " << procTotal << "\t record total " << recordArr.size() << "\t total " << totalT << "\t data " << totalDataT << "\t exact" << totalExactT << "\t sorting " << totalSortT << "\t approx " << totalApproxT << "\t block " << totalBlockT << "\t scc " << totalSCCT << endl;


		cout << "procID: " << procID << " exp: " << loopExp << " count: " << loopCount << " edit: " << editCount << " exact: " << totalExactT << " blk: " << totalBlockT << " gen: " << totalGenListT << endl;

	Finalize();

	return 0;
}


vector<int> arrangeBlockListSorted(vector<int> dataArr)
{
	//cout << "arrangeBlockList" << endl;

	int i, valMax, dataTotal, exp;
	dataTotal	= dataArr.size();
	int dataTempArr[dataTotal][2], tempArr[dataTotal][2];

	for(i = 0; i < dataTotal; ++i)
	{
		dataTempArr[i][0]	= i;
		dataTempArr[i][1]	= dataArr[i];
	}

	valMax		= dataArr[0];
	for(i = 0; i < dataTotal; ++i)
		if(dataArr[i] > valMax)
			valMax	= dataArr[i];

	exp			= 1;
	while((valMax / exp) > 0)
	{
		int bucketArr[10]	= {0};

		for(i = 0; i < dataTotal; ++i)
			bucketArr[(dataTempArr[i][1] / exp) % 10]++;

		for(i = 1; i < 10; ++i)
			bucketArr[i]	+= bucketArr[i - 1];

		for(i = dataTotal - 1; i >= 0; --i)
		{
			--bucketArr[dataTempArr[i][1] / exp % 10];
			tempArr[bucketArr[(dataTempArr[i][1] / exp) % 10]][0]	= dataTempArr[i][0];
			tempArr[bucketArr[(dataTempArr[i][1] / exp) % 10]][1]	= dataTempArr[i][1];
		}

		for(i = 0; i < dataTotal; ++i)
		{
			dataTempArr[i][0]	= tempArr[i][0];
			dataTempArr[i][1]	= tempArr[i][1];
		}

		exp		*= 10;
	}

	//cout << "Sorted min " << dataTempArr[0][1] << " max " << dataTempArr[dataTotal - 1][1] << endl;

	vector<int> dataShuffleArr;
	for(i = 0; i < dataTotal; ++i)
	{
		dataShuffleArr.push_back(dataTempArr[dataTotal - i - 1][0]);

	}
	//cout << "data return " << dataShuffleArr.size() << endl;
	return dataShuffleArr;

}


vector<int> arrangeBlockListRand(int dataTotal)
{
	int i;
	vector<int> dataShuffleArr;

	for(i = 0; i < dataTotal; ++i)
		dataShuffleArr.push_back(i);

	random_shuffle( dataShuffleArr.begin(), dataShuffleArr.end() );

	return dataShuffleArr;

}

// sort blocks using their content total and rearrange them such that smaller and bigger one stay together
vector<int> arrangeBlockList(vector<int> dataArr)
{
	//cout << "arrangeBlockList" << endl;

	int i, valMax, dataTotal, exp;
	dataTotal	= dataArr.size();
	int dataTempArr[dataTotal][2], tempArr[dataTotal][2];

	for(i = 0; i < dataTotal; ++i)
	{
		dataTempArr[i][0]	= i;
		dataTempArr[i][1]	= dataArr[i];
	}

	valMax		= dataArr[0];
	for(i = 0; i < dataTotal; ++i)
		if(dataArr[i] > valMax)
			valMax	= dataArr[i];

	exp			= 1;
	while((valMax / exp) > 0)
	{
		int bucketArr[10]	= {0};

		for(i = 0; i < dataTotal; ++i)
			bucketArr[(dataTempArr[i][1] / exp) % 10]++;

		for(i = 1; i < 10; ++i)
			bucketArr[i]	+= bucketArr[i - 1];

		for(i = dataTotal - 1; i >= 0; --i)
		{
			--bucketArr[dataTempArr[i][1] / exp % 10];
			tempArr[bucketArr[(dataTempArr[i][1] / exp) % 10]][0]	= dataTempArr[i][0];
			tempArr[bucketArr[(dataTempArr[i][1] / exp) % 10]][1]	= dataTempArr[i][1];
		}

		for(i = 0; i < dataTotal; ++i)
		{
			dataTempArr[i][0]	= tempArr[i][0];
			dataTempArr[i][1]	= tempArr[i][1];
		}

		exp		*= 10;
	}

	//cout << "Sorted min " << dataTempArr[0][1] << " max " << dataTempArr[dataTotal - 1][1] << endl;

	vector<int> dataShuffleArr;
	int loopTotal = floor((double)dataTotal / 2.0);
	for(i = 0; i < loopTotal; ++i)
	{
		dataShuffleArr.push_back(dataTempArr[i][0]);
		dataShuffleArr.push_back(dataTempArr[dataTotal - i - 1][0]);
	}

	if((dataTotal % 2) == 1)
		dataShuffleArr.push_back(dataTempArr[i][0]);

	//cout << "data return " << dataShuffleArr.size() << endl;
	return dataShuffleArr;

}

vector<int> generateRandSeq()
{
	vector<int> blockIndArr;
	int i, randInd, temp, blockTotal;

	blockTotal	= pow(26.0, L_MER);
	blockIndArr.resize(blockTotal);

	for(i = 0; i < blockTotal; ++i)
		blockIndArr[i]	= i;

	srand(clock());
	for(i = 0; i < (blockTotal - 1); ++i)
	{
		randInd	= rand() % (blockTotal - i);
		temp	= blockIndArr[i];
		blockIndArr[i]				= blockIndArr[i + randInd];
		blockIndArr[i + randInd]	= temp;
	}

//	for(i = 0; i < blockTotal; ++i)
//		cout << blockIndArr[i] << "\t";
//	cout << endl;

	return blockIndArr;
}


void outputCluster(vector<vector<int> > clusterExactIndMasterArr, vector<int> rootArr)
{
	int i, j, clusterExactTotal, clusterInd;
	vector<int> clusterRowArr;
	vector<vector<int> > clusterArr;

	clusterExactTotal	= clusterExactIndMasterArr.size();
	clusterArr.resize(clusterExactTotal);

	for(i = 0; i < clusterExactTotal; ++i)
	{
		clusterInd		= rootArr.at(i);

		for(j = 0; j < clusterExactIndMasterArr.at(i).size(); ++j)
		{
			clusterArr.at(clusterInd).push_back(clusterExactIndMasterArr.at(i).at(j));
		}
	}

//	for(i = 0; i < clusterExactTotal; ++i)
//	{
//		cout << "cluster " << i << ": ";
//		for(j = 0; j < clusterArr.at(i).size(); ++j)
//		{
//			cout << clusterArr.at(i).at(j) << "\t";
//		}
//		cout << endl;
//	}
//	cout  << " done " << endl;

	clusterInd	= 0;
//	ofstream outFile;
//
//	outFile.open("output.txt", ofstream::out);

	for(i = 0; i < clusterExactTotal; ++i)
	{
		if(clusterArr.at(i).size() > 0)
		{
//			outFile << "cluster " << clusterInd << ": ";
//			for(j = 0; j < clusterArr.at(i).size(); ++j)
//			{
//				outFile << clusterArr.at(i).at(j) << "\t";
//			}
//			outFile << endl;
			++clusterInd;
		}
	}
//	outFile  << endl << "Total Cluster: " << clusterInd << endl;
//
//	outFile.close();

	cout  << endl << "Total Cluster: " << clusterInd << endl;
}


vector<int> findConnComp(vector<int> indArr, int pointTotal)
{
	int i, j, rootU, rootV, edgeTotal;
	vector<int> parentArr, weightArr;

	for(i = 0; i < pointTotal; ++i)
	{
		parentArr.push_back(i);
		weightArr.push_back(1);
	}

	edgeTotal	= indArr.size();

	for(i = 0; i < edgeTotal; i = i + 2)
	{
		rootU	= findRoot(indArr.at(i), parentArr);
		rootV	= findRoot(indArr.at(i + 1), parentArr);

		if(rootU != rootV)
		{
			//cout << i << " e " << rootU << ":" << rootV << endl;
			makeEquivalent(rootU, rootV, parentArr, weightArr);
		}
	}

//	for(i = 0; i < pointTotal; ++i)
//	{
//		cout << i << ":" << parentArr.at(i) << "\t";
//	}
//	cout << endl;

	return parentArr;
}

// collapse find operation
int findRoot(int pointID, vector<int> &parentArr)
{
	int i, j, root;

	i	= pointID;
	while(parentArr[i] != i)
		i	= parentArr[i];
	root	= i;

	i	= pointID;
	while(parentArr[i] != i)
	{
		j				= parentArr[i];
		parentArr[i]	= root;
		i				= j;
	}

	return root;
}

// unio two trees rooted at rootU and rootV
void makeEquivalent(int rootU, int rootV, vector<int> &parentArr, vector<int> &weightArr)
{
	if(weightArr[rootU] > weightArr[rootV])
	{
		parentArr[rootV]	= rootU;
		weightArr[rootU]	+= weightArr[rootV];
	}
	else
	{
		parentArr[rootU]	= rootV;
		weightArr[rootV]	+= weightArr[rootU];
	}
}


/*
 * start
 * find exact clusters
 */


// convert string to int value for partitioning of records among processors
int convertStrToInd(string str)
{
	if(str.length() < 2)
		return 0;

	return (((int)str.at(0) - 97) + 26 * ((int)str.at(1) - 97));
}


// change to lowercase
string convertToLower(string str)
{
	for(int i = 0; i < str.length(); ++i)
		str[i]	= tolower(str[i]);
	return str;
}


// sort data using radix sort to find duplicates
void sortData()
{
	//if(procID == ROOT_ID)


	int i, j, k, recordLimit, recordPerLimit, recordInit, recordLast, recordTotal, lenMax, lenDiff, indVal, recordTotalThis, lastNameInd;
	vector<StrPacket> tempArr;
	vector<int> countArr;
	string strTemp, strSample;

	lastNameInd		= getLNInd();
	strSample		= "0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000";
	lenMax			= 0;
	recordTotal		= recordArr.size();
	recordLimit		= pow(26, 2.0);
	recordPerLimit	= ceil((double)recordLimit / procTotal);
	recordInit		= procID * recordPerLimit;
	recordLast		= (recordInit + recordPerLimit) >= recordLimit? (recordLimit - 1) : (recordInit + recordPerLimit - 1);


	for(i = 0; i < recordTotal; ++i)
	{
		//strTemp		= recordArr[i][lastNameInd] + "#" + recordArr[i][1] + "#" + recordArr[i][2] + "#" + recordArr[i][3] + "#" + recordArr[i][4];
		strTemp		= "";
		for(j = 0; j < attrEditArr.size(); ++j)
			strTemp	+= recordArr[i][attrEditArr.at(j).at(0)];


		indVal		= convertStrToInd(strTemp);
		if((indVal < recordInit) || (indVal > recordLast))
			continue;

		recordStrArr.push_back((StrPacket){strTemp, i});

		if(strTemp.length() > lenMax)
			lenMax	= strTemp.length();
	}

	recordTotalThis	= recordStrArr.size();

	//cout << procID << " fdsortData" << endl;

	for(i = 0; i < recordTotalThis; ++i)
	{
		lenDiff		= lenMax - recordStrArr[i].str.length();
		if(lenDiff > 0)
			recordStrArr[i].str	= recordStrArr[i].str + strSample.substr(0, lenDiff);
	}

	//cout << procID << ":" << recordArr.size() << " this " << recordTotalThis << ":" << lenMax << endl;
	if(recordTotalThis > 0)
	{
		tempArr.resize(recordTotalThis);
		countArr.resize(ALPHABET_SIZE, 0);

		for(i = lenMax - 1; i >= 0; --i)
		{

			for(j = 0; j < recordTotalThis; ++j)
			{
				//cout << i << ":" << j << "::" << recordStrArr[j].str.at(i) << " :v: " << (recordStrArr[j].str.at(i) + 1) << endl;
				countArr[recordStrArr[j].str.at(i) + 1]++;
			}

			for(k = 1; k < ALPHABET_SIZE; ++k)
			{
				countArr[k]			+= countArr[k - 1];
			}

			for(j = 0; j < recordTotalThis; ++j)
			{
				//cout << i << ":" << j << ":tesdfsdfmp:" << recordStrArr[j].str << ":" << countArr[recordStrArr[j].str.at(i)] << " :ds" << endl;
				tempArr[countArr[recordStrArr[j].str.at(i)]++]	= recordStrArr[j];
			}

			for(j = 0; j < recordTotalThis; ++j)
			{
				recordStrArr[j]		= tempArr[j];
				//cout << i << ":" << j << endl;
			}
			countArr.assign(ALPHABET_SIZE, 0);
		}
		//cout << procID <<  " sortData" << endl;
//		if(procID == ROOT_ID)
//		{
//			for(i = 0; i < recordTotalThis; ++i)
//				cout << recordStrArr[i].str << " index " << recordStrArr[i].ind << endl;
//		}


	}
}

// find exact clusters using sorted records
void findExactCluster(vector<int> &clusterExactIndArr)
{
	if(recordStrArr.size() <= 0)
		return;

	//if(procID == ROOT_ID)
		//cout << "findExactCluster" << endl;

	int i, j, recordTotalThis;

	vector<int> clusterSingleArr;

	clusterExactTotalThis	= 1;
	recordTotalThis			= recordStrArr.size();

	clusterSingleArr.push_back(recordStrArr[0].ind);

	for(i = 1; i < recordTotalThis; ++i)
	{
		if(recordStrArr[i].str.compare(recordStrArr[i - 1].str) == 0)
			clusterSingleArr.push_back(recordStrArr[i].ind);
		else
		{
			++clusterExactTotalThis;
			clusterExactIndArr.push_back(clusterSingleArr.size());
			for(j = 0; j < clusterSingleArr.size(); ++j)
				clusterExactIndArr.push_back(clusterSingleArr.at(j));
			clusterSingleArr.clear();
			clusterSingleArr.push_back(recordStrArr[i].ind);
		}
	}
	clusterExactIndArr.push_back(clusterSingleArr.size());
	for(i = 0; i < clusterSingleArr.size(); ++i)
		clusterExactIndArr.push_back(clusterSingleArr.at(i));

	//cout << clusterExactTotalThis<<  " findExactCluster: " << clusterExactIndArr.size() << endl;

	recordStrArr.clear();
	//if(procID == ROOT_ID)
//	{
//		cout << "proc id: " << procID << endl;
//		for(i = 0; i < clusterExactTotalThis; ++i)
//		{
//			cout << clusterExactIndArr.at(i) << " ";
//			//cout << endl;
//		}
//	}
}

/*
 * end
 * find exact clusters
 */

/*
 * start
 * create edgelist within blocks
 */



// split blocks into processors
vector<vector<int> > splitIntoBlock(vector<int> blockMasterArr, vector<int> blockIndArr, vector<int> indArr, vector<int> fracArr, vector<int> blockCountArr)
{
	//cout << procID << " splitIntoBlock" << endl;
	vector<vector<int> > blockArr, blockStartIndArr;
	int i, j, k, t, dataSize, blockTotal, blockTotalThis, blockStart, blockEnd, blockIndThis, val, blockInd, count, ind;

	blockTotal	= pow(26.0, L_MER);
	dataSize	= blockMasterArr.size();
	if(procID == 0)
	{
		blockTotalThis	= indArr[procID] + 1;
		blockStart		= 0;
		fracStart		= 0;
		fracLast		= fracArr[0];
	}
	else
	{
		blockTotalThis	= indArr[procID] - indArr[procID - 1] + 1;
		blockStart		= indArr[procID - 1];
		//fracStart		= (double)blockCountArr[blockIndArr[blockStart]] * blockCountArr[blockIndArr[blockStart]] - (double)fracArr[procID - 1];
		fracStart		= fracArr[procID - 1];
		fracLast		= fracArr[procID];
	}
	blockEnd	= indArr[procID];

	blockcountTest	= blockTotalThis;
//	blockTotalThis	= ceil((double)blockTotal / procTotal);
//	if(procID == (procTotal - 1))
//		blockTotalThis	= blockTotal - procID * blockTotalThis;
//
//	blockStart		= procID * blockTotalThis;
//	blockEnd		= blockStart + blockTotalThis - 1;

	blockStartIndArr.resize(blockTotal);
	blockArr.resize(blockTotalThis);
	count		= 0;
	ind			= 0;

	if(procID == 5)
		cout << procID << " splitIntoBlock: " << blockTotalThis << ":" << blockStart << ":" << fracStart << " :l: " << fracLast << " s:" << (blockCountArr[blockIndArr[blockStart]] * blockCountArr[blockIndArr[blockStart]]) << endl;

	//cout << procID << " bs " << blockStart << " be " << blockEnd << " bt " << blockTotalThis << " fs " << fracStart << " fe " << fracLast << " data " << blockCountArr[blockIndArr[blockStart]] << endl;

	for(i = 0; i < dataSize; )
	{
		blockStartIndArr[ind].push_back(i);

		val	= blockMasterArr[i];
		i	= i + val + 1;

		++ind;
		if(ind >= blockTotal)
			ind	= 0;
	}

	//cout << procID << " t " << blockStartIndArr[17574].size() << " u " << blockStartIndArr[17575].size() << endl;

	for(i = blockStart; i <= blockEnd; ++i)
	{
		blockIndThis	= blockIndArr[i];
		blockInd		= i - blockStart;

		for(j = 0; j < blockStartIndArr[blockIndThis].size(); ++j)
		{
			ind			= blockStartIndArr[blockIndThis][j];
			for(k = 1; k <= blockMasterArr[ind]; ++k)
			{
				++count;
				blockArr[blockInd].push_back(blockMasterArr[ind + k]);
			}
		}
	}

	//cout << procID << " bass " << blockArr.size() << " count: " << count << " t " << blockStartIndArr[17574].size() << " u " << blockStartIndArr[17575].size() << endl;

	blockStartIndArr.clear();
	//cout << procID << " finish splitIntoBlock" << endl;
	return blockArr;
}



// create blocks using l-mer
vector<vector<int> > createBlock()
{
	//if(procID == ROOT_ID)
		//cout << "createBlock" << endl;

	int i, j, k, blockTotal, recordTotalThis, recordStart, recordEnd, strLen, strLenDiff, blockID, lastNameInd, recordExactIndTotal;
	vector<vector<int> > blockArr;
	vector<int> tempArr, codePointArr;
	string lastNameStr, strSample;


	recordExactIndTotal	= recordExactIndArr.size();

	strSample		= "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa";
	lastNameInd		= getLNInd();

	blockTotal		= pow(26.0, L_MER);
	recordTotalThis	= ceil((double)recordExactIndTotal / procTotal);
	recordStart		= recordTotalThis * procID;
	recordEnd		= (recordStart + recordTotalThis) >= recordExactIndTotal? (recordExactIndTotal - 1) : (recordStart + recordTotalThis - 1);
	recordTotalThis	= recordEnd - recordStart + 1;

	blockArr.resize(blockTotal);

	//cout << procID << " block " << blockStart << ":" << blockEnd << " tot " << recordExactIndTotal << endl;

	for(i = recordStart; i <= recordEnd; ++i)
	{
		lastNameStr	= recordArr[recordExactIndArr[i].ind][lastNameInd];
		//if(procID == 1)
			//		cout << i << " for " << clusterExactIndArr.at(i).at(0) << "--" << lastNameStr << " lsadast " << endl;
		strLen		= lastNameStr.length();

		if(strLen < L_MER)
			lastNameStr	= strSample.substr(0, L_MER - strLen) + lastNameStr;
		strLen		= lastNameStr.length();
		codePointArr.clear();
		codePointArr.resize(strLen);

		for(j = 0; j < lastNameStr.length(); ++j)
		{
			if( ((int)lastNameStr.at(j) < 97) || ((int)lastNameStr.at(j) > 122) )
			{
				lastNameStr		= lastNameStr.substr(0, j) + lastNameStr.substr(j + 1);
				--j;
			}
			else
				codePointArr[j]	= (int)lastNameStr[j];
		}

		strLenDiff	= strLen - L_MER + 1;
		for(j = 0; j < strLenDiff; ++j)
		{
			blockID	= 0;
			for(k = 0; k < L_MER; ++k)
				blockID		+= (codePointArr[j + k] - 97) * (int)pow(26.0, L_MER - k - 1);
			if((blockID >= 0) && (blockID < blockTotal))
			{
				//++checkTemp;
				blockArr[blockID].push_back(i);
			}
		}
	}

	//if(procID == ROOT_ID)
	//cout << endl << procID << " easdx " << blockTotal << " asd " << blockArr.size() << endl;

	codePointArr.clear();
	tempArr.clear();
	return blockArr;
}

// create edgelist within all blocks
void createClusterEdgeList(vector<vector<int> > blockArr)
{
	int i, j, blockTotal, blockItemTotal;
	clock_t startT;
	double edgeListT;

	//cout << procID << "createClusterEdgeList " << checkTemp << endl;
	checkTemp	= 0;
	checkTemp2	= 0;
	checkTemp3	= 0.0;
	checkTemp4	= 0.0;

	startT		= clock();
	//if(procID == 1)
	//cout << procID << " blockTotal " << blockTotal << endl;
	blockTotal	= blockArr.size();

	for (i = 0; i < blockTotal; i++)
	{
		if (blockArr.at(i).size() != 0)
			generateEdgilist(blockArr.at(i), i, blockTotal);
	}
	edgeListT	= (double)(clock() - startT) / CLOCKS_PER_SEC;
	//cout << procID <<  " edgeArr " << edgeArr.size() << " listT " << edgeListT << " sqdata " << checkTemp << " sqIdeal " << checkTemp4 << " blockdata " << checkTemp2 << " avg iterate " << (checkTemp3 / blockTotal) << " blocktotal " << blockTotal << " fst " << fracStart << " fls " << fracLast <<  endl;

	//dataArr.clear();
	blockArr.clear();
}



// generate edgelist within a block
void generateEdgilist(vector<int> blockRowArr, int blockInd, int blockTotal)
{
	//if(procID == ROOT_ID)
		//cout << "generateVector " << procID << endl;

	int i, j, n, temp, blockItemTotal, count, indStart, indLast, indInStart, indInLast;
	vector<int> vectorArr, tempArr;
	vector<vector<string> > dataArr;

	count			= 0;
	n				= 0;
	blockItemTotal	= blockRowArr.size();
	indInStart		= 0;
	indStart		= 0;
	indInLast		= blockItemTotal - 1;
	indLast			= blockItemTotal - 1;

	if(blockInd == 0)
	{
		indStart	= floor((double)fracStart / blockItemTotal);
		//indInStart	= indStart;
		//cout << procID << " indstart " << indStart << " frl " << fracStart << " item " << blockItemTotal << endl;
	}
	else if(blockInd == blockTotal - 1)
	{
		indLast		= ceil((double)fracLast / blockItemTotal);
		//indInLast	= indLast;
		//cout << procID << " indlast " << indLast << " frl " << fracLast << " item " << blockItemTotal << endl;
	}


	vectorArr.resize(blockItemTotal, 0);
	for(i = 0; i < blockItemTotal; ++i)
		dataArr.push_back(recordArr[recordExactIndArr[blockRowArr[i]].ind]);

	for(i = indStart; i <= indLast; ++i)
	{
		if(vectorArr[i] == 0)
		{
			++n;
			tempArr.push_back(i);
			vectorArr[i]	= n;

			while(tempArr.size() > 0)
			{
				temp	= tempArr.at(0);
				tempArr.erase(tempArr.begin());
				for(j = indInStart; j <= indInLast; ++j)
				{
					if(vectorArr[j] == 0)
					{
						++count;
						++checkTemp;
//						if(procID == 1 && checkTemp == 3164)
//						{
//							cout << "generateVector " << procID << " thresh " << n << endl;
//							cout << clusterExactIndArr.at(blockRowArr.at(temp)).at(0) << "---" << clusterExactIndArr.at(blockRowArr.at(j)).at(0) << endl;
//						}

						//if(calculateDistTotal(recordArr.at(recordExactIndArr.at(blockRowArr.at(temp)).ind), recordArr.at(recordExactIndArr.at(blockRowArr.at(j)).ind)) <= threshold)
						//if(calculateDistTotal(recordArr[recordExactIndArr[blockRowArr[temp]].ind], recordArr[recordExactIndArr[blockRowArr[j]].ind]) <= threshold)
						if(calculateDistTotal(dataArr[temp], dataArr[j]) <= threshold)
						{
							editCount += 1;
//							if(procID == 1 && checkTemp == 3164)
//								cout << j << ":" << clusterExactIndArr.at(blockRowArr.at(temp)).at(0) << "---" << clusterExactIndArr.at(blockRowArr.at(j)).at(0) << endl;
							//cout << procID << " match " << recordExactIndArr.at(blockRowArr.at(temp)).ind << ":" << recordExactIndArr.at(blockRowArr.at(j)).ind << "\t" << recordArr.at(recordExactIndArr.at(blockRowArr.at(temp)).ind).at(2) << ":" << recordArr.at(recordExactIndArr.at(blockRowArr.at(j)).ind).at(2) << endl;
							tempArr.push_back(j);
							vectorArr[j]	= n;
							edgeArr.push_back(blockRowArr[temp]);
							edgeArr.push_back(blockRowArr[j]);
						}
						//if(procID == 1 && checkTemp == 3164)
						//cout << "generateVector " << procID << " after " << n << endl;
					}
				}
			}
		}
	}

	loopExp		+= blockItemTotal * (indLast - indInStart);

	loopCount	+= count;
	checkTemp2	+= blockItemTotal;
	checkTemp3	+= ((double)count / (blockItemTotal * blockItemTotal));
	checkTemp4	+= (blockItemTotal * blockItemTotal);
	//if(procID == ROOT_ID)
	//cout << "final " << procID << " thresh " << n << endl;

	dataArr.clear();
	vectorArr.clear();
	tempArr.clear();
}



int calculateDistTotal(vector<string> a, vector<string> b)
{
	int w;

	w	= calculateEditDist(a, b, attrEditArr);
//	if(procID == ROOT_ID)
//		cout << w << "--" << a.at(2) << "--" << b.at(2) << endl;
	return w;
}


int calculateEditDist(vector<string> a, vector<string> b, vector<vector<int> > attrEditArr)
{
	int w, i, aInd, bInd, a_f, b_f, temp, threshRem;
	vector<int> attrArr;
	string str1, str2;

	w			= 0;
	a_f			= 0;
	b_f			= 0;
	threshRem	= threshold;

	for(i = 0; i < attrEditArr.size(); ++i)
	{
		//++checkTemp;
		attrArr		= attrEditArr.at(i);
		aInd		= attrArr[a_f];
		bInd		= attrArr[b_f];
		str1		= a.at(aInd);
		str2		= b.at(bInd);
		temp		= calculateBasicED(str1, str2, threshRem);
		w			+= temp;
		threshRem	= threshRem - temp;
	}
	return w;
}


int calculateBasicED(string str1, string str2, int threshRem)
{
	int dist;

	dist	= threshRem;

	//return (threshold + 1);
	//if(procID == 1 && checkTemp == 3164)
		//			cout << str1 << " -- " << str2 << " dist " << dist << endl;
	if(abs((int)(str1.length() - str2.length())) > dist)
		return (threshold + 1);
	else if(str1.compare(str2) == 0)
		return 0;
	else if(dist == 0)
		return (threshold + 1);
	else if((2 * dist + 1) >= max(str1.length(), str2.length()))
		return calculateBasicED2(str1, str2, dist);
	else
	{
		string s1, s2;
		int row, col, diagonal, i, j;
		vector<vector<int> > matArr;

		if (str1.length() > str2.length())
		{
			s1 = str2;
			s2 = str1;
		}
		else
		{
			s1 = str1;
			s2 = str2;
		}

		row	 		= s1.length() + 1;
		col 		= 2 * dist + 1;
		diagonal 	= dist + s2.length() - s1.length();

		matArr.resize(row);
		for(i = 0; i < row; ++i)
			matArr.at(i).resize(col, 0);

		//if(procID == 1 && checkTemp == 3164)
			//	cout << str1 << " -- " << str2 << " rt " << dist << endl;

		for(i = 0; i < dist + 1; i++)
		{
			for(j = dist - i; j < col; j++)
			{
				if (i == 0)
					matArr[i][j]	= j - dist;
				else if(j == (dist - i))
					matArr[i][j] 	= matArr[i - 1][j + 1] + 1;
				else if(j != (col - 1))
				{
					if((int)s1.at(i - 1) == (int)s2.at(j - (dist - i) - 1))
						matArr[i][j]	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1.at(i - 1) == (int)s2.at(j - (dist - i) - 1))
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}

				if((j == diagonal) && matArr[i][j] > dist)
					return (threshold + 1);
			}
		}

		for(i = dist + 1; i < s2.length() - dist + 1; i++)
		{
			for(j = 0; j < col; j++)
			{
				if(j == 0)
				{
					if((int)s1.at(i - 1) == (int)s2.at(j + (i - dist) - 1))
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else if(j != (col - 1))
				{
					if((int)s1.at(i - 1) == (int)s2.at(j + (i - dist) - 1))
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				else
				{
					if((int)s1.at(i - 1) == (int)s2.at(j + (i - dist) - 1))
						matArr[i][j] 	= min(matArr[i - 1][j], matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return (threshold + 1);
			}
		}

		for(i = s2.length() - dist + 1; i < row; i++)
		{
			for(j = 0; j < col - i + s2.length() - dist; j++)
			{
				if(j == 0)
				{
					if((int)s1.at(i - 1) == (int)s2.at(j + (i - dist) - 1))
						matArr[i][j]	= min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1);
					else
						matArr[i][j] 	= min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1);
				}
				else
				{
					if((int)s1.at(i - 1) == (int)s2.at(j + (i - dist) - 1))
						matArr[i][j] 	= min(min(matArr[i - 1][j], matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
					else
						matArr[i][j] 	= min(min(matArr[i - 1][j] + 1, matArr[i - 1][j + 1] + 1), matArr[i][j - 1] + 1);
				}
				if((j == diagonal) && (matArr[i][j] > dist))
					return (threshold + 1);
			}
		}

		//if(procID == 1 && checkTemp == 3164)
			//cout << str1 << " -- " << str2 << " hukjhk " << matArr[row - 1][diagonal] << endl;

		return matArr[row - 1][diagonal];

	}
}


int calculateBasicED2(string str1, string str2, int threshRem)
{
	int row, col, i, j;
	vector<vector<int> > matArr;

	row		= str1.length() + 1;
	col 	= str2.length() + 1;

	matArr.resize(row);
	for(i = 0; i < row; ++i)
		matArr.at(i).resize(col, 0);

	for(i = 0; i < row; i++)
	{
		for(j = 0; j < col; j++)
		{
			if (i == 0)
				matArr[i][j] = j;
			else if (j == 0)
				matArr[i][j] = i;
			else
			{
				if((int)str1.at(i-1) == (int)str2.at(j-1))
					matArr[i][j]	= min(min(matArr[i - 1][j - 1], matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
				else
					matArr[i][j] 	= min(min(matArr[i - 1][j - 1] + 1, matArr[i - 1][j] + 1), matArr[i][j - 1] + 1);
			}

			if((row - col) == (i - j) && (matArr[i][j] > threshRem))
			{
				return (threshold + 1);
			}
		}
	}

	return (matArr[row-1][col-1]);
}



/*
 * end
 * create edgelist within blocks
 */


/*
 * start
 * parameter read from external xml file
 * data read from text files
 */

void readDataFromFile(vector<string> fileNameArr)
{
	int i, j, pointPerFile, count;
	ifstream inFile;
	string lineStr;
	vector<string> rowStrArr;

	regex re("\\s+");

	count		= 0;
	if(FILE_SAME)
		pointPerFile	= RECORD_TOTAL;
	else
		pointPerFile	= (int)ceil((double)RECORD_TOTAL / fileNameArr.size());

	for(i = 0; i < fileNameArr.size(); ++i)
	{
		if(count >= RECORD_TOTAL) break;

		if(!FILE_SAME)
			count	= 0;

		inFile.open(fileNameArr.at(i).c_str(), ifstream::in);
		//if(procID == ROOT_ID)
			//cout << fileNameArr.at(i).c_str() << endl;
		//cout << fileNameArr[i].c_str() << endl;
		if(inFile.good())
		{
			while(getline(inFile, lineStr))
			{
				if(lineStr.length() > 10)
				{
					//if(procID == ROOT_ID)
						//cout << lineStr << " ::: " << recordArr.size() << endl;

					rowStrArr.clear();

					for(sregex_token_iterator it = sregex_token_iterator(lineStr.begin(), lineStr.end(), re, -1); it != sregex_token_iterator(); ++it)
						rowStrArr.push_back(convertToLower((string)*it));

					recordArr.push_back(rowStrArr);

					//cout << rowStrArr.at(2) << ":" << rowStrArr.at(3) << " len " << rowStrArr.size() <<endl;
					//if(recordArr.size() > 400000) break;
					++count;
					if(count >= pointPerFile) break;
				}
			}
		}

		inFile.close();
	}

	if(procID == ROOT_ID)
		cout << pointPerFile << ":" << RECORD_TOTAL << endl;
}

void getInputComparisonPara(vector<vector<int> > &attrEditArr)
{
	xmlDoc *doc;
	xmlNode *nodeCurr, *nodeChild, *nodeRoot, *nodeGChild;
	xmlChar *attrVal;
	xmlXPathContext *xpathCtx;
	xmlXPathObject *xpathObj;
	int i;
	vector<int> attrArr;

	doc	= xmlReadFile(FILE_CONFIG, NULL, 0);

	nodeRoot	= xmlDocGetRootElement(doc);
	xpathCtx	= xmlXPathNewContext(doc);
	xpathObj	= xmlXPathEvalExpression((xmlChar *)"/febrl-config/version-config-param/comparison", xpathCtx);
	xmlXPathFreeContext(xpathCtx);

	for(i = 0; i < xpathObj->nodesetval->nodeNr; ++i)
	{
		//cout << i << endl;
		nodeCurr	= xpathObj->nodesetval->nodeTab[i];
		nodeChild	= nodeCurr->children;
		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"dist_calc_method"))
			{
				nodeGChild	= nodeChild->children;
				while(nodeGChild)
				{
					if(!xmlStrcmp(nodeGChild->name, (xmlChar *)"value"))
					{
						attrVal		= xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1);
						//cout << attrVal << "attr" << endl;
						break;
					}
					nodeGChild	= nodeGChild->next;
				}
			}

			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"comparing_attribute_indices"))
			{
				attrArr.clear();

				nodeGChild	= nodeChild->children;
				while(nodeGChild)
				{
					if(!xmlStrcmp(nodeGChild->name, (xmlChar *)"value"))
					{

						if(atoi((char *)attrVal) == 1)
						{
							attrArr.push_back(atoi((char *)xmlNodeListGetString(nodeCurr->doc, nodeGChild->children, 1)));
							attrEditArr.push_back(attrArr);
						}
						break;
					}
					nodeGChild	= nodeGChild->next;
				}
			}
			nodeChild	= nodeChild->next;
		}

	}

	xmlFreeDoc(doc);
	xmlCleanupParser();
}

void getInputFileNameList(vector<string> &fileNameArr)
{
	xmlDoc *doc;
	xmlNode *nodeCurr, *nodeChild, *nodeRoot;
	xmlChar *attrVal;
	xmlXPathContext *xpathCtx;
	xmlXPathObject *xpathObj;
	int i;

	doc	= xmlReadFile(FILE_CONFIG, NULL, 0);

	nodeRoot	= xmlDocGetRootElement(doc);
	xpathCtx	= xmlXPathNewContext(doc);
	xpathObj	= xmlXPathEvalExpression((xmlChar *)"/febrl-config/dataset", xpathCtx);
	xmlXPathFreeContext(xpathCtx);

	for(i = 0; i < xpathObj->nodesetval->nodeNr; ++i)
	{
		nodeCurr	= xpathObj->nodesetval->nodeTab[i];
		nodeChild	= nodeCurr->children;
		while(nodeChild)
		{
			if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
			{
				fileNameArr.push_back((char *)xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1));
				break;
			}
			nodeChild	= nodeChild->next;
		}

	}

	xmlFreeDoc(doc);
	xmlCleanupParser();
}

int getInputThreshold()
{
	xmlDoc *doc;
	xmlNode *nodeCurr, *nodeChild, *nodeRoot;
	xmlChar *attrVal;
	xmlXPathContext *xpathCtx;
	xmlXPathObject *xpathObj;

	doc	= xmlReadFile(FILE_CONFIG, NULL, 0);

	nodeRoot	= xmlDocGetRootElement(doc);
	xpathCtx	= xmlXPathNewContext(doc);
	xpathObj	= xmlXPathEvalExpression((xmlChar *)"/febrl-config/version-config-param/threshold", xpathCtx);
	xmlXPathFreeContext(xpathCtx);

	nodeCurr	= xpathObj->nodesetval->nodeTab[0];
	nodeChild	= nodeCurr->children;
	while(nodeChild)
	{
		if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
		{
			attrVal	= xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
			break;
		}
		nodeChild	= nodeChild->next;
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return atoi((char *)attrVal);
}

int getLNInd()
{
	xmlDoc *doc;
	xmlNode *nodeCurr, *nodeChild, *nodeRoot;
	xmlChar *attrVal;
	xmlXPathContext *xpathCtx;
	xmlXPathObject *xpathObj;

	doc	= xmlReadFile(FILE_CONFIG, NULL, 0);

	nodeRoot	= xmlDocGetRootElement(doc);
	xpathCtx	= xmlXPathNewContext(doc);
	xpathObj	= xmlXPathEvalExpression((xmlChar *)"/febrl-config/version-config-param/last_name_index", xpathCtx);
	xmlXPathFreeContext(xpathCtx);

	nodeCurr	= xpathObj->nodesetval->nodeTab[0];
	nodeChild	= nodeCurr->children;
	while(nodeChild)
	{
		if(!xmlStrcmp(nodeChild->name, (xmlChar *)"value"))
		{
			attrVal	= xmlNodeListGetString(nodeCurr->doc, nodeChild->children, 1);
			break;
		}
		nodeChild	= nodeChild->next;
	}

	xmlFreeDoc(doc);
	xmlCleanupParser();

	return atoi((char *)attrVal);
}


/*
 * end
 * parameter read from external xml file
 * data read from text files
 */
