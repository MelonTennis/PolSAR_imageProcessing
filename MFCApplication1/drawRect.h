#pragma once

// declare headfile 
	#include <tchar.h>
	#include <assert.h>
	#include "highgui.h"
	#include "cv.h"
	#include <iostream>
	#include <malloc.h>  
	#include <string.h>	
	using namespace cv; 


class drawRect
{
public:
	drawRect(void);
	~drawRect(void);

// declare function 

// function used 
	// main function

	void main_Rect(CString cutpath, CString suspected, CString sample);
	//void main_Rect(void);
	int *findTarget(char *contourPath, char *spanImg, char **countSave, char **traitpath, char **traitD, char **rectCate, int threContour, int traitNum, int traitNumD, int catenum);
	// save targets into bmps
	void drawCategoryRect(char *rawImage,  char *spanImg, char *rectnum, char *carMeanPath, char *vehiclesMean, int meamNum, int *category, int thresholdContour, char **savepath,  int catenum, char **addHead, int Rwidth, int Rheight);
	// given atemplateImg, finding it in image
	void fitTemplate(char * image, char *templateImg);
	// draw inner cars in a contour of cluster of vehicles
	IplImage *drawCars(char *spanImg, IplImage *cut, IplImage *innerCut, CvRect rectB, char *meanCar, char *vehiclesMean, int meanNum, int thresholdContour, char *addHead, int catenum, int Rwidth, int Rheight);
	// given a part of vehicle, return the entire part of it, using min distance
	RotatedRect enlargeRect(char *spanImg, CvBox2D rectCar, char *meanCar, int meanNum, float *minDis, int Rwidth, int Rheight);
	// divide cars of inner contour in a contour of cluster of vehicles
	void innerRect(char *spanImg, IplImage *cut, CvRect rectB, CvRect rectCar, char *meanCar, int meanNum, char *addHead, int Rwidth, int Rheight);
	// judge whether a car is vehicles(can be used to judge other classes)
	bool isVehicles(char *spanImg, CvBox2D rectCar, char *meanVehicles, int meanNum, float minDis);
	// finding rects in cluster of vehicles, fit area max
	RotatedRect enlargeRect2(char *spanImg, IplImage *cut, IplImage *innerCut, CvRect currentRect, CvBox2D rect, int Rwidth, int Rheight);
	// using enlargeRect2 to find inner cars
	void innerRect2(char *spanImg, IplImage *cut, IplImage *Tcut, CvRect rectB, CvRect rectCar, char *meanCar, char *addHead, int Rwidth, int Rheight);
	// find mean minDistance of each catenum's sample to centre model
	float *meanMinDistance(char **countSave, char **traitpath, char *cateArea, int catenum, int traitNum);
	// calculate features of contours
	void calculateFeature(char *contourPath, char *spanImg, char *specimen, char **countSave, char **geometry,  char **texture, char **saveImage, int threContour);

// function.h

	float **BinToMatrix(char *filename, int width, int height);
	// convert bin files to float vector
	float *BinToVector(char *filename, int length);
	// convert bin files to int vector
	int *intBinToVector(char *filename, int length);
	// convert bmp images to unsigned char matrix 
	unsigned char **BmpToMatrix(char *filename, int width, int height);
	// allocate memory for a m*n float matrix
	float **fAllocateMemory(int m, int n);
	// free memory for a m lines float matrix
	void fFreeMemory(float **matrix, int m);
	// allocate memory for a m*n int matrix
	int **intAllocateMemory(int m, int n);
	// allocate memory for a m*n double matrix
	double **dAllocateMemory(int m, int n);
	// allocate memory for a m*n char matrix
	char **charAllocateMemory(int m, int n);
	// free memory for a m lines unsigned char matrix
	void UcFreeMemory(unsigned char **matrix, int m);
	// free memory for a m lines int matrix
	void intFreeMemory(int **matrix, int m);
	// free memory for a m lines double matrix
	void dFreeMemory(double **matrix, int m);
	// free memory for a m lines char matrix
	void charFreeMemory(char **matrix, int m);
	// normalize a float array
	
	// normalize arrayA to [0 ,1]
	void fNormalize(float *arrayA, int arraysize);
	// save contour's bounding rect 
	void saveBoundingRect(IplImage* image, char *address, CvRect rect);

// MinDist
	
	// find the num of min value in vector 返回最小值所在位置
	int findMin(float *vector, int arraysize);
	// find the num of max value in vector 返回最大值所在位置
	int findMaxinA(float *vector, int arraysize);
	// calculate the euclidean distance of sample and mean 计算样本与样本均值间欧氏距离
	float euclidean(float *sample, float *mean, int dim);
	// find mean of trait and whichCate 算whichCate类的traitpath特征均值
	float findMean(char *maskPath, char* traitPath, int contourCount, int whichCate);
	// calculate mean vector of chosen traits 计算每类每个特征均值并以矩阵形式返回
	float **findCenter(char** traitpath, int dim, int num, int catenum, char *maskpath);
	// calculate mean of traits for whichCate rects

// BDistance
	float cateMean(char *maskPath, char* traitPath, int contourCount, int whichCate);
	// find max number in the array
	float findMax(float *array, int arraysize);
	 // function used in chosenmean
	int findminnum(float *array, int lastmin, int arraysize);
	float BDistanceCate(char *mask, char *trait, char *contourC, int cate1, int cate2);

// rect1
	void savefloat(float *data, char* name, int length);  
	// Write length int numbers to file
	void saveint(int *data, char* name, int length);  
	// find k of line(x1, y1)->(x2, y2)
	float findk(int x1, int x2, int y1, int y2);
	// find b of line(x1, y1)->(x2, y2)
	float findb(int x1, int x2, int y1, int y2);
	// calculate GLCM by using matrix span
	void calculateGraylevelmatrix(unsigned char **span, float **GLCMatrix, int x, int y, int width, int height);
	// change data in rawImage to graylevel
	void changeLevel(unsigned char **rawImage, int width, int height);
	// calculate GLDV by using matrix span
	void calculateGLDV(int **GLDVsta, float *GLDVector, int x, int y);
	// Write length float numbers to txt file
	void savetxt(float*data, char* name, int length) ;
	// find minnum in arrayA that is not zero
	int findMinNotZero(float *arrayA, int arraysize);
	float cateDev(char *maskPath, char *traitPath, int contourCount, int cate);
	CString GetModuleDir();
	char* connectChar(char *addHead, char *addTail);
public:
		unsigned char** cmatrix(int w,int h);

	void OnPath();
	void StartDir(const CString &strfile1);
	void RunDir(const CString &strfile2);
	void drawRect::drawAllTemplate(CString image);
};

