#include "stdafx.h"
#include "drawRect.h"
#define GrayLayerNum 16
#define pi 3.1415926

int tempi;//��ʱ�ļ�����ͼ�������
int counts;//������ͼ�������
bool dir;//���ü������־λ
CString m_strPath;//������·��
CString* templaPath[300];//��������ͼ��·��
//char path[300][80];

drawRect::drawRect(void)
{
}


drawRect::~drawRect(void)
{

}

// retunre CString path of exe file
CString drawRect::GetModuleDir()
{
	HMODULE module = GetModuleHandle(0);
	char pFileName[MAX_PATH];
	GetModuleFileName(module, pFileName, MAX_PATH);
 
	CString csFullPath(pFileName);
	int nPos = csFullPath.ReverseFind( _T('\\') );
	if( nPos < 0 )
		return CString("");
	else
		return csFullPath.Left( nPos );
}

// return char* = addHead + addTail
char* drawRect::connectChar(char *addHead, char *addTail)
{
	char *address = (char *)calloc(sizeof(char), 80);
	strcpy(address, addHead);
	strcat(address, addTail);
	return address;
} // connectChar

/* 
 input: (��ťѡ��)
    CString cutpath: ����Ŀ�������ַ
    CString rawImage: ԭͼ��ַ
    CString sample: ������ַ
 output: ���Զ�������exe�����ļ��У�
	roi: Ŀ����Ƭ
*/
void drawRect::main_Rect(CString cutpath, CString rawImage, CString sample)
{
	int Rwidth = 35, Rheight = 70; // һ�ѳ�ճ��һ��ʱ����ÿ�����εĴ�С
	// �浽exe����
	USES_CONVERSION;   
	CString addLocation =  GetModuleDir();
	char *addLoHead=T2A(addLocation.GetBuffer(0));
	CString newAddress = addLocation+"\\newFloder";
	CString ROIfolder = addLocation + "\\roi";
	addLocation.ReleaseBuffer();
	CreateDirectory(newAddress, NULL); 
	CreateDirectory(ROIfolder, NULL); 
	char *addressHead = connectChar(addLoHead, "\\newFloder\\"); 
	char *addROI = connectChar(addLoHead, "\\roi\\");

	const int contourThreshold = 200;   // ��ͨ�����С������������ͨ��ȥ��
	const int catenum = 5; // building, vehicle, vehicles, aircraft, other   ��������Ŀ
	const int traitNum = 13; // geometry + texture = 12 traits ��ȡ��������
	//printf("contourThre = %d, catenum = %d\n", contourThreshold, catenum);
	// input 2* .bmp 
	//char *span = "E:/PlayWithC/data/MiniSAR3_15.bmp";
	//char *cutpath = "E:/PlayWithC/data/MiniSAR3target.bmp";
	//char *specimen = "E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/feature/specimen.bmp";
	// some creating files, save path
	char *contourCate = connectChar(addressHead, "contourCate.bin"); // ��ͨ���ǽ����0-������1-catenum: ����
	char **countPath = charAllocateMemory(2, 80); // ��ͨ����Ŀ����С���о���(MER)��Ŀ
	countPath[0] =  connectChar(addressHead, "contourCount.bin");
	countPath[1] =  connectChar(addressHead, "rectCount.bin");
	char **geometry = charAllocateMemory(8, 80);
	geometry[0] = connectChar(addressHead, "rectCate.bin"); // ��ͨ��MER�ı�ǽ����0-������1-catenum: ����
	geometry[1] = connectChar(addressHead, "rectAngle.bin");// ���νǶ�
	geometry[2] = connectChar(addressHead, "rectArea.bin"); // �������
	geometry[3] = connectChar(addressHead, "rectFit.bin");  // ������϶�
	geometry[4] = connectChar(addressHead, "LWratio.bin");  // �����
	geometry[5] = connectChar(addressHead, "formCoef.bin"); // ��״ϵ��
	geometry[6] = connectChar(addressHead, "graymean.bin"); // �Ҷ�ƽ��ֵ
	char **texture = charAllocateMemory(8, 80);
	texture[0] = connectChar(addressHead, "GLCMappclust.bin");  // ��������
	texture[1] = connectChar(addressHead, "GLCMcontrast.bin");
	texture[2] = connectChar(addressHead, "GLCMentropy.bin");
	texture[3] = connectChar(addressHead, "GLCMmaxProba.bin");
	texture[4] = connectChar(addressHead, "GLDVASM.bin");
	texture[5] = connectChar(addressHead, "GLDVcontrast.bin");
	texture[6] = connectChar(addressHead, "GLDVentropy.bin");
	texture[7] = connectChar(addressHead, "GLDVmean.bin");
	char **saveImage = charAllocateMemory(3, 80);
	saveImage[0] = connectChar(addressHead, "withRect.bmp");  // ��ͼ������
	saveImage[1] = connectChar(addressHead, "drawRect.bmp");
	saveImage[2] = connectChar(addressHead, "approxpoly.bmp");
	// �������̣�����-����-���Ͼ���-С���MER-�������MER(if any)��С����-���ͼ
	/*----------------------------------------------------------------------------------------------------------*/
	/* �������Ҫ�����ͨ��������� */
	// ���㲢�洢С���ε�����
	printf("\nProcess of calculating features.\n");
	// input 3 bmp
	char *cutImg= T2A(cutpath.GetBuffer(0));
	char *spanImg = T2A(rawImage.GetBuffer(0));
	char *speciment = T2A(sample.GetBuffer(0));
	calculateFeature(cutImg, spanImg, speciment, countPath, geometry, texture, saveImage, contourThreshold);
	printf("\nEnd of calculating features.\n");
	// ��ÿ���������ֵ����, ������Ͼ���
	printf("\nProcess of calculating means and bd.\n");
	float **meanCate = fAllocateMemory(catenum, traitNum);  // ��ֵ
	float **bdv = fAllocateMemory(catenum-1, traitNum);	    // ���Ͼ���
	char **traitpath = charAllocateMemory(traitNum, 80);    // ����·��
	traitpath[0] = geometry[2];
	for(int k = 0; k < 4; k++) {
		traitpath[k+1] = geometry[3+k];
		traitpath[k+5] = texture[k];
		traitpath[k+9] = texture[k+4];
	}
	for(int d = 0; d < traitNum; d++ ) {
		for(int ca = 0; ca < catenum; ca++) 
			meanCate[ca][d] = cateMean(geometry[0], traitpath[d], intBinToVector(countPath[0], 1)[0], ca+1);  // ca + 1!! !
		for(int i = 0; i < catenum-1; i++) 
			bdv[i][d] = BDistanceCate(geometry[0], traitpath[d], countPath[0], i+1, catenum); // i + 1!!!!
	} // for d
	char *carMeanPath = connectChar(addressHead, "carMean.bin");
	char *vehiclesMeanPath = connectChar(addressHead, "vehicelsMean.bin");
	int meanNum = 6;	
	float *carMean = (float *)calloc(sizeof(float), 6);
	carMean[0] = meanCate[1][4];
	carMean[1] = meanCate[1][5];
	carMean[2] = meanCate[1][6];
	carMean[3] = meanCate[1][8];
	carMean[4] = meanCate[1][9];
	carMean[5] = meanCate[1][11];
	float *VehiclesMean = (float *)calloc(sizeof(float), 6);
	VehiclesMean[0] = meanCate[2][4];
	VehiclesMean[1] = meanCate[2][5];
	VehiclesMean[2] = meanCate[2][6];
	VehiclesMean[3] = meanCate[2][8];
	VehiclesMean[4] = meanCate[2][9];
	VehiclesMean[5] = meanCate[2][11];
	savefloat(carMean, carMeanPath, meanNum);
	savefloat(VehiclesMean, vehiclesMeanPath, meanNum);
	// ����������ģ���ƽ������
	float *meanDistance = (float *)calloc(sizeof(float), catenum-1);
	printf("\nEnd of calculating means and bd.\n");
	/*----------------------------------------------------------------------------------------------------------*/
	// ŷ�Ͼ��� 
	printf("\nProcess of calculting Euclidistance.\n");
	const int bdim = 8, bdim2 = 3; // ѡbdֵ����bdim�� bidm2��������һ�ѳ��ͷɻ�
	char **traitE = charAllocateMemory(bdim, 80);
	traitE[0] = geometry[4];
	traitE[1] = geometry[6];
	traitE[2] = geometry[3];
	traitE[3] = texture[0];
	traitE[4] = texture[1];
	traitE[5] = texture[3];
	traitE[6] = texture[4];
	traitE[7] = texture[6];
	char **traitD = charAllocateMemory(bdim2, 80);  // for bidm2
	traitD[0] = texture[1];
	traitD[1] = texture[0];
	traitD[2] = texture[6];
	char **usepath = charAllocateMemory(2, 80); 
	usepath[0] =  geometry[0];  // ���ο����
	usepath[1] =  geometry[2];  // ���ο����
	char **cateImage = charAllocateMemory(2, 80); // �������ͼ
	cateImage[0] =connectChar(addressHead, "cateImage.bmp");
	cateImage[1] = connectChar(addressHead, "cateImage2.bmp");
	//char *cateCount = connectChar(addressHead, "cateNumber.txt";
	int *category = findTarget(cutImg, spanImg, countPath, traitE, traitD, usepath, contourThreshold, bdim, bdim2, catenum);
	saveint(category, connectChar(addressHead, "category.bin"), intBinToVector(countPath[1], 1)[0]); // ������
	float *meanMin = meanMinDistance(countPath, traitE, geometry[0], catenum, bdim);
	//savetxt(meanMin, "E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/meanMin.txt", catenum);
	char **addHead = charAllocateMemory(2 , 80);    
	addHead[0] = connectChar(addLoHead, "\\roi\\");
	addHead[1] = connectChar(addLoHead, "\\roi\\");
	// ����������+��Ƭ
	drawCategoryRect(cutImg, spanImg, countPath[1],carMeanPath, vehiclesMeanPath, meanNum, category, contourThreshold, cateImage ,catenum, addHead, Rwidth, Rheight);
	printf("\nEnd of calculting center and Euclidistance.\n");
    /*----------------------------------------------------------------------------------------------------------*/
} // main

/*
	�����������threContour��ÿ����ͨ���geometry��texture���������浽countSave,geomertry,texture�ļ�
	input:
		countourPath������Ŀ������
		spanImg: ԭͼ�Ҷ�ͼ
		speciemen: ����
	output:
		��Ƭ
*/
void drawRect:: calculateFeature(char *contourPath, char *spanImg, char *specimen, char **countSave, char **geometry,  char **texture, char **saveImage, int threContour)
{
	Mat srcmat, spanmat, colormat, spanmat1;
	IplImage *src, *gray, *span, *span1, *bg, *color;  
    CvMemStorage* storage = cvCreateMemStorage(0);  // ��ͨ��
	CvMemStorage* storage1 = cvCreateMemStorage(0);  // ��������
    CvSeq *contour = 0, *vertex = 0;  

	srcmat = imread(contourPath);
	spanmat = imread(spanImg);
	spanmat1 = imread(spanImg);
	colormat = imread(specimen);
	src = &IplImage(srcmat);
	span = &IplImage(spanmat);
	span1 = &IplImage(spanmat1);
	color = &IplImage(colormat);
	//src = cvLoadImage(contourPath,1);   // building area
	//span = cvLoadImage(spanImg,1);   // span.bmp
	//span1 = cvLoadImage(spanImg,1);   // span.bmp
	//color = cvLoadImage(specimen,1); // specimen.bmp
	bg = cvCloneImage(span);
	gray=cvCreateImage(cvSize(src->width,src->height),src->depth,1);  
    cvCvtColor(src,gray,CV_BGR2GRAY); 
	int height=src->height;   
	int width=src->width;
	//printf("height =%d, width = %d\n", height, width);
	int aircraft = 0, building = 0, vehicle = 0, other = 0, vehicles = 0;
	int **colorM = intAllocateMemory(height, width);

	CvScalar ss;   // mark specimens
	for(int i = 0; i < height; i++)
		for(int j = 0; j < width; j++) {
			ss = cvGet2D(color, i, j);   // bgr
			if(ss.val[0]<50 && ss.val[1]<50 && ss.val[2]>200) {  // red - bgr
				colorM[i][j] = 1; // building
				assert(colorM[i][j] >0);
			} // if red
			else if(ss.val[0]<50 && ss.val[1]>200 && ss.val[2]<50) {  // green - bgr
				colorM[i][j] = 5; // other
				assert(colorM[i][j] >0);
			} // if green
			else if(ss.val[0]>200 && ss.val[1]<50 && ss.val[2]<50) {  // blue - bgr
				colorM[i][j] = 2; // vehicle
				assert(colorM[i][j] >0);
			} // if blue
			else if(ss.val[0]>200 && ss.val[1]>200 && ss.val[2]<50) {  // light blue - bgr
				colorM[i][j] = 3; // vehicles
				assert(colorM[i][j] >0);
			} // if blue
			else if(ss.val[0]>200 && ss.val[1]<50 && ss.val[2]>200) {  // pink - bgr
				colorM[i][j] = 4; // aircraft
				assert(colorM[i][j] >0);
			} // if pink
			else 
				colorM[i][j] = 10; // û����
		} // for i,j

	// for texture traits
	unsigned char **spanMatrix = BmpToMatrix(spanImg, width, height);    // data in span [0, 255]
	unsigned char **grayValue = BmpToMatrix(spanImg, width, height);  // gray value in spanImg.bmp
	changeLevel(spanMatrix, width, height);
	int **GLDVsta = intAllocateMemory(height, width);  // for calculating GLDV
	for(int p = 0; p < height - 1; p++)
		for(int q =0; q < width - 1; q++) {
			GLDVsta[p][q]=spanMatrix[p][q]-spanMatrix[p+1][q+1];
			if(GLDVsta[p][q] > 15)
				printf("error: GLDV > 15. \n");
		} // for 
	
    int contourCount = cvFindContours(gray,storage,&contour,sizeof(CvContour),CV_RETR_TREE,CV_CHAIN_APPROX_SIMPLE);  // number of contours
 
	int rectCount = 0; // number of MER
	int smallMER = 0;  // no use
	int maxX = 0, maxY = 0, minX = 0, minY = 0; // coordinate of each MER

	int *rectCate = (int *)calloc(sizeof(int ), contourCount);  // 0-������, 1~catenum��
	float *rectArea = (float *)calloc(sizeof(float), contourCount);    // sava area of rects ��������
	float *contourArea = (float *)calloc(sizeof(float), contourCount);    // save data of contour ����ͨ�����
	float *contourPeri = (float *)calloc(sizeof(float), contourCount);    // save data of contour ����ͨ���ܳ�
	int *rectL = (int *)calloc(sizeof(int), contourCount);  // length of rects ����γ�
	int *rectW = (int *)calloc(sizeof(int), contourCount);  // width of rects ����ο�
	float *centerX = (float *)calloc(sizeof(float), contourCount);  // x of gravity center ����X����
	float *centerY = (float *)calloc(sizeof(float), contourCount);  // y of gravity center ����Y����
	float *rectAngle = (float *)calloc(sizeof(float), contourCount);    // angle of rects ˮƽ����ʱ����ת���������ĵ�һ���߼н�, ���ȱ�ʾ
	float *rectFit = (float *)calloc(sizeof(float), contourCount);  // similarity degree of contour and rect ������϶�
	float *formCoef = (float *)calloc(sizeof(float), contourCount);  // coefficient of form ��ͨ����״ϵ��
	float *LWratio = (float *)calloc(sizeof(float), contourCount);  // rect.Length/rect.Width ���γ����
	float *grayMean = (float *)calloc(sizeof(float), contourCount);  // mean of gray value ��������ҶȾ�ֵ
	//int *rightAngle = (int *)calloc(sizeof(int ), contourCount);  // 0-no, 1-yes

	float *GLDVASM = (float *)calloc(sizeof(float), contourCount); // GLDVASM
	float *GLDVentropy = (float *)calloc(sizeof(float), contourCount); // GLDVentropy
	float *GLDVcontrast = (float *)calloc(sizeof(float), contourCount); //GlDVcontrast
	float *GLDVmean = (float *)calloc(sizeof(float), contourCount); //GlDVmean
	float *GLCMentropy = (float *)calloc(sizeof(float), contourCount);  //GLCMentropy
	float *GLCMcontrast = (float *)calloc(sizeof(float), contourCount);  //GLCMcontrast
	float *GLCMappclust = (float *)calloc(sizeof(float), contourCount);  //GLCMappclust
	float *GLCMmaxProba = (float *)calloc(sizeof(float), contourCount); //GLCMmaxProba
	
	// mѭ��-��ÿ����ͨ����
	for(int m = 0;contour!=0 ;contour=contour->h_next, m++)  {
	vertex = cvApproxPoly(contour, sizeof(CvContour), storage1, CV_POLY_APPROX_DP, cvContourPerimeter(contour)*0.02,0); // ��������
	int countR = 0; // points counted in this rect ���ο��еļ����
	//printf("m = %d\n", m);
	CvBox2D rect=cvMinAreaRect2(contour,storage);  // building box
    CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rect, rect_pts0);
	contourArea[m] = fabs(cvContourArea(contour));
	contourPeri[m] = cvContourPerimeter(contour);
	//printf("��ͨ�������%d\n", contourArea[m]);

	//��ΪcvPolyLineҪ��㼯������������CvPoint**  
	//����Ҫ�� CvPoint2D32f �͵� rect_pts0 ת��Ϊ CvPoint �͵� rect_pts  
	//������һ����Ӧ��ָ�� *pt  
	int npts = 4,k=0;   
	CvPoint rect_pts[4], *pt = rect_pts;  
  
	//printf("��ͨ������С��Ӿ��ζ�������ֱ�Ϊ:\n");  
	for (int i=0; i<4; i++)  {  
		rect_pts[i]= cvPointFrom32f(rect_pts0[i]);
		//printf("%d %d\n",rect_pts[i].x,rect_pts[i].y);  
		rectL[m]=(int)sqrt((pow((rect_pts[0].x-rect_pts[1].x),2)+pow((rect_pts[0].y-rect_pts[1].y),2)));  
		rectW[m]=(int)sqrt((pow((rect_pts[0].x-rect_pts[3].x),2)+pow((rect_pts[0].y-rect_pts[3].y),2)));  
        if(rectL[m] < rectW[m]) { 
             k = rectL[m];  
             rectL[m] = rectW[m];  
             rectW[m] = k;  
           } // if
		}  // for i [0, 3]
	minX = rect_pts[1].x;
	minY = rect_pts[2].y;
	maxX = rect_pts[3].x;
	maxY = rect_pts[0].y;
	if(maxY >= height)
		maxY = height -1;
	if(maxX >= width)
		maxX = width - 1;
	if(minX < 0)
		minX = 0;
	if(minY < 0)
		minY =0;
	//printf("minX = %d, minY = %d, maxX = %d, maxY = %d\n", minX, minY, maxX, maxY);
	rectArea[m] = (float)rectW[m] * rectL[m];
	rectAngle[m] = rect.angle;
	centerX[m] = rect.center.x;
	centerY[m] = rect.center.y;
	if(rectArea[m] > 0 && contourArea[m] > 0) {
		rectFit[m] = (float)contourArea[m]/rectArea[m];
		formCoef[m] = (float)contourPeri[m]/(4 * sqrt(contourArea[m]));
		LWratio[m] = (float)rectL[m]/rectW[m];
	} // if
	else if(contourArea[m] = 0) {
		rectFit[m] = 0;
		formCoef[m] = 0;
		LWratio[m] = 0;
	} // else if
	//printf("��С��Ӿ��εĳ���%d����%d�������%f, �Ƕ�: %f, ���ģ�(%.2f, %.2f)\n",rectL[m], rectW[m], rectArea[m],rectAngle[m],rect.center.x,rect.center.y); 
	//printf("������϶ȣ�%.4f, ��״ϵ��: %.4f\n",rectFit[m], formCoef[m]);
	
	// Draw rects whose area no less than Threshold   
	if(contourArea[m] >= threContour ) { // ȥ��С�������
		rectCate[m] = colorM[(int)rect.center.y][(int)rect.center.x];
		if(rectCate[m] == 1)
			building++;
		else if(rectCate[m] == 2)
			vehicle++;
		else if(rectCate[m] == 3)
			vehicles++;
		else if(rectCate[m] == 4)
			aircraft++;
		else if(rectCate[m] == 5)
			other++;
		//cvFillPoly(bg, &pt, &npts, 1, CV_RGB(255,255,255));  //���黭�ɰ׵�
		//cvPolyLine(src, &pt, &npts, 1, 1, CV_RGB(255,0,0), 2); // �������
		//cvDrawContours (span1, vertex, CV_RGB(255,0,0),CV_RGB(0,0,100),1,2,8,cvPoint(0,0)); // ��������
		//float vcos = mincos(vertex);
		//if(vcos < 0.1)
		//	rightAngle[m] = 1;

		rectCount++;
		printf("areaR = %f\n", rectArea[m]);
		//system("pause");
		smallMER++;
		//cvFillPoly(bg, &pt, &npts, 1, CV_RGB(255,0,0));   // С���������

	// Extra textural features in rects ��ȡ��������������
	float **GLCMatrix = fAllocateMemory(GrayLayerNum, GrayLayerNum);    // GLCM
	float *GLDVector = (float *)calloc(sizeof(float), 2*GrayLayerNum - 1);  // GLDV
	for(int x = minX; x <= maxX; x++)  // ��������ĵ�
		for(int y = minY; y <= maxY; y++) {    
			if(y <= (findk(rect_pts[0].x, rect_pts[1].x, rect_pts[0].y, rect_pts[1].y)*x + findb(rect_pts[0].x, rect_pts[1].x, rect_pts[0].y, rect_pts[1].y))
				&& y <= (findk(rect_pts[0].x, rect_pts[3].x, rect_pts[0].y, rect_pts[3].y)*x + findb(rect_pts[0].x, rect_pts[3].x, rect_pts[0].y, rect_pts[3].y))
				&& y >= (findk(rect_pts[1].x, rect_pts[2].x, rect_pts[1].y, rect_pts[2].y)*x + findb(rect_pts[1].x, rect_pts[2].x, rect_pts[1].y, rect_pts[2].y))
				&& y >= (findk(rect_pts[2].x, rect_pts[3].x, rect_pts[2].y, rect_pts[3].y)*x + findb(rect_pts[2].x, rect_pts[3].x, rect_pts[2].y, rect_pts[3].y))) {
				countR++;
				// calculate GLCM ����GLCM����
				calculateGraylevelmatrix(spanMatrix, GLCMatrix,x,y, width, height);
				// calculate GLDV ����GLDVʸ��
				calculateGLDV(GLDVsta, GLDVector, x, y);	
				// calculae gray value ����ҶȾ�ֵ
				grayMean[m] += grayValue[y][x];
		} // if
	}  // for y
	assert( grayMean[m]>0);
	// normalization of GLCM GLCM�����һ��
	float sumGLCM = 0;
	for(int mm = 0; mm < GrayLayerNum; mm++)
		for(int n =0; n < GrayLayerNum; n++)
			sumGLCM += GLCMatrix[mm][n];
	for(int mm = 0; mm < GrayLayerNum; mm++)
		for(int n = 0; n < GrayLayerNum; n++)
			GLCMatrix[mm][n] = GLCMatrix[mm][n]/sumGLCM;
	//for(int k = 0; k < GrayLayerNum; k++)
	//	for(int f =0; f<GrayLayerNum;f++)
	//		printf("%f\n", GLCMatrix[k][f]);
	// calculate textural feature of GLCM ����GLCM����
	float *ui = (float *)calloc(sizeof(float), GrayLayerNum);
	float *uj = (float *)calloc(sizeof(float), GrayLayerNum);
	float contrasttemp=0;
	float entropytemp=0;
	float appclustemp=0;
	float maxProbatemp=0;
	for(int gi = 0 ; gi < GrayLayerNum; gi++) {  
		for( int gj = 0; gj < GrayLayerNum; gj++)
			{
				uj[gi]+=GLCMatrix[gi][gj];
			}
				uj[gi]=uj[gi]/GrayLayerNum;
		} // for
		for( int gj = 0 ; gj < GrayLayerNum; gj++)
			{    
				for( int gi = 0; gi < GrayLayerNum; gi++)
				{
					ui[gj]+=GLCMatrix[gi][gj];
				}
			ui[gj]=ui[gj]/GrayLayerNum;
			} // for
	for( int k = 0 ; k < GrayLayerNum; k++)
		{    
			for( int z = 0; z< GrayLayerNum; z++)
			{
				contrasttemp += GLCMatrix[k][z]*abs((k-z)*(k-z));

				if(GLCMatrix[k][z]<1&&GLCMatrix[k][z])
				entropytemp -= GLCMatrix[k][z]*log10(GLCMatrix[k][z]);
				if(GLCMatrix[k][z]-maxProbatemp>=0)
					maxProbatemp=GLCMatrix[k][z];				
			appclustemp += (k+z-ui[k]-uj[z])*(k+z-ui[k]-uj[z])*GLCMatrix[k][z];
				}
		} // for
	GLCMcontrast[m]=contrasttemp;			
	GLCMentropy[m]=entropytemp;
	GLCMappclust[m]=appclustemp;
	GLCMmaxProba[m]=maxProbatemp;
	free(ui);
	free(uj);

	// normalization of GLDV GLDVʸ����һ��
	float sumGLDV = 0;
	for(int i= 0; i < 2*GrayLayerNum - 1; i++) 
		sumGLDV += GLDVector[i];
	for(int i= 0; i < 2*GrayLayerNum - 1; i++) 
		GLDVector[i] = GLDVector[i]/sumGLDV;
	//for(int f =0; f<2*GrayLayerNum-1;f++)
	//		printf("%f\n", GLDVector[f]);
	// calculate textural feature of GLDV ����GLDV����
	float GLDVcontrasttemp=0;
	float GLDVentropytemp=0;
	float GLDVASMtemp=0;
	float GLDVmeantemp=0;
	for (int k = 0; k < (2*GrayLayerNum-1); k++) {
		GLDVcontrasttemp = GLDVcontrasttemp+k*k*GLDVector[k];
		if(GLDVector[k])
			GLDVentropytemp-=GLDVector[k]*log(GLDVector[k]);
		GLDVASMtemp+=GLDVector[k]*GLDVector[k];
		GLDVmeantemp+=k*GLDVector[k];
		} // for
	GLDVcontrast[m]=GLDVcontrasttemp;
	GLDVentropy[m]=GLDVentropytemp;
	GLDVASM[m]=GLDVASMtemp;
	GLDVmean[m]=GLDVmeantemp;
	// claculate mean of halpha ����Halpha��ֵ
	if(countR == 0)
		countR = 1;
	// calculate gray mean ����Ҷ�ֵ��ֵ
	grayMean[m] = grayMean[m]/countR;
	fFreeMemory(GLCMatrix, GrayLayerNum);
	free((void *)GLDVector);
 		} // if conoturArea > thresholdContour
		countR = 0;
	} // for m
	printf("building = %d, vehicle = %d, vehicles = %d, aircraft = %d, other = %d\n", building, vehicle, vehicles, aircraft, other);
	printf("������Ŀ�� %d; ��������Ŀ�� %d\n", rectCount, contourCount);

	//int *rightAngle1 = (int *)calloc(sizeof(int), rectCount);
	//int pp =0;
	//for(int p = 0; p < contourCount;p++) {
	//	if(GLCMappclust[p] != 0 ) {
	//		rightAngle1[pp] = rightAngle[p];
	//		pp++;
	//	}
	//}
	//assert(pp == rectCount);

	// �漸��������bin�ļ�
	saveint(&contourCount, countSave[0], 1);
	saveint(&rectCount, countSave[1], 1);
	saveint(rectCate, geometry[0], contourCount);
	savefloat(rectAngle, geometry[1], contourCount);
	savefloat(rectArea, geometry[2], contourCount);
	savefloat(rectFit, geometry[3], contourCount);
	savefloat(LWratio,geometry[4], contourCount);
	savefloat(formCoef, geometry[5], contourCount);
	savefloat(grayMean, geometry[6], contourCount);
	//saveint(rightAngle1, geometry[7], rectCount);

	// ������������bin�ļ�
	savefloat(GLCMappclust, texture[0], contourCount);
	savefloat(GLCMcontrast, texture[1], contourCount);
	savefloat(GLCMentropy, texture[2], contourCount);
	savefloat(GLCMmaxProba, texture[3], contourCount);
	savefloat(GLDVASM, texture[4], contourCount);
	savefloat(GLDVcontrast, texture[5], contourCount);
	savefloat(GLDVentropy, texture[6], contourCount);
	savefloat(GLDVmean, texture[7], contourCount);
	
	// ��ͼ��bmp�ļ� 
	//cvSaveImage(saveImage[0], src,0);
	//cvSaveImage(saveImage[1], bg,0);
	//cvSaveImage(saveImage[2], span1,0);

	// �ͷ��ڴ�
	//cvReleaseImage(&src);    
	cvReleaseImage(&gray);
	//cvReleaseImage(&span);
	//cvReleaseImage(&span1);
	cvReleaseImage(&bg);
	//cvReleaseImage(&color);
	src = NULL;
	span = NULL;
	span1 = NULL;
	color = NULL;
	srcmat.release();
	spanmat.release();
	spanmat1.release();
	colormat.release();
	cvReleaseMemStorage(&storage);
	cvReleaseMemStorage(&storage1);
	free(rectArea);
	free(contourArea);
	free(rectAngle);
	//free(rightAngle1);
	free(rectCate);
	free(rectFit);
	free(rectL);
	free(rectW);
	free(formCoef);
	free(grayMean);
	free(centerX);
	free(centerY);
	free(contourPeri);
	free(LWratio);
	//free(rightAngle);
	free(GLCMappclust);
	free(GLCMcontrast);
	free(GLCMentropy);
	free(GLCMmaxProba);
	free(GLDVASM);
	free(GLDVcontrast);
	free(GLDVentropy);
	free(GLDVmean);
	intFreeMemory(colorM,height);
	UcFreeMemory(spanMatrix, height);
	UcFreeMemory(grayValue, height);
	intFreeMemory(GLDVsta, height);
} // calculateRect

/*
	ŷ�Ͼ����MER���з���
	input:
		contourPath: ����Ŀ������
		spanImg: ԭͼ�Ҷ�ͼ
		countSave: ��ͨ��, MER
		traitPath: �����Ĵ洢·��
		traitD: ����aircraft and vehicles�������洢·��
		rectCate: [0]���ο���� [1]���ο����
		threContour: ��ͨ����ֵ��С���򲻴���
		traitNum: ŷ�Ͼ���ʹ�õ�������Ŀ
		traitNumD������aircraft��vehiclesʹ�õ�������Ŀ
		catenum: �����Ŀ��Ŀǰ5��
	output:
		����ŷ�Ͼ�������int*
*/
int *drawRect::findTarget(char *contourPath, char *spanImg, char **countSave, char **traitpath, char **traitD, char **rectCate, int threContour, int traitNum, int traitNumD, int catenum)
{
	 int carT = 1000, airT = 5000; // delete, minmum area
	 int num = intBinToVector(countSave[0], 1)[0]; // number of contours ��ͨ����Ŀ  
	 //printf("num = %d\n",num);
	 int rectnum = intBinToVector(countSave[1], 1)[0]; // number of speciments ������Ŀ
	 //printf("rectnum = %d\n", rectnum);
	 int *smallMER = intBinToVector(rectCate[0], num); // ��0-С����
	 int *rectArea =  intBinToVector(rectCate[1], num); // area of the rects
	 float **sample = fAllocateMemory(traitNum, num);  // ��������
	 float **sample1 = fAllocateMemory(traitNum, num);  // ��������
	 for(int j = 0; j < traitNum; j++) {
		 sample[j] = BinToVector(traitpath[j], num);
	 } // if j
	 for(int j = 0; j < traitNumD; j++) {
		 sample1[j] = BinToVector(traitD[j], num);
	 } // if j
	 float **centre = findCenter(traitpath, traitNum, num, catenum, rectCate[0]); 
	 float **centre1 = findCenter(traitD, traitNumD, num, catenum, rectCate[0]);   
	 int *EucliCategory = (int *)calloc(sizeof(int), rectnum);   // ŷ�Ͼ��������
	 float *distanceE = (float *)calloc(sizeof(float), catenum); // ŷ�Ͼ���
	 float **tempSample = fAllocateMemory(rectnum, traitNum);
	 int *tempArea = (int *)calloc(sizeof(int), rectnum);
	 float **tempSample1 = fAllocateMemory(rectnum, traitNumD);
	 int aa = 0, nn =0;
	 for(int a = 0; a < num; a++) {
		 if(smallMER[a] != 0) { 
			 tempArea[aa] = rectArea[a]; 
			 for(nn = 0; nn < traitNum; nn++) {
				 tempSample[aa][nn] = sample[nn][a];  // ��������
				 if(nn < traitNumD)
					tempSample1[aa][nn] = sample1[nn][a];  
			 } // for nn
		 aa++;		
		 } // if
	 } // for a 
	 assert(aa == rectnum);
	 for(int b = 0; b < rectnum; b++) {  // ���� [0, catenum-1]
		for(int m = 0; m < catenum; m++) {
			distanceE[m] = euclidean(tempSample[b], centre[m], traitNum);
			//assert(distanceE[m] > 0);
				} // for m
		EucliCategory[b] = findMinNotZero(distanceE, catenum) + 1;  // the address of min value
	 } // for b
	 // if category == 3 or 4, it need to be distingulished
	 for(int cc = 0; cc < rectnum; cc++) {
		 if(EucliCategory[cc] == 3 || EucliCategory[cc] == 4) {
			 float newDistance1 = euclidean(tempSample1[cc], centre1[2], traitNumD);
			 float newDistance2 = euclidean(tempSample1[cc], centre1[3], traitNumD);
			 newDistance1 = newDistance1>0?newDistance1:newDistance2;
			 newDistance2 = newDistance2>0?newDistance2:newDistance1;
			 EucliCategory[cc] = newDistance1 <= newDistance2 ? 3 : 4;
		 } // if 3 || 4
		 if(EucliCategory[cc] == 3 ) {
			 if(tempArea[cc] < carT)
				 EucliCategory[cc] = 5;
		 } // if 3
		 if(EucliCategory[cc] == 4 ) {
			 if(tempArea[cc] < airT)
				 EucliCategory[cc] = 5;
		 } // if 4
	 } // for cc
	 free((void *)smallMER);
	 free((void *)rectArea);
	 free(tempArea);
	 fFreeMemory(sample, traitNum);
	 fFreeMemory(centre, catenum);
	 fFreeMemory(sample1, traitNum);
	 fFreeMemory(centre1, catenum);
	 fFreeMemory(tempSample, rectnum);
	 fFreeMemory(tempSample1, rectnum);
	 free(distanceE);
	 return EucliCategory;
} // findtarget

// find and put templatImg in image, save in exe path 
void drawRect::fitTemplate(char * image, char *templateImg)
{
	// exe path
	USES_CONVERSION;
	CString addLocation =  GetModuleDir();
	char *addLoHead=T2A(addLocation.GetBuffer(0));
	CString addr = addLocation + "\\result.bmp";
	char *addresult=T2A(addr.GetBuffer(0));
	addr.ReleaseBuffer();
	IplImage *src = cvLoadImage(image, 1);  
    IplImage *srcResult = cvLoadImage(image, 3);  //������ʾ  
    IplImage *templat = cvLoadImage(templateImg, 1);  
    IplImage *result;  // ������Ž�� 
    int srcW, srcH, templatW, templatH, resultH, resultW;  
    srcW = src->width;  
    srcH = src->height;  
    templatW = templat->width;  
    templatH = templat->height; 
	//IplImage *tempSrc = cvCreateImage(cvSize(srcW, srcH),32, 1);
	//cvCopyImage(src, tempSrc);
	//int echo = (int)srcW*srcH/templatW/templatH;
	if(srcW <= templatW || srcH <= templatH )
		return;
    resultW = srcW - templatW + 1;  
    resultH = srcH - templatH + 1;  
    result = cvCreateImage(cvSize(resultW, resultH), 32, 1);    //  ƥ�䷽������Ľ����СֵΪfloat  
    cvMatchTemplate(src, templat, result, CV_TM_SQDIFF);     
    double minValue, maxValue;  
    CvPoint minLoc, maxLoc;  
    cvMinMaxLoc(result, &minValue, &maxValue, &minLoc, &maxLoc);  // minmum
	//assert(minLoc.x + templatW < srcW && minLoc.y+ templatH < srcH);
    CvScalar ss = cvGet2D(templat, 0, 0);  // bgr - boundary
	cvRectangle(srcResult, minLoc, cvPoint(minLoc.x + templatW, minLoc.y+ templatH), cvScalar(ss.val[2],ss.val[1],ss.val[0]),2);  // draw rect
	cvSaveImage(image, srcResult,0);
    cvReleaseImage(&result);  
    cvReleaseImage(&templat);  
    cvReleaseImage(&srcResult);  
    cvReleaseImage(&src);   
	return;
} // fitTemplate


/*
	���������ŷ�Ͼ������������ʹ���Ƭ
	input��
		rawImage: ����Ŀ����
		spanImg: ԭ�Ҷ�ͼ
		rectnum: MER��Ŀ
		carMeanPath��car ��������ֵ
		vehiclesMean�� һ�ѳ���ֵ
		meanNum: ��ֵ����������Ŀ
		category: ŷ�Ͼ�����
		threshouldContour����ͨ����ֵ��С�ڲ�����
		savepath: ��ͼ û����
		catenum:����Ŀ������
		addHead: ����Ƭ�ĵ�ַǰ��Ĳ��� [0]������[1]һ�ѳ����һ����
		Rwidth��Rheight: һ�ѳ�ʱÿ����Ŀ�͸�
*/
void drawRect::drawCategoryRect(char *rawImage, char *spanImg, char *rectnum, char *carMeanPath, char *vehiclesMean, int meanNum, int *category, int thresholdContour, char **savepath, int catenum, char **addHead, int Rwidth, int Rheight)
{
	/* src: contour of building and vegetation
	   gray: gray scale of span image
	   cspan, spancopy: copy of span
	   cspan: background of boundingrects, same as span
	   spancopy: show result of enlarge rects of catenum2: single vehicle
	   cutImg: threshold span, used to extract samll contours
	   cutImg1: eroded and opened process of src
	 */
	// cvLoadImage and cvCloneImage may lead to memory leak!!!
	Mat srcmat, spanmat;  // ��˵mat��ʽ���Ա����ڴ�й©
	IplImage *src, *gray, *span, *cspan, *cutImg, *cutImg1, *spancopy;  
    CvMemStorage* storage = cvCreateMemStorage(0);  // Creat new memory
    CvSeq* contour = 0;
	const int threDraw = 6000, threCar = 2000; // < a car >= many cars, a car or a part of it
	const int ostuValue = 170; // threshold image 
	const int threDis = 1.3; // if minDis > threDis, we consider it a mistaked classify
	const int miniflight = 300; // mostly flight has a relativly fixed size
	int addNum = 0;  // ��Ƭ���
	char *addTail = ".bmp";
	char *addHeadRect = (char *)calloc(sizeof(char), 60); 
	char *address = (char *)calloc(sizeof(char), 80); // address of ROI

	// loadimage 
	srcmat = imread(rawImage);
	src = &IplImage(srcmat);
	spanmat = imread(spanImg);
	span = &IplImage(spanmat);
	cspan = cvCreateImage(cvSize(span->width, span->height),span->depth, span->nChannels);
	cvCopy(span, cspan, NULL);
	spancopy = cvCreateImage(cvSize(span->width, span->height),span->depth, span->nChannels);
	cvCopy(span, spancopy, NULL);
	//src = cvLoadImage(rawImage,1);   // building area
	//span = cvLoadImage(spanImg,1);   // span.bmp
	//cspan = cvLoadImage(spanImg,1);
	//spancopy = cvLoadImage(spanImg,1);

	// ������ �ֿ�
	Mat imagemat = imread(spanImg);   // �Ҷ�ͼ
	cvtColor(imagemat, imagemat, CV_BGR2GRAY);  // ��ֵͼ
	threshold(imagemat, imagemat, ostuValue, 255, THRESH_BINARY);  // ��ֵ�ָ�
	Mat imagemat1 = imread(rawImage); // ������
	Mat se(3,3,CV_8U,Scalar(1)); 
	Mat opened, openerode, opened1, openerode1; 
	morphologyEx(imagemat, opened, MORPH_OPEN, se, Point(-1,-1), 4);  // ������
	//erode(opened, openerode, se, Point(-1,-1), 2);
	morphologyEx(imagemat1, opened1, MORPH_OPEN, se, Point(-1,-1), 3); 
	erode(opened1, openerode1, se, Point(-1,-1), 3);  // ��ʴ
	cutImg = &IplImage(opened);  
	cutImg1 = &IplImage(openerode1);
	//cvSaveImage("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/cutImg2D1.bmp", cutImg1,0);
	//MER("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/cutImg2D1.bmp", spanImg, "E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/cutImg1.bmp");
	//cvSaveImage("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/cutImg2D.bmp", cutImg,0);
	//MER("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/cutImg2D.bmp", spanImg, "E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/cutImg.bmp");

	int num = intBinToVector(rectnum, 1)[0]; // number of sitments ������Ŀ 
	float *countT = (float *)calloc(sizeof(float), 2); // countT, countUT Ŀ����Ŀ����Ŀ����Ŀ
	int height=src->height;   
	int width=src->width;
	gray=cvCreateImage(cvSize(src->width,src->height),src->depth,1);  // ��ֵͼ
	cvCvtColor(src,gray,CV_BGR2GRAY);

	// ��ֵ�ָ�
	//cutImg = cvCreateImage(cvSize(span->width, span->height),span->depth,1);  
	//cvThreshold(spancopy, spancopy, ostuValue, 255, CV_THRESH_BINARY);
	//cvCvtColor(spancopy,cutImg,CV_BGR2GRAY); 

    int contourCount = cvFindContours(gray,storage,&contour,sizeof(CvContour),CV_RETR_TREE,CV_CHAIN_APPROX_SIMPLE);  // ��ͨ����Ŀ
	int rectCount = 0; // MER��Ŀ
	float *rectArea = (float *)calloc(sizeof(float), contourCount);    // sava area of rects ��������
	float *contourArea = (float *)calloc(sizeof(float), contourCount);    // save data of contour ����ͨ�����
	int *rectL = (int *)calloc(sizeof(int), contourCount);  // length of rects ����γ�
	int *rectW = (int *)calloc(sizeof(int), contourCount);  // width of rects ����ο�

	// ��ÿ����ͨ����д���
	for(int m = 0;contour!=0 ;contour=contour->h_next, m++)  {
	CvBox2D rect=cvMinAreaRect2(contour,storage);  // building box
    CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rect, rect_pts0);
	contourArea[m] = fabs(cvContourArea(contour)); // ���
  
	//��ΪcvPolyLineҪ��㼯������������CvPoint**  
	//����Ҫ�� CvPoint2D32f �͵� rect_pts0 ת��Ϊ CvPoint �͵� rect_pts  
	//������һ����Ӧ��ָ�� *pt  
	int npts = 4,k=0;   
	CvPoint rect_pts[4], *pt = rect_pts;  
  
	//printf("��ͨ������С��Ӿ��ζ�������ֱ�Ϊ:\n");  
	for (int i=0; i<4; i++)  {  
		rect_pts[i]= cvPointFrom32f(rect_pts0[i]);
/*		printf("%d %d\n",rect_pts[i].x,rect_pts[i].y); */ 
		rectL[m]=(int)sqrt((pow((rect_pts[0].x-rect_pts[1].x),2)+pow((rect_pts[0].y-rect_pts[1].y),2)));  
		rectW[m]=(int)sqrt((pow((rect_pts[0].x-rect_pts[3].x),2)+pow((rect_pts[0].y-rect_pts[3].y),2)));  
        if(rectL[m] < rectW[m]) { 
             k = rectL[m];  
             rectL[m] = rectW[m];  
             rectW[m] = k;  
           } // if
		}  // for i [0, 3]
	rectArea[m] = (float)rectW[m] * rectL[m];

	// �������threshold�ľ��ο򻭳�����
	if(contourArea[m] >= thresholdContour) {
		if(catenum == 2) {  // û����
		if(category[rectCount] == 0) {
			// address of ROI
			addNum++;
			char *addcount = (char *)calloc(sizeof(char), 5);
			_itoa(addNum, addcount,  10);
			address = strcpy(address, addHead[0]);
			address = strcat(address, addcount);
			address = strcat(address, addTail);
			//printf("address = %s\n", address);
			free(addcount);
			// save ROI
			CvRect rectbound =  cvBoundingRect(contour, 0);
			saveBoundingRect(span, address, rectbound);
			RotatedRect rRectDraw = RotatedRect(Point2f(rect.center.x, rect.center.y), Size2f(rect.size.width, rect.size.height), rect.angle);
			CvRect rRect = rRectDraw.boundingRect();
			// draw rects on rawImage
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(255,0,0), 2); // Ŀ����-��ɫ��
			countT[0]++;
		} // if
		else if(category[rectCount] == 1) {
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,255,0), 2); // ��Ŀ����-��ɫ��
			countT[1]++;
		} // else if
		} // catenum = 2
		if(catenum == 3) { // û��
		if(category[rectCount] == 0) 
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(255,0,0), 2); // Ŀ����-��ɫ��
		else if(category[rectCount] == 1)
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,0,255), 2); // Ŀ����1-��ɫ��
		else if(category[rectCount] == 2)
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,255,0), 2); // ��Ŀ����1-��ɫ��
		} // catenum = 3
		if(catenum == 5) { // now
			if(category[rectCount] < catenum ) {
			// address of ROI
			addNum++; // roi ����Ŀ
			char *addcount = (char *)calloc(sizeof(char), 5);
			_itoa(addNum, addcount,  10);
			address = strcpy(address, addHead[0]);
			address = strcat(address, addcount);
			address = strcat(address, addTail);
			addHeadRect = strcpy(addHeadRect, addHead[1]); 
			addHeadRect = strcat(addHeadRect, addcount);  // ���ĵ�ַ
			//printf("address = %s\n", addHeadRect);
			//printf("address = %s\n", address);
			free(addcount);
			if(category[rectCount] != 3 && !(category[rectCount] == 1 && rectArea[m] < threDraw)) { // ����һ�ѳ�����Ŀ��
				RotatedRect rRectDraw;
				float minDis = 0;
				if(category[rectCount] == 2 && rectArea[m] < threCar) { 
					rRectDraw = enlargeRect(spanImg, rect, carMeanPath, meanNum, &minDis, Rwidth, Rheight); // һ���ֳ����һ������
					CvPoint2D32f rect_ptsen[4];   // ����
					cvBoxPoints(rRectDraw, rect_ptsen); 
					CvPoint rect_ptse[4], *pten = rect_ptse;  
					for (int ii=0; ii<4; ii++)   
						rect_ptse[ii]= cvPointFrom32f(rect_ptsen[ii]); 
					cvPolyLine(spancopy, &pten, &npts, 1, 1, CV_RGB(128,0,255), 2); // ������
				} // if
				else if(category[rectCount] == 4) { // �ɻ�
					if(rect.size.width < miniflight) // ��ֹ�ɻ�̫С
						rect.size.width = miniflight;
					if(rect.size.height < miniflight)
						rect.size.height = miniflight;
					rRectDraw = RotatedRect(Point2f(rect.center.x, rect.center.y), Size2f(rect.size.width, rect.size.height), rect.angle);
					}
				else 
					rRectDraw = RotatedRect(Point2f(rect.center.x, rect.center.y), Size2f(rect.size.width, rect.size.height), rect.angle);
				CvRect rRect = rRectDraw.boundingRect(); // ˮƽ��С����
				if(minDis <= threDis)
					saveBoundingRect(cspan, address, rRect); // ��roi��Ƭ
			} // category != 3
			if(category[rectCount] == 3) {  // һ�ѳ�
				CvRect rectbound =  cvBoundingRect(contour, 0);
				//saveBoundingRect(cspan, address, rectbound);
				RotatedRect rRectDraw = RotatedRect(Point2f(rect.center.x, rect.center.y), Size2f(rect.size.width, rect.size.height), rect.angle); 
				CvRect rRect = rRectDraw.boundingRect(); // �����
				//saveBoundingRect(cspan, address, rRect);
				if(rectArea[m] >= threDraw) { // many vehicles
					// draw if the result of draw little vehicles in catenum3:vehicles
					IplImage *draw = drawCars(spanImg, cutImg1, cutImg, rectbound, carMeanPath, vehiclesMean, meanNum, thresholdContour, addHeadRect, catenum, Rwidth, Rheight); // find cars
					//cvSaveImage(address, draw, 0);
					cvReleaseImage(&draw);
				}
				else {// only one vehilce in this rect
					if(isVehicles(spanImg, rect, carMeanPath, meanNum, threDis))
						saveBoundingRect(cspan, address, rRect);
				}
			} // category = 3
			} // category[]<5
		// �������
		if(category[rectCount] == 1 && rectArea[m] > threDraw)  
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(255,0,0), 2); // building-��ɫ��
		else if(category[rectCount] == 2)
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,0,255), 2); // vehicle-��ɫ��
		else if(category[rectCount] == 3) {
			if(rectArea[m] >= threDraw)
				cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,255,255), 2); // vehicles-ǳ��ɫ��
			else if(rectArea[m] < threDraw)
				cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(255,255, 0), 2); // vehicles-��ɫ��
		} // else if
		else if(category[rectCount] == 4)
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(255,0,255), 2); // aircraft-��ɫ��
		else if(category[rectCount] == 5)
			cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,255,0), 2); // other-��ɫ��
		} // if catenum == 5
				rectCount++;
		} // if contourArea[m] >= thresholdContour
	} // for m
	assert(rectCount == num);
	//cvSaveImage(savepath[0], span,0);
	//cvSaveImage(savepath[1], spancopy, 0);
	//savetxt(countT, saveNumber, catenum);
	//cvReleaseImage(&src);      
	//cvReleaseImage(&gray);
	src = NULL;
	span = NULL;
	cutImg = NULL;
	cutImg1 = NULL;
	spanmat.release();
	srcmat.release();
	opened.release();
	opened1.release();
	openerode.release();
	openerode1.release();
	imagemat.release();
	imagemat1.release();
	cvReleaseImage(&cspan);
	cvReleaseImage(&spancopy);
	cvReleaseMemStorage(&storage);
	free(countT);
	free(rectArea);
	free(contourArea);
	free(rectL);
	free(rectW);
	free(address);
	return;
} // drawRect

/*
	��һ��һƬ�������򻭳�һ��һ���ĳ�
	input��
		spanImg: �Ҷ�ͼԭͼ
		cut: ����Ŀ��������ͨ���ж�ΪһƬ������Ӿ���
		innerCut���ֿ��˵�cut
		rectB: cut���������span�����λ��
		meanCar: ����ģ��
		vehiclesMean: һ�ѳ���ģ��
		meanNum���õ���������Ŀ
		thresholdContour����ͨ����ֵ��С�ڴ˲�����
		addHead���洢·����ǰ�벿��
		catenum: ����Ŀ
		Rwidth, Rheight: ����Ĵ�С
	output: ����Ƭ����addHeadָ���ļ���
			drawRect��rectB��С�Ĵ�����ͼ
*/
IplImage *drawRect::drawCars(char *spanImg, IplImage *cut, IplImage *innerCut, CvRect rectB, char *meanCar, char *vehiclesMean, int meanNum, int thresholdContour, char *addHead, int catenum, int Rwidth, int Rheight)
{
	/* 
	   (carMin, carMax)->considered cars
	   (carMin, threCar)->part of a car -> enlargeCars 
	   (threCar, carThre)->single car
	   (carThre, carMax)->finding inner cars in contour
	   gray: gray scale of rectB part of span
	   cutcopy: cutImage, as large as spanImg, eroded building area 
	   cutR: rectB parts of cutcopy - eroded building area
	   spanR: rectB parts of span, backgrounf of saveBoundingRect
	   spanRcopy: same as spanR, drawing results of cars
	   innerCutR: rectB parts of innerCut - eroded threshold span(finding inner cars) 
	*/
	// feature used: graymean, GLCMappclust, GLCMcontrast, GLCMmaxProba, GLDVASM, GLDVentropy
	const int carMax = 40000, carMin = 200, carThre = 6000, threCar = 2000;
	const float minDis = 1.2; // if distance > minDis -> rect not car
	Mat spanmat;
	IplImage *gray, *span, *cutcopy, *cutR, *spanR, *spanRcopy, *innerCutR;  
    CvMemStorage* storage = cvCreateMemStorage(0);  // Creat new memory
    CvSeq* contour = 0;  
	int addNum = 0;  
	char *addTail = ".bmp";
	char *address = (char *)calloc(sizeof(char), 80); // address of ROI

	// mat won't cause memory leak
	spanmat = imread(spanImg);
	span = &IplImage(spanmat);
	//span = cvLoadImage(spanImg,1);   // span.bmp
	cutcopy = cvCloneImage(cut);
	cvSetImageROI(span, rectB);
	cvSetImageROI(cutcopy, rectB);
	cutR = cvCreateImage(cvSize(rectB.width, rectB.height),cutcopy->depth, 3);
	spanR = cvCreateImage(cvSize(rectB.width, rectB.height),span->depth, 3);
	cvCopy(span, spanR, 0);
	cvCopy(cutcopy, cutR, 0);
	spanRcopy = cvCloneImage(spanR);
	gray=cvCreateImage(cvSize(rectB.width,rectB.height),cutcopy->depth,1);
	cvCvtColor(cutR, gray, CV_RGB2GRAY);
	innerCutR = cvCreateImage(cvSize(rectB.width, rectB.height),innerCut->depth, innerCut->nChannels);
	cvSetImageROI(innerCut, rectB);
	cvCopy(innerCut, innerCutR, 0);
	cvResetImageROI(cutcopy);
	cvResetImageROI(span);
	cvResetImageROI(innerCut);
	
	int height=span->height;   
	int width=span->width;
	
    int contourCount = cvFindContours(gray, storage, &contour, sizeof(CvContour), CV_RETR_TREE,CV_CHAIN_APPROX_SIMPLE);  
	int rectCount = 0;
	float *contourArea = (float *)calloc(sizeof(float), contourCount);    // save data of contour ����ͨ�����

	// ��rectB��С�Ŀ��������ͨ����д���
	for(int m = 0;contour!=0 ;contour=contour->h_next, m++)  {
	CvBox2D rect=cvMinAreaRect2(contour,storage);  // building box
    CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rect, rect_pts0);
	contourArea[m] = fabs(cvContourArea(contour)); // ��ͨ�����
  
	//��ΪcvPolyLineҪ��㼯������������CvPoint**  
	//����Ҫ�� CvPoint2D32f �͵� rect_pts0 ת��Ϊ CvPoint �͵� rect_pts  
	//������һ����Ӧ��ָ�� *pt  
	int npts = 4,k=0;   
	CvPoint rect_pts[4], *pt = rect_pts;  
  
	//printf("��ͨ������С��Ӿ��ζ�������ֱ�Ϊ:\n");  
	for (int i=0; i<4; i++)  {  
		rect_pts[i]= cvPointFrom32f(rect_pts0[i]);
		}  // for i [0, 3]

	if(contourArea[m] >= thresholdContour) {  
		int rectArea = rect.size.height*rect.size.width;
		if(catenum == 5) {
			// address of ROI
			addNum++;
			char *addcount = (char *)calloc(sizeof(char), 5);
			_itoa(addNum, addcount,  10);
			address = strcpy(address, addHead);
			address = strcat(address, addcount);
			char *addin = (char *)calloc(sizeof(char), 80);
			addin = strcpy(addin, address);
			//printf("addin = %s\n", addin);
			address = strcat(address, addTail);
			//printf("address = %s\n", address);
			free(addcount);
			if(rectArea >= carMin && rectArea <= carMax) { // car
				// is rect vehicle? 
				// rectM change to spanImg plane
				RotatedRect rectM = RotatedRect(Point2f(rect.center.x + rectB.x, rect.center.y+rectB.y), Size2f(rect.size.width, rect.size.height), rect.angle);
				bool isCar = (isVehicles(spanImg, rectM, vehiclesMean, meanNum, minDis) || isVehicles(spanImg, rectM, meanCar, meanNum, minDis));  
				//CvRect rectbound =  cvBoundingRect(contour, 0);
				//saveBoundingRect(cspan, address, rectbound);
				RotatedRect rRectDraw = RotatedRect(Point2f(rect.center.x, rect.center.y), Size2f(rect.size.width, rect.size.height), rect.angle);
				if(rectArea > carThre && contourCount > 1 && isCar) { // cars
					innerRect(spanImg, innerCutR, rectB, rRectDraw.boundingRect(), meanCar, meanNum, addin,Rwidth, Rheight);
					//innerRect2(spanImg, innerCutR, cutR, rectB, rRectDraw.boundingRect(), meanCar, addin,Rwidth, Rheight);
					cvPolyLine(spanRcopy, &pt, &npts, 1, 1, CV_RGB(255,255,0), 2); 
				} // if > carThre
				else if(rectArea <= carThre && isCar) { // car
					CvRect rRect;
					if(rectArea > threCar) {
						rRect = rRectDraw.boundingRect();
						saveBoundingRect(spanR, address, rRect);
					}
					else { // part of car
						float minD = 0;
						RotatedRect rRectRo = enlargeRect(spanImg, rRectDraw, meanCar, meanNum, &minD, Rwidth, Rheight);
						if (minD <= minDis) {  // car
							rRect = rRectRo.boundingRect();
							saveBoundingRect(spanR, address, rRect);
							CvPoint2D32f vertices[4];   // ����*4
							cvBoxPoints( rRectRo, vertices);
							CvPoint rRect_pts[4], *rRectpt = rRect_pts;  
							for (int i=0; i<4; i++)  // ת�������ʽ
								rRect_pts[i]= cvPointFrom32f(vertices[i]);
							cvPolyLine(spanRcopy, &rRectpt, &npts, 1, 1, CV_RGB(0,0,128), 2); 
						} // if minD
					} // else
					cvPolyLine(spanRcopy, &pt, &npts, 1, 1, CV_RGB(0,255,0), 2); 
				} // else if 
				else { // not car
					CvRect rRect = rRectDraw.boundingRect();
					cvPolyLine(spanRcopy, &pt, &npts, 1, 1, CV_RGB(128,128,0), 2); 
				}
			} // if rectArea
		} // if catenum == 5
				rectCount++;
		} // if contourArea[m] >= thresholdContour
	} // for m  
	//cvNamedWindow("Rect",CV_WINDOW_AUTOSIZE); 
	//cvShowImage("Rect",spanRcopy);
	//cvWaitKey(1); 
	//cvDestroyWindow("Rect");
	cvReleaseImage(&gray);
	//cvReleaseImage(&span);
	span = NULL;
	spanmat.release();
	cvReleaseImage(&cutcopy);
	cvReleaseImage(&cutR);
	cvReleaseImage(&spanR);
	cvReleaseImage(&innerCutR);
	cvReleaseMemStorage(&storage);
	free(contourArea);
	free(address);
	return spanRcopy;
} // drawCar

/*
	һ���ֳ�->��������ǰλ����ת12*30��ȡ��������meanCar���Ǹ�λ��
	input: 
		spanImg��ԭͼ�Ҷ�ͼ
		rectCar: һ���ֳ�
		meanCar: ��������ֵ
		meanNum����ֵ������Ŀ
		minDis: �����ģ��ľ���
		Rwidth, Rheight: ����size
	output:
		return RotetedRect���͵�����
*/
RotatedRect drawRect::enlargeRect(char *spanImg, CvBox2D rectCar, char *meanCar, int meanNum, float *minDis, int Rwidth, int Rheight)
{
	/*
		gray span: as large as the whole image
		rectCar: small rects, a part of a car
		meanCar: feature vector
		Rwidth, Rheight: size fo model
		rotateNum: angle increase by 360/rotateNum
		minDis: result of minDistance, address
	*/
	// feature used: graymean, GLCMappclust, GLCMcontrast, GLCMmaxProba, GLDVASM, GLDVentropy
	Mat spanmat;
	IplImage *gray, *span;  
    CvMemStorage* storage = cvCreateMemStorage(0);  
    CvSeq *contour = 0;  
	/*const int Rwidth = 35, Rheight = 70;*/
	const int rotateNum = 12;   // ��ת����
	// < threPercent can be considered cars, more than precent of neighbouring rects are cars
	// can't devide them _(:�٩f��)_ i draw a rect hoping to include the target
	//float threPercent = 0.7, percent = 0.65;  

	spanmat = imread(spanImg);
	span = &IplImage(spanmat);
	//span = cvLoadImage(spanImg,1);   // span.bmp
	gray=cvCreateImage(cvSize(span->width,span->height),span->depth,1);  
    cvCvtColor(span,gray,CV_BGR2GRAY); 
	int height=span->height;   
	int width=span->width;

	// ������
	unsigned char **spanMatrix = BmpToMatrix(spanImg, width, height);    // data in span [0, 255]
	unsigned char **grayValue = BmpToMatrix(spanImg, width, height);  // gray value in spanImg.bmp
	changeLevel(spanMatrix, width, height);
	int **GLDVsta = intAllocateMemory(height, width);  // for calculating GLDV
	for(int p = 0; p < height - 1; p++)
		for(int q =0; q < width - 1; q++) {
			GLDVsta[p][q]=spanMatrix[p][q]-spanMatrix[p+1][q+1];
			if(GLDVsta[p][q] > 15)
				printf("error: GLDV > 15. \n");
		} // for 

	//float maxX = 0, maxY = 0, minX = 0, minY = 0;
	// ��������
	float grayMean = 0;
	float GLDVASM = 0;
	float GLDVentropy = 0;
	float GLDVcontrast = 0;
	float GLDVmean = 0;
	float GLCMentropy = 0;
	float GLCMcontrast = 0;
	float GLCMappclust = 0;
	float GLCMmaxProba = 0;
	// ����Ŀ����� graymean, GLCMappclust, GLCMcontrast, GLCMmaxProba, GLDVASM, GLDVentropy
	float *feature = (float *)calloc(sizeof(float), 6);
	// ����ο���������
	int angle = 0; // ��ת�Ƕ�
	float *distance = (float *)calloc(sizeof(float), rotateNum);
	// distance between current small rect and result rect
	float dif = sqrt(Rwidth*Rwidth + Rheight*Rheight)*0.5 - 0.4*(rectCar.size.height+rectCar.size.width);
	//float dif = sqrt(Rwidth*Rwidth + Rheight*Rheight)*0.5 ;
	float **GLCMatrix = fAllocateMemory(GrayLayerNum, GrayLayerNum);    // GLCM
	float *GLDVector = (float *)calloc(sizeof(float), 2*GrayLayerNum - 1);  // GLDV
	int dog = 0;
	for(angle = 0; angle < 360; angle += (360/rotateNum)) {
		float pointx = rectCar.center.x + dif*cos((float)angle*pi/180);
		float pointy = rectCar.center.y + dif*sin((float)angle*pi/180);
		//pointx = rectCar.center.x ;
		//pointy = rectCar.center.y ;
		RotatedRect rRect = RotatedRect(Point2f(pointx, pointy), Size2f(Rheight, Rwidth),angle);
		Point2f vertices[4];  // ����*4
		rRect.points(vertices); // ��ȡ����
		CvPoint rRect_pts[4], *rRectpt = rRect_pts;  
		for (int i=0; i<4; i++)  // ת�������ʽ
			rRect_pts[i]= vertices[i];
			// Extra textural features in rects ��ȡ��������������

		int minxx = 0, minyy = 0, maxxx = 0, maxyy = 0;
		float k1,k2,k3,k4,b1,b2,b3,b4;
	    if(angle >= 0 && angle < 90) { // С����������
		minxx = rRect_pts[0].x; 
		minyy = rRect_pts[1].y;
		maxxx = rRect_pts[2].x;
		maxyy = rRect_pts[3].y;
		k1 = findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);     
		b1 = findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k2 = findk(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		b2 = findb(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		k3 = findk(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		b3 = findb(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		k4 = findk(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		b4 = findb(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		} // if angle
		else if(angle >= 90 && angle < 180) { // С����������
		minxx = rRect_pts[3].x; 
		minyy = rRect_pts[0].y;
		maxxx = rRect_pts[1].x;
		maxyy = rRect_pts[2].y;
		k1 = findk(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);   
		b1 = findb(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		k2 = findk(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		b2 = findb(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		k3 = findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		b3 = findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k4 = findk(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		b4 = findb(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		} // if angle
		else if(angle >= 180 && angle < 270) { // С����������
		minxx = rRect_pts[2].x; 
		minyy = rRect_pts[3].y;
		maxxx = rRect_pts[0].x;
		maxyy = rRect_pts[1].y;
		k1 = findk(rRect_pts[1].x, rRect_pts[2].x, rRect_pts[1].y, rRect_pts[2].y);   
		b1 = findb(rRect_pts[1].x, rRect_pts[2].x, rRect_pts[1].y, rRect_pts[2].y);
		k2 = findk(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		b2 = findb(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		k3 = findk(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		b3 = findb(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		k4 = findk(rRect_pts[3].x, rRect_pts[0].x, rRect_pts[3].y, rRect_pts[0].y);
		b4 = findb(rRect_pts[3].x, rRect_pts[0].x, rRect_pts[3].y, rRect_pts[0].y);
		} // if angle
		else if(angle >= 270 && angle < 360) { // С����������
		minxx = rRect_pts[1].x; 
		minyy = rRect_pts[2].y;
		maxxx = rRect_pts[3].x;
		maxyy = rRect_pts[0].y;
		k1 = findk(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);  
		b1 = findb(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		k2 = findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		b2 = findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k3 = findk(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		b3 = findb(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		k4 = findk(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		b4 = findb(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		} // if angle
		int countR = 0; // С������ĵ�����
		int xx =0, yy = 0;
		if(maxyy >= height)
			maxyy = height-1;
		if(maxxx >= width)
			maxxx = width-1;
		if(minxx < 0)
			minxx = 0;
		if(minyy < 0)
			minyy =0;
		for( xx = minxx; xx <= maxxx; xx++)
		for( yy = minyy; yy <= maxyy; yy++) {    
		if(yy <= (k1*xx + b1) && yy <= (k2*xx + b2)
		      && yy >= (k3*xx + b3) && yy >= (k4*xx + b4)
			  && xx >= 0 && xx < width
			  && yy >= 0 && yy < height) {
	countR++;
	// calculate GLCM ����GLCM����
	calculateGraylevelmatrix(spanMatrix, GLCMatrix, xx, yy, width, height);
	// calculate GLDV ����GLDVʸ��
	calculateGLDV(GLDVsta, GLDVector, xx, yy);	
	// calculae gray value ����ҶȾ�ֵ
	grayMean += grayValue[yy][xx];
		} // if xx yy belongs to rect
	} // for xx yy 
	//assert( grayMean>0);
		//if(grayMean <= 0)
		//	system("pause");
	// normalization of GLCM GLCM�����һ��
	float sumGLCM = 0;
	for(int mm = 0; mm < GrayLayerNum; mm++)
		for(int n =0; n < GrayLayerNum; n++)
			sumGLCM += GLCMatrix[mm][n];
	for(int mm = 0; mm < GrayLayerNum; mm++)
		for(int n = 0; n < GrayLayerNum; n++)
			GLCMatrix[mm][n] = GLCMatrix[mm][n]/sumGLCM;
	//for(int k = 0; k < GrayLayerNum; k++)
	//	for(int f =0; f<GrayLayerNum;f++)
	//		printf("%f\n", GLCMatrix[k][f]);
	// calculate textural feature of GLCM ����GLCM����
	float *ui = (float *)calloc(sizeof(float), GrayLayerNum);
	float *uj = (float *)calloc(sizeof(float), GrayLayerNum);
	float contrasttemp=0;
	float entropytemp=0;
	float appclustemp=0;
	float maxProbatemp=0;
	for(int gi = 0 ; gi < GrayLayerNum; gi++) {  
		for( int gj = 0; gj < GrayLayerNum; gj++)
			{
				uj[gi]+=GLCMatrix[gi][gj];
			}
				uj[gi]=uj[gi]/GrayLayerNum;
		} // for
		for( int gj = 0 ; gj < GrayLayerNum; gj++)
			{    
				for( int gi = 0; gi < GrayLayerNum; gi++)
				{
					ui[gj]+=GLCMatrix[gi][gj];
				}
			ui[gj]=ui[gj]/GrayLayerNum;
			} // for
	for( int k = 0 ; k < GrayLayerNum; k++)
		{    
			for( int z = 0; z< GrayLayerNum; z++)
			{
				contrasttemp += GLCMatrix[k][z]*abs((k-z)*(k-z));

				if(GLCMatrix[k][z]<1&&GLCMatrix[k][z])
				entropytemp -= GLCMatrix[k][z]*log10(GLCMatrix[k][z]);
				if(GLCMatrix[k][z]-maxProbatemp>=0)
					maxProbatemp=GLCMatrix[k][z];				
			appclustemp += (k+z-ui[k]-uj[z])*(k+z-ui[k]-uj[z])*GLCMatrix[k][z];
				}
		} // for
	GLCMcontrast=contrasttemp;			
	GLCMentropy=entropytemp;
	GLCMappclust=appclustemp;
	GLCMmaxProba=maxProbatemp;
	free(ui);
	free(uj);

	// normalization of GLDV GLDVʸ����һ��
	float sumGLDV = 0;
	for(int i= 0; i < 2*GrayLayerNum - 1; i++) 
		sumGLDV += GLDVector[i];
	for(int i= 0; i < 2*GrayLayerNum - 1; i++) 
		GLDVector[i] = GLDVector[i]/sumGLDV;
	// calculate textural feature of GLDV ����GLDV����
	float GLDVcontrasttemp=0;
	float GLDVentropytemp=0;
	float GLDVASMtemp=0;
	float GLDVmeantemp=0;
	for (int k = 0; k < (2*GrayLayerNum-1); k++) {
		GLDVcontrasttemp = GLDVcontrasttemp+k*k*GLDVector[k];
		if(GLDVector[k])
			GLDVentropytemp-=GLDVector[k]*log(GLDVector[k]);
		GLDVASMtemp+=GLDVector[k]*GLDVector[k];
		GLDVmeantemp+=k*GLDVector[k];
		} // for
	GLDVcontrast=GLDVcontrasttemp;
	GLDVentropy=GLDVentropytemp;
	GLDVASM=GLDVASMtemp;
	GLDVmean=GLDVmeantemp;
	if(countR == 0)
		countR = 1;
	// calculate gray mean ����Ҷ�ֵ��ֵ
	grayMean = grayMean/countR;
	//����ŷ�Ͼ���
	feature[0] = grayMean;
	feature[1] = GLCMappclust;
	feature[2] = GLCMcontrast;
	feature[3] = GLCMmaxProba;
	feature[4] = GLDVASM;
	feature[5] = GLDVentropy;

	distance[dog] = euclidean(feature, BinToVector(meanCar, meanNum), meanNum); 
	dog++;
	countR = 0;
	// draw
	CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rRect, rect_pts0);
	int npts = 4;   
	CvPoint rect_pts[4], *pt = rect_pts;  
	for (int ii=0; ii<4; ii++)  {  
		rect_pts[ii]= cvPointFrom32f(rect_pts0[ii]);
		}  // if ii
	cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,0,255), 1); //  ÿ����
	} // if angle
	for(int cat = 0; cat < rotateNum; cat++) {
		if(distance[cat] <= 0)
			distance[cat] = findMax(distance, rotateNum)+1; // ���ǳ���
	} // if cat
	int min = findMin(distance, rotateNum); // ��С����

	if(rectCar.size.height*rectCar.size.width >= Rwidth*Rheight) 
		*minDis = 1000;  // cannot be target
	*minDis = distance[min];

	float pointx = rectCar.center.x + dif*cos((float)min*2*pi/rotateNum);  // ���������x
	float pointy = rectCar.center.y + dif*sin((float)min*2*pi/rotateNum);  // ���������y
	// result rect
	RotatedRect enRect = RotatedRect(Point2f(pointx, pointy), Size2f(Rheight, Rwidth), min*360/rotateNum); // ���
	//int similiar = 0; // num of rects considered cars
	//for(int i = 0; i < rotateNum; i++) {
	//	if(distance[i] <= threPercent)
	//		similiar++;
	//} // for i
	//if(similiar > (float)rotateNum * percent)
	//	enRect = RotatedRect(Point2f(rectCar.center.x, rectCar.center.y), Size2f((float)(Rheight+Rwidth)*0.5, (float)(Rwidth+Rheight)*0.5), 0);
	CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(enRect, rect_pts0);
	int npts = 4;   
	CvPoint rect_pts[4], *pt = rect_pts;  
	for (int ii=0; ii<4; ii++)  {  
		rect_pts[ii]= cvPointFrom32f(rect_pts0[ii]);
		}  // if ii
	cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,255,0), 2);  // ���������
	CvPoint2D32f arect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rectCar, arect_pts0);
	CvPoint arect_pts[4], *apt = arect_pts;  
	for (int ii=0; ii<4; ii++)  {  
		arect_pts[ii]= cvPointFrom32f(arect_pts0[ii]);
		}  // if ii
	cvPolyLine(span, &apt, &npts, 1, 1, CV_RGB(255,0,0), 2); // ����С��
	if(*minDis > 1.3) 
		cvPolyLine(span, &apt, &npts, 1, 1, CV_RGB(255,0,255), 2); 
	//cvNamedWindow("Rect",CV_WINDOW_AUTOSIZE); 
	//cvShowImage("Rect",span);
	//cvWaitKey(0);
	//cvDestroyWindow("Rect");
	//cvSaveImage("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/rectenlarge.bmp", span);
	//system("pause");
	// �ͷ��ڴ�
	free(distance); 
	free(feature);
	cvReleaseImage(&gray);
	//cvReleaseImage(&span);
	span = NULL;
	spanmat.release();
	cvReleaseMemStorage(&storage);
	UcFreeMemory(spanMatrix, height);
	UcFreeMemory(grayValue, height);
	intFreeMemory(GLDVsta, height); 
	fFreeMemory(GLCMatrix, GrayLayerNum);
	free(GLDVector);
	return enRect;
} // enlargeRect


/*
	��һ�ѳ��ֿ�
	input:
		spanImg: �Ҷ�ͼԭͼ
		cut: ����Ŀ��������ͨ��
		rectB: ��ͨ��MER��Сλ��
		rectCar: rectB������һ��MER
		meanNum: ģ��������Ŀ
		addHead: ��roi��ַǰ��
		Rwidth��Rheight: ���Ĵ�С
	output: ����һ��roi
*/
void drawRect::innerRect(char *spanImg, IplImage *cut, CvRect rectB, CvRect rectCar, char *meanCar, int meanNum, char *addHead, int Rwidth, int Rheight)
{
	/*
		rectB: the entire contour: connected cars
		rectCar: one connected car contour in rectB
		cut: small contours to divide cars that connected, as large as rectB
		gray: find contours, same as cutR
		cutcopy: = cut
		spanWhole: spanImg, whole image
		span: rectB part of spanImg
		cutR, spanR: middle files, as large as currentRect
		threDis: if dis>threDis, consider it a wrong target
	*/
	Mat spanWholemat;
	IplImage *gray, *cutcopy, *spanWhole, *cutR, *spanR, *span, *anotherSpan;  
    CvMemStorage* storage = cvCreateMemStorage(0);  // Creat new memory
    CvSeq* contour = 0;  
	int addNum = 0;  
	char *addTail = ".bmp";
	char *address = (char *)calloc(sizeof(char), 80); // address of ROI
	const float threDis = 1.3;
	int inCount = 0; // count of cars in currentRect
	float dis = 0; // minDistance of meanCar and newRect

	spanWholemat = imread(spanImg);
	spanWhole = &IplImage(spanWholemat);
	span = cvCreateImage(cvSize(rectB.width, rectB.height),spanWhole->depth, spanWhole->nChannels);
	anotherSpan = cvCreateImage(cvSize(spanWhole->width, spanWhole->height),spanWhole->depth, spanWhole->nChannels);
	cvCopy(spanWhole, anotherSpan);
	cvSetImageROI(spanWhole, rectB);
	cvCopy(spanWhole, span);
	cutcopy = cvCloneImage(cut);

	float pointx = 0, pointy = 0; // current rect coordinate
	CvRect currentRect; // current rect on entire plane
	if(rectCar.x < 0) {
		rectCar.width = rectCar.width+rectCar.x;
		rectCar.x = 0;
	} // if rectCar.x < 0
	if(rectCar.y < 0) {
		rectCar.height = rectCar.height+rectCar.y;
		rectCar.y = 0;
	} // if rectCar.y < 0
	if(rectB.width <= rectCar.x+rectCar.width-1 || rectB.height <= rectCar.y+rectCar.height-1) {   // rectB �� rectCar ����
		cutR = cvCreateImage(cvSize(rectB.width, rectB.height),cutcopy->depth, cutcopy->nChannels);
		spanR = cvCreateImage(cvSize(rectB.width, rectB.height),span->depth, span->nChannels);
		// span, spanR, cutR, cutcopy = rectB
		cvCopy(cutcopy, cutR);
		cvCopy(span, spanR);
		//gray=cvCreateImage(cvSize(rectB.width,rectB.height),cutcopy->depth,1);
		pointx = rectB.x;
		pointy = rectB.y;
		currentRect = rectB;
	} // if small rectB
	else { // rectCar �� rectB����ϵ����
	cvSetImageROI(span, rectCar);
	cvSetImageROI(cutcopy, rectCar);
	cutR = cvCreateImage(cvSize(rectCar.width, rectCar.height),cutcopy->depth, cutcopy->nChannels);
	spanR = cvCreateImage(cvSize(rectCar.width, rectCar.height),span->depth, span->nChannels);
	// span, spanR, cutR, cutcopy on rectB background
	cvCopy(cutcopy, cutR);
	cvCopy(span, spanR);
	//gray=cvCreateImage(cvSize(rectCar.width,rectCar.height),cutcopy->depth,1);
	pointx = rectCar.x+rectB.x;
	pointy = rectCar.y+rectB.y;
	currentRect = cvRect(pointx, pointy, rectCar.width, rectCar.height);
	} // else large rectB
	gray = cvCloneImage(cutR);
	cvResetImageROI(cutcopy);
	cvResetImageROI(span);
	cvResetImageROI(spanWhole);
	cvRectangle(anotherSpan, Point(pointx, pointy), Point(pointx+currentRect.width, pointy+currentRect.height),CV_RGB(255,0,0),2);
	int height=span->height;   
	int width=span->width;
	
	// �Ե�ǰ��rectCar����
    int contourCount = cvFindContours(gray, storage, &contour, sizeof(CvContour), CV_RETR_TREE,CV_CHAIN_APPROX_SIMPLE);  
	for(int m = 0;contour!=0 ;contour=contour->h_next, m++)  {
			// rect on currentRect plane
			CvBox2D rect=cvMinAreaRect2(contour,storage);  
			// address of ROI
			addNum++;
			char *addcount = (char *)calloc(sizeof(char), 5);
			_itoa(addNum, addcount,  10);
			address = strcpy(address, addHead);
			address = strcat(address, addcount);
			address = strcat(address, addTail);
			//printf("address = %s\n", address);
			free(addcount);
			// newRect- spanImg 
			// change rect frtom currentRect to spanImg
			RotatedRect newRect = RotatedRect(Point2f(currentRect.x+rect.center.x, currentRect.y+rect.center.y), Size2f(rect.size.width, rect.size.height),rect.angle);  // һ���ֳ�
			// rRectDraw - spanImg, resule of enlarge
			 
			RotatedRect rRectDraw = enlargeRect(spanImg, newRect, meanCar, meanNum, &dis, Rwidth, Rheight); // ����
			CvRect rRect = rRectDraw.boundingRect();
			CvPoint2D32f rect_pts0[4];   // find 4 vertexs
			cvBoxPoints(rRectDraw, rect_pts0);
			int npts = 4;   
			CvPoint rect_pts[4], *pt = rect_pts;  
			for (int i=0; i<4; i++)  {  
				rect_pts[i]= cvPointFrom32f(rect_pts0[i]);   // ����
			}  // for i [0, 3]			
			if(rRectDraw.center.x < currentRect.x + currentRect.width && rRectDraw.center.y < currentRect.y + currentRect.height) {
				inCount++;
				cvPolyLine(anotherSpan, &pt, &npts, 1, 1, CV_RGB(0,0,255), 2); 
			// rRect - spanImg
			if(dis <= threDis) // ������roi
				saveBoundingRect(spanWhole, address, rRect);
			} // if center in currentRect
			if (inCount == 0) {
				cvPolyLine(anotherSpan, &pt, &npts, 1, 1, CV_RGB(0,255,255), 2); 
				saveBoundingRect(spanWhole, address, currentRect); // ������
			}
			// draw rRect
			//cvRectangle(spanWhole, Point(rRect.x, rRect.y), Point(rRect.x+rRect.width, rRect.y+rRect.height),CV_RGB(0,255,0),2);
	} // for m
	//cvSaveImage("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/innerRect.bmp", anotherSpan);
	cvReleaseImage(&gray);
	cvReleaseImage(&cutcopy);
	spanWhole = NULL;
	spanWholemat.release();
	cvReleaseImage(&cutR);
	cvReleaseImage(&spanR);
	cvReleaseImage(&span);
	cvReleaseImage(&anotherSpan);
	cvReleaseMemStorage(&storage);
	free(address);
	return;
} // innerRect

// �ж�rectCar�ǲ��ǳ�
bool drawRect::isVehicles(char *spanImg, CvBox2D rectCar, char *meanVehicles, int meanNum, float minDis)
{
	// feature used: graymean, GLCMappclust, GLCMcontrast, GLCMmaxProba, GLDVASM, GLDVentropy
	/*
		gray span: as large as the whole image
		rectCar: rect wait to be examed
		meanCar: feature vector of vehicle
		minDis: if distance{to meanCar} > minDis->not car, else is
	*/
	Mat spanmat;
	IplImage *gray, *span;  
    CvMemStorage* storage = cvCreateMemStorage(0);  
    CvSeq *contour = 0;  
	spanmat = imread(spanImg);
	span = &IplImage(spanmat);
	//span = cvLoadImage(spanImg,1);   // span.bmp
	gray=cvCreateImage(cvSize(span->width,span->height),span->depth,1);  
    cvCvtColor(span,gray,CV_BGR2GRAY); 
	int height=span->height;   
	int width=span->width;
	bool result = true;

	unsigned char **spanMatrix = BmpToMatrix(spanImg, width, height);    // data in span [0, 255]
	unsigned char **grayValue = BmpToMatrix(spanImg, width, height);  // gray value in spanImg.bmp
	changeLevel(spanMatrix, width, height);
	int **GLDVsta = intAllocateMemory(height, width);  // for calculating GLDV
	for(int p = 0; p < height - 1; p++)
		for(int q =0; q < width - 1; q++) {
			GLDVsta[p][q]=spanMatrix[p][q]-spanMatrix[p+1][q+1];
			if(GLDVsta[p][q] > 15)
				printf("error: GLDV > 15. \n");
		} // for 

	// ��������
	float grayMean = 0;
	float GLDVASM = 0;
	float GLDVentropy = 0;
	float GLDVcontrast = 0;
	float GLDVmean = 0;
	float GLCMentropy = 0;
	float GLCMcontrast = 0;
	float GLCMappclust = 0;
	float GLCMmaxProba = 0;
	// ����Ŀ����� graymean, GLCMappclust, GLCMcontrast, GLCMmaxProba, GLDVASM, GLDVentropy
	float *feature = (float *)calloc(sizeof(float), meanNum);
	// ����ο���������
	float **GLCMatrix = fAllocateMemory(GrayLayerNum, GrayLayerNum);    // GLCM
	float *GLDVector = (float *)calloc(sizeof(float), 2*GrayLayerNum - 1);  // GLDV
	float pointx = rectCar.center.x ;
	float pointy = rectCar.center.y ;
	CvPoint2D32f vertices[4];   // ����*4
	cvBoxPoints(rectCar, vertices);
	CvPoint rRect_pts[4], *rRectpt = rRect_pts;  
	for (int i=0; i<4; i++)  // ת�������ʽ
		rRect_pts[i]= cvPointFrom32f(vertices[i]);
	// Extra textural features in rects ��ȡ��������������
	int angle = rectCar.angle;
	int minxx = 0, minyy = 0, maxxx = 0, maxyy = 0;
	float k1,k2,k3,k4,b1,b2,b3,b4;
		minxx = rRect_pts[1].x; 
		minyy = rRect_pts[2].y;
		maxxx = rRect_pts[3].x;
		maxyy = rRect_pts[0].y;
		k1 = findk(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);     
		b1 = findb(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		k2 = findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		b2 = findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k3 = findk(rRect_pts[1].x, rRect_pts[2].x, rRect_pts[1].y, rRect_pts[2].y);
		b3 = findb(rRect_pts[1].x, rRect_pts[2].x, rRect_pts[1].y, rRect_pts[2].y);
		k4 = findk(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		b4 = findb(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);

		int countR = 0; // С������ĵ�����
		int xx =0, yy = 0;
		if(maxyy >= height)
			maxyy = height-1;
		if(maxxx >= width)
			maxxx = width-1;
		if(minxx < 0)
			minxx = 0;
		if(minyy < 0)
			minyy =0;
		for( xx = minxx; xx <= maxxx; xx++)
		for( yy = minyy; yy <= maxyy; yy++) {    
		if(yy <= (k1*xx + b1) && yy <= (k2*xx + b2)
		      && yy >= (k3*xx + b3) && yy >= (k4*xx + b4)
			  && xx >= 0 && xx < width
			  && yy >= 0 && yy < height) {
	countR++;
	// calculate GLCM ����GLCM����
	calculateGraylevelmatrix(spanMatrix, GLCMatrix, xx, yy, width, height);
	// calculate GLDV ����GLDVʸ��
	calculateGLDV(GLDVsta, GLDVector, xx, yy);	
	// calculae gray value ����ҶȾ�ֵ
	grayMean += grayValue[yy][xx];
		} // if xx yy belongs to rect
	} // for xx yy 
	assert( grayMean>0);
	// normalization of GLCM GLCM�����һ��
	float sumGLCM = 0;
	for(int mm = 0; mm < GrayLayerNum; mm++)
		for(int n =0; n < GrayLayerNum; n++)
			sumGLCM += GLCMatrix[mm][n];
	for(int mm = 0; mm < GrayLayerNum; mm++)
		for(int n = 0; n < GrayLayerNum; n++)
			GLCMatrix[mm][n] = GLCMatrix[mm][n]/sumGLCM;
	//for(int k = 0; k < GrayLayerNum; k++)
	//	for(int f =0; f<GrayLayerNum;f++)
	//		printf("%f\n", GLCMatrix[k][f]);
	// calculate textural feature of GLCM ����GLCM����
	float *ui = (float *)calloc(sizeof(float), GrayLayerNum);
	float *uj = (float *)calloc(sizeof(float), GrayLayerNum);
	float contrasttemp=0;
	float entropytemp=0;
	float appclustemp=0;
	float maxProbatemp=0;
	for(int gi = 0 ; gi < GrayLayerNum; gi++) {  
		for( int gj = 0; gj < GrayLayerNum; gj++)
			{
				uj[gi]+=GLCMatrix[gi][gj];
			}
				uj[gi]=uj[gi]/GrayLayerNum;
		} // for
		for( int gj = 0 ; gj < GrayLayerNum; gj++)
			{    
				for( int gi = 0; gi < GrayLayerNum; gi++)
				{
					ui[gj]+=GLCMatrix[gi][gj];
				}
			ui[gj]=ui[gj]/GrayLayerNum;
			} // for
	for( int k = 0 ; k < GrayLayerNum; k++)
		{    
			for( int z = 0; z< GrayLayerNum; z++)
			{
				contrasttemp += GLCMatrix[k][z]*abs((k-z)*(k-z));

				if(GLCMatrix[k][z]<1&&GLCMatrix[k][z])
				entropytemp -= GLCMatrix[k][z]*log10(GLCMatrix[k][z]);
				if(GLCMatrix[k][z]-maxProbatemp>=0)
					maxProbatemp=GLCMatrix[k][z];				
			appclustemp += (k+z-ui[k]-uj[z])*(k+z-ui[k]-uj[z])*GLCMatrix[k][z];
				}
		} // for
	GLCMcontrast=contrasttemp;			
	GLCMentropy=entropytemp;
	GLCMappclust=appclustemp;
	GLCMmaxProba=maxProbatemp;
	free(ui);
	free(uj);

	// normalization of GLDV GLDVʸ����һ��
	float sumGLDV = 0;
	for(int i= 0; i < 2*GrayLayerNum - 1; i++) 
		sumGLDV += GLDVector[i];
	for(int i= 0; i < 2*GrayLayerNum - 1; i++) 
		GLDVector[i] = GLDVector[i]/sumGLDV;
	// calculate textural feature of GLDV ����GLDV����
	float GLDVcontrasttemp=0;
	float GLDVentropytemp=0;
	float GLDVASMtemp=0;
	float GLDVmeantemp=0;
	for (int k = 0; k < (2*GrayLayerNum-1); k++) {
		GLDVcontrasttemp = GLDVcontrasttemp+k*k*GLDVector[k];
		if(GLDVector[k])
			GLDVentropytemp-=GLDVector[k]*log(GLDVector[k]);
		GLDVASMtemp+=GLDVector[k]*GLDVector[k];
		GLDVmeantemp+=k*GLDVector[k];
		} // for
	GLDVcontrast=GLDVcontrasttemp;
	GLDVentropy=GLDVentropytemp;
	GLDVASM=GLDVASMtemp;
	GLDVmean=GLDVmeantemp;
	if(countR == 0)
		countR = 1;
	// calculate gray mean ����Ҷ�ֵ��ֵ
	grayMean = grayMean/countR;
	//����ŷ�Ͼ���
	feature[0] = grayMean;
	feature[1] = GLCMappclust;
	feature[2] = GLCMcontrast;
	feature[3] = GLCMmaxProba;
	feature[4] = GLDVASM;
	feature[5] = GLDVentropy;

	float distance = euclidean(feature, BinToVector(meanVehicles, meanNum), meanNum); 
	if(distance > minDis)
		result = false;
	// �ͷ��ڴ�
	free(feature);
	cvReleaseImage(&gray);
	span = NULL;
	spanmat.release();
	cvReleaseMemStorage(&storage);
	UcFreeMemory(spanMatrix, height);
	UcFreeMemory(grayValue, height);
	intFreeMemory(GLDVsta, height); 
	fFreeMemory(GLCMatrix, GrayLayerNum);
	free(GLDVector);
	return result;
} // isVehicles

/*
	��һ�ѳ��ֿ�2...��û���ϣ���ӦenlargeCar2
	input:
		spanImg: �Ҷ�ͼԭͼ
		cut: ����Ŀ��������ͨ��
		rectB: ��ͨ��MER��Сλ��
		rectCar: rectB������һ��MER
		meanNum: ģ��������Ŀ
		addHead: ��roi��ַǰ��
		Rwidth��Rheight: ���Ĵ�С
	output: ����һ��roi
*/
void drawRect::innerRect2(char *spanImg, IplImage *cut, IplImage *Tcut, CvRect rectB, CvRect rectCar, char *meanCar,char *addHead, int Rwidth, int Rheight)
{
	/*
		rectB: the entire contour: connected cars
		rectCar: one connected car contour in rectB
		cut: small contours to divide cars that connected, as large as rectB
		innerCut: rectB size, points
		gray: find contours, same as cutR
		cutcopy: = cut
		spanWhole: spanImg, whole image
		span: rectB part of spanImg
		cutR, spanR: middle files, as large as currentRect
		threDis: if dis>threDis, consider it a wrong target
	*/
	Mat spanWholemat;
	IplImage *gray, *cutcopy, *spanWhole, *cutR, *spanR, *span, *cutT, *Tcutcopy, *anotherSpan;  
    CvMemStorage* storage = cvCreateMemStorage(0);  // Creat new memory
    CvSeq* contour = 0;  
	int addNum = 0;  
	char *addTail = ".bmp";
	char *address = (char *)calloc(sizeof(char), 80); // address of ROI
	float threDis = 1.3;
	//float threDis = findMax(meanMin, catenum); //
	////change here

	int inCount = 0; // count of cars in currentRect
	float dis = 0; // minDistance of meanCar and newRect

	spanWholemat = imread(spanImg);
	spanWhole = &IplImage(spanWholemat);
	span = cvCreateImage(cvSize(rectB.width, rectB.height),spanWhole->depth, spanWhole->nChannels);
	anotherSpan = cvCreateImage(cvSize(spanWhole->width, spanWhole->height),spanWhole->depth, spanWhole->nChannels);
	cvCopy(spanWhole, anotherSpan);
	cvSetImageROI(spanWhole, rectB);
	cvCopy(spanWhole, span);
	cutcopy = cvCloneImage(cut);
	Tcutcopy = cvCloneImage(Tcut);
	float pointx = 0, pointy = 0; // current rect coordinate
	CvRect currentRect; // current rect on entire plane
	if(rectCar.x < 0) {
		rectCar.width = rectCar.width+rectCar.x;
		rectCar.x = 0;
	} // if rectCar.x < 0
	if(rectCar.y < 0) {
		rectCar.height = rectCar.height+rectCar.y;
		rectCar.y = 0;
	} // if rectCar.y < 0
	if(rectB.width <= rectCar.x+rectCar.width-1 || rectB.height <= rectCar.y+rectCar.height-1) {
		cutR = cvCreateImage(cvSize(rectB.width, rectB.height),cutcopy->depth, cutcopy->nChannels);
		cutT = cvCreateImage(cvSize(rectB.width, rectB.height),Tcut->depth, Tcut->nChannels);
		spanR = cvCreateImage(cvSize(rectB.width, rectB.height),span->depth, span->nChannels);
		// span, spanR, cutR, cutcopy = rectB
		cvCopy(cutcopy, cutR);
		cvCopy(Tcutcopy, cutT);
		cvCopy(span, spanR);
		//gray=cvCreateImage(cvSize(rectB.width,rectB.height),cutcopy->depth,1);
		pointx = rectB.x;
		pointy = rectB.y;
		currentRect = rectB;
	} // if small rectB
	else {
	cvSetImageROI(span, rectCar);
	cvSetImageROI(cutcopy, rectCar);
	cvSetImageROI(Tcutcopy, rectCar);
	cutR = cvCreateImage(cvSize(rectCar.width, rectCar.height),cutcopy->depth, cutcopy->nChannels);
	cutT = cvCreateImage(cvSize(rectCar.width, rectCar.height),Tcut->depth, Tcut->nChannels);
	spanR = cvCreateImage(cvSize(rectCar.width, rectCar.height),span->depth, span->nChannels);
	// span, spanR, cutR, cutcopy on rectB background
	cvCopy(cutcopy, cutR);
	cvCopy(Tcutcopy, cutT);
	cvCopy(span, spanR);
	//gray=cvCreateImage(cvSize(rectCar.width,rectCar.height),cutcopy->depth,1);
	pointx = rectCar.x+rectB.x;
	pointy = rectCar.y+rectB.y;
	currentRect = cvRect(pointx, pointy, rectCar.width, rectCar.height);
	} // else large rectB
	gray = cvCloneImage(cutR);
	cvResetImageROI(cutcopy);
	cvResetImageROI(span);
	cvResetImageROI(spanWhole);
	cvResetImageROI(Tcutcopy);
	cvRectangle(anotherSpan, Point(pointx, pointy), Point(pointx+currentRect.width, pointy+currentRect.height),CV_RGB(255,0,0),2);
	int height=span->height;   
	int width=span->width;
	
    int contourCount = cvFindContours(gray, storage, &contour, sizeof(CvContour), CV_RETR_TREE,CV_CHAIN_APPROX_SIMPLE);  
	for(int m = 0;contour!=0 ;contour=contour->h_next, m++)  {
			// rect on currentRect plane
			CvBox2D rect=cvMinAreaRect2(contour,storage);  
			// address of ROI
			addNum++;
			char *addcount = (char *)calloc(sizeof(char), 5);
			_itoa(addNum, addcount,  10);
			address = strcpy(address, addHead);
			address = strcat(address, addcount);
			address = strcat(address, addTail);
			//printf("address = %s\n", address);
			free(addcount);
			// newRect- spanImg 
			// change rect frtom currentRect to spanImg
			RotatedRect newRect = RotatedRect(Point2f(currentRect.x+rect.center.x, currentRect.y+rect.center.y), Size2f(rect.size.width, rect.size.height),rect.angle);
			// rRectDraw - spanImg, resule of enlarge
			IplImage *cutRectB = cvCreateImage(cvSize(currentRect.width, currentRect.height),cut->depth, cutT->nChannels);
			IplImage *innerCutRectB = cvCreateImage(cvSize(currentRect.width, currentRect.height),cut->depth, cut->nChannels);
			//cvCvtColor(cutT,cutRectB,CV_BGR2GRAY); 
			cvCopy(cutR, innerCutRectB);
			cvCopy(cutT, cutRectB);
			RotatedRect rRectDraw = enlargeRect2(spanImg, cutRectB, innerCutRectB, currentRect, rect, 45, 80);
			cvReleaseImage(&cutRectB);
			cvReleaseImage(&innerCutRectB);

			CvRect rRect = rRectDraw.boundingRect();
			CvPoint2D32f rect_pts0[4];   // find 4 vertexs
			cvBoxPoints(rRectDraw, rect_pts0);
			int npts = 4;   
			CvPoint rect_pts[4], *pt = rect_pts;  
			for (int i=0; i<4; i++)  {  
				rect_pts[i]= cvPointFrom32f(rect_pts0[i]);
			}  // for i [0, 3]			
			if(rRectDraw.center.x < currentRect.x + currentRect.width && rRectDraw.center.y < currentRect.y + currentRect.height) {
				inCount++;
				cvPolyLine(anotherSpan, &pt, &npts, 1, 1, CV_RGB(0,0,255), 2); 
			// rRect - spanImg
			if(dis <= threDis)
				saveBoundingRect(spanWhole, address, rRect);
			} // if center in currentRect
			if (inCount == 0) {
				cvPolyLine(anotherSpan, &pt, &npts, 1, 1, CV_RGB(0,255,255), 2); 
				saveBoundingRect(spanWhole, address, currentRect);
			}
			// draw rRect
			cvRectangle(anotherSpan, Point(rRect.x, rRect.y), Point(rRect.x+rRect.width, rRect.y+rRect.height),CV_RGB(0,255,0),2);
	} // for m
	//cvSaveImage("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/innerRect.bmp",anotherSpan);
	cvReleaseImage(&gray);
	cvReleaseImage(&cutcopy);
	spanWhole = NULL;
	spanWholemat.release();
	cvReleaseImage(&cutR);
	cvReleaseImage(&spanR);
	cvReleaseImage(&span);
	cvReleaseImage(&cutT);
	cvReleaseImage(&Tcutcopy);
	cvReleaseImage(&anotherSpan);
	cvReleaseMemStorage(&storage);
	free(address);
	return;
} // innerRect2

// allocate memory for w*h char matrix
unsigned char** drawRect::cmatrix(int w,int h)
{
	unsigned char **m=new unsigned char*[h];
	for(int i=0;i<h;i++)
		m[i]=new unsigned char[w];
	return m;
}

// ���ֳ�->����2 ��ӦinnerRect2 ����û����
RotatedRect drawRect::enlargeRect2(char *spanImg, IplImage *cut, IplImage *innerCut, CvRect currentRect, CvBox2D rect, int Rwidth, int Rheight)
{
	/*
		span: as large as the whole image
		rectCar: contours' bounding, on currentRect plane 
		cut: target cut, on currentRect plane
		innerCut: points cut, on currentRect plane
		rotateNum
	*/
	Mat spanmat;
	IplImage *gray, *span;  
    CvMemStorage* storage = cvCreateMemStorage(0);  
    CvSeq *contour = 0;  
	//const int Rwidth = 35, Rheight = 70;
	const int rotateNum = 12;  
	// > threPercent can be considered cars, more than precent of neighbouring rects are cars
	// can't devide them _(:�٩f��)_ i draw a rect hoping to include the target
	// float threPercent = 0.6, percent = 0.65;  
	assert(cut->width == innerCut->width && cut->height == innerCut->height);

	spanmat = imread(spanImg);
	span = &IplImage(spanmat);
	gray=cvCreateImage(cvSize(innerCut->width,innerCut->height),innerCut->depth,1);  
	cvCopy(innerCut, gray);
	int height=innerCut->height;   
	int width=innerCut->width;

	// matrix of 2 cut
	unsigned char **cutmatrix =cmatrix(cut->width, cut->height);
	unsigned char **innerCutmatrix =cmatrix(innerCut->width, innerCut->height);
	for (int i = 0; i < cut->height; i++)
	{
		for (int j = 0; j < cut->width; j++)
		{
			CvScalar s = cvGet2D(cut,i,j);  
			CvScalar ss = cvGet2D(innerCut,i,j); 
			cutmatrix[i][j]=s.val[0];
			innerCutmatrix[i][j]=ss.val[0];
		} // for
	} // for

	// fit of each angle
	float *fit = (float *)calloc(sizeof(float), rotateNum);

	float dif = sqrt(Rwidth*Rwidth + Rheight*Rheight)*0.5 - 0.4*(rect.size.height+rect.size.width);
	int dog = 0, angle = 0;
	int countCut = 0, countInnerCut = 0;
	for(angle = 0; angle < 360; angle += (360/rotateNum)) {
		float pointx = rect.center.x + dif*cos((float)angle*pi/180);
		float pointy = rect.center.y + dif*sin((float)angle*pi/180);
		//pointx = rectCar.center.x ;
		//pointy = rectCar.center.y ;
		RotatedRect rRect = RotatedRect(Point2f(pointx, pointy), Size2f(Rheight, Rwidth),angle);
		Point2f vertices[4];  // ����*4
		rRect.points(vertices); // ��ȡ����
		CvPoint rRect_pts[4], *rRectpt = rRect_pts;  
		for (int i=0; i<4; i++)  // ת�������ʽ
			rRect_pts[i]= vertices[i];
			
		// fit each rect
		int minxx = 0, minyy = 0, maxxx = 0, maxyy = 0;
		float k1,k2,k3,k4,b1,b2,b3,b4;
	    if(angle >= 0 && angle < 90) { // С����������
		minxx = rRect_pts[0].x; 
		minyy = rRect_pts[1].y;
		maxxx = rRect_pts[2].x;
		maxyy = rRect_pts[3].y;
		k1 =findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);     
		b1 =findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k2 =findk(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		b2 =findb(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		k3 =findk(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		b3 =findb(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		k4 =findk(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		b4 =findb(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		} // if angle
		else if(angle >= 90 && angle < 180) { // С����������
		minxx = rRect_pts[3].x; 
		minyy = rRect_pts[0].y;
		maxxx = rRect_pts[1].x;
		maxyy = rRect_pts[2].y;
		k1 =findk(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);   
		b1 =findb(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		k2 =findk(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		b2 =findb(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		k3 =findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		b3 =findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k4 =findk(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		b4 =findb(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		} // if angle
		else if(angle >= 180 && angle < 270) { // С����������
		minxx = rRect_pts[2].x; 
		minyy = rRect_pts[3].y;
		maxxx = rRect_pts[0].x;
		maxyy = rRect_pts[1].y;
		k1 =findk(rRect_pts[1].x, rRect_pts[2].x, rRect_pts[1].y, rRect_pts[2].y);   
		b1 =findb(rRect_pts[1].x, rRect_pts[2].x, rRect_pts[1].y, rRect_pts[2].y);
		k2 =findk(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		b2 =findb(rRect_pts[1].x, rRect_pts[0].x, rRect_pts[1].y, rRect_pts[0].y);
		k3 =findk(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		b3 =findb(rRect_pts[3].x, rRect_pts[2].x, rRect_pts[3].y, rRect_pts[2].y);
		k4 =findk(rRect_pts[3].x, rRect_pts[0].x, rRect_pts[3].y, rRect_pts[0].y);
		b4 =findb(rRect_pts[3].x, rRect_pts[0].x, rRect_pts[3].y, rRect_pts[0].y);
		} // if angle
		else if(angle >= 270 && angle < 360) { // С����������
		minxx = rRect_pts[1].x; 
		minyy = rRect_pts[2].y;
		maxxx = rRect_pts[3].x;
		maxyy = rRect_pts[0].y;
		k1 =findk(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);  
		b1 =findb(rRect_pts[0].x, rRect_pts[1].x, rRect_pts[0].y, rRect_pts[1].y);
		k2 =findk(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		b2 =findb(rRect_pts[0].x, rRect_pts[3].x, rRect_pts[0].y, rRect_pts[3].y);
		k3 =findk(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		b3 =findb(rRect_pts[2].x, rRect_pts[1].x, rRect_pts[2].y, rRect_pts[1].y);
		k4 =findk(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		b4 =findb(rRect_pts[2].x, rRect_pts[3].x, rRect_pts[2].y, rRect_pts[3].y);
		} // if angle

		int xx =0, yy = 0;
		if(maxyy >= height)
			maxyy = height-1;
		if(maxxx >= width)
			maxxx = width-1;
		if(minxx < 0)
			minxx = 0;
		if(minyy < 0)
			minyy =0;
		for( xx = minxx; xx <= maxxx; xx++)
		for( yy = minyy; yy <= maxyy; yy++) {    
		if(yy <= (k1*xx + b1) && yy <= (k2*xx + b2)
		      && yy >= (k3*xx + b3) && yy >= (k4*xx + b4)
			  && xx >= 0 && xx < width
			  && yy >= 0 && yy < height) {
			if(cutmatrix[yy][xx] > 128)
				countCut++;
			if(innerCutmatrix[yy][xx] > 128)
				countInnerCut++;
		} // if xx yy belongs to rect
	} // for xx yy 
		fit[dog] = (float)countCut/(Rheight*Rwidth);
	dog++;
	// draw
	RotatedRect rRecton = RotatedRect(Point2f(pointx+currentRect.x, pointy+currentRect.y), Size2f(Rheight, Rwidth),angle);
	CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rRecton, rect_pts0);
	int npts = 4;   
	CvPoint rect_pts[4], *pt = rect_pts;  
	for (int ii=0; ii<4; ii++)  {  
		rect_pts[ii]= cvPointFrom32f(rect_pts0[ii]);
		}  // if ii
	cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,0,255), 1); 
	} // if angle
	// largest && in contour
	int max =findMaxinA(fit, rotateNum);
	float pointx = rect.center.x + dif*cos((float)max*2*pi/rotateNum);
	float pointy = rect.center.y + dif*sin((float)max*2*pi/rotateNum);
	// result rect
	RotatedRect enRect = RotatedRect(Point2f(pointx+currentRect.x, pointy+currentRect.y), Size2f(Rheight, Rwidth), max*360/rotateNum);
	//int similiar = 0; // num of rects considered cars
	//for(int i = 0; i < rotateNum; i++) {
	//	if(fit[i] >= threPercent)
	//		similiar++;
	//} // for i
	//if(similiar > (float)rotateNum * percent)
	//	enRect = RotatedRect(Point2f(rect.center.x+currentRect.x, rect.center.y+currentRect.y), Size2f((float)(Rheight+Rwidth)*0.5, (float)(Rwidth+Rheight)*0.5), 0);
	CvPoint2D32f rect_pts0[4];   // find 4 vertexs
    cvBoxPoints(enRect, rect_pts0);
	int npts = 4;   
	CvPoint rect_pts[4], *pt = rect_pts;  
	for (int ii=0; ii<4; ii++)  {  
		rect_pts[ii]= cvPointFrom32f(rect_pts0[ii]);
		}  // if ii
	cvPolyLine(span, &pt, &npts, 1, 1, CV_RGB(0,255,0), 2); 
	RotatedRect rectCaron = RotatedRect(Point2f(currentRect.x, currentRect.y), Size2f(currentRect.width, currentRect.height),0);
	CvPoint2D32f arect_pts0[4];   // find 4 vertexs
    cvBoxPoints(rectCaron, arect_pts0);
	CvPoint arect_pts[4], *apt = arect_pts;  
	for (int ii=0; ii<4; ii++)  {  
		arect_pts[ii]= cvPointFrom32f(arect_pts0[ii]);
		}  // if ii
	//cvSaveImage("E:/PlayWithC/GradDesign/GeoTrait/GeoTrait/miniSAR/mini3/bdv/rectenlarge.bmp", span);
	// �ͷ��ڴ�
	free(fit); 
	cvReleaseImage(&gray);
	span = NULL;
	spanmat.release();
	cvReleaseMemStorage(&storage);
	UcFreeMemory(cutmatrix, height);
	UcFreeMemory(innerCutmatrix, height);
	return enRect;
} // enlargeRect

// mean minDistance of specimen and mean vector
float *drawRect::meanMinDistance(char **countSave, char **traitpath, char *cateArea, int catenum, int traitNum)
{
	int echo = 20; 
	float *meanMin = (float *)calloc(sizeof(float), catenum); 
	int num = intBinToVector(countSave[0], 1)[0]; // number of contours ��ͨ����Ŀ  
	 printf("num = %d\n",num);
	 int rectnum = intBinToVector(countSave[1], 1)[0]; // number of speciments ������Ŀ
	 printf("rectnum = %d\n", rectnum);
	 int *smallMER = intBinToVector(cateArea, num); // ��0-С����
	 float **sample = fAllocateMemory(traitNum, num);  // ��������
	 for(int j = 0; j < traitNum; j++) {
		 sample[j] = BinToVector(traitpath[j], num);
	 } // if j
	 float **centre = findCenter(traitpath, traitNum, num, catenum, cateArea);  // center for categray  
	 int *EucliCategory = (int *)calloc(sizeof(int), rectnum);   // ŷ�Ͼ��������
	 float *distanceE = (float *)calloc(sizeof(float), catenum); // ŷ�Ͼ���
	 float **tempSample = fAllocateMemory(rectnum, traitNum);
	 int aa = 0, nn =0;
	 for(int a = 0; a < num; a++) {
		 if(smallMER[a] != 0) { 
			 for(nn = 0; nn < traitNum; nn++) {
				 tempSample[aa][nn] = sample[nn][a];  // ��������
			 } // for nn
		 aa++;		
		 } // if
	 } // for a 
	 assert(aa == rectnum);
	 for(int b = 0; b < rectnum; b++) {  // ���� [0, catenum-1]
		for(int m = 0; m < catenum; m++) {
			distanceE[m] = euclidean(tempSample[b], centre[m], traitNum);
				} // for m
		EucliCategory[b] = findMinNotZero(distanceE, catenum) + 1;   // the address of min value
	 } // for b
	 for(int cat = 0; cat<catenum; cat++) {
		 int iterate = 0;
		 for(int e = 0; e < rectnum; e++) {
			if( EucliCategory[e] == cat+1) { // +1!!!!!!!!!!!
				meanMin[cat] += euclidean(tempSample[e], centre[cat], traitNum);
				iterate++;
			} // if [cat]
			if(iterate >= echo)
				break;
		} // for e
	 } // for cat
	 for(int cc = 0; cc < catenum; cc++)
		 meanMin[cc] = meanMin[cc]/echo;
	 free((void *)smallMER);
	 fFreeMemory(sample, traitNum);
	 fFreeMemory(centre, catenum);
	 fFreeMemory(tempSample, rectnum);
	 free(EucliCategory);
	 free(distanceE);
	 return meanMin;
} // meanMinDistance

// function
unsigned char **drawRect:: BmpToMatrix(char *filename, int width, int height)
{
	IplImage* src = cvLoadImage(filename,0);    // input image
	int width1 = src->width;    // Ncol
	int height1 = src->height;  // Nrow
	assert(width == width1 && height == height1);
	unsigned char **matrix = cmatrix(width, height);
	for (int i = 0; i < height; i++)
	{
		for (int j = 0; j < width; j++)
		{
			CvScalar s = cvGet2D(src,i,j);         
			matrix[i][j]=s.val[0];
		} // for
	} // for
	cvReleaseImage(&src);
	return matrix;
} // BmpToMatrix

float *drawRect::BinToVector(char *filename, int length)
{
	FILE  *inFile;
	inFile = fopen(filename, "rb");
	if(inFile == NULL){
		printf("Failed to open the file.\n");
		exit(1);
	}
	float *tempMatrix = (float *)calloc(sizeof(float), length);
	fread(tempMatrix, sizeof(float), length, inFile);
	fclose(inFile);
	return tempMatrix;
}// BinToVector

int *drawRect::intBinToVector(char *filename, int length)
{
	FILE  *inFile;
	inFile = fopen(filename, "rb");
	if(inFile == NULL){
		printf("Failed to open the file.\n");
		exit(1);
	}
	int *tempMatrix = (int *)calloc(sizeof(int), length);
	fread(tempMatrix, sizeof(int), length, inFile);
	fclose(inFile);
	return tempMatrix;
}// intBinToVector

float **drawRect::fAllocateMemory(int m, int n)
{
	int k = 0;
	float **matrix = (float **)calloc(sizeof(float *), m);                          
	for(k = 0; k < m; k++){
		matrix[k] = (float *)calloc(sizeof(float), n);
	}	
	return matrix;
} // AllocateMemory float

 int **drawRect::  intAllocateMemory(int m, int n)
{
	int k = 0;
	int **matrix = (int **)calloc(sizeof(int *), m);                          
	for(k = 0; k < m; k++){
		matrix[k] = (int *)calloc(sizeof(int), n);
	}	
	return matrix;
} // AllocateMemory int

char **drawRect::charAllocateMemory(int m, int n)
{
	int k = 0;
	char **matrix = (char **)calloc(sizeof(char *), m);                          
	for(k = 0; k < m; k++){
		matrix[k] = (char *)calloc(sizeof(char), n);
	}	
	return matrix;
} // AllocateMemory char

void drawRect::fFreeMemory(float **matrix, int m)
{
	int k=0;
	for (k = 0; k < m; k++)           
		free((void *)matrix[k]);
	free(matrix);
	return;
} // FreeMemory

void drawRect::UcFreeMemory(unsigned char **matrix, int m)
{
	int k=0;
	for (k = 0; k < m; k++)           
		free((void *)matrix[k]);
	free(matrix);
	return;
} // UcFreeMemory

void drawRect::intFreeMemory(int **matrix, int m)
{
	int k=0;
	for (k = 0; k < m; k++)           
		free((void *)matrix[k]);
	free(matrix);
	return;
} // intFreeMemory

void drawRect::charFreeMemory(char **matrix, int m)
{
	int k=0;
	for (k = 0; k < m; k++)           
		free((void *)matrix[k]);
	free(matrix);
	return;
} // charFreeMemory

void drawRect::fNormalize(float *arrayA, int arraysize)
{
	float min = arrayA[findMin(arrayA, arraysize)];
	float max = findMax(arrayA, arraysize);
	if(max != min) {
	for(int i = 0; i < arraysize; i++) 
		arrayA[i] = (arrayA[i] - min)/(max - min);
	} // if
	return;
} // fNormalize

//  ��image�ϵ�rect��address
void drawRect::saveBoundingRect(IplImage* image, char *address, CvRect rect)
{
	IplImage *image2 = cvCloneImage(image);
	if(rect.x <= 0 || rect.y <= 0)
		return;
	if(rect.x >= image->width)
		rect.x = image->width-1;
	if(rect.y >= image->height)
		rect.y = image->height-1;
	if(rect.x + rect.width >= image->width)
		rect.width = image->width - rect.x-1;
	if(rect.y + rect.height >= image->height)
		rect.height = image->height - rect.y;
	cvSetImageROI(image2, rect);
	//cvNamedWindow("boundingRect",CV_WINDOW_AUTOSIZE);  
	//cvShowImage("boundingRect",image2);
	cvSaveImage(address, image2, 0);
	cvReleaseImage(&image2);
} // saveeBoundingRect



// BDistance
/*
	��whichCate����traitPath�����ľ�ֵʸ��
	input:
		maskPath: ��ͨ���Ƿ������ĸ�����
		traitPath������
		contourCount�� ��ͨ����Ŀ
		whichcate: ��һ��
	output: return ��ֵ
*/
float drawRect::cateMean(char *maskPath, char* traitPath, int contourCount, int whichCate)
{
	float mean = 0;
	int count = 0;
	float *trait = BinToVector(traitPath, contourCount);
	int *mask = intBinToVector(maskPath, contourCount);
	for(int i=0; i< contourCount; i++) {
		if(mask[i] == whichCate) {
			mean += trait[i];
		count++;
		} // if
		if(count > 20)
			break;
	} // for
	//printf("target mean = %f\n", mean/count);
	return mean/count;
} // targetMean

// array�����ֵ
float drawRect::findMax(float *array, int arraysize)
{
	float max = 0;
	for(int i = 0; i< arraysize; i++) {
		if(max < array[i])
			max = array[i];
	} // for i
	//printf("max = %f\n", max);
	return max;
} // findMax

// ����array ��lastmin����С����λ��
int drawRect::findminnum(float *array, int lastmin, int arraysize)
{
	int minnum = 0;
	float min = findMax(array, arraysize);
	for(int i = 0; i < arraysize; i++) {
		if(i != lastmin && array[i] >= array[lastmin]) {
			if(min > array[i]) {
				min = array[i];
				minnum = i;
			} // if
		} // if
	} // for
	return minnum;
} // findmin

// minDist
/* 
	���ؾ�ֵʸ��
	input: 
		traitpath: ����
		dim: ����ά��
		num: ������Ŀ
		maskpath: �������ķ���
		catenum: ��Ŀ
	output������ÿ����Ŀ��������ֵ����
*/
float **drawRect::findCenter(char** traitpath, int dim, int num, int catenum, char *maskpath) // ���ؾ�ֵ����
{
	float **mean = fAllocateMemory(catenum, dim);
	for(int i = 0; i < catenum; i++)
		for(int j = 0; j < dim; j++) 
			mean[i][j] = findMean(maskpath, traitpath[j], num, i+1); // i+1!!
	return mean;
} // findCentre

// ��whichcate�࣬ ��contourCount���ݣ�����traitpath�����ľ�ֵ��������������maskPath
float drawRect::findMean(char *maskPath, char* traitPath, int contourCount, int whichCate)
{
	float mean = 0;
	int count = 0;
	float *trait = BinToVector(traitPath, contourCount);
	int *mask = intBinToVector(maskPath, contourCount);
	for(int i=0; i< contourCount; i++) {
		if(mask[i] == whichCate) {
			mean += trait[i];
		count++;
		} // if
		if(count > 20)
			break;
	} // for
	//printf("target mean = %f\n", mean/count);
	return mean/count;
} // findMean

// ����arraysize��С��vector����Сֵ����λ��
int drawRect::findMin(float *vector, int arraysize)
{
	int minnum = 0;
	float min = vector[0];
	for(int i = 0; i < arraysize; i++) {
			if(min > vector[i]) {
				min = vector[i];
				minnum = i;
		} // if
	} // for
	return minnum;
} // findMin

// ����arraysize��С��vector�����ֵ����λ��
int drawRect::findMaxinA(float *vector, int arraysize)
{
	int maxnum = 0;
	float max = vector[0];
	for(int i = 0; i < arraysize; i++) {
			if(max < vector[i]) {
				max = vector[i];
				maxnum = i;
		} // if
	} // for
	return maxnum;
} // findMax

// ��sample��mean����dimά������ŷ�Ͼ���
float drawRect::euclidean(float *sample, float *mean, int dim) // ŷ�Ͼ���
{
	float result = 0;
	for(int i = 0; i < dim; i++) {
		//assert(sample[i] >= 0 && mean[i] >= 0);
		result +=(sample[i] - mean[i])*(sample[i] - mean[i])/mean[i]/mean[i];
	} // for i
	result = sqrt(result);
	//printf("ED = %f\n", result);
	return result;
} // euclidean

// rect1
// ��length����float* data ��name.bin
void drawRect::savefloat(float *data, char* name, int length)  
{
	FILE *fp = NULL;
	fp = fopen(name,"wb");
	fwrite(data,sizeof(float),length,fp);
	//printf()
	assert(fp != NULL);
	fclose(fp);
	fp = NULL;
} // savefloat

// ��length����int* data ��name.bin
void drawRect::saveint(int *data, char* name, int length)  
{
	FILE *fp = NULL;
	fp = fopen(name,"wb");
	fwrite(data,sizeof(int),length,fp);
	assert(fp != NULL);
	fclose(fp);
	fp = NULL;
} // saveint

// ��length����float* data ��name.txt
void drawRect::savetxt(float *data, char* name, int length)  
{
	FILE *fp = NULL;
	fp = fopen(name,"w");
	for (int i = 0; i < length; i++)
		fprintf(fp,"%f ",data[i]);
	assert(fp != NULL);
	fclose(fp);
	fp = NULL;
} // savetxt

// ����(x1,y1)(x2,y2)��б��
float drawRect::findk(int x1, int x2, int y1, int y2)
{
	if(y1 != y2 && x1 != x2)
		return (float)(y1 - y2)/(x1 - x2);
	else
		return 0;
} // findk

// ����(x1,y1)(x2,y2)�Ľؾ�
float drawRect::findb(int x1, int x2, int y1, int y2)
{
	if(findk(x1, x2, y1, y2) != 0)
		return (float)y2 - findk(x1, x2, y1, y2)*x2;
	else
		return (float)y1;
} // findb

// GLCM 
void drawRect::calculateGraylevelmatrix(unsigned char **span, float **GLCMatrix, int x, int y, int width, int height)
{
	int dis = 1;   // ����
	//printf("x = %d, y = %d\n", x,y);
	//����0�ȻҶȹ�������
	if(x + dis < width && y < height && x >= 0 && y >= 0) {
	GLCMatrix[span[y][x]][span[y][x+dis]] += 1;
	GLCMatrix[span[y][x+dis]][span[y][x]] += 1;	
	}
	//����90�ȵĻҶȹ�����
	if(x < width && y + dis < height && x >= 0 && y >= 0) { 
	GLCMatrix[span[y][x]][span[y+dis][x]] += 1;
	GLCMatrix[span[y+dis][x]][span[y][x]] += 1;
	}
	//����135�ȵĻҶȹ�����
	int newx = x + dis;
	int newy = y + dis;
	if(newx < width && newy < height && x >=0 && y >= 0) {
	GLCMatrix[span[y][x]][span[newy][newx]] += 1;
	GLCMatrix[span[newy][newx]][span[y][x]] += 1;
	}
	//����45�ȵĻҶȹ�����
	newx = x - dis;
	newy = y + dis;
	if(x < width && newy < height && newx >=0 && y >= 0) {
	GLCMatrix[span[y][x]][span[newy][newx]] += 1;
	GLCMatrix[span[newy][newx]][span[y][x]] += 1;
	}
} // Graylevelmatrix

// 256->16
void drawRect::changeLevel(unsigned char **rawImage, int width, int height)
{
	for(int i = 0; i < height; i++)
		for(int j = 0; j < width; j++) {
			rawImage[i][j] = rawImage[i][j]/(256/GrayLayerNum);
			assert(rawImage[i][j] < 16 && rawImage[i][j] >=0);
		}
}// changeLevel

// GLDV
void drawRect::calculateGLDV(int **GLDVsta, float *GLDVector, int x, int y)
{
	//assert(y < 1638 && x < 2510);
	assert((GLDVsta[y][x]+GrayLayerNum-1) < 2*GrayLayerNum-1);
	GLDVector[(GLDVsta[y][x]+GrayLayerNum-1)] += 1;
} // calculateGLDV

// ����arraysize���arrayA�з�����Сֵ����λ��
int drawRect::findMinNotZero(float *arrayA, int arraysize)
{
	float min = 0;
	int minnum = 0;
	for(int i = 0; i < arraysize; i ++) {
		if(arrayA[i] > 0) {
			min = arrayA[i];
			minnum = i;
			break;
		} // if
	} // for i
	for(int j = minnum; j < arraysize; j++) {
		if(arrayA[j] > 0 && arrayA[j] < min ) {
			min = arrayA[j];
			minnum = j;
		} // if
	} // for j
	//assert(min != 0);
	return minnum;
} // findMinNotZero

// ����cate1/cate2���trait������Ӧ���Ͼ��룬������contourC, �Ƿ����¼��mask
float drawRect::BDistanceCate(char *mask, char *trait, char *contourC, int cate1, int cate2)
{
	float bd = 0; // BDistance 
	int contourCount = intBinToVector(contourC, 1)[0];
	float mean1 = cateMean(mask, trait, contourCount, cate1);  
	float mean2 = cateMean(mask, trait, contourCount, cate2);  
	float devi1 = cateDev(mask, trait, contourCount, cate1); 
	float devi2 = cateDev(mask, trait, contourCount, cate2);  
	bd = (0.25)*(mean1-mean2)*(mean1-mean2)/(devi1*devi1+devi2*devi2)+(0.5)*log10((devi1*devi1+devi2*devi2)/(2*devi1*devi2));
	//printf("bd = %f\n", bd);
	return bd;
} // BDistance

// ������ͨ����cate��traitpath�����ı�׼�contourCount������maskPath��¼�Ƿ���
float drawRect::cateDev(char *maskPath, char *traitPath, int contourCount, int cate)
{
	float deviation = 0;
	int N = 0;
	float *trait = BinToVector(traitPath, contourCount);
	int *mask = intBinToVector(maskPath, contourCount);
	for(int i=0; i< contourCount; i++) {
		if(mask[i] == cate) {
		N++;
		} // if
		if(N > 20)
			break;
	} // for
	float mean = cateMean(maskPath, traitPath, contourCount, cate);
	for(int j = 0; j < N; j++) {
		deviation += (trait[j] - mean)*(trait[j] - mean);
	} // for
	deviation = sqrt(deviation/N);
	//printf("target deviation = %f\n", deviation);
	return deviation;
} // targetDev 

// _(:�٩f��)__(:�٩f��)__(:�٩f��)_
void drawRect::OnPath() 
{
	//���ü�����·��
	// TODO: Add your control notification handler code here
	//��ͨ�öԻ���BROWSEINFO�ṹ�а������û�ѡ��Ŀ¼����Ҫ��Ϣ
	BROWSEINFO browse;
	ZeroMemory(&browse,sizeof(browse));//fills a block of memory with zeros.
	browse.hwndOwner = NULL;
	browse.pszDisplayName = m_strPath.GetBuffer(MAX_PATH);
	browse.lpszTitle = "��ѡ��һ��ͼ��Ŀ¼";
	//SHBrowseForFolder��������һ��ITEMIDLIST�ṹ��ָ�룬�������û�ѡ���ļ��е���Ϣ
	LPITEMIDLIST lpItem = SHBrowseForFolder(&browse);
	if(lpItem == NULL) return ;
	m_strPath.ReleaseBuffer();
	//SHGetPathFromIDList����Ŀ��־���б�ת��Ϊ�ĵ�ϵͳ·��
	if(SHGetPathFromIDList(lpItem,m_strPath.GetBuffer(MAX_PATH)) == false) return;
	m_strPath.ReleaseBuffer();
	dir=true; //��־λ����Ϊtrue����ʾ������ͼ��������
	AfxMessageBox("��ѡ���Ŀ¼Ϊ:"+m_strPath,MB_ICONINFORMATION|MB_OK);	
	//ɨ������⣬�õ�ͼ��Ŀ¼���ļ���·��
	CString tempath;
	CString temps;
	tempath=m_strPath;
	tempath.TrimRight();tempath.TrimLeft(); //ȥ��ǰ�����
	CString strfilepath=tempath;
	tempi=0;
	counts=0;//����������
	//��������ͼ���������counts�У���·������temp[300]��
	StartDir(strfilepath);
    temps.Format("��Ŀ¼�¹���%d��ͼ��!",counts);
    AfxMessageBox(temps,MB_ICONINFORMATION|MB_OK);
}

// ɨ��������е�ͼ�����
void drawRect::StartDir(const CString &strfile1)
{
	BOOL yesno; 
	CFileFind find;	
	char tempFileFind[300]; 
	sprintf_s(tempFileFind,"%s\\*.*",strfile1);	
	RunDir(strfile1);//�ڵ�ǰĿ¼�в���ͼ�񣬲�������Ŀ¼
	yesno = (BOOL)find.FindFile(tempFileFind); 
	//���ҵ�ǰĿ¼����Ŀ¼�е��ļ�
	while(yesno) 
	{ 
		yesno = find.FindNextFile(); 
		if (find.IsDots() != TRUE)//����ȱʡĿ¼
		{
			char foundFileName[200];
			strcpy_s(foundFileName,find.GetFileName().GetBuffer(200));			
			if((find.IsDirectory() == TRUE)) //�ж��Ƿ�ΪĿ¼
			{ 
				char tempDir[200];
				sprintf_s(tempDir,"%s\\%s",strfile1,foundFileName);//�����Ŀ¼·��
				StartDir(tempDir); 	// �ݹ���ã�������Ŀ¼��ͼ��
			} 			
		}
	}   
	find.Close(); 
	return; 
}

// ɨ�赱ǰĿ¼�е�ͼ�����
void drawRect::RunDir(const CString &strfile2)
{
	BOOL yesno; 	
	CFileFind find; 
	char tempFileFind[200]; 
	sprintf_s(tempFileFind,"%s\\*.bmp",strfile2); 
	yesno = find.FindFile(tempFileFind); //�ڵ�ǰĿ¼�²���BMP�ļ�
	while(yesno) 
	{ 
		yesno = find.FindNextFile(); 
		char foundFileName[300];//��ʱ�洢���ҵ���ͼ����
		strcpy_s(foundFileName,find.GetFileName().GetBuffer(200));//��ȡͼ����
		if(!find.IsDots()) //����ȱʡĿ¼
		{ 
			char tempFileName[200];
			sprintf_s(tempFileName,"%s\\%s",strfile2,foundFileName);
            CString strfilepath1;
			strfilepath1.Format("%s",tempFileName);//��ȡͼ������·��
			counts++;			
			templaPath[tempi] = new CString(strfilepath1);//����ͼ������·��			
			//USES_CONVERSION;
			//char *temppath=T2A(strfilepath1.GetBuffer(0));
			//path[tempi] = 
			tempi++;
		} 
	} 
	find.Close(); 
	return; 
}

// ������С���λ���image�ϣ��������������������õĵ�ȫ�ֱ���
void drawRect::drawAllTemplate(CString image)
{
	USES_CONVERSION;
	char *image_char=T2A(image.GetBuffer(0));
	image.ReleaseBuffer();
	drawRect dog;
	for(int t = 0; t < tempi; t++) {  // tempi ���о�����Ŀ
		char *path = T2A((*templaPath[t]).GetBuffer(0)); // pathÿ������·��
		dog.fitTemplate(image_char,path); // ÿ������ԭͼ��
	}
} // drawAllTemplate