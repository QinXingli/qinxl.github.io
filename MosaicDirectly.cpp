// MosaicDirectly.cpp : Defines the entry point for the console application.
//
/************************************************* 
Copyright: QinXingli 
Author: QinXingli
Date:2017-11-02
Description: Mosaic Images
**************************************************/  
#include "stdafx.h"
#include "gdal_priv.h"
#include "ogrsf_frmts.h"
#include <math.h>
using namespace std;
typedef unsigned char BYTE;

#define DBL_MAX         1.7976931348623158e+308 /* max value */
#define FLT_MAX         3.402823466e+38F        /* max value */
//冒泡法排序
void Bubble_Sort(double *num, int n)
{
	int i, j;
	for(i = 0; i < n; i++)
	{
		for(j = 0; i + j < n - 1; j++)
		{
			if(num[j] > num[j + 1])
			{
				double temp = num[j];
				num[j] = num[j + 1];
				num[j + 1] = temp;
			}
		}
	}
	return;
}
/************************************************* 
Function: MosaicDirectly
Description: 影像直接拼接
Calls:       无
Input:       nImageCount（待拼接影像数）、 *InputImages（待拼接影像全路径，数组）、OutputImage（输出影像全路径）
             
Output:      输出拼接后的影像
Return:      返回1，拼接成功，返回0，拼接失败 
Others:      
*************************************************/ 
int MosaicDirectly(int nImageCount,string *InputImages,string OutputImage)
{
	//注册GDAL
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF-8","NO");

	//待拼接影像数
	int nImageNum = nImageCount;
	if (nImageNum == 0)
	{
		return 0;
	}
	//获取输入影像路径
	string *strInputImages;
	strInputImages = new string[nImageNum];
	//strInputImages[0] = "E:\\数据\\数据制备示例数据集\\4影像自动拼接\\GF2_PMS2_E116.1_N29.4_20150531_L1A0000836045-MSS2_rectify_byte.tif";
	//strInputImages[1] = "E:\\数据\\数据制备示例数据集\\4影像自动拼接\\GF2_PMS2_E116.1_N29.5_20150531_L1A0000836044-MSS2_rectify_byte.tif";
	for (int i=0;i<nImageNum;i++)
	{
		strInputImages[i] = InputImages[i];
	}

	//输出影像路径
	string strOutputImage = OutputImage;
	//每幅输入影像的最大最小坐标
	double *nMaxLng = new double[nImageNum];
	double *nMaxLat = new double[nImageNum];
	double *nMinLng = new double[nImageNum];
	double *nMinLat = new double[nImageNum];

	//每幅输入影像的像元大小
	double *nSrcCellsizeX = new double[nImageNum];
	double *nSrcCellSizeY = new double[nImageNum];

	//每幅输入影像的波段数、影像大小、数据类型
	int *nSrcBandnum = new int[nImageNum];
	int *nImgCols = new int[nImageNum];
	int *nImgRows = new int[nImageNum];
	GDALDataType *nDataType = new GDALDataType[nImageNum];
	//设置投影
	const char *pszSrcWKT = NULL;

	/* 循环每一幅输入影像 */
	/* 获取输入影像参数：坐标范围，像元大小 */
	for(int i=0; i<nImageNum; i++)
	{
		//打开影像
		GDALDataset *poDataIn;
		poDataIn = (GDALDataset *)GDALOpen(strInputImages[i].c_str(),GA_ReadOnly);
		if (!poDataIn)
		{
			printf("%s 打开失败\n",strInputImages[i]);
			return 0;
		}
		//获取影像参数
		int nRows = poDataIn->GetRasterYSize();//行数
		int nCols = poDataIn->GetRasterXSize();//列数
		nImgCols[i] = nCols;
		nImgRows[i] = nRows;
		nSrcBandnum[i] = poDataIn->GetRasterCount();//波段数
		double dGeoTrans[6];//仿射变换参数
		poDataIn->GetGeoTransform(dGeoTrans);
		nDataType[i] = poDataIn->GetRasterBand(1)->GetRasterDataType();//数据类型
		pszSrcWKT = poDataIn->GetProjectionRef();//获取投影
		//左上角坐标
		nMinLng[i] = dGeoTrans[0];
		nMaxLat[i] = dGeoTrans[3];
		//右下角坐标
		nMaxLng[i] = dGeoTrans[0] + nCols*dGeoTrans[1] + nRows*dGeoTrans[2];
		nMinLat[i] = dGeoTrans[3] + nCols*dGeoTrans[4] + nRows*dGeoTrans[5];
		//确定最大最小值
		if(nMinLat[i]>nMaxLat[i]){double tmp=nMinLat[i]; nMinLat[i]=nMaxLat[i]; nMaxLat[i]=tmp;}
		if(nMinLng[i]>nMaxLng[i]){double tmp=nMinLng[i]; nMinLng[i]=nMaxLng[i]; nMaxLng[i]=tmp;}
		//获取像元大小
		nSrcCellsizeX[i] = dGeoTrans[1];
		nSrcCellSizeY[i] = dGeoTrans[5];
		GDALClose(poDataIn);
	}



	/*获取结果影像参数*/
	//获取结果影像范围:通过输入影像坐标的最大最小值
	double nDstMaxLng = -200,nDstMinLng = 200, nDstMaxLat = -200, nDstMinLat = 200;
	for (int i=0; i<nImageNum; i++)
	{
		nDstMaxLng = max(nDstMaxLng,nMaxLng[i]);
		nDstMinLng = min(nDstMinLng,nMinLng[i]);
		nDstMaxLat = max(nDstMaxLat,nMaxLat[i]);
		nDstMinLat = min(nDstMinLat,nMinLat[i]);
	}
	//获取输出影像的分辨率：暂时由输入影像分辨率取平均
	double nDstCellsizeX = 0.0, nDstCellsizeY = 0.0;
	for (int i=0; i<nImageNum; i++)
	{
		nDstCellsizeX += nSrcCellsizeX[i];
		nDstCellsizeY += nSrcCellSizeY[i];
	}
	nDstCellsizeX = nDstCellsizeX/nImageNum;
	nDstCellsizeY = nDstCellsizeY/nImageNum;
	//获取结果影像仿射变换参数
	double nDstGeotrans[6];
	nDstGeotrans[0] = nDstMinLng;
	nDstGeotrans[3] = nDstMaxLat;
	nDstGeotrans[1] = nDstCellsizeX;
	nDstGeotrans[2] = 0;
	nDstGeotrans[4] = 0;
	nDstGeotrans[5] = nDstCellsizeY;
	//获取结果影像波段数
	int nBandNum = nSrcBandnum[0];
	//获取结果影像大小
	int  nDstCols, nDstRows;
	nDstCols = (int)((nDstMaxLng-nDstMinLng)/nDstCellsizeX+0.5);
	nDstRows = (int)((nDstMaxLat-nDstMinLat)/(-nDstCellsizeY)+0.5);
	//获取结果影像数据类型
	GDALDataType nDstDataType = nDataType[0];

	/*创建结果影像*/
	const char* pzsFormat = "GTiff";//文件格式
	char **papszOptions = NULL; 
	GDALDriver *poDrive = GetGDALDriverManager()->GetDriverByName(pzsFormat);
	GDALDataset *poDataOut = poDrive->Create(strOutputImage.c_str(),nDstCols,nDstRows,nBandNum,nDstDataType,papszOptions);
	poDataOut->SetGeoTransform(nDstGeotrans);
	char *pszWKT = NULL;
	const char* pszOutputWKT;
	//如果没有投影，则人为设置一个投影
	if(pszSrcWKT == NULL)  
	{  
		OGRSpatialReference oSRS;   
		oSRS.SetUTM(50,true); //北半球  东经120度   
		oSRS.SetWellKnownGeogCS("WGS84");   
		oSRS.exportToWkt(&pszWKT); 
		pszOutputWKT = pszWKT;
	}
	else
		pszOutputWKT = pszSrcWKT;
	poDataOut->SetProjection(pszOutputWKT);
	//通过坐标赋值
	int nStripRow = 1024;//逐条带处理，一个条带里有nStripRow行
	int nStripNum=(nDstRows+nStripRow-1)/nStripRow;//条带数
	int nBPP = nBandNum*nDstDataType;

	double maxX = nDstMaxLng;
	double minX = nDstMinLng;
	double maxY = nDstMaxLat;
	double minY = nDstMinLat;
	//判断数据类型
	switch(nDstDataType)
	{
	case GDT_Byte:
		
		//逐波段处理
		for (int nBand=0; nBand<nBandNum; nBand++)
		{
			BYTE *pBuffer = new BYTE[int((nStripRow+1)*nDstCols)];
			//逐条带处理
			for(int nStripIndex=0;nStripIndex<nStripNum;nStripIndex++)
			{
				int nRow0=nStripIndex*nStripRow;
				int nRow1=nRow0+nStripRow;
				double stripMinY,stripMaxY;//条带的Y坐标范围
				if(nRow1>nDstRows)
				{
					nRow1=nDstRows;
				}
				stripMaxY = maxY + nRow0*nDstCellsizeY;//nDstCellsizeY是负数
				stripMinY = maxY + nRow1*nDstCellsizeY;
				int nBackgroundColor = 0;
				switch(nBackgroundColor)//选择背景色
				{
				case 0:
					memset(pBuffer,0,nStripRow*nDstCols*nDstDataType);
					break;
				case 1:
					memset(pBuffer,255,nStripRow*nDstCols*nDstDataType);
					break;
				default:
					memset(pBuffer,255,nStripRow*nDstCols*nDstDataType);
					break;
				}
				//逐影像处理
				for( int i=0;i<nImageNum;i++)
				{
					int tempRow,tempCol;
					tempRow = nImgRows[i];
					tempCol = nImgCols[i];
					double tempX,tempY;
					tempX = nMinLng[i];
					tempY = nMaxLat[i];
					//判断当前影像是否和当前条带有交集
					if (nMaxLat[i] < stripMinY || nMinLat[i] > stripMaxY)
					{
						//当前影像不在条带范围内
						continue;
					}
					//计算当前影像有多少行在当前条带内
					//判断方法：取条带的Y坐标范围和影像Y坐标范围的交集
					double nTempLat[4];
					nTempLat[0] = stripMaxY;
					nTempLat[1] = stripMinY;
					nTempLat[2] = nMaxLat[i];
					nTempLat[3] = nMinLat[i];
					Bubble_Sort(nTempLat,4);
					//当前影像位于条带内的坐标范围
					double curYMin,curYMax;
					curYMin = nTempLat[1];
					curYMax = nTempLat[2];
					//当前影像位于条带内的当前影像的行号
					int curRow0,curRow1;
					curRow0 = (int)((nMaxLat[i] - curYMax)/(-nSrcCellSizeY[i])+0.5);
					curRow1 = (int)((nMaxLat[i] - curYMin)/(-nSrcCellSizeY[i])+0.5);
					//当前影像位于条带内的行数
					int curRowNum = curRow1 - curRow0;
					//分配内存
					BYTE* tempBuffer=new BYTE[curRowNum*tempCol];

					//读取当前影像到内存中
					GDALRasterBand *poBand;
					GDALDataset *poDataIn;
					poDataIn = (GDALDataset *)GDALOpen(strInputImages[i].c_str(),GA_ReadOnly);
					poBand = poDataIn->GetRasterBand(nBand+1);
					poBand->RasterIO(GF_Read,0,curRow0,tempCol,curRowNum,tempBuffer,tempCol,curRowNum,nDstDataType,NULL,NULL);
					//获取nodatavalue值
					double* pnodatavalue = new double[nBandNum];
					memset(pnodatavalue, 0, sizeof(double)*nBandNum);
					int *pbHaveNoData = new int[nBandNum];
					memset(pbHaveNoData, 0, sizeof(int)*nBandNum);

					pnodatavalue[nBand] = GDALGetRasterNoDataValue(poBand, &pbHaveNoData[nBand]);

					//逐像素处理
					for (int row=0;row<curRowNum;row++)
					{
						for (int col=0;col<tempCol;col++)
						{
							//获取当前像素值
							double nCurPixel = tempBuffer[row*tempCol+col];
							//获取当前像素的坐标
							double curX,curY;
							curX = nMinLng[i] + col*nSrcCellsizeX[i];
							curY = curYMax + row*nSrcCellSizeY[i];
							//计算当前像素在结果影像上的像素坐标
							int dstCol,dstRow;
							dstCol = (int)((curX-nDstMinLng)/nDstCellsizeX+0.5);
							dstRow = (int)((nDstMaxLat-curY)/(-nDstCellsizeY)+0.5);
							//转换成条带中的位置
							int stripCol,stripRow;
							stripCol = dstCol;
							stripRow = dstRow-nRow0;

							//判断该点是否在当前条带上
							if (dstRow<nRow0 || dstRow > nRow1)
							{
								continue;
							}

							bool ballzero = true;

							double BackgroundValue = nBackgroundColor;
							//判断当前波段是否有Nodata
							if (pbHaveNoData[nBand])
							{
								//有Nodata
								BackgroundValue = pnodatavalue[nBand];
							}

							if ((nCurPixel != BackgroundValue) && nCurPixel != 0)//background
							{
								//当前值既不等于背景，也不是NaN，也不等于极大值	
								ballzero = false;
							}
							//将当前像素值赋给结果影像
							if (!ballzero)
							{
								pBuffer[stripRow*nDstCols + stripCol] = (BYTE)nCurPixel;
							}

						}
					}
					//释放内存
					delete []tempBuffer;
					GDALClose(poDataIn);
				}
				printf("band=%d,nRow0=%d, nRow1=%d\n",nBand,nRow0,nRow1);
				GDALRasterBand *poBandOut;
				poBandOut = poDataOut->GetRasterBand(nBand+1);
				poBandOut->RasterIO(GF_Write,0,nRow0,nDstCols,(nRow1-nRow0),pBuffer,nDstCols,(nRow1-nRow0),nDstDataType,NULL,NULL);
			}
			delete []pBuffer;
		}	
		break;
		case GDT_UInt16:
			
			//逐波段处理
			for (int nBand=0; nBand<nBandNum; nBand++)
			{
				unsigned short* pBuffer = new unsigned short[int((nStripRow+1)*nDstCols)];
				//逐条带处理
				for(int nStripIndex=0;nStripIndex<nStripNum;nStripIndex++)
				{
					int nRow0=nStripIndex*nStripRow;
					int nRow1=nRow0+nStripRow;
					double stripMinY,stripMaxY;//条带的Y坐标范围
					if(nRow1>nDstRows)
					{
						nRow1=nDstRows;
					}
					stripMaxY = maxY + nRow0*nDstCellsizeY;//nDstCellsizeY是负数
					stripMinY = maxY + nRow1*nDstCellsizeY;
					int nBackgroundColor = 0;
					switch(nBackgroundColor)//选择背景色
					{
					case 0:
						memset(pBuffer,0,nStripRow*nDstCols*nDstDataType);
						break;
					case 1:
						memset(pBuffer,255,nStripRow*nDstCols*nDstDataType);
						break;
					default:
						memset(pBuffer,255,nStripRow*nDstCols*nDstDataType);
						break;
					}

					//逐影像处理
					for( int i=0;i<nImageNum;i++)
					{
						int tempRow,tempCol;
						tempRow = nImgRows[i];
						tempCol = nImgCols[i];
						double tempX,tempY;
						tempX = nMinLng[i];
						tempY = nMaxLat[i];
						if (nMaxLat[i] < stripMinY || nMinLat[i] > stripMaxY)
						{
							//当前影像不在条带范围内
							continue;
						}
						//计算当前影像有多少行在当前条带内
						//判断方法：取条带的Y坐标范围和影像Y坐标范围的交集
						double nTempLat[4];
						nTempLat[0] = stripMaxY;
						nTempLat[1] = stripMinY;
						nTempLat[2] = nMaxLat[i];
						nTempLat[3] = nMinLat[i];
						Bubble_Sort(nTempLat,4);
						//当前影像位于条带内的坐标范围
						double curYMin,curYMax;
						curYMin = nTempLat[1];
						curYMax = nTempLat[2];
						//当前影像位于条带内的当前影像的行号
						int curRow0,curRow1;
						curRow0 = (int)((nMaxLat[i] - curYMax)/(-nSrcCellSizeY[i])+0.5);
						curRow1 = (int)((nMaxLat[i] - curYMin)/(-nSrcCellSizeY[i])+0.5);
						//当前影像位于条带内的行数
						int curRowNum = curRow1 - curRow0;

						//分配内存
						unsigned short* tempBuffer=new unsigned short[curRowNum*tempCol];


						//读取当前影像到内存中
						GDALRasterBand *poBand;
						GDALDataset *poDataIn;
						poDataIn = (GDALDataset *)GDALOpen(strInputImages[i].c_str(),GA_ReadOnly);
						poBand = poDataIn->GetRasterBand(nBand+1);
						poBand->RasterIO(GF_Read,0,curRow0,tempCol,curRowNum,tempBuffer,tempCol,curRowNum,nDstDataType,NULL,NULL);
						//获取nodatavalue值
						double* pnodatavalue = new double[nBandNum];
						memset(pnodatavalue, 0, sizeof(double)*nBandNum);
						int *pbHaveNoData = new int[nBandNum];
						memset(pbHaveNoData, 0, sizeof(int)*nBandNum);

						pnodatavalue[nBand] = GDALGetRasterNoDataValue(poBand, &pbHaveNoData[nBand]);

						//逐像素处理
						for (int row=0;row<curRowNum;row++)
						{
							for (int col=0;col<tempCol;col++)
							{
								//获取当前像素值
								double nCurPixel = tempBuffer[row*tempCol+col];
								//获取当前像素的坐标
								double curX,curY;
								curX = nMinLng[i] + col*nSrcCellsizeX[i];
								curY = curYMax + row*nSrcCellSizeY[i];
								//计算当前像素在结果影像上的像素坐标
								int dstCol,dstRow;
								dstCol = (int)((curX-nDstMinLng)/nDstCellsizeX+0.5);
								dstRow = (int)((nDstMaxLat-curY)/(-nDstCellsizeY)+0.5);
								//转换成条带中的位置
								int stripCol,stripRow;
								stripCol = dstCol;
								stripRow = dstRow-nRow0;

								//判断该点是否在当前条带上
								if (dstRow<nRow0 || dstRow > nRow1)
								{
									continue;
								}

								bool ballzero = true;

								double BackgroundValue = nBackgroundColor;
								//判断当前波段是否有Nodata
								if (pbHaveNoData[nBand])
								{
									//有Nodata
									BackgroundValue = pnodatavalue[nBand];
								}

								if ((nCurPixel != BackgroundValue) && nCurPixel != 0)//background
								{
									//当前值既不等于背景，也不是NaN，也不等于极大值	
									ballzero = false;
								}
								//将当前像素值赋给结果影像
								if (!ballzero)
								{
									pBuffer[stripRow*nDstCols + stripCol] = (unsigned short)nCurPixel;
								}
							}
						}
						//释放内存
						delete []tempBuffer;
						GDALClose(poDataIn);
					}
					printf("band=%d,nRow0=%d, nRow1=%d\n",nBand,nRow0,nRow1);
					GDALRasterBand *poBandOut;
					poBandOut = poDataOut->GetRasterBand(nBand+1);
					poBandOut->RasterIO(GF_Write,0,nRow0,nDstCols,(nRow1-nRow0),pBuffer,nDstCols,(nRow1-nRow0),nDstDataType,NULL,NULL);
				}
				delete []pBuffer;
			}
			
			break;
	}
	
	
	GDALClose(poDataOut);
	return 1;
}