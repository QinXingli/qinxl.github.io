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
//ð�ݷ�����
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
Description: Ӱ��ֱ��ƴ��
Calls:       ��
Input:       nImageCount����ƴ��Ӱ�������� *InputImages����ƴ��Ӱ��ȫ·�������飩��OutputImage�����Ӱ��ȫ·����
             
Output:      ���ƴ�Ӻ��Ӱ��
Return:      ����1��ƴ�ӳɹ�������0��ƴ��ʧ�� 
Others:      
*************************************************/ 
int MosaicDirectly(int nImageCount,string *InputImages,string OutputImage)
{
	//ע��GDAL
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF-8","NO");

	//��ƴ��Ӱ����
	int nImageNum = nImageCount;
	if (nImageNum == 0)
	{
		return 0;
	}
	//��ȡ����Ӱ��·��
	string *strInputImages;
	strInputImages = new string[nImageNum];
	//strInputImages[0] = "E:\\����\\�����Ʊ�ʾ�����ݼ�\\4Ӱ���Զ�ƴ��\\GF2_PMS2_E116.1_N29.4_20150531_L1A0000836045-MSS2_rectify_byte.tif";
	//strInputImages[1] = "E:\\����\\�����Ʊ�ʾ�����ݼ�\\4Ӱ���Զ�ƴ��\\GF2_PMS2_E116.1_N29.5_20150531_L1A0000836044-MSS2_rectify_byte.tif";
	for (int i=0;i<nImageNum;i++)
	{
		strInputImages[i] = InputImages[i];
	}

	//���Ӱ��·��
	string strOutputImage = OutputImage;
	//ÿ������Ӱ��������С����
	double *nMaxLng = new double[nImageNum];
	double *nMaxLat = new double[nImageNum];
	double *nMinLng = new double[nImageNum];
	double *nMinLat = new double[nImageNum];

	//ÿ������Ӱ�����Ԫ��С
	double *nSrcCellsizeX = new double[nImageNum];
	double *nSrcCellSizeY = new double[nImageNum];

	//ÿ������Ӱ��Ĳ�������Ӱ���С����������
	int *nSrcBandnum = new int[nImageNum];
	int *nImgCols = new int[nImageNum];
	int *nImgRows = new int[nImageNum];
	GDALDataType *nDataType = new GDALDataType[nImageNum];
	//����ͶӰ
	const char *pszSrcWKT = NULL;

	/* ѭ��ÿһ������Ӱ�� */
	/* ��ȡ����Ӱ����������귶Χ����Ԫ��С */
	for(int i=0; i<nImageNum; i++)
	{
		//��Ӱ��
		GDALDataset *poDataIn;
		poDataIn = (GDALDataset *)GDALOpen(strInputImages[i].c_str(),GA_ReadOnly);
		if (!poDataIn)
		{
			printf("%s ��ʧ��\n",strInputImages[i]);
			return 0;
		}
		//��ȡӰ�����
		int nRows = poDataIn->GetRasterYSize();//����
		int nCols = poDataIn->GetRasterXSize();//����
		nImgCols[i] = nCols;
		nImgRows[i] = nRows;
		nSrcBandnum[i] = poDataIn->GetRasterCount();//������
		double dGeoTrans[6];//����任����
		poDataIn->GetGeoTransform(dGeoTrans);
		nDataType[i] = poDataIn->GetRasterBand(1)->GetRasterDataType();//��������
		pszSrcWKT = poDataIn->GetProjectionRef();//��ȡͶӰ
		//���Ͻ�����
		nMinLng[i] = dGeoTrans[0];
		nMaxLat[i] = dGeoTrans[3];
		//���½�����
		nMaxLng[i] = dGeoTrans[0] + nCols*dGeoTrans[1] + nRows*dGeoTrans[2];
		nMinLat[i] = dGeoTrans[3] + nCols*dGeoTrans[4] + nRows*dGeoTrans[5];
		//ȷ�������Сֵ
		if(nMinLat[i]>nMaxLat[i]){double tmp=nMinLat[i]; nMinLat[i]=nMaxLat[i]; nMaxLat[i]=tmp;}
		if(nMinLng[i]>nMaxLng[i]){double tmp=nMinLng[i]; nMinLng[i]=nMaxLng[i]; nMaxLng[i]=tmp;}
		//��ȡ��Ԫ��С
		nSrcCellsizeX[i] = dGeoTrans[1];
		nSrcCellSizeY[i] = dGeoTrans[5];
		GDALClose(poDataIn);
	}



	/*��ȡ���Ӱ�����*/
	//��ȡ���Ӱ��Χ:ͨ������Ӱ������������Сֵ
	double nDstMaxLng = -200,nDstMinLng = 200, nDstMaxLat = -200, nDstMinLat = 200;
	for (int i=0; i<nImageNum; i++)
	{
		nDstMaxLng = max(nDstMaxLng,nMaxLng[i]);
		nDstMinLng = min(nDstMinLng,nMinLng[i]);
		nDstMaxLat = max(nDstMaxLat,nMaxLat[i]);
		nDstMinLat = min(nDstMinLat,nMinLat[i]);
	}
	//��ȡ���Ӱ��ķֱ��ʣ���ʱ������Ӱ��ֱ���ȡƽ��
	double nDstCellsizeX = 0.0, nDstCellsizeY = 0.0;
	for (int i=0; i<nImageNum; i++)
	{
		nDstCellsizeX += nSrcCellsizeX[i];
		nDstCellsizeY += nSrcCellSizeY[i];
	}
	nDstCellsizeX = nDstCellsizeX/nImageNum;
	nDstCellsizeY = nDstCellsizeY/nImageNum;
	//��ȡ���Ӱ�����任����
	double nDstGeotrans[6];
	nDstGeotrans[0] = nDstMinLng;
	nDstGeotrans[3] = nDstMaxLat;
	nDstGeotrans[1] = nDstCellsizeX;
	nDstGeotrans[2] = 0;
	nDstGeotrans[4] = 0;
	nDstGeotrans[5] = nDstCellsizeY;
	//��ȡ���Ӱ�񲨶���
	int nBandNum = nSrcBandnum[0];
	//��ȡ���Ӱ���С
	int  nDstCols, nDstRows;
	nDstCols = (int)((nDstMaxLng-nDstMinLng)/nDstCellsizeX+0.5);
	nDstRows = (int)((nDstMaxLat-nDstMinLat)/(-nDstCellsizeY)+0.5);
	//��ȡ���Ӱ����������
	GDALDataType nDstDataType = nDataType[0];

	/*�������Ӱ��*/
	const char* pzsFormat = "GTiff";//�ļ���ʽ
	char **papszOptions = NULL; 
	GDALDriver *poDrive = GetGDALDriverManager()->GetDriverByName(pzsFormat);
	GDALDataset *poDataOut = poDrive->Create(strOutputImage.c_str(),nDstCols,nDstRows,nBandNum,nDstDataType,papszOptions);
	poDataOut->SetGeoTransform(nDstGeotrans);
	char *pszWKT = NULL;
	const char* pszOutputWKT;
	//���û��ͶӰ������Ϊ����һ��ͶӰ
	if(pszSrcWKT == NULL)  
	{  
		OGRSpatialReference oSRS;   
		oSRS.SetUTM(50,true); //������  ����120��   
		oSRS.SetWellKnownGeogCS("WGS84");   
		oSRS.exportToWkt(&pszWKT); 
		pszOutputWKT = pszWKT;
	}
	else
		pszOutputWKT = pszSrcWKT;
	poDataOut->SetProjection(pszOutputWKT);
	//ͨ�����긳ֵ
	int nStripRow = 1024;//����������һ����������nStripRow��
	int nStripNum=(nDstRows+nStripRow-1)/nStripRow;//������
	int nBPP = nBandNum*nDstDataType;

	double maxX = nDstMaxLng;
	double minX = nDstMinLng;
	double maxY = nDstMaxLat;
	double minY = nDstMinLat;
	//�ж���������
	switch(nDstDataType)
	{
	case GDT_Byte:
		
		//�𲨶δ���
		for (int nBand=0; nBand<nBandNum; nBand++)
		{
			BYTE *pBuffer = new BYTE[int((nStripRow+1)*nDstCols)];
			//����������
			for(int nStripIndex=0;nStripIndex<nStripNum;nStripIndex++)
			{
				int nRow0=nStripIndex*nStripRow;
				int nRow1=nRow0+nStripRow;
				double stripMinY,stripMaxY;//������Y���귶Χ
				if(nRow1>nDstRows)
				{
					nRow1=nDstRows;
				}
				stripMaxY = maxY + nRow0*nDstCellsizeY;//nDstCellsizeY�Ǹ���
				stripMinY = maxY + nRow1*nDstCellsizeY;
				int nBackgroundColor = 0;
				switch(nBackgroundColor)//ѡ�񱳾�ɫ
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
				//��Ӱ����
				for( int i=0;i<nImageNum;i++)
				{
					int tempRow,tempCol;
					tempRow = nImgRows[i];
					tempCol = nImgCols[i];
					double tempX,tempY;
					tempX = nMinLng[i];
					tempY = nMaxLat[i];
					//�жϵ�ǰӰ���Ƿ�͵�ǰ�����н���
					if (nMaxLat[i] < stripMinY || nMinLat[i] > stripMaxY)
					{
						//��ǰӰ����������Χ��
						continue;
					}
					//���㵱ǰӰ���ж������ڵ�ǰ������
					//�жϷ�����ȡ������Y���귶Χ��Ӱ��Y���귶Χ�Ľ���
					double nTempLat[4];
					nTempLat[0] = stripMaxY;
					nTempLat[1] = stripMinY;
					nTempLat[2] = nMaxLat[i];
					nTempLat[3] = nMinLat[i];
					Bubble_Sort(nTempLat,4);
					//��ǰӰ��λ�������ڵ����귶Χ
					double curYMin,curYMax;
					curYMin = nTempLat[1];
					curYMax = nTempLat[2];
					//��ǰӰ��λ�������ڵĵ�ǰӰ����к�
					int curRow0,curRow1;
					curRow0 = (int)((nMaxLat[i] - curYMax)/(-nSrcCellSizeY[i])+0.5);
					curRow1 = (int)((nMaxLat[i] - curYMin)/(-nSrcCellSizeY[i])+0.5);
					//��ǰӰ��λ�������ڵ�����
					int curRowNum = curRow1 - curRow0;
					//�����ڴ�
					BYTE* tempBuffer=new BYTE[curRowNum*tempCol];

					//��ȡ��ǰӰ���ڴ���
					GDALRasterBand *poBand;
					GDALDataset *poDataIn;
					poDataIn = (GDALDataset *)GDALOpen(strInputImages[i].c_str(),GA_ReadOnly);
					poBand = poDataIn->GetRasterBand(nBand+1);
					poBand->RasterIO(GF_Read,0,curRow0,tempCol,curRowNum,tempBuffer,tempCol,curRowNum,nDstDataType,NULL,NULL);
					//��ȡnodatavalueֵ
					double* pnodatavalue = new double[nBandNum];
					memset(pnodatavalue, 0, sizeof(double)*nBandNum);
					int *pbHaveNoData = new int[nBandNum];
					memset(pbHaveNoData, 0, sizeof(int)*nBandNum);

					pnodatavalue[nBand] = GDALGetRasterNoDataValue(poBand, &pbHaveNoData[nBand]);

					//�����ش���
					for (int row=0;row<curRowNum;row++)
					{
						for (int col=0;col<tempCol;col++)
						{
							//��ȡ��ǰ����ֵ
							double nCurPixel = tempBuffer[row*tempCol+col];
							//��ȡ��ǰ���ص�����
							double curX,curY;
							curX = nMinLng[i] + col*nSrcCellsizeX[i];
							curY = curYMax + row*nSrcCellSizeY[i];
							//���㵱ǰ�����ڽ��Ӱ���ϵ���������
							int dstCol,dstRow;
							dstCol = (int)((curX-nDstMinLng)/nDstCellsizeX+0.5);
							dstRow = (int)((nDstMaxLat-curY)/(-nDstCellsizeY)+0.5);
							//ת���������е�λ��
							int stripCol,stripRow;
							stripCol = dstCol;
							stripRow = dstRow-nRow0;

							//�жϸõ��Ƿ��ڵ�ǰ������
							if (dstRow<nRow0 || dstRow > nRow1)
							{
								continue;
							}

							bool ballzero = true;

							double BackgroundValue = nBackgroundColor;
							//�жϵ�ǰ�����Ƿ���Nodata
							if (pbHaveNoData[nBand])
							{
								//��Nodata
								BackgroundValue = pnodatavalue[nBand];
							}

							if ((nCurPixel != BackgroundValue) && nCurPixel != 0)//background
							{
								//��ǰֵ�Ȳ����ڱ�����Ҳ����NaN��Ҳ�����ڼ���ֵ	
								ballzero = false;
							}
							//����ǰ����ֵ�������Ӱ��
							if (!ballzero)
							{
								pBuffer[stripRow*nDstCols + stripCol] = (BYTE)nCurPixel;
							}

						}
					}
					//�ͷ��ڴ�
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
			
			//�𲨶δ���
			for (int nBand=0; nBand<nBandNum; nBand++)
			{
				unsigned short* pBuffer = new unsigned short[int((nStripRow+1)*nDstCols)];
				//����������
				for(int nStripIndex=0;nStripIndex<nStripNum;nStripIndex++)
				{
					int nRow0=nStripIndex*nStripRow;
					int nRow1=nRow0+nStripRow;
					double stripMinY,stripMaxY;//������Y���귶Χ
					if(nRow1>nDstRows)
					{
						nRow1=nDstRows;
					}
					stripMaxY = maxY + nRow0*nDstCellsizeY;//nDstCellsizeY�Ǹ���
					stripMinY = maxY + nRow1*nDstCellsizeY;
					int nBackgroundColor = 0;
					switch(nBackgroundColor)//ѡ�񱳾�ɫ
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

					//��Ӱ����
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
							//��ǰӰ����������Χ��
							continue;
						}
						//���㵱ǰӰ���ж������ڵ�ǰ������
						//�жϷ�����ȡ������Y���귶Χ��Ӱ��Y���귶Χ�Ľ���
						double nTempLat[4];
						nTempLat[0] = stripMaxY;
						nTempLat[1] = stripMinY;
						nTempLat[2] = nMaxLat[i];
						nTempLat[3] = nMinLat[i];
						Bubble_Sort(nTempLat,4);
						//��ǰӰ��λ�������ڵ����귶Χ
						double curYMin,curYMax;
						curYMin = nTempLat[1];
						curYMax = nTempLat[2];
						//��ǰӰ��λ�������ڵĵ�ǰӰ����к�
						int curRow0,curRow1;
						curRow0 = (int)((nMaxLat[i] - curYMax)/(-nSrcCellSizeY[i])+0.5);
						curRow1 = (int)((nMaxLat[i] - curYMin)/(-nSrcCellSizeY[i])+0.5);
						//��ǰӰ��λ�������ڵ�����
						int curRowNum = curRow1 - curRow0;

						//�����ڴ�
						unsigned short* tempBuffer=new unsigned short[curRowNum*tempCol];


						//��ȡ��ǰӰ���ڴ���
						GDALRasterBand *poBand;
						GDALDataset *poDataIn;
						poDataIn = (GDALDataset *)GDALOpen(strInputImages[i].c_str(),GA_ReadOnly);
						poBand = poDataIn->GetRasterBand(nBand+1);
						poBand->RasterIO(GF_Read,0,curRow0,tempCol,curRowNum,tempBuffer,tempCol,curRowNum,nDstDataType,NULL,NULL);
						//��ȡnodatavalueֵ
						double* pnodatavalue = new double[nBandNum];
						memset(pnodatavalue, 0, sizeof(double)*nBandNum);
						int *pbHaveNoData = new int[nBandNum];
						memset(pbHaveNoData, 0, sizeof(int)*nBandNum);

						pnodatavalue[nBand] = GDALGetRasterNoDataValue(poBand, &pbHaveNoData[nBand]);

						//�����ش���
						for (int row=0;row<curRowNum;row++)
						{
							for (int col=0;col<tempCol;col++)
							{
								//��ȡ��ǰ����ֵ
								double nCurPixel = tempBuffer[row*tempCol+col];
								//��ȡ��ǰ���ص�����
								double curX,curY;
								curX = nMinLng[i] + col*nSrcCellsizeX[i];
								curY = curYMax + row*nSrcCellSizeY[i];
								//���㵱ǰ�����ڽ��Ӱ���ϵ���������
								int dstCol,dstRow;
								dstCol = (int)((curX-nDstMinLng)/nDstCellsizeX+0.5);
								dstRow = (int)((nDstMaxLat-curY)/(-nDstCellsizeY)+0.5);
								//ת���������е�λ��
								int stripCol,stripRow;
								stripCol = dstCol;
								stripRow = dstRow-nRow0;

								//�жϸõ��Ƿ��ڵ�ǰ������
								if (dstRow<nRow0 || dstRow > nRow1)
								{
									continue;
								}

								bool ballzero = true;

								double BackgroundValue = nBackgroundColor;
								//�жϵ�ǰ�����Ƿ���Nodata
								if (pbHaveNoData[nBand])
								{
									//��Nodata
									BackgroundValue = pnodatavalue[nBand];
								}

								if ((nCurPixel != BackgroundValue) && nCurPixel != 0)//background
								{
									//��ǰֵ�Ȳ����ڱ�����Ҳ����NaN��Ҳ�����ڼ���ֵ	
									ballzero = false;
								}
								//����ǰ����ֵ�������Ӱ��
								if (!ballzero)
								{
									pBuffer[stripRow*nDstCols + stripCol] = (unsigned short)nCurPixel;
								}
							}
						}
						//�ͷ��ڴ�
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