package PIV;

import ij.*;
import ij.process.*;

import java.awt.image.BufferedImage;
import java.nio.FloatBuffer;

import static com.googlecode.javacv.cpp.opencv_core.*;
import static com.googlecode.javacv.cpp.opencv_imgproc.*;


/*
this class implement the template match method from OpenCV library. The interface between OpenCV and Java is 
using Samuel Audet's JavaCV code from: http://code.google.com/p/javacv/
This is a short version of the original cvMatch_Template plugin, which could be found at:
https://sites.google.com/site/qingzongtseng/template-matching-ij-plugin
By TSENG Qingzong (qztseng /at/ gmail.com)
2013/05/14 Add support for 32bit grayscale
2014/07/23 recompiled for better library support
*/
public class cvMatchTemplate {


    public static FloatProcessor doMatch(ImageProcessor src, ImageProcessor tpl, int method, boolean showR) {

        BufferedImage bi = null, bi2 = null;
        FloatProcessor resultFp = null;
        int srcW = src.getWidth();
        int srcH = src.getHeight();
        int tplW = tpl.getWidth();
        int tplH = tpl.getHeight();
        IplImage temp, temp2,res;
        IplImage iplSrc = null;
        IplImage iplTpl = null;
        
        switch (src.getBitDepth()) {

            case 32:
                //convert source imageProcessor into iplImage
                CvMat srcMat = CvMat.create(srcH, srcW, CV_32FC1);
                double[] dArr1 = float2DtoDouble1DArray(src.getFloatArray(), srcW, srcH);
                srcMat.put(0, dArr1, 0, dArr1.length);
                iplSrc = srcMat.asIplImage();
                //iplSrc = temp.clone();
                
                //convert template imageProcessor into iplImage
                CvMat tplMat = CvMat.create(tplH, tplW, CV_32FC1);
                double[] dArr2 = float2DtoDouble1DArray(tpl.getFloatArray(), tplW, tplH);
                tplMat.put(0, dArr2, 0, dArr2.length);
                iplTpl = tplMat.asIplImage();
                //iplTpl = temp2.clone();
                
                break;
            case 16:
                //since cvMatchTemplate don't accept 16bit image, we have to convert it to 32bit
                iplSrc = cvCreateImage(cvSize(srcW, srcH), IPL_DEPTH_32F, 1);
                bi = ((ShortProcessor) src).get16BitBufferedImage();
                temp = IplImage.createFrom(bi);
                cvConvertScale(temp, iplSrc, 1 / 65535.0, 0);

                iplTpl = cvCreateImage(cvSize(tplW, tplH), IPL_DEPTH_32F, 1);
                bi2 = ((ShortProcessor) tpl).get16BitBufferedImage();
                temp2 = IplImage.createFrom(bi2);
                cvConvertScale(temp2, iplTpl, 1 / 65535.0, 0);

                temp.release();
                temp2.release();

                break;
            case 8:

                bi = src.getBufferedImage();
                iplSrc = IplImage.createFrom(bi);

                bi2 = tpl.getBufferedImage();
                iplTpl = IplImage.createFrom(bi2);

                break;
            default:
                IJ.error("Unsupported image type");
                break;
        }

        res = cvCreateImage(cvSize(srcW - tplW + 1, srcH - tplH + 1), IPL_DEPTH_32F, 1);

        /*
        CV_TM_SQDIFF        = 0,
        CV_TM_SQDIFF_NORMED = 1,
        CV_TM_CCORR         = 2,
        CV_TM_CCORR_NORMED  = 3,
        CV_TM_CCOEFF        = 4,
        CV_TM_CCOEFF_NORMED = 5;
         */

        cvMatchTemplate(iplSrc, iplTpl, res, method);
        FloatBuffer fb = res.getFloatBuffer();
        float[] f = new float[res.width() * res.height()];
        fb.get(f, 0, f.length);
        resultFp = new FloatProcessor(res.width(), res.height(), f, null);
        cvReleaseImage(res);
        
        switch (src.getBitDepth()) {
            case 32:
                iplSrc.release();
                iplTpl.release();
                break;
            case 16:
                cvReleaseImage(iplSrc);
                cvReleaseImage(iplTpl);
                break;
            case 8:
                iplSrc.release();
                iplTpl.release();
                break;
            default:
                break;    
        }

        return resultFp;
    }

	private static double[] float2DtoDouble1DArray(float[][] arr2d, int column, int row) {


		double[] arr1d = new double[column * row];
		for (int y = 0; y < row; y++) {
		    for (int x = 0; x < column; x++) {
		        arr1d[y * column + x] = (double) arr2d[x][y];
		    }
		}

		 return arr1d;
    }

    
}
