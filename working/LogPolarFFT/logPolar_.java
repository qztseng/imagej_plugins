import ij.*;
import ij.process.*;
import ij.plugin.filter.PlugInFilter;
import java.nio.FloatBuffer;
import java.lang.Math;
import java.nio.FloatBuffer;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import org.bytedeco.javacv.Java2DFrameUtils;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;

public class logPolar_ implements PlugInFilter {
	
	String title;
	
	public int setup(String arg, ImagePlus imp) {
			if (IJ.versionLessThan("1.5"))
				return DONE;
			else	
				title = imp.getTitle();
				return DOES_8G+DOES_16+DOES_32;
	}
	
	public void run(ImageProcessor ip) {
	
		int srcW = ip.getWidth(), 
			srcH = ip.getHeight(),
			M, N; 
		IplImage ipl_src, ipl_LP;
		double maxR, mag;
		FloatBuffer fb;
		float[] f;
		FloatProcessor m = null;
		
		Loader.load(opencv_core.class);         //call the loader before calling opencv core functions
		
/* 		//ipl_src = Java2DFrameUtils.toIplImage(ip.getBufferedImage());
		CvMat srcMat = new CvMat(srcH, srcW, CV_32FC1);
		//Mat srcMat = new Mat(srcH, srcW, CV_32FC1);
		double[] dArr1 = float2DtoDouble1DArray(ip.getFloatArray(), srcW, srcH);
		srcMat.put(0, dArr1, 0, dArr1.length);
		//srcMat.put(srcH, srcW, ip.getFloatArray()); */
		
		CvMat srcMat = img32toMat(ip, srcW, srcH);
		ipl_src = srcMat.asIplImage();
		
		/* 		
		src = Java2DFrameUtils.toMat(ip.getBufferedImage());
		ipl_src = Java2DFrameUtils.toIplImage(src);
		*/
		
		ipl_LP = cvCreateImage(cvSize(srcW, srcH), 32, 1); 
		
		//maxR = ((srcW>srcH)? srcH/2 : srcW/2);
		maxR = 0.5* Math.sqrt(srcW*srcW + srcH*srcH);
		mag = srcW/Math.log(maxR);
		IJ.log("maxR:"+maxR+", mag:"+mag);
			
		cvLogPolar(ipl_src, ipl_LP, cvPoint2D32f(srcW/2, srcH/2), mag, CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS);
		
		fb = FloatBuffer.allocate(ipl_LP.width() * ipl_LP.height());
		fb = ipl_LP.createBuffer();
		f = new float[fb.capacity()];
		fb.get(f);
		m = new FloatProcessor(ipl_LP.width(), ipl_LP.height(), f, null); 
		new ImagePlus("LPT of "+title, m).show();

	}
	
	/* Convert 32bit input image to CvMat */
	public static CvMat img32toMat(ImageProcessor ip32, int W, int H) {
		CvMat mat = CvMat.create(H, W, CV_32FC1);
		double[] dArr1 = float2DtoDouble1DArray(ip32.getFloatArray(), W, H);
		mat.put(0, dArr1, 0, dArr1.length);
		
		return mat;
	}
	
	
	/*
    *  The 2D array from ImageJ is [x][y], while the cvMat is in [y][x]
    */
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