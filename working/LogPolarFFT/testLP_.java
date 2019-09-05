import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.*;
import ij.plugin.filter.PlugInFilter;
import ij.measure.ResultsTable;
import java.awt.image.BufferedImage;
import java.nio.FloatBuffer;
import java.nio.*;
import java.util.*;
import java.util.concurrent.TimeUnit;
import java.lang.Math;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import org.bytedeco.javacv.Java2DFrameUtils;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;


public class testLP_ implements PlugInFilter {

	public int setup(String arg, ImagePlus imp) {
		if (IJ.versionLessThan("1.5"))
			return DONE;
		else	
			return DOES_8G;
	}

	public void run(ImageProcessor ip) {
		
		BufferedImage bi = null, Res = null;
        FloatProcessor m1 = null, m2 = null;
        int srcW, srcH; 
		Mat Src = new Mat();
		IplImage ipl_magI, ipl_LP;
		double maxR, Mag;
		FloatBuffer fb1, fb2;
		float[] f1, f2;
		float ff;

		//bi = ip.getBufferedImage();
		//ipl_Src = Java2DFrameUtils.toIplImage(bi);
		Src = Java2DFrameUtils.toMat(ip.getBufferedImage());
		srcW = Src.cols();
		srcH = Src.rows();
		
		
		
		int M = getOptimalDFTSize(srcH);
		int N = getOptimalDFTSize(srcW);
		
		Mat padded = new Mat();
		copyMakeBorder(Src, padded, 0, M - Src.rows(), 0, N - Src.cols(), BORDER_CONSTANT, Scalar.all(0));
		
        padded.convertTo(padded, CV_32F);
		
		MatVector planes = new MatVector(2);
		planes.put(0, padded);
        planes.put(1, Mat.zeros(padded.size(), CV_32F).asMat());
        
		Mat complexI = new Mat();
        merge(planes, complexI);         // Add to the expanded another plane with zeros
        dft(complexI, complexI);         // this way the result may fit in the source matrix
        // compute the magnitude and switch to logarithmic scale
        // => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
        split(complexI, planes);                               // planes.get(0) = Re(DFT(I)
                                                                    // planes.get(1) = Im(DFT(I))
        magnitude(planes.get(0), planes.get(1), planes.get(0));// planes.get(0) = magnitude
        Mat magI = planes.get(0);
		
		fb1 = magI.<FloatBuffer>createBuffer();
		IJ.log("fb1 status: " + fb1.toString());
		IJ.log("fb1 has array? " + fb1.hasArray());
		f1 = new float[magI.cols() * magI.rows()];
        fb1.get(f1);
        m1 = new FloatProcessor(magI.cols(), magI.rows(), f1, null); 
		new ImagePlus("m1", m1).show();	
		
		ipl_magI = Java2DFrameUtils.toIplImage(magI);
		
		IJ.log("magI depth: " + ipl_magI.depth());
		IJ.log("magI width: " + ipl_magI.width());
		IJ.log("magI height: " + ipl_magI.height());
		
		//ipl_LP = Java2DFrameUtils.toIplImage(Mat.zeros(ipl_magI.height(), ipl_magI.width(), CV_32F).asMat());
		//ipl_LP = new IplImage(Mat.zeros(ipl_magI.height(), ipl_magI.width(), CV_32F).asMat());
		ipl_LP = cvCreateImage(cvSize(ipl_magI.width(), ipl_magI.height()), 32, 1);
		
		maxR = ((srcW>srcH)? srcH/2 : srcW/2);
		Mag = srcW/Math.log(maxR);
		
		IJ.log("maxR:" + maxR);
		IJ.log("Mag: " + Mag);
		
		cvLogPolar(ipl_magI, ipl_LP, cvPoint2D32f(srcW/2, srcH/2), Mag);
		IJ.log("LP depth: " + ipl_LP.depth());
		IJ.log("LP width: " + ipl_LP.width());
		IJ.log("LP height: " + ipl_LP.height());
		IJ.log("magI to string: " + ipl_magI.toString());
		IJ.log("LP to string: " + ipl_LP.toString());
		
		fb2 = FloatBuffer.allocate(ipl_LP.width() * ipl_LP.height());
 		fb2 = ipl_LP.createBuffer();
		//fb2 = ipl_LP.getFloatBuffer();
		IJ.log("fb2 status: " + fb2.toString());
		IJ.log("fb2 has array? " + fb2.hasArray());
		f2 = new float[fb2.capacity()];
		
/* 		for(int i=0; i<fb2.capacity();i++){
			ff = fb2.get(i); 
			IJ.log("f2 pos: " + i + "value: " + ff );
			try
			{
				Thread.sleep(10);
			}
			catch(InterruptedException ex)
			{
				Thread.currentThread().interrupt();
			} 
		} */
		
		
		
		//f2 = new float[ipl_LP.width() * ipl_LP.height()];
        fb2.get(f2);
		
        m2 = new FloatProcessor(ipl_LP.width(), ipl_LP.height(), f2, null); 
		new ImagePlus("m2", m2).show();	 
		//Res = Java2DFrameUtils.toBufferedImage(padded);
		//ByteProcessor resultBp = new ByteProcessor(Res);
		//new ImagePlus("result", resultBp).show();
	
		return;
	}
	
	 public void swapQuadrants(ImageProcessor ip) {
        ImageProcessor t1, t2;
        int size = ip.getWidth()/2;
        ip.setRoi(size,0,size,size);
        t1 = ip.crop();
        ip.setRoi(0,size,size,size);
        t2 = ip.crop();
        ip.insert(t1,0,size);
        ip.insert(t2,size,0);
        ip.setRoi(0,0,size,size);
        t1 = ip.crop();
        ip.setRoi(size,size,size,size);
        t2 = ip.crop();
        ip.insert(t1,size,size);
        ip.insert(t2,0,0);
        ip.resetRoi();
    }
	
	

}
