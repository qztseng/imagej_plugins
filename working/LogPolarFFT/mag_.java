import ij.*;
import ij.process.*;
import ij.plugin.filter.PlugInFilter;
import java.nio.FloatBuffer;
import org.bytedeco.javacv.Java2DFrameUtils;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;

public class mag_ implements PlugInFilter {
	
	String title;
	
	public int setup(String arg, ImagePlus imp) {
			if (IJ.versionLessThan("1.5"))
				return DONE;
			else	
				title = imp.getTitle();
				return DOES_8G;
	}
	
	public void run(ImageProcessor ip) {
	
	Mat src = new Mat(),
		padded = new Mat(),
		complexI = new Mat(),
		magI = new Mat();
	int srcW, srcH, M, N; 
	FloatBuffer fb1;
	float[] f1;
	FloatProcessor m1 = null;
	
	src = Java2DFrameUtils.toMat(ip.getBufferedImage());
	srcW = src.cols();
	srcH = src.rows();
	
	M = getOptimalDFTSize(srcH);
	N = getOptimalDFTSize(srcW);
	
	copyMakeBorder(src, padded, 0, M - src.rows(), 0, N - src.cols(), BORDER_CONSTANT, Scalar.all(0));
	padded.convertTo(padded, CV_32F);
	
	MatVector planes = new MatVector(2);
	planes.put(0, padded);
    planes.put(1, Mat.zeros(padded.size(), CV_32F).asMat());
	merge(planes, complexI);         // Add to the expanded another plane with zeros
    dft(complexI, complexI);         // this way the result may fit in the source matrix
	
	/* compute the magnitude and switch to logarithmic scale
	=> log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2)) */
	split(complexI, planes);                               // planes.get(0) = Re(DFT(I)
																// planes.get(1) = Im(DFT(I))
	magnitude(planes.get(0), planes.get(1), planes.get(0));// planes.get(0) = magnitude
	magI = planes.get(0);
	
	fftShift(magI);
		
	fb1 = magI.<FloatBuffer>createBuffer();
	IJ.log("fb1 status: " + fb1.toString());
	IJ.log("fb1 has array? " + fb1.hasArray());
	f1 = new float[magI.cols() * magI.rows()];
	fb1.get(f1);
	m1 = new FloatProcessor(magI.cols(), magI.rows(), f1, null); 
	new ImagePlus("mag of "+title, m1).show();	
	
	}
	
	public void fftShift(Mat m){
		
		int xMid = m.cols()/2;
		int yMid = m.rows()/2;
				
		Mat q0 = new Mat(m, new Rect(0, 0, xMid, yMid));
		Mat q1 = new Mat(m, new Rect(xMid, 0, xMid, yMid));
		Mat q2 = new Mat(m, new Rect(0, yMid, xMid, yMid));
		Mat q3 = new Mat(m, new Rect(xMid, yMid, xMid, yMid));
		
		Mat tmp = new Mat();
		q0.copyTo(tmp);
		q3.copyTo(q0);
		tmp.copyTo(q3);

		q1.copyTo(tmp);
		q2.copyTo(q1);
		tmp.copyTo(q2);
	}
	
	

}
