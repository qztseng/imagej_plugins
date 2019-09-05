import ij.*;
import ij.process.*;
import ij.plugin.*;
import ij.gui.*;
import java.lang.Math;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import org.bytedeco.javacv.Java2DFrameUtils;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;


public class Phase_Correlation implements PlugIn {
	int[] wList;
	int wCount;
	String title1 = "",
		   title2 = "";
	String[] titles;
	ImagePlus img1, img2;
	Mat mat1, mat2, hann;
	Point2d result;

	@Override
	public void run(String arg) {
	
		 /* check how many existing images */
        wList = WindowManager.getIDList();
        if (wList == null) {
            IJ.noImage();
            return;
        }
        wCount = WindowManager.getImageCount();
        if (wCount < 2) {
            IJ.error("we need two images to do template matching");
            return;
        }

        /* show Dialog and get user input parameters */
        if (!getParams()) {
            return;
        }
        /* check parameter compatibility */

        if (img1.getBitDepth() != img2.getBitDepth()) {
            IJ.error("Images need to have the same type (bit depth)");
            return;
        }

        /* start matching */

		Loader.load(opencv_core.class);         //call the loader before calling opencv core functions
		
		int img1H = img1.getHeight();
		int img1W = img1.getWidth();
		int img2H = img2.getHeight();
		int img2W = img2.getWidth();
		
/* 		mat1 = cvarrToMat(img32toMat(img1.getProcessor(), img1W, img1H));
		mat2 = cvarrToMat(img32toMat(img2.getProcessor(), img2W, img2H)); */
		
		mat1 = ipToMat(img1.getProcessor());
		mat1.convertTo(mat1, CV_64FC1);
		mat2 = ipToMat(img2.getProcessor());
		mat2.convertTo(mat2, CV_64FC1);
		
		Point2d pt = phaseCorrelate(mat1, mat2);
		//C = phC(mat1,mat2);

		double maxR = ((img1W>img1H)? img1H/2 : img1W/2);
		double mag = img1W/Math.log(maxR);
		IJ.log("maxR:"+maxR+", mag:"+mag);

		double scale = Math.exp( pt.x() / mag);
		double rotate = pt.y()*360/img1.getWidth();

		IJ.log("pt.x:"+pt.x()+", pt.y:"+pt.y()+", scale:"+scale+", rotate:"+rotate);
		
	}
	
	
	/*phase correlate function to output correlation map instead of point of maximum)*/
	
/* 	private Mat phC(Mat m1, Mat m2){
		
		Mat result;
		int M = getOptimalDFTSize(m1.rows());
		int N = getOptimalDFTSize(m1.cols());
		
		Mat padded1 = new Mat();
		copyMakeBorder(m1, padded, 0, M - m1.rows(), 0, N - m1.cols(), BORDER_CONSTANT, Scalar.all(0));
		padded.convertTo(padded1, CV_32F);
		Mat padded2 = new Mat();
		copyMakeBorder(m2, padded, 0, M - m1.rows(), 0, N - m1.cols(), BORDER_CONSTANT, Scalar.all(0));
		padded.convertTo(padded2, CV_32F);
		
		Mat dft1, dft2;
		dft(padded1, dft1, DFT_REAL_OUTPUT);
		dft(padded2, dft2, DFT_REAL_OUTPUT);
		
		mulSpectrums(FFT1, FFT2, P, 0, true);
		magSpectrums(P, Pm);
		divSpectrums(P, Pm, C, 0, false);
		idft(C, C);
		fftShift(C); 
		
		
		
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
		
		
		
		return result;
	}
	 */

	/* Referenced from IJ-OpenCV's ImagePlusMatConverter.java
	** https://github.com/joheras/IJ-OpenCV
	*/
	
	private Mat ipToMat(ImageProcessor ip) {
        
		Mat mat = null;
        int w = ip.getWidth();
        int h = ip.getHeight();
				
		if (ip instanceof ByteProcessor) {
            byte[] pixels = (byte[]) ip.getPixels();
			mat = new Mat(h, w, opencv_core.CV_8UC1, new BytePointer(pixels));
        } else if (ip instanceof ShortProcessor) {
			short[] pixels = (short[]) ip.getPixels();
            mat = new Mat(h, w, opencv_core.CV_16UC1, new ShortPointer(pixels));
        } else if (ip instanceof FloatProcessor) {
            float[] pixels = (float[]) ip.getPixels();
			mat = new Mat(h, w, opencv_core.CV_32FC1, new FloatPointer(pixels));
        } else {
            throw new IllegalArgumentException("Unsupported pixel format: " + ip);
        }
        return mat;
	}
	
	private boolean getParams() {

        titles = new String[wList.length];
        for (int i = 0; i < wList.length; i++) {
            ImagePlus imp = WindowManager.getImage(wList[i]);
            if (imp != null) {
                titles[i] = imp.getTitle();
            } else {
                titles[i] = "";
            }
        }

        GenericDialog gd = new GenericDialog("PhaseCorrelation", IJ.getInstance());
        String defaultItem;

        if (title1.equals("")) {
            defaultItem = titles[0];
        } else {
            defaultItem = title1;
        }
        gd.addChoice("Image:", titles, defaultItem);
        if (title2.equals("")) {
            defaultItem = titles[1];
        } else {
            defaultItem = title2;
        }
        gd.addChoice("Template:", titles, defaultItem);
        /*gd.addNumericField("Magnitude for LogPolar transform", magnitude, 0);
		gd.addCheckbox("Output in resultTable ?", showRT);
        gd.addCheckbox("log result?", log); */
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        int index1 = gd.getNextChoiceIndex();
        title1 = titles[index1];
        int index2 = gd.getNextChoiceIndex();
        title2 = titles[index2];
        img1 = WindowManager.getImage(wList[index1]);
        img2 = WindowManager.getImage(wList[index2]);
        /* magnitude = gd.getNextNumber();
		showRT = gd.getNextBoolean();
        log = gd.getNextBoolean(); */

        return true;
    }
		
	/* Convert 32bit input image to CvMat */
/* 	public static CvMat img32toMat(ImageProcessor ip32, int W, int H) {
		CvMat mat = CvMat.create(H, W, CV_32FC1);
		double[] dArr1 = float2DtoDouble1DArray(ip32.getFloatArray(), W, H);
		mat.put(0, dArr1, 0, dArr1.length);
		
		return mat;
	} */

	/*
    *  The 2D array from ImageJ is [x][y], while the cvMat is in [y][x]
    */
	/*     private static double[] float2DtoDouble1DArray(float[][] arr2d, int column, int row) {

			double[] arr1d = new double[column * row];
			for (int y = 0; y < row; y++) {
				for (int x = 0; x < column; x++) {
					arr1d[y * column + x] = (double) arr2d[x][y];
				}
			}

			return arr1d;
		} */
	
}