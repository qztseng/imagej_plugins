package FourierMellinRegistration;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Color;
import ij.plugin.*;
import ij.plugin.filter.MaximumFinder;
import ij.measure.*; //image stats.mean
import ij.measure.ResultsTable;
import java.awt.image.BufferedImage;
import java.nio.FloatBuffer;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import org.bytedeco.javacv.Java2DFrameUtils;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;


/*

 */
public class FM_Registration implements PlugIn {

    static String title1 = "";
    static String title2 = "";
    String[] titles;
	static int cutoffHigh=20, cutoffLow=30;
    ImagePlus img1, img2;
    static boolean showR = true, 
					subpixel = true;
	//private int size;
    int[] dxdy = new int[2];
    int[] wList;
	static float magLP;
	static double rot;
	static float scale;
	static int dx, dy;

    @Override
    public void run(String arg) {

        if (arg.equals("about")) {
            showAbout();
            return;
        }
        /*
         ** check how many existing images
         */
        wList = WindowManager.getIDList();
        if (wList == null) {
            IJ.noImage();
            return;
        }
        int wCount = WindowManager.getImageCount();
        if (wCount < 2) {
            IJ.error("we need two images to do registration");
            return;
        }
        /*
         ** show Dialog and get user input parameters
         */
        if (!getParams()) {
            return;
        }
        /*
         ** check image compatibility
         */

        if (img1.getBitDepth() != img2.getBitDepth()) {
            IJ.error("Images need to have the same type (bit depth)");
            return;
        }
		
		if (img1.getWidth() != img2.getWidth() | img1.getHeight() != img2.getHeight() ) {
            IJ.error("Images need to have the same dimensions");
            return;
        } 
		
		if (img1.getBitDepth() == 24) {
            IJ.error("Please convert color (RGB) images to grayscale first.");
            return;
        }
		IJ.showProgress(1,5);
		long startTime = System.currentTimeMillis();
		
		/* get ImageProcessor */
		ImageProcessor ip1 = img1.getProcessor();
		ImageProcessor ip2 = img2.getProcessor();
		
		/* do match (which will write rot, scale, dx, dy to the class fields */
		doMatch(ip1,ip2, cutoffHigh, cutoffLow);
						
		IJ.showProgress(4,5);
		long elapsedTime = System.currentTimeMillis() - startTime;

		IJ.log("rotate second image: " + rot + " , move x: "+dx+" , move y: "+dy);
		ImageStack stk = img1.createEmptyStack();
		stk.addSlice("original", ip1);
		ImageProcessor ip2dup = ip2.duplicate();
		ip2dup.rotate(rot);
		ip2dup.translate(dx, dy);
		stk.addSlice("registered", ip2dup);
		new ImagePlus("Aligned Images", stk).show();

		IJ.showProgress(5,5);
        IJ.showStatus("Matching done in " + elapsedTime + "msec");
        IJ.log("Matching done in " + elapsedTime + "msec");

		return;
    }
	
	public static void doMatch(ImageProcessor ip1, ImageProcessor ip2, int cHigh, int cLow){
	
		/*pad each image*/
		ImageProcessor paddedOri = getPaddedIP(ip1);
		ImageProcessor padded1 = getPaddedIP(ip2);	
		
		/*get filtered magnitude map and Log-Polar transformed*/
		FloatProcessor mag_ori = transformFM(paddedOri, cutoffHigh, cutoffLow);
		FloatProcessor mag_trf = transformFM(padded1, cutoffHigh, cutoffLow);
		
		/* Do cross correlation */
		FloatProcessor result =  xcorr(mag_ori, mag_trf);
		int size_result = result.getWidth();
		
		/* find maximum on the correlation map */
		int[] dxdy = findMax(result);
	
		double rot1 = dxdy[1]*360f/size_result;
		double rot2 = dxdy[3]*360f/size_result;		
		float scale1 = (float)Math.exp(dxdy[0]/magLP);
		float scale2 = (float)Math.exp(dxdy[2]/magLP);
		      
		//IJ.log("rotation1: " + rot1 + " ; rotation2: "+rot2);
		//IJ.log("scale1: " + scale1 + " ; scale2: "+scale2);
		/* select the correct rotation (due to 180degree ambiguity introduced by using the magnitude spectra) */
		ImageProcessor padded2 = padded1.duplicate();
		/* get image mean for rotation fill */
		ImageStatistics stats = ip2.getStats();	
		
		padded1.setBackgroundValue(stats.mean);
		padded1.rotate(rot1);
		padded2.setBackgroundValue(stats.mean);
		padded2.rotate(rot2);
		
		//new ImagePlus("ip1", padded1).show();
		//new ImagePlus("ip2", padded2).show();
		
		FloatProcessor result1 = xcorr(paddedOri, padded1);
		FloatProcessor result2 = xcorr(paddedOri, padded2);
		
		//new ImagePlus("result1", result1).show();
		//new ImagePlus("result2", result2).show();
		
		int[] dxdy1 = findMax(result1);
		int[] dxdy2 = findMax(result2);
		double max1 = result1.getMax();
		double max2 = result2.getMax();
		
		dx = dxdy1[0];
		dy = dxdy1[1];
		rot = rot1;
		scale = scale1;
		
		if(max2 > max1){
			dx = dxdy2[0];
			dy = dxdy2[1];
			rot = rot2;
			scale = scale2;
		}
		
/* 		if (showR) {
			ResultsTable rtNew = new ResultsTable();
			rtNew.incrementCounter();
			rtNew.addLabel("Match 1");
			rtNew.addValue("X", dxdy[0]);
			rtNew.addValue("Y", dxdy[1]);
			rtNew.addValue("rotation", rot1);
			rtNew.addValue("scale", scale1);
			rtNew.addValue("xcorrMax", max1);
			rtNew.addValue("dx", dxdy1[0]);
			rtNew.addValue("dy", dxdy1[1]);
			rtNew.incrementCounter();
			rtNew.addLabel("Match 2");
			rtNew.addValue("X", dxdy[2]);
			rtNew.addValue("Y", dxdy[3]);
			rtNew.addValue("rotation", rot2);
			rtNew.addValue("scale", scale2);
			rtNew.addValue("xcorrMax", max2);
			rtNew.addValue("dx", dxdy2[0]);
			rtNew.addValue("dy", dxdy2[1]);
			rtNew.updateResults();
			rtNew.show("FMReg_Results");
		}				
		if(showR){
			new ImagePlus("xcorr result", result).show();
		} */
		
	}
	
	public static FloatProcessor transformFM(ImagePlus img, int cHigh, int cLow){
		
		return transformFM(img.getProcessor(), cHigh, cLow);
	}
	
	public static FloatProcessor transformFM(ImageProcessor ip, int cHigh, int cLow){
		
		/*image padding*/
		ImageProcessor ip2 = getPaddedIP(ip);
		/*get padded size (width = height & img1 = img2 )*/
		int size = ip2.getWidth();
		
		/* pre-processing by high and low pass filter and return a filtered amplitude spectrum */
		FloatProcessor amp_filtered = HLpassFiltering(ip2, cHigh, cLow);
	
		/*call the loader before calling opencv core functions*/
		Loader.load(opencv_core.class);
	
		/*Do log-polar transformation*/
		FloatProcessor m = logPolarT(amp_filtered);
		
		return m;	
		
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

        GenericDialog gd = new GenericDialog("Fourier-Mellin Registration", IJ.getInstance());
        String defaultItem;

        if (title1.equals("")) {
            defaultItem = titles[0];
        } else {
            defaultItem = title1;
        }
        gd.addChoice("Orignal Image:", titles, defaultItem);
		if (title2.equals("")) {
            defaultItem = titles[1];
        } else {
            defaultItem = title2;
        }
        gd.addChoice("Transformed Image:", titles, defaultItem);
        gd.addNumericField("Highpass filter threshold ", cutoffHigh, 0);   
        gd.addNumericField("Lowpass filter threshold ", cutoffLow, 0);  
        gd.addCheckbox("Show correlation map image?", showR);
		gd.addCheckbox("subpixel?", subpixel);
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
        cutoffHigh = (int)gd.getNextNumber();
        cutoffLow = (int)gd.getNextNumber();
		showR = gd.getNextBoolean();
		subpixel = gd.getNextBoolean();

        return true;
    }

	public static FloatProcessor HLpassFiltering(ImageProcessor ip, int cHigh, int cLow){
			
		/*get filter size*/
		int size = Math.max(ip.getWidth(), ip.getHeight());
		
		/* Calculate filter */
		FloatProcessor ipFilter = calcFilter(size,cHigh, cLow);
			
		/*Do filtering*/
		FHT fht = new FHT(ip);
		fht.transform();	// calculates the Fourier transformation			
		fht.swapQuadrants(ipFilter);
		float[] pixels_fht = (float[])fht.getPixels();
		for (int i=0; i<pixels_fht.length; i++) {
			pixels_fht[i] *= ipFilter.getf(i);
		}
		/* get the magnitude spectrum */
		FloatProcessor amp = calculateAmplitude(pixels_fht, size);
		fht.swapQuadrants(amp);
		

		return amp;	
	}
	
	/* Pad image to have power of 2 dimension */
	private static ImageProcessor getPaddedIP(ImagePlus img){
		
		return getPaddedIP(img.getProcessor());
	}
	
	/* Pad image to have power of 2 dimension */
	private static ImageProcessor getPaddedIP(ImageProcessor ip){
		
		int M = ip.getWidth();
		int N = ip.getHeight();
		int maxN = Math.max(M, N);

		int size = 2;
		while (size<maxN) size *= 2;   
		if (size == maxN && M == N) {   //padding not required
            return ip;
        }
		ImageStatistics stats = ip.getStats();
		ImageProcessor ip2 = ip.createProcessor(size, size);  // processor of the padded image
		ip2.setValue(stats.mean);
		ip2.fill();
		ip2.insert(ip, (int)Math.round( (size - M) / 2.0 ), (int)Math.round( (size - N) / 2.0 ));	
		
		return ip2;	
	} 
	
	private static FloatProcessor calcFilter(int N, int cHigh, int cLow){
				
		FloatProcessor filterH = Gaussian(N, cHigh, true);
		FloatProcessor filterL = Gaussian(N, cLow, false);

		float[] filter = new float[N*N];
		for (int i=0; i<N*N; i++) {
			filter[i] = filterL.getf(i) - (1f-filterH.getf(i));
		}

		FloatProcessor ipFilter = new FloatProcessor(N, N, filter);
	
		return ipFilter;
	}
	
	private static FloatProcessor Gaussian(int N, int cutoff, boolean high){
		FloatProcessor fp = new FloatProcessor(N,N);
        double distance = 0;
		double value = 0;
        int xcenter = (N/2);
        int ycenter = (N/2);

        for (int y = 0; y < N; y++) {
            for (int x = 0; x < N; x++) {
                distance = Math.abs(x-xcenter)*Math.abs(x-xcenter)+Math.abs(y-ycenter)*Math.abs(y-ycenter);
                distance = Math.sqrt(distance);
                value = Math.exp((-1*distance*distance)/(2*cutoff*cutoff));
                if(high) fp.putPixelValue(x,y,1f-value);
				else fp.putPixelValue(x,y,value);
            }
        }

        return fp;	
	}
	
	private static FloatProcessor calculateAmplitude(float[] fhtArray, int maxN) {
        float[] amp = new float[maxN*maxN];
        for (int row=0; row<maxN; row++) {
            amplitude(row, maxN, fhtArray, amp);
        }
        FloatProcessor fp = new FloatProcessor(maxN, maxN, amp, null);

        return fp;
	} 
	
	private static void amplitude(int row, int maxN, float[] fhtArray, float[] amplitude) {
        int base = row*maxN;
        int l;
        for (int c=0; c<maxN; c++) {
            l = ((maxN-row)%maxN) * maxN + (maxN-c)%maxN;
			
			//The DFT output Xk has a real part (Hk + HN-k)/2 and an imaginary part (HN-k - Hk)/2. To calculate amplitude/magnitude by sqrt(re^2+im^2), there is a factor of sqrt(2) to take into //consideration
            amplitude[base+c] = (float)(Math.sqrt(fhtArray[base+c]*fhtArray[base+c] + fhtArray[l]*fhtArray[l]) / Math.sqrt(2));  
        }
    } 
	
	/* log polar transformation using javacv library */
	private static FloatProcessor logPolarT(ImageProcessor ip){
		
		int w = ip.getWidth();
		int h = ip.getHeight();
		IplImage ipl_src,ipl_LP; 
		FloatBuffer fb;
		float[] f;
		
		CvMat srcMat = img32toMat(ip, w, h);
		ipl_src = srcMat.asIplImage();
		ipl_LP = cvCreateImage(cvSize(w, h), 32, 1); 
		
		//maxR = ((srcW>srcH)? srcH/2 : srcW/2);
		//maxR = 0.5* Math.sqrt(M*M + N*N);
		float maxR = 0.5f* w* 1.41421f;
		magLP = (float)(w/Math.log(maxR));
		//IJ.log("maxR: "+maxR+", mag: "+mag); 
			
		cvLogPolar(ipl_src, ipl_LP, cvPoint2D32f(w/2, h/2), magLP, CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS);
		fb = FloatBuffer.allocate(ipl_LP.width() * ipl_LP.height());
		fb = ipl_LP.createBuffer();
		f = new float[fb.capacity()];
		fb.get(f);
		FloatProcessor m = new FloatProcessor(ipl_LP.width(), ipl_LP.height(), f, null); 
		
		return m;
	
	}

	public static FloatProcessor xcorr(ImageProcessor ip1, ImageProcessor ip2){
		
		FHT F1, F2, resultF;
		F1 = new FHT(ip1);
		F2 = new FHT(ip2);
		F1.transform();
		F2.transform();

		resultF = F1.conjugateMultiply(F2);
		resultF.inverseTransform();
		resultF.swapQuadrants();
		
		return (FloatProcessor)resultF;
	}
	
	/*
	*  return maximum points at the correlation map (substracted with midpoint)
	*/
	private static int[] findMax(ImageProcessor ip) {

        ResultsTable rt = ResultsTable.getResultsTable();
        rt.reset();

        MaximumFinder fd = new MaximumFinder();

        double sd = ip.getStatistics().stdDev;
        double mean = ip.getStatistics().mean;
		double tolerance = sd;
        int count = 0;
		int mid = ip.getWidth()/2;  //assume xcorr result has equal width and height, default condition when using IJ's FHT

        while (rt.getCounter() < 2 && count <10) {
            rt.reset();
            fd.findMaxima(ip, tolerance, ImageProcessor.NO_THRESHOLD, 4, true, false);   //output type is list, exclude edge, not EDM
            tolerance = tolerance / 2;
            count++;
        }
        int[] coord = new int[4];
        if (rt.getCounter() == 1) {
            coord[0] = (int) rt.getValue("X", 0)-mid;
            coord[1] = (int) rt.getValue("Y", 0)-mid;
            coord[2] = (int) rt.getValue("X", 0)-mid;
            coord[3] = (int) rt.getValue("Y", 0)-mid;
        } else if (rt.getCounter() > 1) {
            coord[0] = (int) rt.getValue("X", 0)-mid;
            coord[1] = (int) rt.getValue("Y", 0)-mid;
            coord[2] = (int) rt.getValue("X", 1)-mid;
            coord[3] = (int) rt.getValue("Y", 1)-mid;
        } else { // no maxima found
            coord[0] = -999;
            coord[1] = -999;
            coord[2] = -999;
            coord[3] = -999;
        }

        return coord;
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
	
	/* Convert 32bit input image to CvMat */
	private static CvMat img32toMat(ImageProcessor ip32, int W, int H) {
		CvMat mat = CvMat.create(H, W, CV_32FC1);
		double[] dArr1 = float2DtoDouble1DArray(ip32.getFloatArray(), W, H);
		mat.put(0, dArr1, 0, dArr1.length);
		
		return mat;
	}
	
    public void showAbout() {
        IJ.showMessage("Registration by Fourier-Mellin Transformation", "FMT");
    }
}
