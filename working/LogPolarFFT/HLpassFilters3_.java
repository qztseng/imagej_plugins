import ij.*;
import ij.process.*;
import ij.gui.*; // for generic dialog class
import ij.plugin.filter.*;   //field constant DOES_8G+DOES_16+DOES_32
import ij.plugin.filter.FFTFilter;
import ij.measure.*; //image stats.mean
import java.lang.Math.*;
//import java.awt.*;  // class rect
import java.nio.FloatBuffer;
import org.bytedeco.javacv.*;
import org.bytedeco.javacpp.*;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_imgproc.*;


public class HLpassFilters3_ implements PlugInFilter {
    static private int cutoffHigh=30, cutoffLow=50;
    private ImagePlus imp;
	private int M,N;
	private boolean filterOnly = false;

    // method from PlugInFilter Interface
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
		
		return DOES_8G+DOES_16+DOES_32;
    }

    // method from PlugInFilter Interface
    public void run(ImageProcessor ip) {
             
		M = ip.getWidth();
		N = ip.getHeight();
		
		if (!showDialog()) {
            return;
        }
        
		/* Pad image */
		ImageProcessor ip2 = getPaddedIP(ip);
		int size = ip2.getWidth(); 
	
		/* Calculate filter */
		FloatProcessor ipFilter = calcFilter(size,cutoffHigh, cutoffLow);
		
		/*Do filtering*/
		FHT fht = new FHT(ip2);
		fht.transform();	// calculates the Fourier transformation			
		fht.swapQuadrants(ipFilter);
		float[] pixels_fht = (float[])fht.getPixels();
		if(!filterOnly){
			for (int i=0; i<pixels_fht.length; i++) {
				pixels_fht[i] *= ipFilter.getf(i);
			}
		}else{
			fht.swapQuadrants(ipFilter);
			new ImagePlus("Filter of "+cutoffHigh+","+cutoffLow, ipFilter).show();	
		}
		
		ImageProcessor ipA = calculateAmplitude(pixels_fht, size);
		fht.swapQuadrants(ipA);
		if(filterOnly) new ImagePlus("mag. of "+imp.getTitle(), ipA).show();
		
		/*call the loader before calling opencv core functions*/
		Loader.load(opencv_core.class);
		
		/*Do log-polar transformation*/
		FloatProcessor m = logPolarT(ipA);
		new ImagePlus("LPT of mag. of "+imp.getTitle(), m).show();

    }

    /*
	Dialog for plugin parameters
	*/
	boolean showDialog() {
		
		int limit = Math.min(M,N)/2;
		GenericDialog gd = new GenericDialog("Filters");
		gd.addMessage(" 0< cut-off freq <"+limit);
		gd.addMessage("highpass cut-off < lowpass cut-off");
		gd.addNumericField("Highpass cut-off frequency:", cutoffHigh, 0);
		gd.addNumericField("Lowpass cut-off frequency:", cutoffLow, 0);
		gd.addCheckbox("only calculate filter?", filterOnly);
		gd.showDialog();
		cutoffHigh = (int)gd.getNextNumber();
		cutoffLow = (int)gd.getNextNumber();
		filterOnly = gd.getNextBoolean();
		
		return true;
	}

	private FloatProcessor Gaussian(int size, int cutoff, boolean high){
		FloatProcessor fp = new FloatProcessor(size,size);
        double distance = 0;
		double value = 0;
        int xcenter = (size/2);
        int ycenter = (size/2);

        for (int y = 0; y < size; y++) {
            for (int x = 0; x < size; x++) {
                distance = Math.abs(x-xcenter)*Math.abs(x-xcenter)+Math.abs(y-ycenter)*Math.abs(y-ycenter);
                distance = Math.sqrt(distance);
                value = Math.exp((-1*distance*distance)/(2*cutoff*cutoff));
                if(high) fp.putPixelValue(x,y,1f-value);
				else fp.putPixelValue(x,y,value);
            }
        }

        return fp;	
	}
		
	private FloatProcessor calcFilter(int size, int cHigh, int cLow){
				
		FloatProcessor filterH = Gaussian(size, cHigh, true);
		FloatProcessor filterL = Gaussian(size, cLow, false);

		float[] filter = new float[size*size];
		for (int i=0; i<size*size; i++) {
			filter[i] = filterL.getf(i) - (1f-filterH.getf(i));
		}

		FloatProcessor ipFilter = new FloatProcessor(size, size, filter);
	
		return ipFilter;
	}
	 
	private ImageProcessor getPaddedIP(ImageProcessor ip){
		
		int maxN = Math.max(M, N);

		int size = 2;
        while (size<maxN) size *= 2;   

		ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, null);	 
		ImageProcessor ip2 = ip.createProcessor(size, size);  // processor of the padded image
		ip2.setValue(stats.mean);
		ip2.fill();
        ip2.insert(ip, (int)Math.round( (size - M) / 2.0 ), (int)Math.round( (size - N) / 2.0 ));	
		
		return ip2;	
		
	} 
	
 	private ImageProcessor calculateAmplitude(float[] fhtArray, int maxN) {
        float[] amp = new float[maxN*maxN];
        for (int row=0; row<maxN; row++) {
            amplitude(row, maxN, fhtArray, amp);
        }
        ImageProcessor ip = new FloatProcessor(maxN, maxN, amp, null);

        return ip;
	} 
	
	private void amplitude(int row, int maxN, float[] fhtArray, float[] amplitude) {
        int base = row*maxN;
        int l;
        for (int c=0; c<maxN; c++) {
            l = ((maxN-row)%maxN) * maxN + (maxN-c)%maxN;
			
			//The DFT output Xk has a real part (Hk + HN-k)/2 and an imaginary part (HN-k - Hk)/2. To calculate amplitude/magnitude by sqrt(re^2+im^2), there is a factor of sqrt(2) to take into //consideration
            amplitude[base+c] = (float)(Math.sqrt(fhtArray[base+c]*fhtArray[base+c] + fhtArray[l]*fhtArray[l]) / Math.sqrt(2));  
        }
    } 
	
	/* Convert 32bit input image to CvMat */
	private static CvMat img32toMat(ImageProcessor ip32, int W, int H) {
		CvMat mat = CvMat.create(H, W, CV_32FC1);
		double[] dArr1 = float2DtoDouble1DArray(ip32.getFloatArray(), W, H);
		mat.put(0, dArr1, 0, dArr1.length);
		
		return mat;
	}
	
	/* The 2D array from ImageJ is [x][y], while the cvMat is in [y][x] */
    private static double[] float2DtoDouble1DArray(float[][] arr2d, int column, int row) {

        double[] arr1d = new double[column * row];
        for (int y = 0; y < row; y++) {
            for (int x = 0; x < column; x++) {
                arr1d[y * column + x] = (double) arr2d[x][y];
            }
        }

        return arr1d;
    }
	
	/* log polar transformation using javacv library */
	private FloatProcessor logPolarT(ImageProcessor ip){
		
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
		double maxR = 0.5* w* 1.41421;
		double mag = w/Math.log(maxR);
		IJ.log("maxR: "+maxR+", mag: "+mag); 
			
		cvLogPolar(ipl_src, ipl_LP, cvPoint2D32f(w/2, h/2), mag, CV_INTER_LINEAR+CV_WARP_FILL_OUTLIERS);
		fb = FloatBuffer.allocate(ipl_LP.width() * ipl_LP.height());
		fb = ipl_LP.createBuffer();
		f = new float[fb.capacity()];
		fb.get(f);
		FloatProcessor m = new FloatProcessor(ipl_LP.width(), ipl_LP.height(), f, null); 
		
		return m;
	
	}
	
 	/* Cross correlation */
	// private FloatProcessor xcorr(ImageProcessor im1, ImageProcessor im2){
		
		// FHT fht1, fht2, resultF;
		// fht1 = new FHT(im1);
		// fht2 = new FHT(im2);
		// fht1.transform();
		// fht2.transform();

		// resultF = fht1.conjugateMultiply(fht2);
		// resultF.inverseTransform();
		// resultF.swapQuadrants();
				
		// return resultF;
		
	// }
	
	// private int[] findMax(ImageProcessor ip) {

        // ResultsTable rt = ResultsTable.getResultsTable();
        // rt.reset();

        // MaximumFinder fd = new MaximumFinder();

        // double sd = ip.getStatistics().stdDev;
        // double mean = ip.getStatistics().mean;
        // double tolerance = sd;
        // int count = 0;


        // while ((rt.getCounter() < 2 && count < 5) || (rt.getCounter() == 0 && count <10)) {
            // rt.reset();
            // fd.findMaxima(ip, tolerance, ImageProcessor.NO_THRESHOLD, 4, true, false);   //output type is list, exclude edge, not EDM
            // tolerance = tolerance / 2;
            // count++;
        // }
        // int[] coord = new int[4];
        // if (rt.getCounter() == 1) {
            // coord[0] = (int) rt.getValue("X", 0);
            // coord[1] = (int) rt.getValue("Y", 0);
            // coord[2] = (int) rt.getValue("X", 0);
            // coord[3] = (int) rt.getValue("Y", 0);
        // } else if (rt.getCounter() > 1) {
            // coord[0] = (int) rt.getValue("X", 0);
            // coord[1] = (int) rt.getValue("Y", 0);
            // coord[2] = (int) rt.getValue("X", 1);
            // coord[3] = (int) rt.getValue("Y", 1);
        // } else { // no maxima found
            // coord[0] = -999;
            // coord[1] = -999;
            // coord[2] = -999;
            // coord[3] = -999;
        // }

        // return coord;
    // }

}

