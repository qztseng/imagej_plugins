import ij.*;
import ij.measure.*;
import ij.process.ImageProcessor.*;
import ij.process.*;
import ij.process.FHT.*;
import ij.plugin.filter.*;
import ij.plugin.*;
import ij.plugin.FFT.*;
import ij.plugin.filter.FFTFilter;
import java.lang.Math.*;
import java.awt.*;
import ij.measure.Calibration;
import ij.gui.*;

public class HLpassFilters2_ implements PlugInFilter {
    static private int cutoffHigh=20, cutoffLow=30;
    private int M, N, size, w, h;
    private ImagePlus imp;
    private FHT fht;
    private ImageProcessor mask, ipFilter;
    private String filter;
    private boolean displayFilter = true;
	private boolean filterOnly = true;

    // method from PlugInFilter Interface
    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_ALL;
    }

    // method from PlugInFilter Interface
    public void run(ImageProcessor ip) {
        ip = imp.getProcessor();

        if (showDialog(ip)) {
            filtering(ip,imp);
        }
        IJ.showProgress(1.0);
    }

    // the following method opens a window for users
	boolean showDialog(ImageProcessor ip) {
		
		M = ip.getWidth();
        N = ip.getHeight();
		int limit = Math.min(M,N)/2;
		GenericDialog gd = new GenericDialog("Filters");
		gd.addMessage(" 0< cut-off freq <"+limit);
		gd.addMessage("highpass cut-off < lowpass cut-off");
		gd.addNumericField("Highpass cut-off frequency:", cutoffHigh, 0);
		gd.addNumericField("Lowpass cut-off frequency:", cutoffLow, 0);
		gd.addCheckbox("only calculate filter?", filterOnly);
		gd.showDialog();
		if (gd.wasCanceled()) {
            return false;
        }
		cutoffHigh = (int)gd.getNextNumber();
		cutoffLow = (int)gd.getNextNumber();
		filterOnly = gd.getNextBoolean();
		
		return true;
	}

    public void filtering(ImageProcessor ip, ImagePlus imp) {
        int maxN = Math.max(M, N);

		size = 2;
        while (size<maxN) size *= 2;   
		Rectangle fitRect = new Rectangle();
        fitRect.x = (int) Math.round( (size - M) / 2.0 );
        fitRect.y = (int) Math.round( (size - N) / 2.0 );
        fitRect.width = M;
        fitRect.height = N;
/* 
        FFTFilter fff = new FFTFilter();
        ip2 = fff.tileMirror(ip2, size, size, fitRect.x, fitRect.y);
		new ImagePlus("tileMirror", ip2).show();
		 */	
		ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, null);	 
		ImageProcessor ip2 = ip.createProcessor(size, size);  // processor of the padded image
		ip2.setValue(stats.mean);
		ip2.fill();
        ip2.insert(ip, fitRect.x, fitRect.y);	 	 
			 
		fht = new FHT(ip2);
        fht.transform();	// calculates the Fourier transformation
		
		FloatProcessor filterH = Gaussian(cutoffHigh, true);
		new ImagePlus("filterHP", filterH).show();
		FloatProcessor filterL = Gaussian(cutoffLow, false);
		new ImagePlus("filterLP", filterL).show();

		//byte[] pixels_HP = (byte[])filterH.getPixels();
		//byte[] pixels_LP = (byte[])filterL.getPixels();
		float[] filter = new float[size*size];
		for (int i=0; i<size*size; i++) {
					//float a = filterH.get(i)/255f;
					//float b = filterL.get(i)/255f;
					//float c = filterH.get(i)*filterL.get(i);
					//float d = c/255f;
/* 					if(i%128==0){
						IJ.log("i:"+i);
						IJ.log("HP:"+filterH.get(i));
						IJ.log("LP:"+filterL.get(i));
						IJ.log("a:"+a);
						IJ.log("b:"+b);
						IJ.log("c:"+c);
						IJ.log("d:"+d);
					} */
					//filter[i] = (float)Math.sqrt((filterH.get(i) * filterL.get(i)));
/* 					if(i%512==0){
						IJ.log("i:"+i);
						IJ.log("L:"+filterL.getf(i));
						IJ.log("H:"+filterH.getf(i));
					} */
/* 					float ih = 1f-filterH.getf(i);
					float ih2 = filterHi.getf(i);
					if(i%512==0){
						IJ.log("iH:"+ih);
						IJ.log("iH2:"+ih2);
					} */
					filter[i] = filterL.getf(i) - (1f-filterH.getf(i));
		}

		FloatProcessor ipFilter = new FloatProcessor(size, size, filter);
        
		//new ImagePlus("filter1", ipFilter).show();
		//filter = (float[])ipFilter.getPixels();
        //if (displayFilter) new ImagePlus(filter, ip2).show();
        //byte[] pixels_filter = (byte[])ipFilter.getPixels();
		
		//float[] fht_freq0 = pixels_fht.clone();
		//new ImagePlus("FHT", new FloatProcessor(size,size,fht_freq0)).show();
		
		float[] pixels_fht = (float[])fht.getPixels();

		if(!filterOnly){
			fht.swapQuadrants(ipFilter);
			for (int i=0; i<pixels_fht.length; i++) {
				pixels_fht[i] *= filter[i];
			}
			fht.swapQuadrants(ipFilter);
		}
		new ImagePlus("Filter", new FloatProcessor(size,size,filter)).show();
		
		//float[] fht_freq = pixels_fht.clone();
		//new ImagePlus("FHT2", new FloatProcessor(size,size,fht_freq)).show();
		  
        //mask = fht.getPowerSpectrum();
        //ImagePlus imp2 = new ImagePlus("inverse FFT of "+imp.getTitle(), mask);
        //imp2.setProperty("FHT", fht);
        //imp2.setCalibration(imp.getCalibration());
        
			ImageProcessor ipA = calculateAmplitude(pixels_fht, size);
			ipA.setRoi(fitRect);
			ipA = ipA.crop();
			new ImagePlus("filtered mag. of "+imp.getTitle(), ipA).show();
		
		//fht.inverseTransform();
		//ip2 = fht;
		//fht.setRoi(fitRect);
        //ip2 = fht.crop();
		//int bitDepth = imp.getBitDepth(); 
        /* switch (bitDepth) {
            case 8: ip2 = ip2.convertToByte(false); break;
            case 16: ip2 = ip2.convertToShort(false); break;
            case 24:
                ip.snapshot();
                ((ColorProcessor)ip).setBrightness((FloatProcessor)ip2);
                break;
            case 32: break;
        } */
		//new ImagePlus("ip2", ip2).show();
        // copy filtered image back into original image
/*         if (bitDepth!=24) {
            ip.snapshot();
            ip.copyBits(ip2, roiRect.x, roiRect.y, Blitter.COPY);
        }
        ip.resetMinAndMax(); */
		
		//doInverseTransform(fht, ip);
		//ip.resetMinAndMax();
/*         if (displayFilter) {
			fht.swapQuadrants(ipFilter);
			new ImagePlus("filter2", ipFilter).show();
		}  */

    }

	private FloatProcessor Gaussian(int cutoff, boolean high){
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
	
	
	/* // creates an ideal lowpass filter
    public ByteProcessor Ideal() {
        ByteProcessor ip = new ByteProcessor(M,N);
        int xcenter = M/2;
        int ycenter = N/2;
		ip.setValue(255);
		ip.fill();

        for (int radius=0; radius<threshold;radius++) {
            for (double counter = 0; counter < 10; counter = counter + 0.001) {
                double x = Math.sin(counter) * radius + xcenter;
                double y = Math.cos(counter) * radius + ycenter;
                ip.putPixel((int)x, (int)y, 0);
            }
        }

        ByteProcessor ip2 = new ByteProcessor(size,size);
        ip2.insert(ip, w, h);
        return ip2;
    }

    // creates a Butterworth lowpass filter
    public ByteProcessor Butterworth(int n) {
        ByteProcessor ip = new ByteProcessor(M,N);
        double value = 0;
        double distance = 0;
        int xcenter = (M/2)+1;
        int ycenter = (N/2)+1;

        for (int y = 0; y < N; y++) {
            for (int x = 0; x < M; x++) {
                distance = Math.abs(x-xcenter)*Math.abs(x-xcenter)+Math.abs(y-ycenter)*Math.abs(y-ycenter);
                distance = Math.sqrt(distance);
                double parz = Math.pow(distance/threshold,2*n);
                value = 255*(1/(1+parz));
                if (high) ip.putPixelValue(x,y,255-value);
				else ip.putPixelValue(x,y,value);
            }
        }

        ByteProcessor ip2 = new ByteProcessor(size,size);
        ip2.fill();
        ip2.insert(ip, w, h);
        return ip2;
    }

    // creates a gaussian highpass filter
    public ByteProcessor Gaussian() {
        ByteProcessor ip = new ByteProcessor(M,N);
        double value = 0;
        double distance = 0;
        int xcenter = (M/2);
        int ycenter = (N/2);

        for (int y = 0; y < N; y++) {
            for (int x = 0; x < M; x++) {
                distance = Math.abs(x-xcenter)*Math.abs(x-xcenter)+Math.abs(y-ycenter)*Math.abs(y-ycenter);
                distance = Math.sqrt(distance);
                value = 255*Math.exp((-1*distance*distance)/(2*threshold*threshold));
                if (high) ip.putPixelValue(x,y,255-value);
                else ip.putPixelValue(x,y,value);
            }
        }

        ByteProcessor ip2 = new ByteProcessor(size,size);
        ip2.fill();
        ip2.insert(ip, w, h);
        return ip2;
    }
 */
    // applies the inverse Fourier transform to the filtered image
     void doInverseTransform(FHT fht, ImageProcessor ip) {
        //fht = fht.getCopy();
        fht.inverseTransform();
        fht.resetMinAndMax();
        ImageProcessor ip2 = fht;
        fht.setRoi(w, h, M, N);
        ip2 = fht.crop();

        int bitDepth = fht.originalBitDepth>0?fht.originalBitDepth:imp.getBitDepth();
        switch (bitDepth) {
        case 8:
            ip2 = ip2.convertToByte(true);
            break;
        case 16:
            ip2 = ip2.convertToShort(true);
            break;
        case 24:
            if (fht.rgb==null || ip2==null) {
                IJ.error("FFT", "Unable to set brightness");
                return;
            }
            ColorProcessor rgb = (ColorProcessor)fht.rgb.duplicate();
            rgb.setBrightness((FloatProcessor)ip2);
            ip2 = rgb;
            fht.rgb = null;
            break;
        case 32:
            break;
        }
        if (bitDepth!=24 && fht.originalColorModel!=null)
            ip2.setColorModel(fht.originalColorModel);
		
		ip.insert(ip2,0,0);
		
/*         String title = imp.getTitle();
        if (title.startsWith("FFT of "))
            title = title.substring(7, title.length());
        ImagePlus imp2 = new ImagePlus("Inverse FFT of "+title, ip2);
        if (imp2.getWidth()==imp.getWidth())
            imp2.setCalibration(imp.getCalibration());
        imp2.show(); */ 
    }
	
 	public ImageProcessor calculateAmplitude(float[] fhtArray, int maxN) {
        float[] amp = new float[maxN*maxN];
        for (int row=0; row<maxN; row++) {
            amplitude(row, maxN, fhtArray, amp);
        }
        ImageProcessor ip = new FloatProcessor(maxN, maxN, amp, null);
        fht.swapQuadrants(ip);
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
	
}

