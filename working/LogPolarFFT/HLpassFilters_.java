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
import java.lang.String.*;
import javax.swing.*;
import ij.measure.Calibration;
import ij.gui.*;

public class HLpassFilters_ implements PlugInFilter {
    private int threshold, order;
    private int M, N, size, w, h;
    private ImagePlus imp;
    private FHT fht;
    private ImageProcessor mask, ipFilter;
    private String filter;
    private boolean displayFilter = true, high = true;
	private int testCase;

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
        int dim = 0;
        M = ip.getWidth();
        N = ip.getHeight();
        if (M!=N) dim = (int)(Math.min(M,N)/2);
        else dim = M/2;
        threshold = 10;
        order = 1;
        String[] choices = {"Ideal","Butterworth","Gaussian"};
		String[] tests = {"0","1","2"};

        GenericDialog gd = new GenericDialog("Filters");
        gd.addChoice("Frequency Filters: ",choices, "Gaussian");
        gd.addCheckbox("do highpass?", high);
		gd.addChoice("testCase: ",tests, "0");
        gd.showDialog();
        if (gd.wasCanceled())
            return false;
        int choiceIndex = gd.getNextChoiceIndex();
        filter = choices[choiceIndex];
        high = gd.getNextBoolean();
		testCase = gd.getNextChoiceIndex();

        GenericDialog gd2 = new GenericDialog("Filter Parameters");
        gd2.addNumericField("cut-off frequency:", threshold, 0);
        if (filter.equals("Butterworth")) gd2.addNumericField("Order:", order, 0);
        gd2.addCheckbox("Display Filter", displayFilter);
        gd2.showDialog();
        if (gd2.wasCanceled())return false;
        if (gd2.invalidNumber()) {
            IJ.error("Error", "Invalid input number");
            return false;
        }
        threshold = (int) gd2.getNextNumber();
        if (filter.equals("Butterworth")) order = (int) gd2.getNextNumber();
        displayFilter = gd2.getNextBoolean();
        if (threshold>=0 && threshold<=dim)
            return true;
        else {
            GenericDialog gd3;
            boolean flag = true;
            while (flag) {
                threshold = 20;
                JOptionPane.showMessageDialog(null,"error, cut-off frequency must belong to [" + 0 + "," + dim + "]");
                gd3 = new GenericDialog(" Threshold ");
                gd3.addNumericField("Threshold Factor:", threshold, 0);
                gd3.showDialog();
                if (gd3.wasCanceled() || gd3.invalidNumber())
                    return false;
                else {
                    threshold = (int) gd3.getNextNumber();
                    if (threshold>=0 && threshold<=dim)
                        flag = false;
                }
            }
        }
        return true;
    }

    // shows the power spectrum and filters the image
    public void filtering(ImageProcessor ip, ImagePlus imp) {
        int maxN = Math.max(M, N);
        ImageProcessor ip1 = ip;
		size = 2;
        while (size<maxN) size *= 2;   
        w = Math.round((size-M)/2);
		h = Math.round((size-N)/2);
		
		ImageStatistics stats = ImageStatistics.getStatistics(ip, Measurements.MEAN, null);
		ImageProcessor ip2 = ip.createProcessor(size, size);  // processor of the padded image
		ip2.setValue(stats.mean);
		ip2.fill();
        ip2.insert(ip, w, h);
        if (ip instanceof ColorProcessor) {
            ImageProcessor bright = ((ColorProcessor)ip).getBrightness();
            fht = new FHT(bright);
            fht.rgb = (ColorProcessor)ip.duplicate(); // get a duplication of brightness in order to add it after filtering
        } else  fht = new FHT(ip2);

        fht.originalColorModel = ip.getColorModel();
        fht.originalBitDepth = imp.getBitDepth();
        fht.transform();	// calculates the Fourier transformation

        if (filter.equals("Ideal"))        ipFilter = Ideal();
        if (filter.equals("Butterworth"))  ipFilter = Butterworth(order);
        if (filter.equals("Gaussian"))     ipFilter = Gaussian();

        fht.swapQuadrants(ipFilter);
        //if (displayFilter) new ImagePlus(filter, ip2).show();
        byte[] pixels_filter = (byte[])ipFilter.getPixels();
        float[] pixels_fht = (float[])fht.getPixels();

		switch(testCase){
			case 0:
				for (int i=0; i<size*size; i++) {
					pixels_fht[i] = (float)(pixels_fht[i]*(pixels_filter[i]&255)/255.0);
				}
				break;
				
			case 1:
				for (int i=0; i<pixels_fht.length; i++) {
					pixels_fht[i] = (float)(pixels_fht[i]*(ipFilter.get(i)/255.0));
				}
				break;
			
			case 2:
				break;
		}
        //mask = fht.getPowerSpectrum();
        //ImagePlus imp2 = new ImagePlus("inverse FFT of "+imp.getTitle(), mask);
        //imp2.setProperty("FHT", fht);
        //imp2.setCalibration(imp.getCalibration());
		
		ImageProcessor ipA = calculateAmplitude(pixels_fht, size);
		ipA.setRoi(w,h,M,N);
		ipA = ipA.crop();
		new ImagePlus("filtered mag of"+imp.getTitle(), ipA).show();
		
		doInverseTransform(fht, ip);
		ip.resetMinAndMax();
        if (displayFilter) {
			fht.swapQuadrants(ipFilter);
			new ImagePlus(filter, ipFilter).show();
		}
			
    }

    // creates an ideal lowpass filter
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
		
/*      String title = imp.getTitle();
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
			
			/* The DFT output Xk has a real part (Hk + HN-k)/2 and an imaginary part (HN-k - Hk)/2. To calculate amplitude/magnitude by sqrt(re^2+im^2), there is a factor of sqrt(2) to take into consideration */
            amplitude[base+c] = (float)(Math.sqrt(fhtArray[base+c]*fhtArray[base+c] + fhtArray[l]*fhtArray[l]) / Math.sqrt(2));  
        }
    }
	
}

