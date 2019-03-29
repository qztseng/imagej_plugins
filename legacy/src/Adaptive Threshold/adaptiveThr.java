
import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.gui.DialogListener;
import java.awt.image.BufferedImage;
import java.awt.*;
import java.awt.image.*;

import static com.googlecode.javacv.cpp.opencv_core.*;
import static com.googlecode.javacv.cpp.opencv_imgproc.*;


 /*
 *this plugin implement the adaptive thresholding method from OpenCV library. The interface between OpenCV and Java is
 *using Samuel Audet's JavaCV code from: http://code.google.com/p/javacv/
 *By Qingzong TSENG 2011/5/19
 *Add preview 2012/5/6
 *Using javacv, javacpp and opencv 2.3  2012/5/10
 */
 
 
 /*  How to compile:
 **  javac -d . -classpath (ImageJ folder)/ij.jar:(javacv folder)javacv.jar adaptiveThr_.java
 **
 */
 
 
public class adaptiveThr_ implements ExtendedPlugInFilter, DialogListener {

    ImagePlus img1;
    ImagePlus imgThr;
    ColorProcessor ipTemp, ipC;
    BufferedImage bi1,bi2;
    static int method = -1;
    String[] methods = {"Mean", "Weighted mean"};
    int[] mm = {CV_ADAPTIVE_THRESH_MEAN_C, CV_ADAPTIVE_THRESH_GAUSSIAN_C};
    int tType = CV_THRESH_BINARY;
    static int bSize = 15;
    double maxValue = 255;
    static double param1 = 0;
    
    private PlugInFilterRunner pfr;
    private int nPasses = 1;            // The number of passes (color channels * stack slices)
    private int pass;
    int flags = KEEP_PREVIEW;
    int[] dimensions;
    boolean ok = false;                 //whether the ok button was pressed in dialog
    boolean cancel = false;                 //whether the cancel button was pressed in dialog
    boolean output,pv;
    int maxB = 3;
    
    public int setup(String arg, ImagePlus imp) {
        this.img1 = imp;
        dimensions = img1.getDimensions();

        if(dimensions[0]>=dimensions[1]){
            maxB = dimensions[0]*2/3;
        }else{
            maxB = dimensions[1]*2/3;
        }
        if(dimensions[0]<=2 || dimensions[1]<=2){
            IJ.error("Image too small!");
            return DONE;
        }
        
        return DOES_8G;
    }

    public int showDialog(ImagePlus imp, String command, PlugInFilterRunner pfr) {
    	
    	GenericDialog gd = new GenericDialog("Adaptive Threshold", IJ.getInstance());
        String defaultItem;
    	
        ipC = (ColorProcessor)img1.getProcessor().convertToRGB();
        
    	gd.addMessage("The local threshold will be calculated by");
    	if (method == -1) {
            defaultItem = methods[0];
        } else {
            defaultItem = methods[method];
        }
    	gd.addChoice("using the", methods, defaultItem);    
    	gd.addSlider("from the pixel block with size =:", 3, maxB, bSize);
    	//gd.addMessage("(Block size has to be an odd number larger than 1. ex: 3,5,7...)");
        gd.addSlider("then substract:", -100.0, 100.0, param1);
        gd.addMessage("(Substract a smaller/negative value will result in a higher \nthreshold thus less pixels are thresholded while substract \na larger value will give a lower threshold)");
        gd.addCheckbox(" Output Mask?", false);
        gd.addPreviewCheckbox(pfr);
        gd.addDialogListener(this);
        
        gd.showDialog();
        if (gd.wasCanceled()) {
            if(imgThr!=null){
                imgThr.close();
            }
            return DONE;
        }
        IJ.register(this.getClass());   //protect static class variables (filter parameters) from garbage collection

        this.pfr = pfr;
        flags = IJ.setupDialog(imp, flags); //ask whether to process all slices of stack (if a stack)

        if(gd.wasOKed()){
            ok = true;
            run(imp.getProcessor());
            flags |= DONE;
            return flags;
        }else{
            ok = false;
        }
        if(gd.wasCanceled()){
            cancel = true;
        }else{
            cancel = false;
        }
        
        return flags;
    
    }
    
     public boolean dialogItemChanged(GenericDialog gd, AWTEvent e) {
        
        method = mm[gd.getNextChoiceIndex()];
        bSize = (int) gd.getNextNumber();
        param1 = (double)gd.getNextNumber();
        output = gd.getNextBoolean();
        pv = gd.getPreviewCheckbox().getState();
        
        if(imgThr!= null && !pv){
            imgThr.close();
            imgThr = null;
        }
        
        //check block size. Block size has to be an odd number larger than 1. ex: 3,5,7...
        if(bSize ==1 ){
        	bSize = 3;
        }else if(bSize%2 != 1){
        	bSize +=1;
        }
        
		if(gd.invalidNumber()){
			return false;
		}
        

        return true;
    }
    
    
    
    
    public void run(ImageProcessor ip) {
        
        if(bSize%2!=1){             //force bSize to be odd number, otherwise take 3 as default value to avoid crash during preview
            bSize = 3;
        }
        
        bi1 = ip.getBufferedImage();
        IplImage ipl = IplImage.createFrom(bi1);
        IplImage res = IplImage.create(ipl.width(), ipl.height(), IPL_DEPTH_8U, 1);
       
        cvAdaptiveThreshold(ipl, res, maxValue, method, tType, bSize, param1);

        bi2 = res.getBufferedImage();
        ByteProcessor resultBb = new ByteProcessor(bi2);
        
        res.release();
        ipl.release();
        
        if(pv){
            if(imgThr==null){
                imgThr = new ImagePlus("Preview",ipC.duplicate());
                imgThr.show();
            }else{
                imgThr.setProcessor(ipC.duplicate());
            }

            ipTemp = (ColorProcessor)imgThr.getProcessor();
            resultBb.setColorModel(getCM("RED"));
            ipTemp.copyBits(resultBb, 0, 0, Blitter.COPY_ZERO_TRANSPARENT);
            imgThr.updateAndDraw();
        }
        
        if(cancel){
            if(pv)imgThr.close();
        }
        if(ok){
            if(pv)imgThr.close();
            resultBb.setColorModel(getCM("WHITE"));
            
            if(output){
                new ImagePlus("Mask", resultBb).show();
            }else{
                ip.copyBits(resultBb, 0, 0, Blitter.COPY);
                if(!pv)img1.setProcessor(ip);
            }
        }
		
	if (IJ.escapePressed()){                                // interrupted by user?
            ip.reset();
    	}

        pass++;
    }
    
    private ColorModel getCM (String color){
        
        ColorModel cm;
        byte[] rLUT, gLUT, bLUT;
        
        rLUT = new byte[256]; gLUT = new byte[256]; bLUT = new byte[256];
        
        for (int i=0; i<254; i++) {
                rLUT[i] = (byte)0;
                gLUT[i] = (byte)0;
                bLUT[i] = (byte)0;
        }
        
        
        if(color.equals("RED")){
            rLUT[255] = (byte)255;
        }else if (color.equals("WHITE")){
            rLUT[255] = (byte)255;
            gLUT[255] = (byte)255;
            bLUT[255] = (byte)255;
        }else{
            rLUT[255] = (byte)0;
            gLUT[255] = (byte)0;
            bLUT[255] = (byte)0;
        }
        
        cm = new IndexColorModel(8, 256, rLUT, gLUT, bLUT);
        
        return cm;
    }
    
    
    
    /** This method is called by ImageJ to set the number of calls to run(ip)
     *  corresponding to 100% of the progress bar */
    public void setNPasses (int nPasses) {
        this.nPasses = nPasses;
        pass = 0;
    }

    private void showProgress(double percent, boolean rgb) {
        int nPasses2 = rgb?nPasses*3:nPasses;
        percent = (double)pass/nPasses2 + percent/nPasses2;
        IJ.showProgress(percent);
    }
    

}
