package PIV;

import ij.*;
import ij.process.*;
import ij.gui.*;
import ij.plugin.filter.*;
import ij.measure.ResultsTable;
import ij.io.*;
import java.util.Arrays;
import java.io.*;
import java.text.*;
import org.opensourcephysics.display2d.GridPointData;
import java.util.Locale;


public class iterative_PIV implements PlugInFilter {

    String arg;
    ImagePlus imp;
    int width, height;
    String title, file, piv0Path;

    private static int winS1, vecS1, sW1, winS2, vecS2, sW2, winS3, vecS3, sW3;
    int nPass = 3;
    private static double cThr;     // correlation peak threshold
    boolean db = false, batch;
    boolean pp = false, dCanceled = false, xc = true , chkPeakA = false, noChkPeak = false;
    String ppMethod;
    double pp1, pp2;
    int dbX = -1, dbY = -1;
    private static String dir = "";

    int nx, ny;			//number of vector in horizontal dir. and verticle dir.
    double[][] PIVdata1, PIVdata, PIVdata0;         //working PIV data
    double[][][] PIV0;          //PIV data from the previous iteration
    
    int action = 3;
    double noiseNMT1 = 0.2,thrNMT1 = 5, c1DMT = 3, c2DMT = 1;   //parameters for PIV postprocessing
    
    double sdR, meanR;          //pixel statistics of the correlation result image
    double max0 = 0.0d;         //estimation of the maximum displacement from either preloaded PIV0 or the first pass PIV1
 

    public int setup(String arg, ImagePlus imp) {
        this.arg = arg;
        this.imp = imp;
        return DOES_8G + DOES_16 + STACK_REQUIRED;
    }

    public void run(ImageProcessor ip) {

        int size = imp.getImageStackSize();
        title = imp.getTitle();
        width = imp.getWidth();
        height = imp.getHeight();
        int[] dim;
        String sf = "_PIV1";
        StringBuffer sb;


        if (size != 2) {
            IJ.error("2 slices stack is required");
        }

        if(arg.equals("Cross-correlation")){
            if (!getParamsC()) {
                imp.changes = false;
                return;
            }
        }else if(arg.equals("Basic")){
            if (!getParamsB()) {
                imp.changes = false;
                return;
            }
        }else if(arg.equals("Debug")){
            if (!getParamsD()) {
                imp.changes = false;
                return;
            }
        }else{
            if (!getParamsA()) {
                imp.changes = false;
                return;
            }
        }

        /*Log PIV parameters*/
        IJ.log("PIV paramters: ");
        IJ.log("pass1: Interrogation window="+winS1+" search window="+sW1+" vector spacing="+vecS1);
        IJ.log("pass2: Interrogation window="+winS2+" search window="+sW2+" vector spacing="+vecS2);
        IJ.log("pass3: Interrogation window="+winS3+" search window="+sW3+" vector spacing="+vecS3);
        if(noChkPeak){
            IJ.log("Peak check disabled");
        }else if(chkPeakA){
            IJ.log("Using emperical parameters for peak check");
        }

        /*Start the PIV iteration*/
        for (int np= 1; np <=nPass ; np++) {

            int winS=winS1,
                vecS=vecS1,
                sW=sW1;

            if (np==1){
                if(piv0Path!=null){
                    try {
                        PIVdata0=plot_.load2DArrayFromFile(piv0Path);
                    } catch (Exception e) {
                        IJ.error(e.getMessage());
                    }
                    plotPIV(PIVdata0, title + "_PIV0", false);
                }
            }else if(np==2){
                PIVdata0 = PIVdata;
                winS=winS2;
                vecS=vecS2;
                sW=sW2;
                sf = "_PIV2";
            }else if(np==3){
                PIVdata0 = PIVdata;
                winS=winS3;
                vecS=vecS3;
                sW=sW3;
                sf = "_PIV3";
            }

            //The PIV0 is the PIV data from the previous iteration except for the first iteration
            if(PIVdata0==null){
                PIV0= new double[1][1][1];
            }else{
                dim = plot_.getDimensions(PIVdata0);
                PIV0=plot_.convert2DPivTo3D(PIVdata0, dim[0], dim[1]);
            }

            /*The main PIV function*/
            PIVdata = doPIV(imp, winS, vecS, sW, PIV0);
            //sb = generatePIVToPrint(PIVdata);
            //write2File(dir, title +"_PIVtemp0.txt", sb.toString());
            //IJ.log("before");
            //logPIV(PIVdata);
            /*interpolate the invalid vector by median*/
            if(!pp){
                PIVdata = replaceByMedian(PIVdata);
                PIVdata = replaceByMedian2(PIVdata);
            }

            /*plot the PIV*/
            //sb = generatePIVToPrint(PIVdata);
            //write2File(dir, title +"_PIVtemp1.txt", sb.toString());
            //IJ.log("after");
            //logPIV(PIVdata);
            plotPIV(PIVdata, title + sf, false);

            if (db) {
                IJ.log(""+title+sf+":");
                logPIV(PIVdata);
            }
            if (batch) {
                sb = generatePIVToPrint(PIVdata);
                write2File(dir, title + sf+"_disp.txt", sb.toString());
            }
        }

        /*Show the post processing dialog for the normalized median test or dynamic mean test*/
        if(!batch){
            PIVdata1 = pivPostProcess(PIVdata);
            //After post processing, close the original PIV plot, and redraw the processed PIV
            ImagePlus vPlot = WindowManager.getImage(title + sf);
            if (vPlot != null) {
                vPlot.close();
            }
            plotPIV(PIVdata1, title + sf, true);
            if (dCanceled) {
                IJ.log(""+title+sf+":");
                logPIV(PIVdata1);
            }
        }else if (ppMethod != "None") {  //In batch mode and post-processing != None
            PIVdata1 = pivPostProcess_batch(PIVdata);
			sb = generatePIVToPrint(PIVdata1);
            write2File(dir, title + sf+"_"+ppMethod+"_disp.txt", sb.toString());
        }else{
			//In batch mode but without post-processing --> Do Nothing
		}
        
        imp.changes = false;
        IJ.freeMemory();

    }


    private double[][] pivPostProcess(double[][] _PIV) {

        double[][] _PIVa = new double[_PIV.length][_PIV[0].length];
        for (int i = 0; i < _PIVa.length; i++) {
            System.arraycopy(_PIV[i], 0, _PIVa[i], 0, _PIV[i].length);
        }


        if (db) {
            // show a waitForUser dialog, so that we can check the vector and the log value
            WaitForUserDialog wd = new WaitForUserDialog("pause");
            wd.show();
        }

        boolean OK = false;
        do {
            //show the post-process dialog
            if (!getParamsP()) {
                dCanceled = true;
                return _PIV;
            }
            ImagePlus vPlot = WindowManager.getImage(title + "_temp");
            switch (action) {
                case 0:
                    _PIVa = normalizedMedianTest(_PIVa, noiseNMT1, thrNMT1);
                    _PIVa = replaceByMedian(_PIVa);
                    break;
                case 1:
                    _PIVa = dynamicMeanTest(_PIVa, c1DMT, c2DMT);
                    _PIVa = replaceByMedian(_PIVa);
                    break;
                case 2:         //restore the original PIV
                    for (int i = 0; i < _PIVa.length; i++) {
                        System.arraycopy(_PIV[i], 0, _PIVa[i], 0, _PIV[i].length);
                    }
                    break;
                case 3:         //accept this PIV
                    //_PIVa = _PIV;
                    OK = true;
            }
            if (vPlot != null) {
                vPlot.close();
            }
            plotPIV(_PIVa, title + "_temp", false);
            if (db) {
                logPIV(_PIVa);
            }

        } while (!OK);

        ImagePlus vPlot = WindowManager.getImage(title + "_temp");
        if (vPlot != null) {
            vPlot.close();
        }

        if (db) {
            IJ.log("PIV post process:");
            logPIV(_PIVa);
        }
        if (action == 3) {
            StringBuffer sb = generatePIVToPrint(_PIVa);
            write2File(dir, file, sb.toString());
        }
        dCanceled = false;
        return _PIVa;

    }

    private double[][] pivPostProcess_batch(double[][] _PIV) {

		double[][] _PIVa = new double[_PIV.length][_PIV[0].length];
        for (int i = 0; i < _PIVa.length; i++) {
            System.arraycopy(_PIV[i], 0, _PIVa[i], 0, _PIV[i].length);
        }

		if(ppMethod=="NMT") {
			_PIVa = normalizedMedianTest(_PIVa, noiseNMT1, thrNMT1);
				_PIVa = replaceByMedian(_PIVa);
		}else{
			_PIVa = dynamicMeanTest(_PIVa, c1DMT, c2DMT);
				_PIVa = replaceByMedian(_PIVa);
		}
		
		return _PIVa;
    }
	

    private void plotPIV(double[][] PIV, String title, boolean scaleGraph) {

        ImageProcessor ip = new ColorProcessor(width, height);
        int[] dim = plot_.getDimensions(PIV);
        double[][] mag = plot_.get2DElement(PIV, dim[0], dim[1], 4);
        double max = plot_.findMax2DArray(mag);
        double sc;
        plot_.colorMax = max;
        // if the max0 was not set (in the case of first pass PIV without preloaded PIV0, or when drawing the plot for preloaded PIV0)
        if(max0==0.0d){
            sc = 24/max;
            max0 = max;
        // otherwise, use the max0 (the maximum displacement from the very first PIV pass) to set the plot scale, and use it for all iteration so that every PIV plot will have the same scale. 
        }else{
            sc = 24/max0;
            plot_.colorMax = max0;
        }
        
        plot_.loadLut("S_Pet");
        plot_.drawVectors(ip, dim, PIV, mag, sc, plot_.colors);
        ImagePlus vp = new ImagePlus(title, ip);
        vp.show();
        if(scaleGraph){
            plot_.makeScaleGraph(sc);
        }

        if (batch) {
            FileSaver fs = new FileSaver(vp);
            fs.saveAsTiff(dir + title + "_vPlot.tif");
        }
        
    }
    
    /*The dialog for PIV Advanced Mode*/
    private boolean getParamsA() {

        GenericDialog gd = new GenericDialog("Iterative PIV (Advanced Mode)");
        gd.addCheckbox("Load file as 0th pass PIV data?", false);
        gd.addMessage("(All sizes are in pixels)");
        gd.addMessage("1st pass PIV parameters:");
        if (winS1 == 0) {
            winS1 = 128;
        }
        gd.addNumericField("PIV1 interrogation window size", winS1, 0);
        if (sW1 == 0) {
            sW1 = 256;
        }
        gd.addMessage("(If search window size=window size, conventional xcorr will be used)");
        gd.addNumericField("SW1 :search window size", sW1, 0);
        if (vecS1 == 0) {
            vecS1 = 64;
        }
        gd.addNumericField("VS1 :Vector spacing", vecS1, 0);
        
        gd.addMessage("-----------------------");
        gd.addMessage("2nd pass PIV parameters: (set window size to zero to do only 1pass PIV)");
        if (winS2 == 0) {
            winS2 = 64;
        }
        gd.addNumericField("PIV2 interrogation window size", winS2, 0);
        if (sW2 == 0) {
            sW2 = 128;
        }
        gd.addNumericField("SW2 :search window size", sW2, 0);
        if (vecS2 == 0) {
            vecS2 = 32;
        }
        gd.addNumericField("VS2 :Vector spacing", vecS2, 0);
        gd.addMessage("-----------------------");
        gd.addMessage("3rd pass PIV parameters: (set window size to zero to do only 2pass PIV)");
        if (winS3 == 0) {
            winS3 = 48;
        }
        gd.addNumericField("PIV3 interrogation window size", winS3, 0);
        if (sW3 == 0) {
            sW3 = 128;
        }
        gd.addNumericField("SW3 :search window size", sW3, 0);
        if (vecS3 == 0) {
            vecS3 = 16;
        }
        gd.addNumericField("VS3 :Vector spacing", vecS3, 0);
        gd.addMessage("-----------------------");
        if (cThr == 0.0D) {
            cThr = 0.60;
        }
        gd.addNumericField("correlation threshold", cThr, 2);
        gd.addCheckbox("Use advanced peak check? (empirical parameters)", false);
        gd.addCheckbox("Disable all peak checking?", false);
        gd.addCheckbox("Don't replace invalid vector by median?", false);
        //gd.addCheckbox("Find maximum correlation at edge?", false);
        gd.addCheckbox("batch mode?", false);
		gd.addChoice("Postprocessing", (new String[]{"None", "NMT", "DMT"}), "None");
		gd.addNumericField("Postprocessing parameter1", 0.2, 2);
        gd.addNumericField("Postprocessing parameter1", 5, 2);
        if (dir.equals("")) {
            dir = "/";
        }
        gd.addStringField("Path to save outputs", dir, 30);
        gd.showDialog();

        boolean load = gd.getNextBoolean();
        winS1 = (int) gd.getNextNumber();
        sW1 = (int) gd.getNextNumber();
        vecS1 = (int) gd.getNextNumber();

        winS2 = (int) gd.getNextNumber();
        sW2 = (int) gd.getNextNumber();
        vecS2 = (int) gd.getNextNumber();

        winS3 = (int) gd.getNextNumber();
        sW3 = (int) gd.getNextNumber();
        vecS3 = (int) gd.getNextNumber();

        cThr = (double) gd.getNextNumber();
        chkPeakA = gd.getNextBoolean();
        noChkPeak = gd.getNextBoolean();
        pp = gd.getNextBoolean();
        //edgeUser = gd.getNextBoolean();
        batch = gd.getNextBoolean();
        ppMethod = gd.getNextChoice();
		pp1 = (double) gd.getNextNumber();
		pp2 = (double) gd.getNextNumber();
		if (ppMethod=="NMT"){
			noiseNMT1 = pp1;
			thrNMT1 = pp2;
		}else{
			c1DMT = pp1;
			c2DMT = pp2;
		}
		dir = gd.getNextString();

        /*determine the number of PIV iteration*/
        if (vecS3 == 0 || sW3 == 0 || winS3 == 0) {
            nPass = 2;
        }
        if (vecS2 == 0 || sW2 == 0 || winS2 == 0) {
            nPass = 1;
        }

        if (!gd.wasCanceled()) {
           /*check parameters*/
            if(!checkParams()){
                IJ.error("Incompatible PIV parameters");
                return false;
            }

            if(load){
                OpenDialog od = new OpenDialog("Select the PIV data", "");
                if (od.getDirectory() == null || od.getFileName() == null){ 
                    return false;
                }
                piv0Path = od.getDirectory();
                piv0Path += od.getFileName();
            }
        }else{
            return false;
        }

        return true;
    }
	
    /*The dialog for Basic PIV mode*/
    private boolean getParamsB() {

        GenericDialog gd = new GenericDialog("Iterative PIV (Basic)");
//        gd.addCheckbox("Load file as 0th pass PIV data?", false);
        gd.addMessage("(All sizes are in pixels)");
        gd.addMessage("1st pass PIV parameters:");
        if (winS1 == 0) {
            winS1 = 128;
        }
        gd.addNumericField("PIV1 interrogation window size", winS1, 0);
        if (sW1 == 0) {
            sW1 = 256;
        }
        gd.addMessage("(If search window size=window size, conventional xcorr will be used)");
        gd.addNumericField("SW1 :search window size", sW1, 0);
        gd.addMessage("-----------------------");
        gd.addMessage("2nd pass PIV parameters: (set window size to zero to do only 1pass PIV)");
        if (winS2 == 0) {
            winS2 = 64;
        }
        gd.addNumericField("PIV2 interrogation window size", winS2, 0);
        if (sW2 == 0) {
            sW2 = 128;
        }
        gd.addNumericField("SW2 :search window size", sW2, 0);
        gd.addMessage("-----------------------");
        gd.addMessage("3rd pass PIV parameters: (set window size to zero to do only 2pass PIV)");
        if (winS3 == 0) {
            winS3 = 32;
        }
        gd.addNumericField("PIV3 interrogation window size", winS3, 0);
        if (sW3 == 0) {
            sW3 = 96;
        }
        gd.addNumericField("SW3 :search window size", sW3, 0);
        gd.addMessage("-----------------------");
        if (cThr == 0.0D) {
            cThr = 0.6;
        }
        gd.addNumericField("correlation threshold", cThr, 2);

        //gd.addCheckbox("debug?", false);
        gd.showDialog();

        //boolean load = gd.getNextBoolean();
        winS1 = (int) gd.getNextNumber();
        sW1 = (int) gd.getNextNumber();
        vecS1 = winS1 / 2;

        winS2 = (int) gd.getNextNumber();
        sW2 = (int) gd.getNextNumber();
        vecS2 = winS2 / 2;

        winS3 = (int) gd.getNextNumber();
        sW3 = (int) gd.getNextNumber();
        vecS3 = winS3 / 2;

        cThr = (double) gd.getNextNumber();
        //db = gd.getNextBoolean();

        if (vecS3 == 0 || sW3 == 0 || winS3 == 0) {
            nPass = 2;
        }
        if (vecS2 == 0 || sW2 == 0 || winS2 == 0) {
            nPass = 1;
        }

        /*check parameters*/
        if(!checkParams()){
            IJ.error("Incompatible PIV parameters");
            return false;
        }

        if (gd.wasCanceled()) {
            return false;
        }

        return true;
    }    
    
    /*The dialog for Cross-correlation PIV */    
    private boolean getParamsC() {

        GenericDialog gd = new GenericDialog("Iterative PIV (Cross-Correlation)");
//        gd.addCheckbox("Load file as 0th pass PIV data?", false);
        gd.addMessage("(All sizes are in pixels)");
        if (winS1 == 0) {
            winS1 = 128;
        }
        gd.addNumericField("PIV1 interrogation window size", winS1, 0);
        if (winS2 == 0) {
            winS2 = 64;
        }
        gd.addMessage("(set PIV2 window size to zero to do only 1 pass PIV)");
        gd.addNumericField("PIV2 interrogation window size", winS2, 0);
        if (winS3 == 0) {
            winS3 = 32;
        }
        gd.addMessage("(set PIV3 window size to zero to do only 2 pass PIV)");
        gd.addNumericField("PIV3 interrogation window size", winS3, 0);
        gd.showDialog();

        //boolean load = gd.getNextBoolean();
        winS1 = (int) gd.getNextNumber();
        sW1 = winS1;
        vecS1 = winS1 / 2;

        winS2 = (int) gd.getNextNumber();
        sW2 = winS2;
        vecS2 = winS2 / 2;

        winS3 = (int) gd.getNextNumber();
        sW3 = winS3;
        vecS3 = winS3 / 2;

        if (vecS3 == 0 || sW3 == 0 || winS3 == 0) {
            nPass = 2;
        }
        if (vecS2 == 0 || sW2 == 0 || winS2 == 0) {
            nPass = 1;
        }

        /*check parameters*/
        if(!checkParams()){
            IJ.error("Incompatible PIV parameters");
            return false;
        }

        if (gd.wasCanceled()) {
            return false;
        }

        return true;
    }
    
    /*The dialog for debugging PIV */
    private boolean getParamsD() {

        GenericDialog gd = new GenericDialog("Iterative PIV (Debug mode)");
        gd.addMessage("(All sizes are in pixels)");
        gd.addMessage("1st pass PIV parameters:");
        if (winS1 == 0) {
            winS1 = 128;
        }
        gd.addNumericField("PIV1 interrogation window size", winS1, 0);
        if (sW1 == 0) {
            sW1 = 256;
        }
        gd.addMessage("(If search window size=window size, conventional xcorr will be used)");
        gd.addNumericField("SW1 :search window size", sW1, 0);
        if (vecS1 == 0) {
            vecS1 = 64;
        }
        gd.addNumericField("VS1 :Vector spacing", vecS1, 0);
        
        gd.addMessage("-----------------------");
        gd.addMessage("2nd pass PIV parameters: (set window size to zero to do only 1pass PIV)");
        if (winS2 == 0) {
            winS2 = 64;
        }
        gd.addNumericField("PIV2 interrogation window size", winS2, 0);
        if (sW2 == 0) {
            sW2 = 128;
        }
        gd.addNumericField("SW2 :search window size", sW2, 0);
        if (vecS2 == 0) {
            vecS2 = 32;
        }
        gd.addNumericField("VS2 :Vector spacing", vecS2, 0);
        gd.addMessage("-----------------------");
        gd.addMessage("3rd pass PIV parameters: (set window size to zero to do only 2pass PIV)");
        if (winS3 == 0) {
            winS3 = 48;
        }
        gd.addNumericField("PIV3 interrogation window size", winS3, 0);
        if (sW3 == 0) {
            sW3 = 128;
        }
        gd.addNumericField("SW3 :search window size", sW3, 0);
        if (vecS3 == 0) {
            vecS3 = 16;
        }
        gd.addNumericField("VS3 :Vector spacing", vecS3, 0);
        gd.addMessage("-----------------------");
        if (cThr == 0.0D) {
            cThr = 0.60;
        }
        gd.addNumericField("correlation threshold", cThr, 2);
        gd.addCheckbox("Use advanced peak check? (empirical parameters)", false);
        gd.addCheckbox("Disable all peak checking?", true);
        gd.addCheckbox("Don't replace invalid vector by median?", true);
        //gd.addCheckbox("Find maximum correlation at edge?", false);
        gd.addMessage("-----------------------");
        gd.addNumericField("debug_X", -1, 0);
        gd.addNumericField("debug_Y", -1, 0);
        gd.showDialog();

        winS1 = (int) gd.getNextNumber();
        sW1 = (int) gd.getNextNumber();
        vecS1 = (int) gd.getNextNumber();

        winS2 = (int) gd.getNextNumber();
        sW2 = (int) gd.getNextNumber();
        vecS2 = (int) gd.getNextNumber();

        winS3 = (int) gd.getNextNumber();
        sW3 = (int) gd.getNextNumber();
        vecS3 = (int) gd.getNextNumber();

        cThr = (double) gd.getNextNumber();
        chkPeakA = gd.getNextBoolean();
        noChkPeak = gd.getNextBoolean();
        pp = gd.getNextBoolean();
        //edgeUser = gd.getNextBoolean();
        db = true;
        dbX = (int) gd.getNextNumber();
        dbY = (int) gd.getNextNumber();
        batch = false;

        /*determine the number of PIV iteration*/
        if (vecS3 == 0 || sW3 == 0 || winS3 == 0) {
            nPass = 2;
        }
        if (vecS2 == 0 || sW2 == 0 || winS2 == 0) {
            nPass = 1;
        }

        if (!gd.wasCanceled()) {
           /*check parameters*/
            if(!checkParams()){
                IJ.error("Incompatible PIV parameters");
                return false;
            }
        }else{
            return false;
        }

        return true;
    }

    /*The dialog for PIV post-processing*/
    private boolean getParamsP() {

        String[] actions = {"Normalized median test and replace invalid by median", "Dynamic mean test and replace invalid by median", "Restore unprocessed PIV", "Accept this PIV and output"};

        GenericDialog gd = new GenericDialog("PIV post-processing");
        gd.addChoice("What to do?", actions, actions[action]);
        gd.addMessage("-----------------------");
        gd.addMessage("Normalized median test parameters:");
        gd.addNumericField("noise for NMT", noiseNMT1, 2);
        gd.addNumericField("Threshold for NMT", thrNMT1, 2);
        gd.addMessage("-----------------------");
        gd.addMessage("Dynamic mean test parameters:");
        gd.addNumericField("C1 for DMT", c1DMT, 2);
        gd.addNumericField("C2 for DMT", c2DMT, 2);
        gd.addMessage("Dynamic threshold = C1+C2*(StdDev within the surrounding 3x3 vectors)");
        gd.showDialog();

        action = gd.getNextChoiceIndex();
        noiseNMT1 = gd.getNextNumber();
        thrNMT1 = gd.getNextNumber();
        c1DMT = gd.getNextNumber();
        c2DMT = gd.getNextNumber();

        if (gd.wasCanceled()) {
            return false;
        }
        if (action == 3) {
            SaveDialog sd = new SaveDialog("Save PIVdata", IJ.getDirectory("home"), "PIV_" + imp.getTitle(), ".txt");
            if (sd.getDirectory() == null || sd.getFileName() == null) {
                return false;
            }
            dir = sd.getDirectory();
            file = sd.getFileName();
        }

        return true;
    }

    private boolean checkParams(){

        //first check if we are doing NCC or CC
        if(winS1 == sW1){
            xc = true;
        }else{
            xc = false;
        }
        if(xc == true && !powerOf2Size(winS1)){
                IJ.error("PIV using conventional cross-correlation need the window size to be power of 2");
                return false;
        }
        if(winS1>sW1){
            IJ.error("Search window must be larger than interrogation window");
            return false;
        }
        if(vecS1>winS1){
            IJ.error("PIV vector spacing must be smaller or equal to interrogation window size");
            return false;
        }
        if (nPass != 1) {
            if(winS2>=winS1){
                IJ.error("Interrogation window of second pass should be smaller than that of first pass");
                return false;
            }
            if(xc == true && !powerOf2Size(winS2)){
                IJ.error("PIV using conventional cross-correlation need the window size to be power of 2");
                return false;
            }
            if(winS2>sW2){
                IJ.error("Search window must be larger than interrogation window");
                return false;
            }else if(vecS2>winS2){
                IJ.error("PIV vector spacing must be smaller or equal to interrogation window size");
                return false;
            }
        }
        if (nPass == 3){
            if(winS3>=winS2){
                IJ.error("Interrogation window of third pass should be smaller than that of second pass");
                return false;
            }
            if(xc == true && !powerOf2Size(winS3)){
                IJ.error("PIV using conventional cross-correlation need the window size to be power of 2");
                return false;
            }
            if(winS3>sW3){
                IJ.error("Search window must be larger than interrogation window");
                return false;
            }else if(vecS3>winS3){
                IJ.error("PIV vector spacing must be smaller or equal to interrogation window size");
                return false;
            }
        }
		
		if(!xc) cvMatchTemplate.init();  //initialize the javacv library

        return true;
    }

    /* Copied from the ImageJ FHT code */
    private boolean powerOf2Size(int size) {
        int i=2;
        while(i<size) i *= 2;
        return i==size;
    }

    double[][] doPIV(ImagePlus imp, int winS, int vecS, int sW, double[][][] PIV1) {

        double[][] vec;
        ImageStack stack = imp.getStack();
        ImageProcessor ref = stack.getProcessor(1),
                       tar = stack.getProcessor(2);
        FloatProcessor result;
        ImageProcessor result16;
        int[] dxdy = new int[6];
        double[] dxdyG = new double[2];
        double[] dxdyG2 = new double[2];
        double dx0 = 0d, dy0 = 0d, dx1 = 0d, dy1 = 0d, dx2 = 0d, dy2 = 0d;
        double mag0 = 0d, mag1 = 0d, mag2 = 0d;
        boolean edge = true;
        boolean firstPass = PIV1.length == 1;	//if the prePIV	has only one entry, we are in the first pass
        int xOri, yOri;
        int shiftX = 0, shiftY = 0;
        int refXStart, refYStart, xStart, yStart;
        int border = winS / 4;        //To avoid erroneous vector at the border
        int invCount = 0;
        int thrCount = 0;
        boolean bkResult = false;

        /*Determine how many vectors in x and y*/
        nx = (int) Math.floor((width - (border * 2) - winS) / vecS) + 1;
        ny = (int) Math.floor((height - (border * 2) - winS) / vecS) + 1;
        vec = new double[nx * ny][16];      //the array to store the piv result

        if (db) {
            IJ.log("nx=" + nx);
            IJ.log("ny=" + ny);
        }

        /*Start looping over all the interrogation windows*/
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                IJ.showProgress((j * nx) + i, nx * ny);

                /*vector position (x,y)*/
                vec[j * nx + i][0] = border + (winS / 2) + vecS * (i);
                vec[j * nx + i][1] = border + (winS / 2) + vecS * (j);

                /*interpolate the preshift from the previous PIV pass*/
                if (!firstPass) {
                    double[] lerpRes = lerpData(vec[j * nx + i][0], vec[j * nx + i][1], PIV1);
                    dx0 = lerpRes[0];
                    dy0 = lerpRes[1];
                    mag0 = Math.sqrt(lerpRes[0] * lerpRes[0] + lerpRes[1] * lerpRes[1]);
                    vec[j * nx + i][12] = dx0;
                    vec[j * nx + i][13] = dy0;
                    vec[j * nx + i][14] = mag0;
                }

                /*set the window pre-shift for conventional PIV  using cross-correlation*/
                if (xc) {         //
                    shiftX = (int) dx0;
                    shiftY = (int) dy0;
                } else {
                    shiftX = 0;
                    shiftY = 0;
                }

                /*define the ROI of searching window*/
                yStart = (border + j * vecS - (sW - winS) / 2) + shiftY;
                xStart = (border + i * vecS - (sW - winS) / 2) + shiftX;
                if (yStart < 0) {
                    yStart = 0;
                }
                if (xStart < 0) {
                    xStart = 0;
                }
                if (yStart + sW > height) {
                    yStart = height - sW;
                }
                if (xStart + sW > width) {
                    xStart = width - sW;
                }

//                if (border + (winS / 2) + vecS * i == dbX && border + (winS / 2) + vecS * j == dbY) {
//                    IJ.log("xStart:"+xStart+" yStart:"+yStart+" sW:"+sW);
//                }

                tar.setRoi(xStart, yStart, sW, sW);

                /*define the ROI of interrogation window*/
                refXStart = border + i * vecS;
                refYStart = border + j * vecS;
                if (refYStart + winS > (height - border)) {
                    refYStart = (height - border) - winS;
                }
                if (refXStart + winS > (width - border)) {
                    refXStart = (width - border) - winS;
                }
                ref.setRoi(refXStart, refYStart, winS, winS);

                /*check if we are matching features at the border*/
                if (refYStart == 0 || refXStart == 0 || refYStart == height - winS || refXStart == width - winS) {
                    edge = false;   //if we are at the border , don't exclude edge maximum when findmax
                } else if (sW - winS < 20) {
                    edge = false;   // if search window minus interrogation window is no larger than 10pixel, the correlation peak is very like to be found at the edge. So we should not exclude edge maximum
                } else {
                    edge = true;
                }
//                    
//                }else{
//                    edge = edgeUser;
//                }

                /*perform the cross-correlation or the normalized cross-correlation between the two window */
                if (xc) { //use conventional cross-correlation instead of template matching
                    FHT F_ref, F_tar, resultF;
                    F_ref = new FHT(ref.crop());
                    F_tar = new FHT(tar.crop());
                    F_ref.transform();
                    F_tar.transform();

                    resultF = F_tar.conjugateMultiply(F_ref);
                    resultF.inverseTransform();
                    resultF.swapQuadrants();
                    result = resultF;

                    xOri = winS / 2;
                    yOri = winS / 2;
                } else {  //use the normalized cross correlation as the matching method
                    result = cvMatchTemplate.doMatch(tar.crop(), ref.crop(), 5, false);
                    xOri = border + i * vecS - xStart;
                    yOri = border + j * vecS - yStart;
                }

                /*for debugging*/
                if (border + (winS / 2) + vecS * i == dbX && border + (winS / 2) + vecS * j == dbY) {
                    IJ.log("position: " + vec[j * nx + i][0] + "," + vec[j * nx + i][1]);
                    IJ.log("edge:" + edge);
                    //result = cv_MatchTemplate.doMatch(tar.crop(), ref.crop(), 5, true);
                    new ImagePlus("Match result", result).show();
                    new ImagePlus("tar_" + dbX + "," + dbY, tar.crop()).show();
                    new ImagePlus("ref_" + dbX + "," + dbY, ref.crop()).show();
                }

                /*if result is a blank image (no net displacement)*/
                if (result.getStatistics().stdDev == 0) {
                    vec[j * nx + i][2] = 0;				// valid displacement
                    vec[j * nx + i][3] = 0;
                    vec[j * nx + i][4] = 0;
                    vec[j * nx + i][5] = 0;
                    vec[j * nx + i][6] = 0;
                    vec[j * nx + i][7] = 0;				// record the second peak as well
                    vec[j * nx + i][8] = 0;
                    vec[j * nx + i][9] = 0;
                    vec[j * nx + i][10] = 0;
                    vec[j * nx + i][11] = 0;
                    vec[j * nx + i][15] = 0;
                } else {


                    /*Convert the 32-bit float correlation result image to 16-bit,
                     *in order to avoid negative pixel value which caused problem for gaussian peak fitting.
                     *The compromise of accuracy of subpixel peak fitting needs to be checked.
                     */
                    
                    int[] xyOri = {xOri, yOri};
                    int[] curPos = {(int) vec[j * nx + i][0], (int) vec[j * nx + i][1]};
                    
                    result16 = result.convertToShort(true);

                    if (!firstPass) {

                        double[] ang1 = new double[2];
                        double[] ang2 = new double[2];
                        double p1, p2;

                        //debug
                        if (db) {
                            IJ.log("position: " + curPos[0] + "," + curPos[1]);
                            IJ.log("dx0, dy0: " + dx0 + "," + dy0);
                            IJ.log("xyOri: " + xyOri[0] + "," + xyOri[1]);
                        }

                        // find the two most probable peaks from the correlation map
                        if (xc) {
                            dxdy = findMaxA(result16, edge);
                            if (dxdy[0] == -999 && dxdy[1] == -999 && dxdy[2] == -999 && dxdy[3] == -999) {
                                IJ.log("no maximum found at:");
                                IJ.log("position: " + curPos[0] + "," + curPos[1]);
                            }
                        } else {
                            double[] dxy = {dx0, dy0};
                            dxdy = findMaxC(result, edge, dxy, xyOri, curPos);
                        }

                        int c = 1;      // varialbes for the switch case

                        if (dxdy[0] == -999) {     //no significant peak found
                            dxdyG[0] = -999;
                            dxdyG[1] = -999;
                            dxdyG2[0] = -999;
                            dxdyG2[1] = -999;
                            c = 2;
                            p1 = 0;
                            p2 = 0;
                        } else {
                            if (dxdy[0] == dxdy[2] && dxdy[1] == dxdy[3]) {  //only one peak found
                                dxdyG = gaussianPeakFit(result16, dxdy[0], dxdy[1]);
                                if (db) {
                                    IJ.log("dxdyG[0]:" + dxdyG[0]);
                                    IJ.log("dxdyG[1]:" + dxdyG[1]);
                                }
                                dxdyG2 = dxdyG;
                                dx1 = dxdyG[0] - xOri;
                                dy1 = dxdyG[1] - yOri;
                               
                               	if (xc){
                               		dx1 += shiftX;
                               		dy1 += shiftY;
                               	}
                               
                                dx2 = dx1;
                                dy2 = dy1;
                                mag1 = Math.sqrt(dx1 * dx1 + dy1 * dy1);
                                mag2 = mag1;
                                ang1 = checkVector(dx1, dy1, mag1, dx0, dy0, mag0);
                                ang2 = ang1;
                                p1 = result.getf(dxdy[0], dxdy[1]);
                                p2 = p1;
                                if (noChkPeak) {
                                    c = checkThr(p1);
                                } else {
                                    c = checkPeakB1(p1, mag1, ang1);
                                }
                            } else {                                      // two peaks found
                                dxdyG = gaussianPeakFit(result16, dxdy[0], dxdy[1]);
                                dxdyG2 = gaussianPeakFit(result16, dxdy[2], dxdy[3]);
                                if (db) {
                                    IJ.log("dxdyG[0]:" + dxdyG[0]);
                                    IJ.log("dxdyG[1]:" + dxdyG[1]);
                                    IJ.log("dxdyG2[0]:" + dxdyG2[0]);
                                    IJ.log("dxdyG2[1]:" + dxdyG2[1]);
                                }
                                dx1 = dxdyG[0] - xOri;
                                dy1 = dxdyG[1] - yOri;
                                dx2 = dxdyG2[0] - xOri;
                                dy2 = dxdyG2[1] - yOri;
                                
                                if(xc){
                                	dx1 += shiftX;
                               		dy1 += shiftY;
                               		dx2 += shiftX;
                               		dy2 += shiftY;
                                }
                                
                                mag1 = Math.sqrt(dx1 * dx1 + dy1 * dy1);
                                mag2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
                                ang1 = checkVector(dx1, dy1, mag1, dx0, dy0, mag0);
                                ang2 = checkVector(dx2, dy2, mag2, dx0, dy0, mag0);
                                p1 = result.getf(dxdy[0], dxdy[1]);
                                p2 = result.getf(dxdy[2], dxdy[3]);
                                if (!xc) {          //peak checking only worked with NCC, since in XC the correlation value is not normalized
                                    if (noChkPeak) {
                                        c = checkThr(p1);
                                    } else if (chkPeakA) {
                                        c = checkPeakA(p1, p2, mag1, mag2, ang1, ang2);
                                    } else {
                                        c = checkPeakB(p1, p2, mag1, mag2, ang1, ang2);
                                    }
                                }
                            }

                            if (db) {
                                IJ.log("dx1: " + dx1);
                                IJ.log("dy1: " + dy1);
                                IJ.log("dx2: " + dx2);
                                IJ.log("dy2: " + dy2);
                                IJ.log("dx0: " + dx0);
                                IJ.log("dy0: " + dy0);
                                IJ.log("mag0: " + mag0);
                            }
                            if (p1 < cThr) {
                                thrCount++;
                            }

                        }
                        if (db) {
                            IJ.log("ang1: " + ang1[0]);
                            IJ.log("ang2: " + ang2[0]);
                            IJ.log("p1: " + p1);
                            IJ.log("p2: " + p2);
                            IJ.log("dL1: " + ang1[1]);
                            IJ.log("dL2: " + ang2[1]);
                            IJ.log("dL2-dL1: " + (Math.abs(ang2[1]) - Math.abs(ang1[1])));
                        	IJ.log("Choice:"+c);
                        }

                        switch (c) {
                            
                            // Use the second highest correlation peak
                            case 0:
                                vec[j * nx + i][2] = dx2;
                                vec[j * nx + i][3] = dy2;
                                vec[j * nx + i][4] = mag2;
                                vec[j * nx + i][5] = ang2[0];
                                vec[j * nx + i][6] = p2;
                                vec[j * nx + i][7] = dx1;
                                vec[j * nx + i][8] = dy1;
                                vec[j * nx + i][9] = mag1;
                                vec[j * nx + i][10] = ang1[0];
                                vec[j * nx + i][11] = p1;
                                vec[j * nx + i][15] = 21;           //flag it as changed to second peak
                                break;
                            // Use the first highest correlation peak
                            case 1:
                                vec[j * nx + i][2] = dx1;
                                vec[j * nx + i][3] = dy1;
                                vec[j * nx + i][4] = mag1;
                                vec[j * nx + i][5] = ang1[0];
                                vec[j * nx + i][6] = p1;
                                vec[j * nx + i][7] = dx2;
                                vec[j * nx + i][8] = dy2;
                                vec[j * nx + i][9] = mag2;
                                vec[j * nx + i][10] = ang2[0];
                                vec[j * nx + i][11] = p2;
                                break;
                            // Flag this vector as invalid. (To be interpolated or not)
                            case 2:
                                vec[j * nx + i][2] = dx1;
                                vec[j * nx + i][3] = dy1;
                                vec[j * nx + i][4] = mag1;
                                vec[j * nx + i][5] = ang1[0];
                                vec[j * nx + i][6] = p1;
                                vec[j * nx + i][7] = dx2;
                                vec[j * nx + i][8] = dy2;
                                vec[j * nx + i][9] = mag2;
                                vec[j * nx + i][10] = ang2[0];
                                vec[j * nx + i][11] = p2;
                                vec[j * nx + i][15] = -1;           //flag it as invalid
                                invCount++;                                
                                break;
                        }

                    } else {   // for the first PIV iteration (without information from previous PIV)
                        
                        //debug
                        if (db) {
                            IJ.log("position: " + curPos[0] + "," + curPos[1]);
                            IJ.log("dx0, dy0: " + dx0 + "," + dy0);
                            IJ.log("xyOri: " + xyOri[0] + "," + xyOri[1]);
                            
                        }
                        
                        dxdy = findMaxA(result, edge);
                       
//                        IJ.log("dxdy[0], dxdy[1]:"+dxdy[0]+","+dxdy[1]);
                        
                        if(dxdy[0]==-999 & dxdy [1] == -999){
                            dxdyG = new double[] {-999,-999};
                            dxdyG2 = new double[] {-999,-999};
                            dx1 = 0;
                            dy1 = 0;
                            dx2 = 0;
                            dy2 = 0;
                            mag1 = 0;
                            mag2 = 0;
                        }else{
                            
                            dxdyG = gaussianPeakFit(result16, dxdy[0], dxdy[1]);
                            dxdyG2 = gaussianPeakFit(result16, dxdy[2], dxdy[3]);
//                        if (Double.isNaN(dxdyG[0])) {
//                            new ImagePlus("Match result", result).show();
//                            new ImagePlus("tar_" + vec[j * nx + i][0] + "," + vec[j * nx + i][1], tar.crop()).show();
//                            new ImagePlus("ref_" + vec[j * nx + i][0] + "," + vec[j * nx + i][1], ref.crop()).show();
//                        }
                            dx1 = dxdyG[0] - xOri;
                            dy1 = dxdyG[1] - yOri;
                            dx2 = dxdyG2[0] - xOri;
                            dy2 = dxdyG2[1] - yOri;
                            mag1 = Math.sqrt(dx1 * dx1 + dy1 * dy1);
                            mag2 = Math.sqrt(dx2 * dx2 + dy2 * dy2);
//
//
//                            IJ.log("dx1, dy1: " + dx1 + "," + dy1);
//                            IJ.log("mag1: " + mag1);
                            
                        }                 

                        //Mark invalid vector by gaussianPeakFit, adopted from JPiv package (http://www.jpiv.vennemann-online.de/)
                        int rW = result.getWidth();
                        int rH = result.getHeight();
                        if (dxdy[0] == -999 && dxdy[1] == -999){    //invalid displacement due to no peak found
                            vec[j * nx + i][2] = dx1;				
                            vec[j * nx + i][3] = dy1;
                            vec[j * nx + i][4] = mag1;
                            vec[j * nx + i][5] = 0;
                            vec[j * nx + i][6] = -1;                //no peak found, set correlation value to -1
                            vec[j * nx + i][7] = dx2;				
                            vec[j * nx + i][8] = dy2;
                            vec[j * nx + i][9] = mag2;
                            vec[j * nx + i][10] = 0;                //no peak found, set correlation value to -1
                            vec[j * nx + i][11] = -1;
                            invCount++;
                        }else if (Math.abs(dxdyG[0] - dxdy[0]) < rW / 2 && Math.abs(dxdyG[1] - dxdy[1]) < rH / 2) {
                            vec[j * nx + i][2] = dx1;				// valid displacement
                            vec[j * nx + i][3] = dy1;
                            vec[j * nx + i][4] = mag1;
                            vec[j * nx + i][5] = 0;
                            vec[j * nx + i][6] = result.getf(dxdy[0], dxdy[1]);
                            vec[j * nx + i][7] = dx2;				// record the second peak as well
                            vec[j * nx + i][8] = dy2;
                            vec[j * nx + i][9] = mag2;
                            vec[j * nx + i][10] = 0;
                            vec[j * nx + i][11] = result.getf(dxdy[2], dxdy[3]);
                        } else {
                            vec[j * nx + i][2] = dx1;				// valid displacement
                            vec[j * nx + i][3] = dy1;
                            vec[j * nx + i][4] = mag1;
                            vec[j * nx + i][5] = 0;
                            vec[j * nx + i][6] = result.getf(dxdy[0], dxdy[1]);
                            vec[j * nx + i][7] = dx2;				// record the second peak as well
                            vec[j * nx + i][8] = dy2;
                            vec[j * nx + i][9] = mag2;
                            vec[j * nx + i][10] = 0;
                            vec[j * nx + i][11] = result.getf(dxdy[2], dxdy[3]);
                            vec[j * nx + i][15] = -1;                               // flag it as invalid
                            invCount++;
                        }
                    }
                }
            }
        }

        IJ.log("#interpolated vector / #total vector = " + invCount + "/" + nx * ny);
        IJ.log("#vector with corr. value lower than threshold / #total vector = " + thrCount + "/" + nx * ny);
        return vec;
    }

    /**
     * Checking the correlation peaks using empirical parameters.
     */
    private int checkPeakA(double p1, double p2, double mag1, double mag2, double[] ang1, double[] ang2) {

        int c = 1;
        double p2p = p1/p2;    //the peak to peak ratio
        double mag0 = mag1-ang1[1];

        if (p2p>1.5 && p1>meanR+(2*sdR)) {
            c = 1;        //if the first peak is distinct enough, keep it.
        }else if (p1>meanR+(3*sdR) && (ang1[0]<20 && ang2[0]<20)){  //two similar peaks with the same small deviation from previous vector, and peak1 is still distinct enough (>3sigma)
            c = 1;
        }else if (p1>meanR+(2*sdR) && (ang1[0]<20 && ang2[0]<20) && Math.abs(ang1[0]-ang2[0])<5 ){
            if (mag1 >= mag0*0.8 && mag1/mag0 < 3) {
                c = 1;
            } else if (mag2 >= mag0*0.8 && mag2/mag0 <3) {
                c = 0;
            } else {
                c = 2;
            }
        }else if (p1 < cThr) {
            c = 2;       //when the correlation value is too low and the highest peak is not distinct enough, mark this vector as invalid (due to too low correlation)
        }else if (p1 - p2 < 0.1 || p1/p2<1.2 ) {  //if we have two similar peaks
            if (ang1[0] - ang2[0] > 90 && (mag2/mag0 <3 && mag2/mag0 >0.33) ) { //peak1 deviate from the previous vector more than 120 degree than the peak2, and vector2 is no more than 3 times longer than the previous vector
                c = 0;
            }else if(ang1[0] - ang2[0] > 30 && (mag2/mag0<1.5 && mag2/mag0 >0.67 )){
                c = 0;
            }
        }else if (ang1[0] - ang2[0] > 50 && ang2[0]<6 && (mag2/mag0 <3 && mag2/mag0 >0.33) ){
            c = 0;
//            }else if (ang1[0] - ang2[0] > 30 &&)
//                // if the first peak deviate from the previous vector more than 30 degree than the 2nd peak
//            else if (Math.abs(ang2[1]) / Math.abs(ang1[1]) < 2 && Math.abs(ang2[1]) / Math.abs(ang1[1]) > 0.5) {
//                    c = 0;
//
//            }
        }
        //experimental criteria to switch peaks
//            } else if (exp && p1 - p2 < 0.1 && ang2[0] < ang1[0] && Math.abs(ang2[1]) - Math.abs(ang1[1]) < 3) { // two peaks have similar height
//                c = 0;         // use second peak
//            } else if (exp && p1 - p2 < 0.1 && ang1[0] - ang2[0] > 30 && vec[j * nx + i][12] > 3) {
//                c = 0;
//            } else if (exp && ang2[0] < 5 && ang1[0] - ang2[0] > 20 && Math.abs(ang2[1]) - Math.abs(ang1[1]) < 3 && p1 - p2 < 0.2) {
//                c = 0;
//            }
        return c;
    }

    /**
     * check peaks only by the threshold correlation value, as well as the deviation from previous vector.
     */
    private int checkPeakB(double p1, double p2, double mag1, double mag2, double[] ang1, double[] ang2) {

        int c = 1;
        //double p2p = p1/p2;    //the peak to peak ratio
        double mag0 = mag1-ang1[1];

        if(p1<cThr){
            c = 2;
        }else if(ang1[0]>30 && mag0>1 ){
            if(ang2[0]<5 && mag2/mag0<2 && mag2/mag0>0.5){
                c = 0;
            }else{
                c = 2;
            }
        }else if ( (mag1/mag0>2 || mag1/mag0<0.5) && mag0>1){
            if(ang2[0]<5 && mag2/mag0<2 && mag2/mag0>0.5){
                c = 0;
            }else{
                c = 2;
            }
        }


        return c;
    }

    /**
     * Check the single peak by the threshold correlation value, as well as the deviation from previous vector.
     */
       private int checkPeakB1(double p1, double mag1, double[] ang1) {

        int c = 1;
        //double p2p = p1/p2;    //the peak to peak ratio
        double mag0 = mag1-ang1[1];

        if (mag1/mag0>5 && mag1>1 && p1<cThr){
            c=2;
        }
        if (ang1[0] >60 && mag1>1 && p1<cThr){
            c=2;
        }


        return c;
    }

    /*Mark all the vectors with peak value lower than the threshold*/
   private int checkThr(double p1) {

        int c=1;

        if(p1<cThr){
            c=2;
        }else{
            c=1;
        }

        return c;
    }

    private void logPIV(double[][] PIVdata) {
        for (int i = 0; i < PIVdata.length; i++) {
            IJ.log("\t" + PIVdata[i][0] + "\t" + PIVdata[i][1] + "\t" + PIVdata[i][2] + "\t" + PIVdata[i][3] + "\t" + PIVdata[i][4] + "\t" + PIVdata[i][5] + "\t" + PIVdata[i][6] + "\t" + PIVdata[i][7] + "\t" + PIVdata[i][8] + "\t" + PIVdata[i][9] + "\t" + PIVdata[i][10] + "\t" + PIVdata[i][11] + "\t" + PIVdata[i][12]+ "\t" + PIVdata[i][13] + "\t" + PIVdata[i][14] + "\t" + PIVdata[i][15]);
        }
    }

    /*The peak finding function for first pass PIV or conventional cross-correlation */
    private int[] findMaxA(ImageProcessor ip, boolean edge) {

        ResultsTable rt = ResultsTable.getResultsTable();
        rt.reset();

        MaximumFinder fd = new MaximumFinder();

        double sd = ip.getStatistics().stdDev;
        double mean = ip.getStatistics().mean;
        double tolerance = sd;
        int count = 0;


        while ((rt.getCounter() < 2 && count < 5) || (rt.getCounter() == 0 && count <10)) {
            rt.reset();
            fd.findMaxima(ip, tolerance, ImageProcessor.NO_THRESHOLD, 4, edge, false);   //output type is list, exclude edge, not EDM
            tolerance = tolerance / 2;
            count++;
        }
        int[] coord = new int[4];
        if (rt.getCounter() == 1) {
            coord[0] = (int) rt.getValue("X", 0);
            coord[1] = (int) rt.getValue("Y", 0);
            coord[2] = (int) rt.getValue("X", 0);
            coord[3] = (int) rt.getValue("Y", 0);
        } else if (rt.getCounter() > 1) {
            coord[0] = (int) rt.getValue("X", 0);
            coord[1] = (int) rt.getValue("Y", 0);
            coord[2] = (int) rt.getValue("X", 1);
            coord[3] = (int) rt.getValue("Y", 1);
        } else { // no maxima found
            coord[0] = -999;
            coord[1] = -999;
            coord[2] = -999;
            coord[3] = -999;
        }

        return (coord);
    }

    /*Peak finding function using the prior knowledge from previous PIV iteration */
    private int[] findMaxC(ImageProcessor ip, boolean edge, double[] dxy, int[] x0y0, int[] curPos) {

        //boolean whole = false;
        int limitX, limitY, c;
        ResultsTable rt = ResultsTable.getResultsTable();
        rt.reset();
        MaximumFinder fd = new MaximumFinder();

        double sd = ip.getStatistics().stdDev;
        double mean = ip.getStatistics().mean;
        sdR = sd;
        meanR = mean;
        double max = ip.getStatistics().max;
        double kurt = ip.getStatistics().kurtosis;
        double skew = ip.getStatistics().skewness;
        double tolerance = sd;
        int[] coord = new int[4];
        int sW = ip.getWidth() / 2;		//searching window size

        //sW = 28;

        if (kurt > 6 || (kurt > 3 && skew > -0.05) || (max - mean) / sd > 4) {
            fd.findMaxima(ip, sd / 2, mean + 2 * sd, 4, edge, false);
            if (rt.getCounter() > 1) {
                coord[0] = (int) rt.getValue("X", 0);
                coord[1] = (int) rt.getValue("Y", 0);
                coord[2] = (int) rt.getValue("X", 1);
                coord[3] = (int) rt.getValue("Y", 1);

                double p1 = ip.getf(coord[0], coord[1]);
                double p2 = ip.getf(coord[2], coord[3]);
                double p2p = p1 - p2;
                double pRatio = p1 / p2;
                double peakH = p2p / sd;

                if (peakH > 2 || pRatio > 2) {
                    if (db) {
                        IJ.log("Z1: significant peak found. curPos = " + curPos[0] + "," + curPos[1]);
                        IJ.log("peakH: " + peakH);
                        IJ.log("p/p: " + pRatio);
                        IJ.log("X: " + coord[0]);
                        IJ.log("Y: " + coord[1]);
                    }
                    coord[2] = coord[0];
                    coord[3] = coord[1];
                    return coord;
                }

            } else if (rt.getCounter() == 1) {
                if (db) {
                    IJ.log("Z3: significant peak found. curPos = " + curPos[0] + "," + curPos[1]);
                    IJ.log("kurt= " + kurt);
                    IJ.log("skew= " + skew);
                    IJ.log("X: " + rt.getValue("X", 0));
                    IJ.log("Y: " + rt.getValue("Y", 0));
                }
                coord[0] = (int) rt.getValue("X", 0);
                coord[1] = (int) rt.getValue("Y", 0);
                coord[2] = (int) rt.getValue("X", 0);
                coord[3] = (int) rt.getValue("Y", 0);
                return coord;
            }
        }
        //first set the window preshift to the first PIV vector
        limitX = (int) (x0y0[0] - Math.round(sW / 2) + Math.round(dxy[0]));
        limitY = (int) (x0y0[1] - Math.round(sW / 2) + Math.round(dxy[1]));

        // define the further shift according to the first vector value
        double biasX, biasY;
        int moment = 4;

        if (dxy[0] < 0) {
            biasX = dxy[0] * (-1) < sW / moment ? dxy[0] : (-sW) / moment;
        } else {
            biasX = dxy[0] < sW / moment ? dxy[0] : sW / moment;
        }
        if (dxy[1] < 0) {
            biasY = dxy[1] * (-1) < sW / moment ? dxy[1] : (-sW) / moment;
        } else {
            biasY = dxy[1] < sW / moment ? dxy[1] : sW / moment;
        }

        limitX += biasX;
        limitY += biasY;

        ip.setRoi(limitX, limitY, sW, sW);
        double sd2 = ip.getStatistics().stdDev;
        rt.reset();
        fd.findMaxima(ip, sd / 20, mean, 4, edge, false);

        if (rt.getCounter() == 0) {
            // search the peak in the whole map
            ip.resetRoi();
            rt.reset();
            fd.findMaxima(ip, sd / 10, mean + 2 * sd, 4, edge, false);
            if (rt.getCounter() == 0) {
                //no significant peak found
                if (db) {
                    IJ.log("no significant peak found. curPos = " + curPos[0] + "," + curPos[1]);
                    IJ.log("limitX: " + limitX);
                    IJ.log("limitY: " + limitY);
                    IJ.log("tolerance: " + sd / 10);
                    IJ.log("threshold: " + (mean + 2 * sd));
                }
                coord[0] = -999;
                coord[1] = -999;
                coord[2] = -999;
                coord[3] = -999;
                // to be replaced by mean value later
            } else {
                // in the whole map, only record the first peak
                if (db) {

                    IJ.log("A1: one peak found in the whole map. curPos = " + curPos[0] + "," + curPos[1]);
                    IJ.log("limitX: " + limitX);
                    IJ.log("limitY: " + limitY);
                    IJ.log("tolerance: " + sd / 10);
                }
                coord[0] = (int) rt.getValue("X", 0);
                coord[1] = (int) rt.getValue("Y", 0);
                coord[2] = (int) rt.getValue("X", 0);
                coord[3] = (int) rt.getValue("Y", 0);
            }
        } else {   //result >0
            if (db) {

                IJ.log("B1: peaks found in the preshift window. curPos = " + curPos[0] + "," + curPos[1]);
                IJ.log("limitX: " + limitX);
                IJ.log("limitY: " + limitY);
                IJ.log("tolerance: " + sd / 20);

                for (int i = 0; i < rt.getCounter(); i++) {
                    if (i >= 2) {
                        break;
                    }
                    IJ.log("X: " + rt.getValue("X", i));
                    IJ.log("Y: " + rt.getValue("Y", i));
                }
            }
            coord[0] = (int) rt.getValue("X", 0);
            coord[1] = (int) rt.getValue("Y", 0);
            if (rt.getCounter() > 1) {
                coord[2] = (int) rt.getValue("X", 1);
                coord[3] = (int) rt.getValue("Y", 1);
            } else {  //only one peak was found

                coord[2] = coord[0];
                coord[3] = coord[1];
            }

            //check the peak to peak ratio
            double p1 = ip.getf(coord[0], coord[1]);
            double p2 = ip.getf(coord[2], coord[3]);
            double p2p = p1 - p2;
            double pRatio = p1 / p2;
            double peakH = p2p / sd;
            double peakH2 = p2p / sd2;
            if (db) {

                IJ.log("peakH= " + peakH);
                IJ.log("peakH2= " + peakH2);
                IJ.log("pRatio= " + pRatio);
            }
            if (p1 > 0.5 && ((peakH + peakH2 > 2 && pRatio > 2) || peakH + peakH2 > 3)) {
                if (db) {
                    IJ.log("peak1 muck significant than peak2");
                }


                coord[2] = coord[0];
                coord[3] = coord[1];
            } else if ((p1 <= 0.5 && pRatio > 2) || peakH + peakH2 > 4) {
                if (db) {
                    IJ.log("peak1 2 times higher than peak2");
                    IJ.log("peak1: " + p1);
                }

                coord[2] = coord[0];
                coord[3] = coord[1];
            }
        }
        return coord;
    }

 /** Adapted from the JPiv package (http://www.jpiv.vennemann-online.de/)
  *  See Raffel, Willert, Kompenhans; Particle Image Velocimetry; 3rd printing; S. 131 for details
  */
    private double[] gaussianPeakFit(ImageProcessor ip, int x, int y) {
        double[] coord = new double[2];
        double a = 0, b = 0, c = 0;
        // border values
        if (x == 0
                || x == ip.getWidth() - 1
                || y == 0
                || y == ip.getHeight() - 1) {
            coord[0] = x;
            coord[1] = y;
        } else {
            if (ip.getPixel(x - 1, y) != 0) {
                a = Math.log(ip.getPixel(x - 1, y));
            }
            
            if (ip.getPixel(x, y) != 0) {
                b = Math.log(ip.getPixel(x, y));
            }
            if (ip.getPixel(x + 1, y) != 0) {
                c = Math.log(ip.getPixel(x + 1, y));
            }
            coord[0] = x + (a - c) / ((2 * a) - (4 * b) + (2 * c));
            if(Double.isNaN(coord[0]) || Double.isInfinite(coord[0])){
                coord[0] = x;
            }

            if (ip.getPixel(x, y - 1) != 0) {
                a = Math.log(ip.getPixel(x , y - 1));
            }

            if (ip.getPixel(x, y + 1) != 0) {
                c = Math.log(ip.getPixel(x , y + 1));
            }
            coord[1] = y + (a - c) / ((2 * a) - (4 * b) + (2 * c));
            if(Double.isNaN(coord[1]) || Double.isInfinite(coord[1])){
                coord[1] = y;
            }

//            coord[0] = x
//                    + (Math.log(ip.getPixel(x - 1, y))
//                    - Math.log(ip.getPixel(x + 1, y)))
//                    / (2 * Math.log(ip.getPixel(x - 1, y))
//                    - 4 * Math.log(ip.getPixel(x, y))
//                    + 2 * Math.log(ip.getPixel(x + 1, y)));
//            coord[1] = y
//                    + (Math.log(ip.getPixel(x, y - 1))
//                    - Math.log(ip.getPixel(x, y + 1)))
//                    / (2 * Math.log(ip.getPixel(x, y - 1))
//                    - 4 * Math.log(ip.getPixel(x, y))
//                    + 2 * Math.log(ip.getPixel(x, y + 1)));
        }
        return (coord);
    }

    /** Modified from the JPiv package (http://www.jpiv.vennemann-online.de/)
     * Mark vectors as invalid using the normalized median test proposed by
     * Jerry Westerweel and Fulvio Scarano (Experiments in Fluids (2005) 39: 1096-1100).
     * The velocity data remains unchanged, but the value in the fifth column
     * (generally the peak height) is set to -1.
     * @param noiseLevel Noise Level of the velocity data in pixel units.
     * @param threshold Data above this threshold will be discarted (the default is 2).
     */
    double[][] normalizedMedianTest(double[][] pivData, double noiseLevel, double threshold) {

        double[] nb;
        double normResidu;
        int flag = 15;
        // loop over data points
        for (int i = 0; i < pivData.length; ++i) {
            // loop over velocity components x and y
            for (int c = 2; c < 4; c++) {
                nb = getNeighbours(pivData, i, c, flag);
                if (nb != null) {
                    normResidu = Math.abs(pivData[i][c] - getMedian(nb))
                            / (getMedian(getResidualsOfMedian(nb)) + noiseLevel);
                    if (normResidu > threshold) {
                        pivData[i][flag] = -1;
                    }
                } else {
                    pivData[i][flag] = -1;
                }
            }
        }

        return pivData;
    }

    /**
     * Implemented as described in "Raffel M. (2007). Particle image velocimetry: a practical guide. 2nd ed. Springer"
     */
    double[][] dynamicMeanTest(double[][] pivData, double c1, double c2) {

        double[] nb;
        double thr;
        double meanNb, stdNb;
        int flag = 15;
        // loop over data points
        for (int i = 0; i < pivData.length; ++i) {
            // loop over velocity components x and y
            for (int c = 2; c < 4; c++) {
                nb = getNeighbours(pivData, i, c, flag);
                if (nb != null) {
                    meanNb = getMean(nb);
                    stdNb = calcStd(nb, meanNb);
                    thr = c1 + c2 * stdNb;
                    if (Math.abs(pivData[i][c] - meanNb) > thr) {
                        pivData[i][flag] = -1;
                    }
                } else {
                    pivData[i][flag] = -1;
                }
            }
        }

        return pivData;
    }

    /** 
     * Modified from the JPiv package (http://www.jpiv.vennemann-online.de/)
     * Gets the neigbhouring values of a vector component.
     * The values are returned in ascending index order.
     * @param cnt The index of the center vector
     * @param col The vector component to return
     * @return An array of up to nine values from column <code>col</code>
     * or <code>null</null> if no neighbours could be found.
     */
    double[] getNeighbours(double[][] pivData, int cnt, int col, int _pos) {
        double[] nb = new double[9];
        int i = cnt / nx;
        int j = cnt - i * nx;
        int numOfNb = 0;
        int n = 0;
        // first row neighbours
        n = cnt - nx - 1;
        if (i - 1 >= 0
                && j - 1 >= 0
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt - nx;
        if (i - 1 >= 0
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt - nx + 1;
        if (i - 1 >= 0
                && j + 1 < nx
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        // secont row neighbours
        n = cnt - 1;
        if (j - 1 >= 0
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        
        n = cnt + 1;
        if (j + 1 < nx
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        // third row neighbours
        n = cnt + nx - 1;
        if (i + 1 < ny
                && j - 1 >= 0
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt + nx;
        if (i + 1 < ny
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt + nx + 1;
        if (i + 1 < ny
                && j + 1 < nx
                && (pivData[n][_pos] != -1 )) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        if (numOfNb > 0) {
            double[] ret = new double[numOfNb];
            System.arraycopy(nb, 0, ret, 0, numOfNb);
            return (ret);
        } else {
            return (null);
        }
    }

    /*Get the values having at least one valid neighbor */
    double[] getNeighbours2(double[][] pivData, int cnt, int col, int _pos) {
        double[] nb = new double[9];
        int i = cnt / nx;
        int j = cnt - i * nx;
        int numOfNb = 0;
        int n = 0;
        // first row neighbours
        n = cnt - nx - 1;
        if (i - 1 >= 0
                && j - 1 >= 0
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt - nx;
        if (i - 1 >= 0
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt - nx + 1;
        if (i - 1 >= 0
                && j + 1 < nx
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        // secont row neighbours
        n = cnt - 1;
        if (j - 1 >= 0
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }

        n = cnt + 1;
        if (j + 1 < nx
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        // third row neighbours
        n = cnt + nx - 1;
        if (i + 1 < ny
                && j - 1 >= 0
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt + nx;
        if (i + 1 < ny
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        n = cnt + nx + 1;
        if (i + 1 < ny
                && j + 1 < nx
                && (pivData[n][_pos] != -2)) {
            numOfNb++;
            nb[numOfNb - 1] = pivData[n][col];
        }
        if (numOfNb > 0) {
            double[] ret = new double[numOfNb];
            System.arraycopy(nb, 0, ret, 0, numOfNb);
            return (ret);
        } else {
            return (null);
        }
    }

    /** Copied from the JPiv package (http://www.jpiv.vennemann-online.de/)
     * Gets the median of an array of double values.
     * @param val The data
     * @return The median value. If the number of data elements is even the
     * artithmetic mean of the two median values is returned.
     */
    double getMedian(double[] val) {
        Arrays.sort(val);
        int mid = val.length / 2;
        // uneven
        if ((val.length % 2) > 0) {
            return (val[mid]);
        } // even
        else {
            return ((val[mid] + val[mid - 1]) / 2);
        }
    }

    double getMean(double[] val) {
        double sum = 0;
        for (int i = 0; i < val.length; i++) {
            sum += val[i];
        }
        return sum / val.length;
    }

    double calcStd(double[] val, double mean) {
        double sum = 0;
        for (int i = 0; i < val.length; i++) {
            sum += (val[i] - mean) * (val[i] - mean);
        }
        return sum / val.length;
    }

    /** Copied from the JPiv package (http://www.jpiv.vennemann-online.de/)
     * Gets the residuals r of the median x_m of an array of
     * double values x_i ( r_i = |x_i - x_m| ).
     * @param val The data x_i.
     * @return The residuals r_i of the median value.
     */
    double[] getResidualsOfMedian(double[] val) {
        double median = getMedian(val);
        for (int i = 0; i < val.length; i++) {
            val[i] = Math.abs(val[i] - median);
        }
        return val;
    }

    /** Modified from the JPiv package (http://www.jpiv.vennemann-online.de/)
     *
     * Replaces invalid vectors by the median of their surrounding vectors.
     * A vector is considered to be invalid if its fifth column (peak height)
     * is negative or null.
     * @param all If all is set <code>true</code> every data-point is replaced by
     * the median of a 3x3 array inclusive center point (median filtering).
     * @param includeInvalid set <code>false</code> to exclude invalid vectors
     * in the median calculation.
     */
    double[][] replaceByMedian(double[][] pivData) {

        //int _p = 6;
        int _flag = 15;
        double[][] newPIV = new double[pivData.length][pivData[0].length];
        // deep-copy of pivData
        for (int i = 0; i < newPIV.length; i++) {
            System.arraycopy(pivData[i], 0, newPIV[i], 0, pivData[i].length);
        }
        // median filtering
        double[] nb;
        for (int i = 0; i < newPIV.length; ++i) {
            if (newPIV[i][_flag] == -1 ) {
                for (int c= 2; c <= 3; c++) {
                    nb = getNeighbours(pivData, i, c, _flag);
                    if (nb != null) {
                        newPIV[i][c] = getMedian(nb);
                        newPIV[i][_flag] = 999;
                    } else {
                        newPIV[i][c] = 0.0d;
                        newPIV[i][_flag] = -2;
                    }
                }
                if(newPIV[i][_flag] != -2){
                    newPIV[i][4] = Math.sqrt(newPIV[i][2]*newPIV[i][2]+newPIV[i][3]*newPIV[i][3]);
                }
            }
        }

        return newPIV;
    }

    double[][] replaceByMedian2(double[][] pivData) {

        int _flag = 15;
        double[][] newPIV = new double[pivData.length][pivData[0].length];
        // deep-copy of pivData
        for (int i = 0; i < newPIV.length; i++) {
            System.arraycopy(pivData[i], 0, newPIV[i], 0, pivData[i].length);
        }
        // median filtering
        double[] nb;
        for (int i = 0; i < newPIV.length; ++i) {
            if (newPIV[i][_flag] == -2) {
                for (int c= 2; c <= 3; c++) {
                    nb = getNeighbours2(pivData, i, c, _flag);
                    if (nb != null) {
                        newPIV[i][c] = getMedian(nb);
                        newPIV[i][_flag] = 9999;
                    } else {
                        newPIV[i][c] = 0.0d;
                        newPIV[i][_flag] = -22;
                    }
                }
                if(newPIV[i][_flag] != -22){
                    newPIV[i][4] = Math.sqrt(newPIV[i][2]*newPIV[i][2]+newPIV[i][3]*newPIV[i][3]);
                }
            }
        }

        return newPIV;
    }

    /**
     *  compare vector 1 (x1,y1) with vector 0 (x0,y0), the length of the vector 0 is already given (L0)
     */
    public static double[] checkVector(double _x, double _y, double _l, double x0, double y0, double L0) {

        double[] resultCV = new double[2];

        double dot = _x * x0 + _y * y0;
        double abs_1 = _l;
        double abs_0 = L0;
        double angleRad = Math.acos(dot / (abs_1 * abs_0));
        if (Double.isNaN(angleRad) || Double.isInfinite(angleRad)) {
            resultCV[0] = 0;
        } else {
            resultCV[0] = angleRad * 180 / Math.PI;    //the difference of angle in degree
        }
        resultCV[1] = abs_1 - abs_0;                 // the difference of vector length (vector1 - vector0)

        return resultCV;
    }

    /**
     * Generate PIV data as ready to print string
     */
    private StringBuffer generatePIVToPrint(double[][] data) {
        NumberFormat nf = NumberFormat.getInstance(Locale.US);
        if (nf instanceof DecimalFormat) {
            ((DecimalFormat) nf).applyPattern("###.##;-###.##");
        }
        nf.setMaximumFractionDigits(12);
        nf.setMinimumFractionDigits(0);

        StringBuffer info = new StringBuffer();

        for (int i = 0; i < data.length; i++) {
            for (int c = 0; c < data[0].length; c++) {
                info.append(nf.format(data[i][c]));
                //info.append(data[i][c]);
                info.append(" ");
            }
            info.append("\n");
        }

        return info;
    }

    /**Copied from the ImageJ Particle Tracker plugin (https://weeman.inf.ethz.ch/ParticleTracker/)
     *
     * Writes the given info to given file information.
     * info will be written to the beginning of the file, overwriting older information
     * If the file doesn't exists it will be created.
     * Any problem creating, writing to or closing the file will generate an ImageJ error
     * @param directory location of the file to write to
     * @param file_name file name to write to
     * @param info info the write to file
     * @see java.io.FileOutputStream#FileOutputStream(java.lang.String)
     */
    public boolean write2File(String directory, String file_name, String info) {

        PrintWriter print_writer = null;

        try {
            FileOutputStream fos = new FileOutputStream(directory + file_name);
            BufferedOutputStream bos = new BufferedOutputStream(fos);
            print_writer = new PrintWriter(bos);
            print_writer.print(info);
            print_writer.close();
            return true;
        } catch (IOException e) {
            IJ.error("" + e);
            return false;
        }

    }

    /*Function to interpolate the displacement at the current position from previous PIV iteration*/
    double[] lerpData(double _x, double _y, double[][][] _data) {

        double[] result = new double[_data[0][0].length - 2];
        GridPointData gd = new GridPointData(_data.length, _data[0].length, _data[0][0].length);
        gd.setData(_data);
        int[] ind = {0, 1, 2};
        result = gd.interpolate(_x, _y, ind, result);

        return result;
    }
 
}
