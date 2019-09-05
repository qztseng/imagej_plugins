package FourierMellinRegistration;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.plugin.frame.Recorder;
import ij.measure.ResultsTable;

public class Align_slices_FMR implements PlugInFilter {

    ImagePlus imp;
    ImageProcessor ref, tar;
    ImageStack stack;
    Rectangle rect;
    Roi r;
    int refSlice;
    int width, height;
	int cutoffHigh=20, cutoffLow=30;
    double rotation, disX, disY;
    ResultsTable rt;
    String arg;
    int windowSizeX, windowSizeY, iniX, iniY;
    boolean showRT = false;
    boolean refine = true;

    public int setup(String arg, ImagePlus imp) {
        this.imp = imp;
        return DOES_8G + DOES_16 + STACK_REQUIRED;
    }

    public void run(ImageProcessor ip) {

        stack = imp.getStack();
        int type = imp.getType();
        int size = stack.getSize();
        String str = Macro.getOptions();
        width = imp.getWidth();
        height = imp.getHeight();

        if (str != null) {
            getMacroParameters(str);
            imp.setSlice(refSlice);
            imp.setRoi(iniX, iniY, windowSizeX, windowSizeY);
            r = imp.getRoi();
            rect = r.getBounds();
        } else {
            if (!getUserParameters()) {
            }
            IJ.setTool("rect");
            new WaitForUserDialog("Align slices", "Select a rectangle region as the landmark\non a reference slice").show();
            refSlice = imp.getCurrentSlice();
            r = imp.getRoi();
            if (r != null && r.isArea()) {
                rect = r.getBounds();
                if (Recorder.record) {
                    Recorder.setCommand("Align slices in stack...");
                    Recorder.recordOption("windowsizex", "" + rect.width);
                    Recorder.recordOption("windowsizey", "" + rect.height);
                    Recorder.recordOption("x0", "" + rect.x);
                    Recorder.recordOption("y0", "" + rect.y);
                    Recorder.recordOption("ref.slice", "" + refSlice);
                    Recorder.recordOption("cutoffHigh", "" + cutoffHigh);
                    Recorder.recordOption("cutoffLow", "" + cutoffLow);
                    Recorder.recordOption("show", "" + showRT);
                    Recorder.recordOption("refine", "" + refine);
                }
            } else {
                IJ.showMessage("Error", "rectangular ROI needed");
                return;
            }
        }

        ref = imp.getProcessor().crop();

        if (showRT) {
            rt = new ResultsTable();
            ResultsTable.getResultsTable();
        }

        for (int i = refSlice - 1; i > 0; i--) {     //align slices before reference slice.
            alignSlices(i);
            IJ.log("Slice:" + i + " rotation:" + rotation +" X displacement:" + disX + " Y displacement:" + disY);
            if (showRT) {
                rt.incrementCounter();
                rt.addValue("Slice", i);
				rt.addValue("rot", rotation);
                rt.addValue("dX", disX);
                rt.addValue("dY", disY);
                rt.updateResults();
            }
        }

        for (int i = refSlice + 1; i < size + 1; i++) {     //align slices after reference slice.
            alignSlices(i);
            IJ.log("Slice:" + i + " rotation:" + rotation +" X displacement:" + disX + " Y displacement:" + disY);
            if (showRT) {
                rt.incrementCounter();
                rt.addValue("Slice", i);
				rt.addValue("rot", rotation);
                rt.addValue("dX", disX);
                rt.addValue("dY", disY);
                rt.updateResults();
            }
        }
        imp.updateAndDraw();
        rt.show("Results");
    }

    private void alignSlices(int slice) {

        int[] dxdy = new int[2];
        boolean edge = false;

        tar = stack.getProcessor(slice);
		tar.setRoi(r);

        //result = cvMatch_Template.doMatch(tar.crop(), ref, method, edge);
		FM_Registration.doMatch(ref, tar.crop(),cutoffHigh, cutoffLow);
		rotation = FM_Registration.rot;
		disX = FM_Registration.dx;
		disY = FM_Registration.dy;
        
		tar.resetRoi();
		tar.rotate(rotation);
        tar.translate(disX, disY);
        
        if(refine){
            tar.setRoi(r);
            
            
            
        }

    }


    private void getMacroParameters(String str) {

        windowSizeX = Integer.parseInt(Macro.getValue(str, "windowsizex", "512"));
        windowSizeY = Integer.parseInt(Macro.getValue(str, "windowsizey", "512"));
        iniX = Integer.parseInt(Macro.getValue(str, "x0", "0"));
        iniY = Integer.parseInt(Macro.getValue(str, "y0", "0"));
        refSlice = Integer.parseInt(Macro.getValue(str, "ref.slice", "1"));
		cutoffHigh = Integer.parseInt(Macro.getValue(str, "cutoffHigh", "20"));
        cutoffLow = Integer.parseInt(Macro.getValue(str, "cutoffLow", "30"));
        showRT = Boolean.parseBoolean(Macro.getValue(str, "show", "TRUE"));
        refine = Boolean.parseBoolean(Macro.getValue(str, "refine", "TRUE"));

    }

    private boolean getUserParameters() {

        GenericDialog gd = new GenericDialog("Align slices by Fourier Mellin Transformation");
        gd.addMessage("Select a rectangle region as the landmark for alignment.");
		gd.addNumericField("Highpass filter threshold ", cutoffHigh, 0);   
        gd.addNumericField("Lowpass filter threshold ", cutoffLow, 0);   
        gd.addCheckbox("show align coordinates in results table?", true);
        gd.addCheckbox("Do a refine translational registration?", true);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        cutoffHigh = (int)gd.getNextNumber();
        cutoffLow = (int)gd.getNextNumber();
		showRT = gd.getNextBoolean();
        refine = gd.getNextBoolean();
        return true;
    }
	
	
}