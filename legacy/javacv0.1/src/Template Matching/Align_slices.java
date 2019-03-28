package TemplateMatching;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.*;
import ij.plugin.filter.*;
import ij.plugin.frame.Recorder;
import ij.measure.ResultsTable;

//import Template_Matching.*;

/* Align slices in a stack using cross-correlation. Developped for aligning fluorescent 
pattern images or beads images in TFM experiments. 
By TSENG Qingzong  2010/9/8,  qztseng /at/ gmail dot com
# modified to add marco compatibility. 2011/3/14 by Zong
 */
public class Align_slices implements PlugInFilter {

    ImagePlus imp;
    ImageProcessor ref, tar;
    ImageStack stack;
    Rectangle rect;
    Roi r;
    int method, refSlice, sArea = 0;
    int width, height;
    double disX, disY;
    double rMax, rMin;
    FloatProcessor result;
    ResultsTable rt;
    String arg;
    int windowSizeX, windowSizeY, iniX, iniY;
    boolean subPixel = false;
    boolean showRT = false;
    int itpMethod = 0;

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
            //imp.repaintWindow();
            IJ.setTool("rect");
            new WaitForUserDialog("Align slices", "Select a rectangle region as the landmark\non a reference slice").show();
            refSlice = imp.getCurrentSlice();
            r = imp.getRoi();
            if (r != null && r.isArea()) {
                rect = r.getBounds();
                if (Recorder.record) {
                    Recorder.setCommand("Align slices in stack...");
                    Recorder.recordOption("method", "" + method);
                    Recorder.recordOption("windowsizex", "" + rect.width);
                    Recorder.recordOption("windowsizey", "" + rect.height);
                    Recorder.recordOption("x0", "" + rect.x);
                    Recorder.recordOption("y0", "" + rect.y);
                    Recorder.recordOption("sWindow", "" + sArea);
                    Recorder.recordOption("subPixel", "" + subPixel);
                    Recorder.recordOption("itpMethod", "" + itpMethod);
                    Recorder.recordOption("ref.slice", "" + refSlice);
                    Recorder.recordOption("show", "" + showRT);
                }
            } else {
                IJ.showMessage("Error", "rectangular ROI needed");
                return;
            }
        }

        ref = imp.getProcessor().crop();


        if (showRT) {
            rt = new ResultsTable();
            rt.getResultsTable();
        }

        for (int i = refSlice - 1; i > 0; i--) {     //align slices before reference slice.
            alignSlices(i);
            IJ.log("Slice:" + i + " X displacement:" + disX + " Y displacement:" + disY);
            if (showRT) {
                rt.incrementCounter();
                rt.addValue("Slice", i);
                rt.addValue("dX", disX);
                rt.addValue("dY", disY);
                rt.updateResults();
            }
        }

        for (int i = refSlice + 1; i < size + 1; i++) {     //align slices after reference slice.
            alignSlices(i);
            IJ.log("Slice:" + i + " X displacement:" + disX + " Y displacement:" + disY);
            if (showRT) {
                rt.incrementCounter();
                rt.addValue("Slice", i);
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

        if (sArea != 0) {

//            IJ.log(""+rect.x);
//            IJ.log(""+rect.y);
//            IJ.log(""+rect.width);
//            IJ.log(""+rect.height);

            int xStart = rect.x - sArea;
            int yStart = rect.y - sArea;
            int sWX = rect.width + 2 * sArea;
            int sWY = rect.height + 2 * sArea;

            if (xStart < 0) {
                xStart = 0;
            }
            if (yStart < 0) {
                yStart = 0;
            }
            if (xStart + sWX > width) {
                xStart = width - sWX;
            }
            if (yStart + sWY > height) {
                yStart = height - sWY;
            }
            tar.setRoi(xStart, yStart, sWX, sWY);
            //edge = true;

        } else {
            tar.resetRoi();
        }
//
//        tar2 = tar.crop();
//        new ImagePlus("tar", tar.crop()).show();
//        new ImagePlus("ref", ref).show();
//        IJ.log("testbbbb");
        result = cvMatch_Template.doMatch(tar.crop(), ref, method, edge);
//        IJ.log("testaaa");
//        new ImagePlus("result", result).show();

        dxdy = findMax(result, 0);
        if (subPixel) {
            double[] dxdyG = new double[2];

            dxdyG = gaussianPeakFit(result, dxdy[0], dxdy[1]);
            if(sArea==0){
                disX = rect.x - dxdyG[0];
                disY = rect.y - dxdyG[1];
            }else{
                disX = sArea - dxdyG[0];
                disY = sArea - dxdyG[1];
            }
            tar.setInterpolationMethod(itpMethod);
        } else {
            if(sArea==0){
                disX = rect.x - dxdy[0];
                disY = rect.y - dxdy[1];
            }else{
                disX = sArea - dxdy[0];
                disY = sArea - dxdy[1];
            }
        }

        tar.resetRoi();
        tar.translate(disX, disY);

        return;
    }

    public static int[] findMax(ImageProcessor ip, int sW) {
        int[] coord = new int[2];
        float max = ip.getPixel(0, 0);
        int sWh, sWw;

        if (sW == 0) {
            sWh = ip.getHeight();
            sWw = ip.getWidth();
        } else {
            sWh = sW;
            sWw = sW;
        }

        for (int j = (int) (ip.getHeight() - sWh) / 2; j < (ip.getHeight() + sWh) / 2; j++) {
            for (int i = (ip.getWidth() - sWw) / 2; i < (ip.getWidth() + sWw) / 2; i++) {
                if (ip.getPixel(i, j) > max) {
                    max = ip.getPixel(i, j);
                    coord[0] = i;
                    coord[1] = j;
                }
            }
        }
        return (coord);
    }

    private double[] gaussianPeakFit(ImageProcessor ip, int x, int y) {
        double[] coord = new double[2];
        // border values
        if (x == 0
                || x == ip.getWidth() - 1
                || y == 0
                || y == ip.getHeight() - 1) {
            coord[0] = x;
            coord[1] = y;
        } else {
            coord[0] = x
                    + (Math.log(ip.getPixel(x - 1, y))
                    - Math.log(ip.getPixel(x + 1, y)))
                    / (2 * Math.log(ip.getPixel(x - 1, y))
                    - 4 * Math.log(ip.getPixel(x, y))
                    + 2 * Math.log(ip.getPixel(x + 1, y)));
            coord[1] = y
                    + (Math.log(ip.getPixel(x, y - 1))
                    - Math.log(ip.getPixel(x, y + 1)))
                    / (2 * Math.log(ip.getPixel(x, y - 1))
                    - 4 * Math.log(ip.getPixel(x, y))
                    + 2 * Math.log(ip.getPixel(x, y + 1)));
        }
        return (coord);
    }

    private void getMacroParameters(String str) {

        method = Integer.parseInt(Macro.getValue(str, "method", "5"));
        windowSizeX = Integer.parseInt(Macro.getValue(str, "windowsizex", "512"));
        windowSizeY = Integer.parseInt(Macro.getValue(str, "windowsizey", "512"));
        iniX = Integer.parseInt(Macro.getValue(str, "x0", "0"));
        iniY = Integer.parseInt(Macro.getValue(str, "y0", "0"));
        sArea = Integer.parseInt(Macro.getValue(str, "sWindow", "0"));
        subPixel = Boolean.parseBoolean(Macro.getValue(str, "subPixel", "TRUE"));
        itpMethod = Integer.parseInt(Macro.getValue(str, "itpMethod", "0"));
        refSlice = Integer.parseInt(Macro.getValue(str, "ref.slice", "1"));
        showRT = Boolean.parseBoolean(Macro.getValue(str, "show", "TRUE"));

    }

    private boolean getUserParameters() {

        String[] methods = {"Square difference", "Normalized square difference", "Cross correlation", "Normalized cross correlation", "Correlation coefficient", "Normalized correlation coefficient"};
        String[] itpMethods = {"Bilinear", "Bicubic"};

        GenericDialog gd = new GenericDialog("Align slices by cvMatchTemplate");
        gd.addMessage("Select a rectangle region as the landmark for alignment.");
        gd.addChoice("Matching method", methods, methods[5]);
        gd.addNumericField("Search area(pixels around ROI) ", 0, 0);
        gd.addMessage("(Template will be searched on the whole image if search area =0)");
        gd.addCheckbox("Subpixel registration?", subPixel);
        gd.addChoice("Interpolation method for subpixel translation", itpMethods, itpMethods[itpMethod]);
        gd.addCheckbox("show align coordinates in results table?", true);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        method = gd.getNextChoiceIndex();
        sArea = (int) gd.getNextNumber();
        subPixel = gd.getNextBoolean();
        itpMethod = gd.getNextChoiceIndex();
        showRT = gd.getNextBoolean();

        return true;



    }
}
