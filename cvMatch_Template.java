package TemplateMatching;

import ij.*;
import ij.process.*;
import ij.gui.*;
import java.awt.Color;
import ij.plugin.*;
import ij.plugin.filter.MaximumFinder;
import ij.measure.ResultsTable;
import java.awt.image.BufferedImage;
import java.nio.FloatBuffer;

import static com.googlecode.javacv.cpp.opencv_core.*;
import static com.googlecode.javacv.cpp.opencv_imgproc.*;


/*
 * This plugin implement the template match method from OpenCV library. The interface between OpenCV and Java is 
 * using Samuel Audet's JavaCV code from: http://code.google.com/p/javacv/
 * The detailed algorithm of each matching method can be found at opencv's documentation page:
 * http://opencv.willowgarage.com/documentation/object_detection.html?highlight=matchtemplate#cvMatchTemplate
 * It supports 8 bit and 16 bit grayscale only. 
 * 
 * 2013/05/14 Add support for 32bit grayscale
 * 
 * By TSENG Qingzong (qztseng /at/ gmail dot com) 
 */
public class cvMatch_Template implements PlugIn {

    private static String title1 = "";
    private static String title2 = "";
    private static String[] methods = {"Square difference", "Normalized square difference", "Cross correlation", "Normalized cross correlation", "Correlation coefficient", "Normalized correlation coefficient"};
    private static String[] titles;
    private static int method = -1;
    private static boolean showR = true;
    private static boolean multiM = false;
    private static boolean showRT = true;
    private static boolean log = false;
    private static double mmt = 0.1;
    private static double mmth = 0.0;
    int[] dxdy = new int[2];
    private static ImagePlus img1, img2;
    private BufferedImage bi, bi2;
    private int[] wList;

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
            IJ.error("we need two images to do template matching");
            return;
        }

        /*
         ** show Dialog and get user input parameters
         */


        if (!getParams()) {
            return;
        }

        /*
         ** check parameter compatibility
         */

        if ((method == 0 || method == 1) && multiM) {
            IJ.error("Multiple match is not compatible with the square difference method");
            return;
        }
        if (img1.getBitDepth() != img2.getBitDepth()) {
            IJ.error("Images need to have the same type (bit depth)");
            return;
        }

        /*
         ** start matching
         */


        long startTime = System.currentTimeMillis();

        FloatProcessor rFp = doMatch(img1, img2, method, showR);

        long elapsedTime = System.currentTimeMillis() - startTime;

        if (multiM) {

            MaximumFinder fd = new MaximumFinder();

            ByteProcessor fdB = fd.findMaxima(rFp, mmt, mmth, 4, false, false);

            ResultsTable rt = ResultsTable.getResultsTable();

            int tmpW = img2.getWidth();
            int tmpH = img2.getHeight();

            int npoints = rt.getCounter();

            int[] xpoints = new int[npoints];
            int[] ypoints = new int[npoints];
            Overlay ol = new Overlay();

            for (int i = 0; i < npoints; i++) {
                xpoints[i] = (int) rt.getValue("X", i);
                ypoints[i] = (int) rt.getValue("Y", i);

                Roi r = new Roi(xpoints[i], ypoints[i], tmpW, tmpH);
                r.setStrokeColor(Color.green);
                ol.add(r);
            }
            img1.setOverlay(ol);

            if (log) {
                IJ.log("Best match found at= " + xpoints[0] + "," + ypoints[0]);
            }

        } else {
            if (method == 0 || method == 1) {
                dxdy = findMin(rFp, 0);
            } else {
                dxdy = findMax(rFp, 0);
            }
            if (showRT) {
                ResultsTable rt = new ResultsTable();
                rt.getResultsTable();
                rt.reset();
                rt.incrementCounter();
                rt.addValue("X", dxdy[0]);
                rt.addValue("Y", dxdy[1]);
                rt.updateResults();
                rt.show("Results");
            }
            if (log) {
                IJ.log("match found at= " + dxdy[0] + "," + dxdy[1]);
            }
            img1.setRoi(dxdy[0], dxdy[1], img2.getWidth(), img2.getHeight());

        }

        if (showR) {
            new ImagePlus("Match result", rFp).show();
        }

        IJ.showStatus("Matching done in " + elapsedTime + "msec");
        //IJ.log("Matching done in "+elapsedTime+"msec");

        return;
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

        GenericDialog gd = new GenericDialog("MatchTemplate", IJ.getInstance());
        String defaultItem;

        if (title1.equals("")) {
            defaultItem = titles[0];
        } else {
            defaultItem = title1;
        }

        gd.addChoice("Image:", titles, defaultItem);

        if (method == -1) {
            defaultItem = methods[3];
        } else {
            defaultItem = methods[method];
        }

        gd.addChoice("Method:", methods, defaultItem);    //default method is set to normalized cross correlation

        if (title2.equals("")) {
            defaultItem = titles[1];
        } else {
            defaultItem = title2;
        }

        gd.addChoice("Template:", titles, defaultItem);
        gd.addCheckbox("Output in resultTable ?", showRT);
        gd.addCheckbox("log result?", log);
        gd.addCheckbox("Show correlation map image?", showR);
        gd.addCheckbox("Multiple match?", multiM);
        gd.addNumericField("tolerence for Multi match ", mmt, 2);   //as in imageJ built-in find maxima
        gd.addNumericField("threshold for Multi match ", mmth, 2);   //minimum value for a maxima
        gd.showDialog();
        if (gd.wasCanceled()) {
            return false;
        }
        int index1 = gd.getNextChoiceIndex();
        title1 = titles[index1];
        method = gd.getNextChoiceIndex();
        int index2 = gd.getNextChoiceIndex();
        title2 = titles[index2];
        img1 = WindowManager.getImage(wList[index1]);
        img2 = WindowManager.getImage(wList[index2]);
        showRT = gd.getNextBoolean();
        log = gd.getNextBoolean();
        showR = gd.getNextBoolean();
        multiM = gd.getNextBoolean();
        mmt = gd.getNextNumber();
        mmth = gd.getNextNumber();

        return true;
    }

    public static FloatProcessor doMatch(ImagePlus src, ImagePlus tpl, int method, boolean showR) {
        
        return doMatch(src.getProcessor(), tpl.getProcessor(), method, showR);

    }

    public static FloatProcessor doMatch(ImageProcessor src, ImageProcessor tpl, int method, boolean showR) {

        BufferedImage bi = null, bi2 = null;
        FloatProcessor resultFp = null;
        int srcW = src.getWidth();
        int srcH = src.getHeight();
        int tplW = tpl.getWidth();
        int tplH = tpl.getHeight();
        IplImage temp, temp2,res;
        IplImage iplSrc = null;
        IplImage iplTpl = null;
        
        switch (src.getBitDepth()) {

            case 32:
                //convert source imageProcessor into iplImage
                CvMat srcMat = CvMat.create(srcH, srcW, CV_32FC1);
                double[] dArr1 = float2DtoDouble1DArray(src.getFloatArray(), srcW, srcH);
                srcMat.put(0, dArr1, 0, dArr1.length);
                iplSrc = srcMat.asIplImage();
                //iplSrc = temp.clone();
                
                //convert template imageProcessor into iplImage
                CvMat tplMat = CvMat.create(tplH, tplW, CV_32FC1);
                double[] dArr2 = float2DtoDouble1DArray(tpl.getFloatArray(), tplW, tplH);
                tplMat.put(0, dArr2, 0, dArr2.length);
                iplTpl = tplMat.asIplImage();
                //iplTpl = temp2.clone();
                
                break;
            case 16:
                //since cvMatchTemplate don't accept 16bit image, we have to convert it to 32bit
                iplSrc = cvCreateImage(cvSize(srcW, srcH), IPL_DEPTH_32F, 1);
                bi = ((ShortProcessor) src).get16BitBufferedImage();
                temp = IplImage.createFrom(bi);
                cvConvertScale(temp, iplSrc, 1 / 65535.0, 0);

                iplTpl = cvCreateImage(cvSize(tplW, tplH), IPL_DEPTH_32F, 1);
                bi2 = ((ShortProcessor) tpl).get16BitBufferedImage();
                temp2 = IplImage.createFrom(bi2);
                cvConvertScale(temp2, iplTpl, 1 / 65535.0, 0);

                temp.release();
                temp2.release();

                break;
            case 8:

                bi = src.getBufferedImage();
                iplSrc = IplImage.createFrom(bi);

                bi2 = tpl.getBufferedImage();
                iplTpl = IplImage.createFrom(bi2);

                break;
            default:
                IJ.error("Unsupported image type");
                break;
        }

        res = cvCreateImage(cvSize(srcW - tplW + 1, srcH - tplH + 1), IPL_DEPTH_32F, 1);

        /*
        CV_TM_SQDIFF        = 0,
        CV_TM_SQDIFF_NORMED = 1,
        CV_TM_CCORR         = 2,
        CV_TM_CCORR_NORMED  = 3,
        CV_TM_CCOEFF        = 4,
        CV_TM_CCOEFF_NORMED = 5;
         */

        cvMatchTemplate(iplSrc, iplTpl, res, method);
        FloatBuffer fb = res.getFloatBuffer();
        float[] f = new float[res.width() * res.height()];
        fb.get(f, 0, f.length);
        resultFp = new FloatProcessor(res.width(), res.height(), f, null);
        cvReleaseImage(res);
        
        switch (src.getBitDepth()) {
            case 32:
                iplSrc.release();
                iplTpl.release();
                break;
            case 16:
                cvReleaseImage(iplSrc);
                cvReleaseImage(iplTpl);
                break;
            case 8:
                iplSrc.release();
                iplTpl.release();
                break;
            default:
                break;    
        }

        return resultFp;
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


        for (int j = (ip.getHeight() - sWh) / 2; j < (ip.getHeight() + sWh) / 2; j++) {
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

    public static int[] findMin(ImageProcessor ip, int sW) {
        int[] coord = new int[2];
        float min = ip.getPixel(0, 0);
        int sWh, sWw;

        if (sW == 0) {
            sWh = ip.getHeight();
            sWw = ip.getWidth();
        } else {
            sWh = sW;
            sWw = sW;
        }


        for (int j = (ip.getHeight() - sWh) / 2; j < (ip.getHeight() + sWh) / 2; j++) {
            for (int i = (ip.getWidth() - sWw) / 2; i < (ip.getWidth() + sWw) / 2; i++) {
                if (ip.getPixel(i, j) < min) {
                    min = ip.getPixel(i, j);
                    coord[0] = i;
                    coord[1] = j;
                }
            }
        }
        return (coord);
    }
    
    /*
     *  The 2D array from ImageJ is [x][y], while the cvMat is in [y][x]
     */
    private static double[] float2DtoDouble1DArray(float[][] arr2d, int column, int row) {



//        IJ.log("src.col: " + column);
//        IJ.log("src.row: "+ row);
//        IJ.log("arr2d[]: "+ arr2d.length);
//        IJ.log("arr2d[][]: "+ arr2d[0].length);


        double[] arr1d = new double[column * row];
        for (int y = 0; y < row; y++) {
            for (int x = 0; x < column; x++) {
//                IJ.log("convert x:"+x+" y: "+y );
//                IJ.log("arr1d x*row+y:"+(y*column+x));
//                IJ.log("arr1d value:"+ arr2d[x][y]);
                arr1d[y * column + x] = (double) arr2d[x][y];
            }
        }

        return arr1d;
    }

    public void showAbout() {
        IJ.showMessage("cvMatch Template", "This plugin implements the tempalte matching function from\n"
                + "the OpenCV library. It will try to find an object (template)\n"
                + "within a given image (image).Six different matching algorithms\n"
                + "(methods)could be used. The matching result could be printed\n"
                + "in the log window or the result table. \n"
                + "You can also decide to display the correlation map. The coordinates\n"
                + "of the maximum (or minimum for the square difference methods)\n"
                + "correspond to the best match.\n"
                + "By checking the multimatch option, not only the best match will\n"
                + "be shown, but also all the similar pattern above the defined\n"
                + "threshold will be shown (Find maximum function on the correlation map)\n"
                + "More details on \nhttps://sites.google.com/site/qingzongtseng/template-matching-ij-plugin");
    }
}
